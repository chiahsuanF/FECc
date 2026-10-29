#!/usr/bin/env python3

import argparse, os, sys, shutil, logging

from fec_login import *
from fec_minisup import *
from fec_amplified import *
from fec_bpmq import *
from fec_refine import *
from fec_junctiongraph import *
from fec_cycled import *
from fec_filterc import *
from fec_cycle import *
from fec_consensus import run_consensus, add_consensus_arguments
from fec_plot import plot_bigcycles_from_summary


# -----------------------

DEPTH = None                        
COPY_TOLERANCE = 0                  
REQUIRE_DEPTH_SUPPORT = False       
APPLY_DEPTH_FILTER = True           
DEPTH_FILTER_FOLD_GAP = 1.4         
MIN_JUNCTION_SUPPORT = 3             
MAX_REGION_EXTENSION = 50000        
JUNCTION_PARTNER_MQ = 0              
PARTNER_MQ_MIN_REGION_LEN = 20000    
MIN_JUNCTION_DEPTH_RATIO = 0         
MIN_PAIR_JUNCTION_RATIO = 0.1        
LOW_MAPPABILITY_DEPTH = False        
MIN_VARIANT_DEPTH_RATIO = 0.05       
LOW_DEPTH_VARIANTS = True           


# -----------------------
# Main
# -----------------------
if __name__=="__main__":
    p=argparse.ArgumentParser(description="Amplified-region -> graph -> cycles -> consensus sequence pipeline")
    p.add_argument('-b',"--bam", required=True, type=str, help='[STR] Input sorted bam, index required')
    p.add_argument('-o',"--output_dir", default='./fec_output/', type=str, help='[STR] Name of output directory')
    p.add_argument('-r',"--reference_genome", required=True, type=str, help='[STR] Import reference chromosome genome')
    p.add_argument('-ms',"--min_support", type=int, help='[INT] Minimum number of support reads for insertion reads')
    p.add_argument('-t',"--threads", type=int, default=4, help="[INT] Number of CPU threads to use")
    p.add_argument('-mc',"--max_cycle_len", type=int, default=20, help="[INT] Maximum number of fragments allowed in each eccDNA cycle")
    p.add_argument('-mq',"--min_mapping_quality", type=int, default=30, help="[INT] Minimum read mapping quality for big-cycle region depth; 0 counts multi-mapping reads and inflates depth in repeat regions")

    p.add_argument("--consensus", dest="consensus", action="store_true", default=True, help="Build each big cycle's sequence from its reads along the FEC path with bsalign POA (default)")
    p.add_argument("--no_consensus", dest="consensus", action="store_false", help="Skip the consensus step; no summary_bigcycle.fa is written")
    add_consensus_arguments(p)

    args=p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logger(args.output_dir)

    # step 9 needs bsalign; fail now rather than after the cycle search
    if args.consensus and shutil.which("bsalign") is None:
        logging.error("bsalign not found on PATH (needed for the consensus step; --no_consensus skips it); exiting"); sys.exit(1)

    sample=os.path.splitext(os.path.basename(args.bam))[0]
    logging.info(f"Starting pipeline for sample {sample}")

    # 1) Depth & MinSupport
    depth, minsupp = estimate_depth_and_minsupp( args.bam, args.min_support, DEPTH, args.output_dir, sample)
    
    # 2) Coarse amplified regions
    regs, amp_thresh, gaps = identify_amplified_regions( args.bam, args.output_dir, sample, threads=args.threads)
    if not regs:
        logging.error("No amplified regions found; exiting"); sys.exit(1)
        
    # 2b) Merge overlapping regions
    merged = merge_regions(regs)
    fnm = os.path.join(args.output_dir, f"02b_{sample}_merged_regions.txt")
    with open(fnm,'w') as f:
        f.write("chrom\tstart\tend\n")
        for c,s,e in merged:
            f.write(f"{c}\t{s}\t{e}\n")
    logging.info(f"Wrote {len(merged)} merged regions → {fnm}")

    # 3) Breakpoints
    bps = extract_breakpoints( args.bam, args.output_dir, sample)
    
    # 3a) Split regions back at bridged depth dips that carry a junction
    merged = split_regions_at_junction_gaps( merged, gaps, bps, args.output_dir, sample, min_reads=MIN_JUNCTION_SUPPORT,
        partner_mq=JUNCTION_PARTNER_MQ)

    # 3b) Extend regions whose edge stops short of an anchored junction
    if MAX_REGION_EXTENSION > 0:
        merged = extend_regions_to_anchored_breakpoints( merged, bps, args.output_dir, sample, min_reads=MIN_JUNCTION_SUPPORT,
            partner_mq=JUNCTION_PARTNER_MQ, max_extension=MAX_REGION_EXTENSION,
            partner_min_len=PARTNER_MQ_MIN_REGION_LEN)

    # 4) Refine on merged regions
    refined = refine_regions_by_highQ_breakpoints( merged, bps, args.output_dir, sample, window=200000, min_cluster_reads=MIN_JUNCTION_SUPPORT) #120,000

    # 5) Build Junction Graph
    graph = build_graph(refined, bps, args.bam, args.output_dir, sample, partner_mq=JUNCTION_PARTNER_MQ,
        partner_min_len=PARTNER_MQ_MIN_REGION_LEN)
    
    # 6) Detect cycles & write raw coords
    raw_cycles = find_cycles( graph, args.output_dir, sample, max_cycle_len=args.max_cycle_len)
    if not raw_cycles:
        logging.error("No cycles detected; exiting")
        sys.exit(1)

    # Log which nodes occur in cycles and which graph nodes are
    log_cycle_diagnostics( raw_cycles, refined, graph)

    # Write individual simple cycles.
    write_cycles_with_coords( raw_cycles, refined, args.output_dir, sample, "raw")

    # Merge overlapping cycles and add directly attached fragments.
    merged_events = write_merged_eccdna_events( cycles=raw_cycles, regions=refined, graph=graph, outdir=args.output_dir, sample=sample,
        label="raw", include_attached=True, max_hops=1, min_edge_support=1)
    
    # 7) Filter & write filtered coords
    filtered = filter_and_write_cycles( raw_cycles, refined, args.bam, args.output_dir, sample)
    write_cycles_with_coords(filtered, refined, args.output_dir, sample, "filtered")
    
    if not filtered:
        logging.error("No cycles pass quality filter; exiting"); sys.exit(1)
        
    # 8) Merge small cycles into big cycles    
    summary_file = os.path.join( args.output_dir, "summary.txt")
    bigcycle_file = os.path.join( args.output_dir, "summary_bigcycle.txt")
    
    merge_small_cycles(
    summary_file=summary_file,
    outfile=bigcycle_file,
    bam_path=args.bam,
    baseline_depth=depth,
    min_mapping_quality=args.min_mapping_quality,
    copy_tolerance=COPY_TOLERANCE,
    require_depth_support=REQUIRE_DEPTH_SUPPORT,
    apply_depth_filter=APPLY_DEPTH_FILTER,
    depth_filter_fold_gap=DEPTH_FILTER_FOLD_GAP,
    breakpoints=bps,
    min_junction_support=MIN_JUNCTION_SUPPORT,
    junction_partner_mapping_quality=JUNCTION_PARTNER_MQ,
    junction_partner_min_region_len=PARTNER_MQ_MIN_REGION_LEN,
    min_pair_junction_ratio=MIN_PAIR_JUNCTION_RATIO,
    min_junction_depth_ratio=MIN_JUNCTION_DEPTH_RATIO,
    low_depth_variants=LOW_DEPTH_VARIANTS,
    min_variant_depth_ratio=MIN_VARIANT_DEPTH_RATIO,
    low_mappability_depth=LOW_MAPPABILITY_DEPTH
)
    bigcycle_plot_dir = os.path.join(
    args.output_dir,
    "bigcycle_plots"
)
    plot_bigcycles_from_summary(
        summary_bigcycle_file=bigcycle_file,
        output_dir=bigcycle_plot_dir,
        dpi=200,
    )

    # 9) Consensus sequence per big cycle, from its reads along the FEC path
    if args.consensus:
        consensus_dir = os.path.join(args.output_dir, "consensus")
        run_consensus(bigcycle_file, args.bam, args.reference_genome, consensus_dir, args)
        shutil.copyfile(os.path.join(consensus_dir, "consensus_bigcycle.fa"), os.path.join(args.output_dir, "summary_bigcycle.fa"))
        logging.info(f"Wrote big-cycle sequences → {os.path.join(args.output_dir, 'summary_bigcycle.fa')}")
    else:
        logging.info("Skipping the consensus step (--no_consensus)")

    logging.info("Pipeline completed successfully")
