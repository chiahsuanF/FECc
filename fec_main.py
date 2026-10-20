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
# Main
# -----------------------
if __name__=="__main__":
    p=argparse.ArgumentParser(description="Amplified-region -> graph -> cycles -> consensus sequence pipeline")
    p.add_argument('-b',"--bam", required=True, type=str, help='[STR] Input sorted bam, index required')
    p.add_argument('-o',"--output_dir", default='./fec_output/', type=str, help='[STR] Name of output directory')
    p.add_argument('-r',"--reference_genome", required=True, type=str, help='[STR] Import reference chromosome genome')
    p.add_argument('-ms',"--min_support", type=int, help='[INT] Minimum number of support reads for insertion reads')
    p.add_argument('-d',"--depth", type=float, help='[FLOAT]Approximately total depth')
    p.add_argument('-t',"--threads", type=int, default=4, help="[INT] Number of CPU threads to use")
    p.add_argument('-mc',"--max_cycle_len", type=int, default=20, help="[INT] Maximum number of fragments allowed in each eccDNA cycle")
    p.add_argument('-ct',"--copy_tolerance", type=int, default=0, help="[INT] Extra region traversals allowed beyond the depth-estimated copy number")
    p.add_argument("--require_depth_support", action="store_true", help="Drop big cycles whose region usage exceeds the depth-estimated copy number, instead of only marking them")
    p.add_argument('-mq',"--min_mapping_quality", type=int, default=30, help="[INT] Minimum read mapping quality for big-cycle region depth; 0 counts multi-mapping reads and inflates depth in repeat regions")
    p.add_argument("--no_depth_filter", action="store_true", help="Merge every cycle of a group, instead of holding back the cycles that visit a region below the group's depth cutoff")
    p.add_argument("--depth_filter_fold_gap", type=float, default=1.4, help="[[FLOAT] Smallest fold gap in sorted region depths that may be used as the cutoff; the widest qualifying gap separates the amplified core from the low-depth backgroun")

    p.add_argument('-mjs',"--min_junction_support", type=int, default=3, help="[INT] Split reads a segment-end junction needs before a big cycle may use it")
    p.add_argument("--max_region_extension", type=int, default=50000, help="[INT] Furthest (bp) a region edge may be extended to a breakpoint that anchored split reads place beyond it; 0 turns extension off")
    p.add_argument('-jpq',"--junction_partner_mq", type=int, default=0, help="[INT] One side of a junction split read must reach MQ 30; the other side only needs this, so junctions into repeats are kept. Applies only between regions of at least --partner_mq_min_region_len; 30 requires both sides everywhere")
    p.add_argument("--partner_mq_min_region_len", type=int, default=20000, help="[INT] A junction read unique on one side only is counted only when both regions it joins are at least this long (bp); short pieces need both sides at MQ 30. 0 allows any length")
    p.add_argument("--min_junction_depth_ratio", type=float, default=0, help="[FLOAT] A big-cycle junction also needs split reads of at least this fraction of the depth of its deeper flanking region; drops few-read chimeric joints in deep amplicons. 0 turns it off")
    p.add_argument("--min_pair_junction_ratio", type=float, default=0.1, help="[FLOAT] A 2-fragment big cycle whose weaker junction has fewer split reads than this fraction of its median region depth goes to summary_bigcycle_weak_pairs.txt instead; 0 keeps them all")
    p.add_argument("--low_mappability_depth", action="store_true", help="Use MQ 0 depth for structure in regions where most reads are multi-mapped (MQ >= 30 depth below 0.8x the MQ 0 depth); off by default")
    p.add_argument("--min_variant_depth_ratio", type=float, default=0.05, help="[FLOAT] A held-back low-depth region is added back as a variant only when its depth is at least this fraction of the median depth of the core regions it joins; 0 adds every one")
    p.add_argument("--no_low_depth_variants", action="store_true", help="Do not add variant big cycles that carry low-depth regions held back by the depth filter (summary_bigcycle_variants.txt)")
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
    depth, minsupp = estimate_depth_and_minsupp( args.bam, args.min_support, args.depth, args.output_dir, sample)
    
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
    merged = split_regions_at_junction_gaps( merged, gaps, bps, args.output_dir, sample, min_reads=args.min_junction_support,
        partner_mq=args.junction_partner_mq)

    # 3b) Extend regions whose edge stops short of an anchored junction
    if args.max_region_extension > 0:
        merged = extend_regions_to_anchored_breakpoints( merged, bps, args.output_dir, sample, min_reads=args.min_junction_support,
            partner_mq=args.junction_partner_mq, max_extension=args.max_region_extension,
            partner_min_len=args.partner_mq_min_region_len)

    # 4) Refine on merged regions
    refined = refine_regions_by_highQ_breakpoints( merged, bps, args.output_dir, sample, window=200000, min_cluster_reads=args.min_junction_support) #120,000

    # 5) Build Junction Graph
    graph = build_graph(refined, bps, args.bam, args.output_dir, sample, partner_mq=args.junction_partner_mq,
        partner_min_len=args.partner_mq_min_region_len)
    
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
    copy_tolerance=args.copy_tolerance,
    require_depth_support=args.require_depth_support,
    apply_depth_filter=not args.no_depth_filter,
    depth_filter_fold_gap=args.depth_filter_fold_gap,
    breakpoints=bps,
    min_junction_support=args.min_junction_support,
    junction_partner_mapping_quality=args.junction_partner_mq,
    junction_partner_min_region_len=args.partner_mq_min_region_len,
    min_pair_junction_ratio=args.min_pair_junction_ratio,
    min_junction_depth_ratio=args.min_junction_depth_ratio,
    low_depth_variants=not args.no_low_depth_variants,
    min_variant_depth_ratio=args.min_variant_depth_ratio,
    low_mappability_depth=args.low_mappability_depth
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
