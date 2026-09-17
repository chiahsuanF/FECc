#!/usr/bin/env python3

import argparse, os, sys,logging

from fec_login import *
from fec_minisup import *
from fec_amplified import *
from fec_bpmq import *
from fec_refine import *
from fec_junctiongraph import *
from fec_cycled import *
from fec_filterc import *
from fec_poa import *
from fec_cycle import *
from fec_fa import *
from fec_plot import plot_bigcycles_from_summary



# -----------------------
# Main
# -----------------------
if __name__=="__main__":
    p=argparse.ArgumentParser(description="Amplified-region -> graph -> cycles -> POA pipeline")
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

    args=p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logger(args.output_dir)

    sample=os.path.splitext(os.path.basename(args.bam))[0]
    logging.info(f"Starting pipeline for sample {sample}")

    # 1) Depth & MinSupport
    depth, minsupp = estimate_depth_and_minsupp( args.bam, args.min_support, args.depth, args.output_dir, sample)
    
    # 2) Coarse amplified regions
    regs, amp_thresh = identify_amplified_regions( args.bam, args.output_dir, sample)
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
    
    # 4) Refine on merged regions
    refined = refine_regions_by_highQ_breakpoints( merged, bps, args.output_dir, sample, window=200000) #120,000

    # 5) Build Junction Graph
    graph = build_graph(refined, bps, args.output_dir, sample, depth)
    
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
        
    # 8) POA & refined connections    
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
    depth_filter_fold_gap=args.depth_filter_fold_gap
)
    bigcycle_plot_dir = os.path.join(
    args.output_dir,
    "bigcycle_plots"
)
    junctions = build_junctions_from_bigcycles(
        bigcycle_file
    )

    results = process_junctions_fast(
        junctions,
        args.bam,
        args.reference_genome,
        os.path.join(args.output_dir, "poa"),
        threads=args.threads
    )

    plot_bigcycles_from_summary(
        summary_bigcycle_file=bigcycle_file,
        output_dir=bigcycle_plot_dir,
        dpi=200,
    )
    
    build_bigcycle_fasta(
        summary_bigcycle_file=bigcycle_file,
        reference_genome=args.reference_genome,
        output_dir=args.output_dir,
        sample=sample
    )    
    logging.info("Pipeline completed successfully")
