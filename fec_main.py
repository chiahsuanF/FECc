#!/usr/bin/env python3

import argparse, os, sys,logging
import fec_login, fec_minisup, fec_amplified, fec_bpmq, fec_refine, fec_junctiongraph, fec_cycled, fec_filterc, fec_poa, fec_fa

# -----------------------
# Main
# -----------------------
if __name__=="__main__":
    p=argparse.ArgumentParser(description="Amplified-region → graph → cycles → POA pipeline")
    p.add_argument('-b',"--bam", required=True, type=str, help='[STR] Input sorted bam, index required')
    p.add_argument('-o',"--output_dir", required=True, type=str, help='[STR] Name of output directory')
    p.add_argument('-r',"--reference_genome", required=True, type=str, help='[STR] Import reference chromosome genome')
    p.add_argument('-ms',"--min_support", type=int, help='[INT] Minimum number of support reads for insertion reads')
    p.add_argument('-d',"--depth", type=float, help='[FLOAT]Approximately total depth')
    p.add_argument('-t',"--threads", type=int, default=4, help="[INT] Number of CPU threads to use")

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
    fnm = os.path.join(args.output_dir, f"10_{sample}_merged_regions.txt")
    with open(fnm,'w') as f:
        f.write("chrom\tstart\tend\n")
        for c,s,e in merged:
            f.write(f"{c}\t{s}\t{e}\n")
    logging.info(f"Wrote {len(merged)} merged regions → {fnm}")

    # 3) Breakpoints
    bps = extract_breakpoints( args.bam, args.output_dir, sample)
    
    # 4) Refine on merged regions
    refined = refine_regions_by_highQ_breakpoints( merged, bps, args.output_dir, sample, window=100000) #1000

    # 5) Build Junction Graph
    graph = build_graph(refined, bps, amp_thresh, args.output_dir, sample)
    
    # 6) Detect cycles & write raw coords
    raw_cycles = find_cycles( graph, args.output_dir, sample)
    write_cycles_with_coords(raw_cycles, refined, args.output_dir, sample, "raw")
    if not raw_cycles:
        logging.error("No cycles detected; exiting"); sys.exit(1)
        
    # 7) Filter & write filtered coords
    filtered = filter_and_write_cycles( raw_cycles, refined, args.bam, args.output_dir, sample)
    write_cycles_with_coords(filtered, refined, args.output_dir, sample, "filtered")
    if not filtered:
        logging.error("No cycles pass quality filter; exiting"); sys.exit(1)
        
    # 8) POA & refined connections
    junctions = build_unique_junctions(filtered, refined)

    results = process_junctions_fast(
        junctions,
        args.bam,
        args.reference_genome,
        os.path.join(args.output_dir, "poa"),
        threads=args.threads)
    
    build_cycle_fasta( filtered, refined, args.reference_genome, args.output_dir, sample)
    
    logging.info("Pipeline completed successfully")

    os.replace('07_merge.v4.1_filtered_cycles_coords.txt', 'summary.txt')
    os.system('rm 0*.txt ')
