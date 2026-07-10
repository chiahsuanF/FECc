#!/usr/bin/env python3

import os, logging,  pysam, csv
from bisect import bisect_left, bisect_right

# -----------------------
# 3. Extract Breakpoints with MapQ  // c1, bp1 connect c2, bp2 
# -----------------------
def extract_breakpoints(bam_path, outdir, sample):
    t0 = tic("Collecting split-read breakpoints")
    bam = pysam.AlignmentFile(bam_path,'rb')
    bps = []  # (qname, c1, bp1, mq1, c2, bp2, mq2)
    for r in bam.fetch(until_eof=True):
        if not r.has_tag('SA'): continue
        blocks = r.get_blocks()
        if not blocks: continue
        bp1 = blocks[-1][1]; mq1 = r.mapping_quality; c1 = r.reference_name
        qname = r.query_name
        for sa in r.get_tag('SA').strip(';').split(';'):
            if not sa: continue
            c2, bp2, *_ = sa.split(',')[:3]
            mq2 = int(sa.split(',')[4])
            bps.append((qname, c1, bp1, mq1, c2, int(bp2), mq2))
    bam.close()

    fn = os.path.join(outdir, f"03_{sample}_breakpoints.txt")
    with open(fn,'w') as f:
        f.write("qname\tc1\tbp1\tmq1\tc2\tbp2\tmq2\n")
        for q,c1,b1,m1,c2,b2,m2 in bps:
            f.write(f"{q}\t{c1}\t{b1}\t{m1}\t{c2}\t{b2}\t{m2}\n")
    logging.info(f"Wrote {len(bps)} breakpoints → {fn}")
    toc(t0, "Collecting split-read breakpoints")
    return bps