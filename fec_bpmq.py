#!/usr/bin/env python3

import os, logging, pysam, re
from fec_login import *

# -----------------------
# 3. Extract Breakpoints with MapQ  // c1, bp1 connect c2, bp2 
# -----------------------
def extract_breakpoints(bam_path, outdir, sample):
    t0 = tic("Collecting split-read breakpoints")
    bam = pysam.AlignmentFile(bam_path,'rb')
    bps = []  # (qname, c1, bp1, mq1, strand1, c2, bp2, mq2, strand2)
    for r in bam.fetch(until_eof=True):
        if not r.has_tag('SA'): continue
        blocks = r.get_blocks()
        if not blocks: continue
        bp1 = blocks[-1][1]
        mq1 = r.mapping_quality
        c1 = r.reference_name
        qname = r.query_name
        
        sa_entries = r.get_tag("SA").strip(";").split(";")
        
        for sa in sa_entries:

            if not sa:continue

            fields = sa.split(",")

            c2 = fields[0]
            pos2 = int(fields[1])
            strand2 = fields[2]
            cigar2 = fields[3]
            mq2 = int(fields[4])

            bp2 = breakpoint_from_sa( pos2, strand2, cigar2)
            
            strand1 = "-" if r.is_reverse else "+"
            bps.append(( qname, c1, bp1, mq1, strand1, c2, bp2, mq2, strand2 ))
    bam.close()

    fn = os.path.join(outdir, f"03_{sample}_breakpoints.txt")
    with open(fn,'w') as f:
        f.write("qname\tc1\tbp1\tmq1\tstrand1\tc2\tbp2\tmq2\tstrand2\n")
        for q,c1,b1,m1,s1,c2,b2,m2,s2 in bps:
            f.write(f"{q}\t{c1}\t{b1}\t{m1}\t{s1}\t{c2}\t{b2}\t{m2}\t{s2}\n")
    logging.info(f"Wrote {len(bps)} breakpoints → {fn}")
    toc(t0, "Collecting split-read breakpoints")
    return bps

def cigar_ref_len(cigar):
    ref_len = 0

    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):

        if op in ('M', 'D', 'N', '=', 'X'):
            ref_len += int(n)

    return ref_len

def breakpoint_from_sa(pos2, strand2, cigar2):

    left_soft = 0
    right_soft = 0

    m = re.match(r'^(\d+)S', cigar2)
    if m:
        left_soft = int(m.group(1))

    m = re.search(r'(\d+)S$', cigar2)
    if m:
        right_soft = int(m.group(1))

    ref_len = cigar_ref_len(cigar2)

    if strand2 == "+":
        return pos2 if left_soft > right_soft else pos2 + ref_len
    else:
        return pos2 + ref_len if left_soft > right_soft else pos2
