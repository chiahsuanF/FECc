#!/usr/bin/env python3

import os, time, logging


# -----------------------
# Logger setup
# -----------------------
def setup_logger(output_dir):
    logfile = os.path.join(output_dir, "login_note.txt")
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(ch)
    logging.info("Logger initialized")

def tic(msg):
    logging.info(f"START: {msg}")
    return time.time()

def toc(start, msg):
    logging.info(f" END: {msg} ({time.time()-start:.1f}s)")

def junction_mq_passes(mq1, mq2, min_mq, partner_mq):
    return max(mq1, mq2) >= min_mq and min(mq1, mq2) >= partner_mq

def partner_side_allowed(mq1, mq2, min_mq, len1, len2, partner_min_len):
    if min(mq1, mq2) >= min_mq or partner_min_len <= 0:
        return True
    return min(len1, len2) >= partner_min_len
