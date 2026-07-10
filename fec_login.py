#!/usr/bin/env python3

import  os, time, logging
from bisect import bisect_left, bisect_right

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