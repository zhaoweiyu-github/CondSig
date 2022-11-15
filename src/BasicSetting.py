#!/usr/bin/env python3


# ------------------
#  python modules
# ------------------

import os,sys
import logging


# --------------
#  constants
# --------------

logging.basicConfig(level = logging.INFO,
                    format = '#%(asctime)s %(levelname)s: %(message)s', 
                    datefmt = '%a,%d %b %Y %H:%M:%S',
                    stream = sys.stderr,
                    filemode = 'w'
                    )

software_version = "v1.2.0"

# ----------------
#  Misc functions
# ----------------

error	= logging.error
warn    = logging.warning
debug   = logging.debug
info    = logging.info