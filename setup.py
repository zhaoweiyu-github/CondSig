#!/usr/bin/env python3
# Time-stamp: <2022-11-08 Zhaowei Yu>

"""
Description: CondSig detector v1.2.0
Copyright (c) 2022 Zhaowei Yu <zhaoweiyu@tongji.edu.cn>

@status: first release 
@version: v1.2.0
@author:  Zhaowei Yu
@contact: zhaoweiyu@tongji.edu.cn
"""

import sys,os
try:
	from setuptools import setup, find_packages
except ImportError:
	print("Could not load setuptools. Please install the setuptools package.")

def main():
	setup(
		name = "CondSig",
		version = "1.2.0",
		author = "Zhaowei Yu",
		author_email = "zhaoweiyu@tongji.edu.cn",
		description = "A computational framework to detect Condensate-like chromatin-associated protein co-occupancy Signatures.",
		package_dir = {"CondSig":"src"},
		package_data = {"CondSig":["resource/*"]},
		packages = ["CondSig"],
		scripts = ["bin/condsig_detector", "bin/biterm_run", "bin/indexDocs.py"],
		include_package_data = True
		)

if __name__ == "__main__":
    main()
