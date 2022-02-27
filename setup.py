#!/usr/bin/env python3
# Time-stamp: <2022-02-26 Zhaowe Yu>

"""
Description: CondSig v1.0.0
Copyright (c) 2022 Zhaowei Yu <zhaoweiyu@tongji.edu.cn>

@status: first release 
@version: v1.0.0
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
		name = "condsig",
		version = "1.0.0",
		package_dir = {"CondSig":"src"},
		packages = ["CondSig"],
		scripts = ["bin/condsig", "bin/biterm_local.sh"]
		)

if __name__ == "__main__":
    main()