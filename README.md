# CondSig

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version: 1.0.0](https://img.shields.io/badge/Version-1.0.0-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)

**CondSig** (**Cond**ensate-like chromatin-associated proteins co-occupancy **Sig**nature) is a comprehensive computational framework to predict component CAPs (chromatin-associated proteins) and genomic loci of potential chromatin-associated biomolecular condensates. The computation framework first detect genome-wide CAP collaborations in specific loci by integrating ChIP-seq datasets and then screen out some may be involved in biomolecular condensation based on known condensation characteristics.

<p align="center">
<img src="./image/Schematic.png"/>
</p>


## Change Log

### v1.0.0
* Release CondSig.

## System requirements
* Linux/Unix

## Install CondSig

### Installing the CondSig workflow through conda

CondSig uses the [Anaconda3](http://conda.pydata.org/miniconda.html) package management system to harmonize all of the software packages. Users can install the full solution of CondSig using the conda environment.

Use the following commands to install Minicoda3：
``` bash
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```
And then users can create an isolated environment for CondSig and install through the following commands:
``` bash
# Create environment for CondSig
$ conda env crate -n CondSig_env python=3.7
# Activate the environment
$ conda activate CondSig_env
# Install CondSig
$ conda install -c yuzhaowei condsig_alpha
```

### Test CondSig

```bash
# Test CondSig command
condsig LearnSig --help
condsig FilterSig --help
```

## Citation

-
