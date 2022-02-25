# CondSig

**CondSig** (**Cond**ensate-like chromatin-associated proteins co-occupancy **Sig**nature) is a comprehensive computational framework to predict component collaborations and genomic loci of potential chromatin-associated biomolecular condensates. The computation framework first detect genome-wide chromatin-associated proteins collaborations in specific loci by integrating ChIP-seq datasets and then screen out some may be involved in biomolecular condensation based on known condensation characteristics.

<p align="center">
<img src="./Image/Schematic.png"/>
</p>



## Documentation

We are hosting MAESTRO documentation, instruction and tutorials at .

## Change Log

### v1.0.0
* Release CondSig.

## System requirements
* Linux/Unix
* Python (>= 3.6) for CondSig workflow

### Install CondSig

#### Installing the CondSig workflow through conda

CondSig uses the [Anaconda3](http://conda.pydata.org/miniconda.html) package management system to harmonize all of the software packages. Users can install the full solution of MAESTRO using the conda environment.

Use the following commands to install Minicoda3：
``` bash
$ wget run.sh
$ bash run.sh
```
And then users can create an isolated environment for CondSig and install through the following commands:
``` bash
$ conda config --add channels bioconda
# To make the installation faster, we recommend using mamba
$ conda install mamba -c conda-forge
$ mamba create -n CondSig maestro=1.5.0 -c bioconda
# Activate the environment
$ conda activate CondSig
```

## Citation
-
