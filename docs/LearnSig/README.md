# CondSigDetector (LearnSig)



## User can type in 'condsig_detector LearnSig -h' for detailed usage.

``` bash
$condsig_detector LearnSig -h
usage: condsig_detector LearnSig [-h] --DataAnnotation DATA_ANNOTATION
                                 --DataDir DATA_DIR [--Mode MODE]
                                 [--Focus FOCUS] [--FocusNumber FOCUS_NUMBER]
                                 [--MinSignatures MIN_SIGNATURES]
                                 [--MaxSignatures MAX_SIGNATURES]
                                 [--Threads THREADS] --GenomeVersion
                                 GENOME_VERSION [--Name NAME]
                                 [--OutDir OUT_DIR] [--Zscore ZSCORE]

optional arguments:
  -h, --help            show this help message and exit

Input files arguments:
  --DataAnnotation DATA_ANNOTATION
                        Annotation file of CAPs in the same cell type. 4
                        columns of annotation file are factor, label, peak
                        file name, uniprot id with tab delimited. Example:
                        CTCFL K562_GSM803401_CTCFL K562_GSM803401_CTCFL.bed
                        Q8NI51. REQUIRED.
  --DataDir DATA_DIR    Input Peak files directory. REQUIRED.
  --Mode MODE           Mode of detecting CAP co-occupancy signatures. 'all'
                        means detecting co-occupancy signatures for each CAP
                        iteatively and 'focus' means detecting co-occupancy
                        signatures for focus CAP. Default:all.
  --Focus FOCUS         Focus CAP, avaiable when focus mode is set.
  --FocusNumber FOCUS_NUMBER
                        The number of focused potential combinatorial CAPs.
                        Default:50.
  --MinSignatures MIN_SIGNATURES
                        The minimum number of signatures. Default:2.
  --MaxSignatures MAX_SIGNATURES
                        The maximum number of signatures. Default:10.
  --Threads THREADS     The processes used of the job. Default:8.
  --GenomeVersion GENOME_VERSION
                        The UCSC genome version,only hg38 and mm10 is avaiable
                        now. Default:hg38.

Output files arguments:
  --Name NAME           Output name. Default:Test
  --OutDir OUT_DIR      Output directory (absolute path required). Default:.
  --Zscore ZSCORE       Z-score threshold for identifying component CAPs of
                        signatures. Default:1.3.
```

## Run LearnSig module

### Demo for mESC
```bash
condsig_detector LearnSig \
--DataAnnotation mESC_dataset_meta.txt \
--DataDir mESC_peaks \
--Mode all \
--FocusNumber 50 \
--MinSignatures 2 \
--MaxSignatures 10 \
--Threads 10 \
--GenomeVersion mm10 \
--Name mESC \
--OutDir mESC \
--Zscore 1.3
```

### Demo for K562
```bash
condsig_detector LearnSig \
--DataAnnotation K562_dataset_meta.txt \
--DataDir K562_peaks \
--Mode all \
--FocusNumber 50 \
--MinSignatures 2 \
--MaxSignatures 10 \
--Threads 10 \
--GenomeVersion hg38 \
--Name K562 \
--OutDir K562 \
--Zscore 1.3
```

## Notes

​The ChIP-seq data of CAPs were collected from Cistrome Data Browser and filtrated using a stringent quality control procedures. Datasets and meta data for mESC([mESC_peaks.tar.gz](https://compbio-zhanglab.org/CondSigDB/data/GitHub/mESC_peaks.tar.gz), [mESC_dataset_meta.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/mESC_dataset_meta.txt)) and K562([K562_peaks.tar.gz](https://compbio-zhanglab.org/CondSigDB/data/GitHub/K562_peaks.tar.gz), [K562_dataset_meta.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/K562_dataset_meta.txt)) were released.

It should be recognized that the presence of random processes in biterm topic model may introduce a slight bias to the outcomes of each trail.

