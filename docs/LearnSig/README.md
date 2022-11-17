# CondSigDetector (LearnSig)



User can type in 'condsig_detector LearnSig -h' for detailed usage.

``` bash
$ condsig_detector LearnSig -h

usage: condsig_detector LearnSig [-h] --DataAnnotation DATA_ANNOTATION
                                 [--Mode MODE] [--Focus FOCUS]
                                 [--FocusNumber FOCUS_NUMBER]
                                 [--MinSignatures MIN_SIGNATURES]
                                 [--MaxSignatures MAX_SIGNATURES]
                                 [--Threads THREADS]
                                 [--GenomeVersion GENOME_VERSION]
                                 [--Name NAME] [--OutDir OUT_DIR]
                                 [--Zscore ZSCORE]

optional arguments:
  -h, --help            show this help message and exit

Input files arguments:
  --DataAnnotation DATA_ANNOTATION
                        Annotation file of CAPs in the same cell type. 4
                        columns of annotation file are factor, label, peak
                        file directory, uniprot id with tab delimited.
                        Example: CTCFL K562_GSM803401_CTCFL
                        K562_GSM803401_CTCFL.bed Q8NI51. REQUIRED.
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
                        The genome version of peak file. Default:hg38.

Output files arguments:
  --Name NAME           Output name. Default:Test
  --OutDir OUT_DIR      Output directory (absolute path required). Default:.
  --Zscore ZSCORE       Z-score threshold for identifying component CAPs of
                        signatures. Default:1.
```

Run LearnSig module

```ba
condsig_detector LearnSig \
--DataAnnotation mESC_selected_dataset_meta_220922.txt \
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

Notes

​	

