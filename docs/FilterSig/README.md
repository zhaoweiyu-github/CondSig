# CondSigDetector (FilterSig)



## User can type in 'condsig_detector FilterSig -h' for detailed usage. FilterSig module must be executed after [LearnSig](../LearnSig/README.md).

``` bash
$condsig_detector FilterSig -h
usage: condsig_detector FilterSig [-h] --DataAnnotation DATA_ANNOTATION
                                  --SigInputPath SIG_INPUTPATH --SigInputName
                                  SIG_INPUTNAME --LLPS LLPS --MLO MLO --PPI
                                  PPI --IDR IDR --RBD RBD [--RBS RBS]
                                  [--Threads THREADS] [--Name NAME]
                                  [--OutDir OUT_DIR]

optional arguments:
  -h, --help            show this help message and exit

Input files arguments:
  --DataAnnotation DATA_ANNOTATION
                        Data annotation file used in LearnSig. REQUIRED.
  --SigInputPath SIG_INPUTPATH
                        Input LearnSig path. REQUIRED.
  --SigInputName SIG_INPUTNAME
                        Input LearnSig name. REQUIRED.
  --LLPS LLPS           Annotation of LLPS proteins, check demo file and
                        format details in GitHub. REQUIRED.
  --MLO MLO             Annotation of MLO component, check demo file and
                        format details in GitHub. REQUIRED.
  --PPI PPI             Annotation of protein-protein interaction, check demo
                        file and format details in GitHub. REQUIRED.
  --IDR IDR             Annotation of protein IDR properties, check demo file
                        and format details in GitHub. REQUIRED.
  --RBD RBD             Annotation of RNA-binding domain content of proteins,
                        check demo file and format details in GitHub.
                        REQUIRED.
  --RBS RBS             Big wiggle track of genome-wide RNA-binding density,
                        check demo file and format details in GitHub. Default:
                        None.
  --Threads THREADS     The processes used of the job. Default:8.

Output files arguments:
  --Name NAME           Output name. Default: Input LearnSig name.
  --OutDir OUT_DIR      Output directory. Default: Input LearnSig path.

```

## Run FilterSig module

### Demo for mESC
```bash
condsig_detector FilterSig \
--DataAnnotation mESC_dataset_meta.txt \
--SigInputPath mESC \
--SigInputName mESC \
--Threads 10 \
--LLPS LLPS_mouseExtended.txt \
--MLO MLO_mouseExtended_merged.txt \
--PPI PPI_mouse_merged.txt \
--IDR mouse_MobiDB_lite.txt \
--RBD mouse_TriPepSVM_RBP.txt
```

### Demo for K562
```bash
condsig_detector FilterSig \
--DataAnnotation K562_dataset_meta.txt \
--SigInputPath K562 \
--SigInputName K562 \
--Threads 10 \
--LLPS LLPS_human_merged.txt \
--MLO MLO_human_merged.txt \
--PPI PPI_human_merged.txt \
--IDR human_MobiDB_lite.txt \
--RBD human_TriPepSVM_RBP.txt \
--RBS K562-D210N-V5ChIP_merged.bw
```

## Notes

### Annotation files
All annotation files and source codes to generate annotation file in demo are listed below:
* mESC
  * Annotation of LLPS proteins ([LLPS_mouseExtended.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/LLPS_mouseExtended.txt), [source code](../../example/generate_LLPS_annotation.ipynb)).
  * Annotation of MLO component ([MLO_mouseExtended_merged.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/MLO_mouseExtended_merged.txt), [source code](../../example/generate_MLO_annotation.ipynb)).
  * Annotation of protein-protein interaction ([PPI_mouse_merged.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/PPI_mouse_merged.txt), [source code](../../example/generate_PPI_annotation.ipynb)).
  * Annotation of protein IDR properties ([mouse_MobiDB_lite.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/mouse_MobiDB_lite.txt), [source code](../../example/generate_IDR_annotation.ipynb)).
  * Annotation of RNA-binding domain content of proteins ([mouse_TriPepSVM_RBP.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/mouse_TriPepSVM_RBP.txt), [source code](../../example/generate_RBP_annotation.ipynb)).
* K562
  * Annotation of LLPS proteins ([LLPS_mouseExtended.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/LLPS_mouseExtended.txt), [source code](../../example/generate_LLPS_annotation.ipynb)).
  * Annotation of MLO component ([MLO_mouseExtended_merged.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/MLO_mouseExtended_merged.txt), [source code](../../example/generate_MLO_annotation.ipynb)).
  * Annotation of protein-protein interaction ([PPI_mouse_merged.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/PPI_mouse_merged.txt), [source code](../../example/generate_PPI_annotation.ipynb)).
  * Annotation of protein IDR properties ([mouse_MobiDB_lite.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/mouse_MobiDB_lite.txt), [source code](../../example/generate_IDR_annotation.ipynb)).
  * Annotation of RNA-binding domain content of proteins ([mouse_TriPepSVM_RBP.txt](https://compbio-zhanglab.org/CondSigDB/data/GitHub/mouse_TriPepSVM_RBP.txt), [source code](../../example/generate_RBP_annotation.ipynb)).
  * Big wiggle track of genome-wide RNA-binding density ([K562-D210N-V5ChIP_merged.bw](https://compbio-zhanglab.org/CondSigDB/data/GitHub/K562-D210N-V5ChIP_merged.bw)).

### Development
Only condensation-related feature annotations for CAPs, for which ChIP-seq data is available in the Cistrome Data Browser, are provided here. Users can generate annotation files for additional CAPs by following the source codes provided above (the correct format of annotaiton files is requried).