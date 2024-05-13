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

### Output
FilterSig for mESC or K562 is expected to run for approximately 8 or 5 hours on a high-performance computing cluster equipped with 10 CPUs (Central Processing Units) and a large amount of RAM (Random Access Memory) using the parameter settings described above. 

All identified CondSigs and associated genomic loci can be found within "FilterSig/Summary" folder. This folder is located under the output directory of "FilterSig" module.CondSigDetector produces two distinct summary TXT files for promoter and non-promoter CondSigs separately. For instance, you can find files named "mESC_promoter_CondSigs.txt" and "mESC_nonpromoter_CondSigs.txt" after running the demo for mESC, and each summary TXT file contains five columns:
  1. The name of CondSig
  2. Component CAPs
  3. Condensate-like features of the CondSig
  4. The count of condensate-like features
  5. Mean AUROC

```
CondSig	component_CAP	qualified_CL_features	qualified_CL_features_count	mean_AUROC
mESC_rep3_promoter_CondSig_1	mESC_GSM930151_KDM2B,mESC_GSM2460999_KDM4C,mESC_GSM1341311_TDG,mESC_GSM1603269_RXRA,mESC_GSM3196078_TET2	IDR,LLPS,MLO,PPI,RBP	5	0.6678
mESC_rep3_promoter_CondSig_2	mESC_GSM930151_KDM2B,mESC_GSM2192644_AEBP2,mESC_GSM1693794_HEXIM1,mESC_GSM2460999_KDM4C,mESC_GSM651192_DPY30,mESC_GSM1372576_ATRX,mESC_GSM2142337_FAM60A,mESC_GSM1399511_SUZ12	IDR,LLPS,MLO,PPI,RBP	5	0.6538
mESC_rep3_promoter_CondSig_3	mESC_GSM1355155_POU5F1,mESC_GSM687282_KDM1A,mESC_GSM2417143_SOX2,mESC_GSM1563242_RAD23B,mESC_GSM2123560_NANOG,mESC_GSM1208218_KLF5,mESC_GSM2588408_KMT2D,mESC_GSM1258240_ASH2L	IDR,LLPS,MLO,PPI	4	0.7625
```

Users can checkout all identified CondSigs in mESC and K562 in our database ([CondSigDB](https://compbio-zhanglab.org/CondSigDB/index.html)).

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
Only condensation-related feature annotations for CAPs, for which ChIP-seq data are available in the Cistrome Data Browser, are provided here. Users can generate annotation files for additional CAPs by following the source codes provided above (the correct format of annotation files is required).