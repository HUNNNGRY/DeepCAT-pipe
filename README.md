# Detect Cancer Associated TCR in Blood

This chapter focuses on peripheral blood mononuclear cell (**PBMC**) or tissue samples' TCR-seq or bulk RNA-seq data for cancer sample classification, including TCR calling, filtering, clustering, prediction et al.

Tumor microenvironment immunity plays a key role in tumorigenesis and progression, and in recent years, with the rapid development of liquid biopsy, various signals from the tumor and the surrounding immune environment can be detected in peripheral blood, including reactive T cells that can recognize tumor **neoantigens**, which makes it possible to directly distinguish **non-invasive diagnosis** from normal people and people with asymptomatic tumors in the early stage based on peripheral blood samples. 

Here is an analytical workflow that integrates multiple published methods is provided to detect potential small amounts of tumor-immune environment signals represented by **caTCR** (cancer-associated TCR) in peripheral blood cells based on TCR-seq or bulk RNA-seq.


## 1) Pipeline
![figs/TCR-pipeline.png](figs/TCR-pipeline.png)

## 2) Data Structure
### 2a) getting software & data
```bash
#download test repo from gh
git clone https://github.com/HUNNNGRY/TCR.git 
```

**data** ([download](https://cloud.tsinghua.edu.cn/f/c255750539384dddb98b/?dl=1))

1.TCR-seq dir
```./TCR/sample_data/PBMC_TCR-seq```

2.lulab RNA-seq dir
```./TCR/sample_data/PBMC_RNA-seq ```


**software/env requirement**

1. anaconda3

```bash
conda create -n TCR python==3.7 biopython tensorflow==1.14 matplotlib scikit-learn -y
# activate TCR env
source activate TCR
```

2. TRUST4

```bash
cd ./TCR/packages/TRUST4/
make
cd -
cd ./TCR
```

environment configuration see 
[https://github.com/liulab-dfci/TRUST4](https://github.com/liulab-dfci/TRUST4)

3. DeepCAT

[https://github.com/HUNNNGRY/DeepCAT](https://github.com/HUNNNGRY/DeepCAT)

original repo: [https://github.com/s175573/DeepCAT](https://github.com/s175573/DeepCAT)

4. iSMART

original repo: [https://github.com/s175573/iSMART](https://github.com/s175573/iSMART)

### 2b) input
* option1: pre-processed tsv TCR-seq file from adaptativebiotech

* option2: raw fastq file from PBMC RNA seq

### 2c) output
* cancer score: DeepCAT CNN model prediction
* visualization: box plot and ROC curve 

## 3) Running Steps
Enter into bash interface, make sure you are in TCR dir and under TCR conda env
Option1: Using PBMC TCR-seq as input

1. view TCR_input data
```bash
# check header
head ./sample_data/PBMC_TCR-seq/HIP09046.tsv 
# count row number
wc -l ./sample_data/PBMC_TCR-seq/control/HIP09046.tsv
```

2. make working directory tree 
```bash
mkdir -p ./test_PBMC_TCR-seq/{01_filter_output,02_cluster_output,03_deepcat_output}
mkdir -p ./test_PBMC_TCR-seq/{01_filter,02_cluster}_output/{disease,control}
```

3. filter/prepare input  
```bash
# filter invalid TCR sequence, using TCRB V region CDR3 sequence only
python ./PrepareAdaptiveFile_corrected.py ./sample_data/PBMC_TCR-seq/disease   ./test_PBMC_TCR-seq/01_filter_output/disease
python ./PrepareAdaptiveFile_corrected.py ./sample_data/PBMC_TCR-seq/control   ./test_PBMC_TCR-seq/01_filter_output/control

# keep an eye on the change of row/record number 
wc -l ./test_PBMC_TCR-seq/01_filter_output/control/TestReal-HIP09046.tsv
```

4. cluster using iSMART

```bash
# cluster similar TCR sequences using iSMART
python ./iSMARTv3.py -d ./test_PBMC_TCR-seq/01_filter_output/disease -o ./test_PBMC_TCR-seq/02_cluster_output/disease

python ./iSMARTv3.py -d ./test_PBMC_TCR-seq/01_filter_output/control -o ./test_PBMC_TCR-seq/02_cluster_output/control

# keep an eye on the change of row/record number 
wc -l ./test_PBMC_TCR-seq/02_cluster_output/control/TestReal-HIP09046.tsv_ClusteredCDR3s_7.5.txt
```

5. predict cancer score(probability) using DeepCAT

```bash
# predict cancer score in two groups
bash  ./Script_DeepCAT.sh -t ./test_PBMC_TCR-seq/02_cluster_output/disease 
bash  ./Script_DeepCAT.sh -t ./test_PBMC_TCR-seq/02_cluster_output/control
# note:
# do not add “/” the in end of cmd

# check output cancer score
head ./Cancer_score_control.txt
head ./Cancer_score_disease.txt

# move output to destination dir
mv ./Cancer_score_{control,disease}.txt ./test_PBMC_TCR-seq/03_deepcat_output
```

6. visualize cancer score result

```bash
# make boxplot and ROC curve using Rscripts
Rscript ./plot.R ./test_PBMC_TCR-seq/03_deepcat_output
# note:
# you might meet some errors, mostly lack specific packages, just install them as instructed
# eg.  install.packages("ROCR")
# the final results are generated under ./test_PBMC_TCR-seq/03_deepcat_output
```

* x-axis means different group (cancer vs. normal), y-axis means cancer scores from all samples
* ROC is used for evaluation of classification performance (near topleft, better effect) 


### Option2: Using PBMC RNA-seq as input
Overall similar with Option1, the main differences:

* need extract TCR sequences from bulk RNA-seq using TRUST4
* bulk RNA-seq is non-targeted sequencing, available TCR records can be much more fewer, thus we will skip clustering step, and predict using DeepCAT directly

1. make working directory tree 

```bash
mkdir -p ./test_PBMC_RNA-seq/{01_TCRcalling_output,02_filter_output,03_deepcat_output}
mkdir -p ./test_PBMC_RNA-seq/{01_TCRcalling,02_filter}_output/{disease,control}
```

2. TCR calling

```bash
# extract TCR records from fastq(.gz)
./packages/TRUST4/run-trust4  -1 controlID_1.fastq.gz -2 controlID_2.fastq.gz -f ./reference/TRUST4/hg38_bcrtcr.fa --ref ./reference/TRUST4/human_IMGT+C.fa -t 4 -o ./test_PBMC_RNA-seq/01_TCRcalling_output/control/controlID

# PE raw files as e.g., we use for loop to get TRUST4 TCR calling records of all samples:
for idx in `cat list.txt`;
do 
./packages/TRUST4/run-trust4  -1 ${idx}_1.fastq.gz -2 ${idx}_2.fastq.gz -f ./reference/TRUST4/hg38_bcrtcr.fa --ref ./reference/TRUST4/human_IMGT+C.fa -t 4 -o ./test_PBMC_RNA-seq/01_TCRcalling_output/${idx} ;
done

# view TRUST4 TCR calling output data
head ./test_PBMC_RNA-seq/01_TCRcalling_output/disease/CRC-2415350_report.tsv
wc -l ./test_PBMC_RNA-seq/01_TCRcalling_output/disease/CRC-2415350_report.tsv
```

> Note: need some time to finish this step, consider to use pre-calculated TRUST4 TCR calling files: move files under ```./sample_data/PBMC_RNA-seq/``` to ```./test_PBMC_RNA-seq/01_TCRcalling_output/```

3. filter/prepare input  

```bash
# filter invalid TCR sequence, using TCRB V region CDR3 sequence only
# For loop to filter cancer and ctrl samples (need prepare sampleID files in advance)
for idx in `cat ./test_PBMC_RNA-seq/01_TCRcalling_output/diseaseID.txt`
do
Rscript ./filter_TRUST4.R ./test_PBMC_RNA-seq/01_TCRcalling_output/disease   ./test_PBMC_RNA-seq/02_filter_output/disease ${idx}
done

for idx in `cat ./test_PBMC_RNA-seq/01_TCRcalling_output/controlID.txt`
do
Rscript ./filter_TRUST4.R ./test_PBMC_RNA-seq/01_TCRcalling_output/control   ./test_PBMC_RNA-seq/02_filter_output/control ${idx}
done

# keep an eye on the change of row/record number 
wc -l ./test_PBMC_RNA-seq/02_filter_output/disease/CRC-2415350_report_filter.tsv
```

4. predict cancer score (probability) using DeepCAT
```bash
# similar with TCR-seq, but without clustering, notice the difference of input dir: ./test_PBMC_RNA-seq/02_filter_output
bash  ./Script_DeepCAT.sh -t ./test_PBMC_RNA-seq/02_filter_output/disease 

bash  ./Script_DeepCAT.sh -t ./test_PBMC_RNA-seq/02_filter_output/control

mv ./Cancer_score_{control,disease}.txt ./test_PBMC_RNA-seq/03_deepcat_output
```

5. visualize cancer score result
```bash
# the same wiht TCR-seq
Rscript ./plot.R ./test_PBMC_RNA-seq/03_deepcat_output
```


## 4) Tips/Utilities

**Wiki**

* Liquid biopsy

A liquid biopsy, also known as fluid biopsy or fluid phase biopsy, is the sampling and analysis of non-solid biological tissue, primarily blood. Like traditional biopsy this type of technique is mainly used as a diagnostic and monitoring tool for diseases such as cancer, with the added benefit of being largely non-invasive.

* PBMC

A peripheral blood mononuclear cell (PBMC) is any peripheral blood cell having a round nucleus. These cells consist of lymphocytes (T cells, B cells, NK cells) and monocytes, whereas erythrocytes and platelets have no nuclei, and granulocytes (neutrophils, basophils, and eosinophils) have multi-lobed nuclei. 

* TCR

The T-cell receptor (TCR) is a protein complex found on the surface of T cells, or T lymphocytes, that is responsible for recognizing fragments of antigen as peptides bound to major histocompatibility complex (MHC) molecules.

* CDR3

Complementarity-determining regions (CDRs) are part of the variable chains in immunoglobulins (antibodies) and T cell receptors, generated by B-cells and T-cells respectively, where these molecules bind to their specific antigen.
There are three CDRs (CDR1, CDR2 and CDR3), arranged non-consecutively, on the amino acid sequence of a variable domain of an antigen receptor.
CDR3 is the most variable.

* Immune repertoire

The immune repertoire encompasses the different sub-types an organism's immune system makes of immunoglobulins or T-cell receptors. These help recognise pathogens in most vertebrates. The sub-types, all differing slightly from each other, can amount to tens of thousands, or millions in a given organism. Such a wide variety increases the odds of having a sub-type that recognises one of the many pathogens an organism may encounter.

* neoantigen

Neoantigens are mutated antigens specifically expressed by tumor tissue and are not expressed on the surface of normal cells. Development of sequencing technology has improved the accuracy of identification and localization of neoantigens.

* RNA-seq

RNA-Seq (named as an abbreviation of "RNA sequencing") is a particular technology-based sequencing technique which uses next-generation sequencing (NGS) to reveal the presence and quantity of RNA in a biological sample at a given moment, analyzing the continuously changing cellular transcriptome.

* TCR-seq

TCR-seq using multiple-PCR, 5' RACE or target enrichment methods, followed by deep sequencing and data analysis.
![TCR-seq.png](figs/TCR-seq.png)



**Reference**
* TRUST4: [Ultrasensitive detection of TCR hypervariable-region sequences in solid-tissue RNA–seq data](https://www.nature.com/articles/ng.3820)
* iSMART: [Investigation of Antigen-Specific T-Cell Receptor Clusters in Human Cancers](https://clincancerres.aacrjournals.org/content/26/6/1359)
* DeepCAT: [De novo prediction of cancer-associated T cell receptors for noninvasive cancer detection](https://stm.sciencemag.org/content/12/557/eaaz3738)


> if you want to exercise with raw fastq input with out RNA-seq pipeline, please download raw data in [here](https://cloud.tsinghua.edu.cn/library/14156f8d-93f5-496d-8837-90a8e0d24e4e/shared_data/) (Tsinghua cloud)
or ```/BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/pico/02.rawdata_PBMC```



## 5) Homework

1. Finish above 2 pipelines using provided input data

2. Compare output between TCR-seq and RNA-seq, which sequencing method performs better in cancer diagnosis and why, please summarize (dis)advantages of both sequencing methods.
