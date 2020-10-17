# Detect Cancer Associated TCR in Blood

本章主要介绍公共数据库中PBMC或tissue样本的TCR-seq，bulk RNA-seq数据和用于肿瘤分类的方法，包括 peak calling，filter，cluster, predict等步骤。

肿瘤微环境免疫在肿瘤发生和发展中起到关键作用，近年来随着液体活检的快速发展，人们能够在外周血中检测到各种来自肿瘤和周边免疫环境的信号，包括能够识别肿瘤新抗原 (neoantigen)的反应性T细胞，这使得直接基于外周血样本区分正常人群和无症状肿瘤早期人群的无创诊断 (non-invasive diagnosis)成为可能。
这里提供了一个整合了多种已发表方法的分析流程，基于PBMC的TCR-seq或bulk RNA-seq，检测外周血细胞中潜在的少量以caTCR(cancer-associated TCR)为代表的肿瘤免疫环境的信号。


## 1) Pipeline
![7aec70ff7adc8f0329b14d0e288ba0e8.png](evernotecid://6BB70527-BB50-40E3-8023-A07526CA92F8/appyinxiangcom/28382199/ENResource/p33)

## 2) Data Structure
### 2a) getting software & data
方法1: 使用docker
not distributed yet...

方法2: 自行下载和安装
raw data
1.TCR-seq
2.lulab RNA-seq

software

1. anaconda3

2. TRUST4

environment configuration see 
[https://github.com/HUNNNGRY/TRUST4](https://github.com/HUNNNGRY/TRUST4)

3. DeepCAT (fork AND correct)
[https://github.com/HUNNNGRY/DeepCAT](https://github.com/HUNNNGRY/DeepCAT)

4. iSMART
[https://github.com/HUNNNGRY/iSMART](https://github.com/HUNNNGRY/iSMART)

### 2b) input
* pre-processed tsv TCR-seq file from adaptativebiotech
or
* raw fastq file from PBMC RNA seq

### 2c) output
* cancer score: DeepCAT CNN model prediction
* visualization: box plot and ROC curve 

## 3) Running Steps
Enter into bash interface

```bash
git clone https://github.com/HUNNNGRY/TCR.git
cd ./TCR
```
### using PBMC TCR-seq as input

1. view TCR_input data
head -n3 
wc -l 

2. make working directory tree 
mkdir -p ./test_PBMC_TCR-seq/{01_filter_output,02_cluster_output,03_deepcat_output}

3. filter/prepare input  
```bash
# filter invalid TCR sequence
# using TCRB V region CDR3 sequence only
python ./PrepareAdaptiveFile_corrected.py ./sample_data/PBMC_TCR-seq/input  ./test_PBMC_TCR-seq/01_filter_output

# keep an eye on the change of row/record number;留意查看filter前后行数（TCR）的数目变化 
wc -l ./test_PBMC_TCR-seq/01_filter_output
```
4. cluster using iSMART
```bash
# cluster similar TCR sequences using iSMART;使用iSMART聚类
python ./iSMARTv3.py -d ./test_PBMC_TCR-seq/01-filter-output/ -o ./test_PBMC_TCR-seq/02-cluster-output/
#keep an eye on the change of row/record number;留意查看cluster前后行数（TCR）的数目变化 
wc -l ./test_PBMC_TCR-seq/02-cluster-output/
```
5. predict cancer score(probability) using DeepCAT
```bash
# move files of into specific diretory 之前都是统一处理，这里开始区分肿瘤和对照组，同类型文件移动到相同文件夹
mkdir -p ./02-cluster-output/{disease,control}
mv ./02-cluster-output/*BR*ClusteredCDR3s_7.5.txt ./02-cluster-output/disease/
mv ./02-cluster-output/*HIP*ClusteredCDR3s_7.5.txt ./02-cluster-output/control/
# 预测肿瘤样本组和对照样本组
bash  ./Script_DeepCAT.sh -t ./02-cluster-output/disease/ 
bash  ./Script_DeepCAT.sh -t ./02-cluster-output/control/
# 成功后当前文件夹会产生两个txt文件，分布记录肿瘤和对照组的不同样本的预测分数
# 查看样本肿瘤预测得分
head ./Cancer_score_control.txt
head ./Cancer_score_disease.txt
# 把预测结果放入输出目录中
mv ./Cancer_score_{control,disease}.txt ./test_PBMC_TCR-seq/03_deepcat_output
```
6. visualize cancer score result
```bash
# make boxplot and ROC curve using Rscripts;利用已有脚本画boxplot和ROC curve
Rscript ./scripts/plot.R ./test_PBMC_TCR-seq
# 生成的结果在./test_PBMC_TCR-seq/03_deepcat_output路径下
```
* boxplot横坐标两列分别代表肿瘤和对照组，纵坐标是一个样本（每个点）的cancer score
* ROC曲线用于分类效果的评估，线下面积越接近1（越接近左上角）说明分类效果越好

#### using PBMC RNA-seq as input
总体上PBMC RNA-seq作为input和上一步TCR-seq相似，主要差别在于
* 需要先用TRUST4从RNA-seq中得到包含TCR序列的文件
* 由于RNA-seq是非靶向测序，得到的TCR记录可能会很少，这里统一不再经过聚类步骤而直接用DeepCAT预测
1. make working directory tree 
```bash
# 确认当前路径在DeepCAT目录中
mkdir -p ./test_PBMC_RNA-seq/{01_TCRcalling_output,02_filter_output,03_deepcat_output}
```
2. TCR calling
```bash
# 从fastq(.gz)原始文件得到TCR序列信息
run-trust4  -1 sample_1.fastq.gz -2 sample_2.fastq.gz -f ./reference/TRUST4/hg38_bcrtcr.fa --ref ./reference/TRUST4/human_IMGT+C.fa -t 4 -o ./test_PBMC_RNA-seq/01_TCRcalling_output/sample

# 本步骤运行以常见的双端测序PE的原始文件作为输入文件，得到的只是单一样本的TRUST4 TCR calling输出文件，对于其他样本可以用for loop循环执行,如：
for idx in `cat list.txt`;
do 
run-trust4  -1 ${idx}_1.fastq.gz -2 ${idx}_2.fastq.gz -f ./reference/TRUST4/hg38_bcrtcr.fa --ref ./reference/TRUST4/human_IMGT+C.fa -t 4 -o ./test_PBMC_RNA-seq/01_TCRcalling_output/${idx} ;
done

# view TRUST4 TCR calling output data;留意查看TRUST4后的行数（TCR）的数目
head sample....  
wc -l sample....
```
**由于本步骤需要一定时间，可以考虑直接从已有TRUST4 TCR calling输出文件开始下一步
将```./sample_data/PBMC_TCR-seq/input```目录下的文件移动到```./test_PBMC_RNA-seq/01_TCRcalling_output/```目录下即可**

3. filter/prepare input  
```bash
# filter invalid TCR sequence
# using TCRB V region CDR3 sequence only
Rscript ./scripts/filter_TRUST4.R ./test_PBMC_RNA-seq/01_TCRcalling_output   ./test_PBMC_RNA-seq/02_filter_output

# keep an eye on the change of row/record number;留意查看filter前后行数（TCR）的数目变化 
wc -l ./test_PBMC_TCR-seq/01_filter_output/sample...
```
**注意这里使用的是R脚本而不是之前的python脚本PrepareAdaptiveFile_corrected.py**

4. predict cancer score(probability) using DeepCAT
与TCR-seq相同，只是注意因为RNA-seq没有聚类，这里DeepCAT的input目录是./test_PBMC_RNA-seq/02_filter_output，而不是之前的./test_PBMC_RNA-seq/02_cluster_output

5. visualize cancer score result
与TCR-seq相同

### 


## 4) Tips/Utilities
wiki
* Non-invasive dianosis
* PBMC
* TCR
* CDR3
* Immune repertoire
* RNA-seq
* TCR-seq

reference
* TRUST4: summary + paper link
* iSMART: summary + paper link
* DeepCAT: summary + paper link

> if you want to exercise with raw fastq input with out RNA-seq pipeline, please download raw data in [here](https://cloud.tsinghua.edu.cn/library/14156f8d-93f5-496d-8837-90a8e0d24e4e/shared_data/) (Tsinghua cloud)

## 5) Homework
1. finish above 2 pipeline using provided input
2. compare output between TCR-seq and RNA-seq, which sequencing methods perform better in cancer diagnosis and why, please summarize (dis)advantages of both sequencing methods.
