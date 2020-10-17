# TCR
lulab TCR pipeline using available tools

# Detect Cancer Associated TCR in Blood

本章主要介绍公共数据库中PBMC或tissue样本的TCR-seq，bulk RNA-seq数据和用于肿瘤分类的方法，包括 peak calling，filter，cluster, predict等步骤。

肿瘤微环境免疫在肿瘤发生和发展中起到关键作用，近年来随着液体活检的快速发展，人们能够在外周血中检测到各种来自肿瘤和周边免疫环境的信号，包括能够识别肿瘤新抗原 (neoantigen)的反应性T细胞，这使得直接基于外周血样本区分正常人群和无症状肿瘤早期人群的无创诊断 (non-invasive diagnosis)成为可能。
这里提供了一个整合了多种已发表方法的分析流程，基于PBMC的TCR-seq或bulk RNA-seq，检测外周血细胞中潜在的少量以caTCR(cancer-associated TCR)为代表的肿瘤免疫环境的信号。


## 1) Pipeline
![7aec70ff7adc8f0329b14d0e288ba0e8.png](evernotecid://6BB70527-BB50-40E3-8023-A07526CA92F8/appyinxiangcom/28382199/ENResource/p33)

## 2) Data Structure
2a) getting software & data
方法1: 使用docker

方法2: 自行下载和安装
raw data
1.TCR-seq
2.lulab RNA-seq

software

1. anaconda3
2. TRUST4
[https://github.com/HUNNNGRY/TRUST4](https://github.com/HUNNNGRY/TRUST4)
3. DeepCAT (fork AND correct)

environment configuration see [https://github.com/HUNNNGRY/DeepCAT/blob/master/README.md](https://github.com/HUNNNGRY/DeepCAT/blob/master/README.md)

4. iSMART

2b) input

2c) output


## 3) Running Steps

## 4) Tips/Utilities

## 5) Homework
