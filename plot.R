#!/usr/bin/Rscript

#FOR visualization of DeepCAT without iSMART

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("usage: Rscript inputdir")
} else if (length(args)==1) {
  print(paste("cancer score txt dir:",args[1]))
}

inputdir=args[1]
cancername <- paste(inputdir,"/Cancer_score_disease.txt",sep = "")
normalname <- paste(inputdir,"/Cancer_score_control.txt",sep = "")

# import and merge 2 groups
Cancer <- read.delim2(cancername,header = F,stringsAsFactors = F) 
Control <- read.delim2(normalname,header = F,stringsAsFactors = F) 
Cancer$disease <- "Cancer"
Control$disease <- "Control"
score <- rbind(Cancer,Control)

colnames(score) <- c("patients","DeepCAT_Score","disease")
score$DeepCAT_Score <- as.numeric(score$DeepCAT_Score)

score$labels <- "0"
score$labels[score$disease == "Cancer" ] <- "1"

# boxplot
if(!require(devtools)) install.packages("ggpubr")
if(!require(devtools)) install.packages("ggplot2")
library(ggpubr)
library(ggplot2,quietly = T)
ggplot(score,aes(x=disease,y=as.numeric(DeepCAT_Score))) + 
  geom_boxplot(fill = c("firebrick","steelblue"),outlier.shape = NA,show.legend = T) +
  stat_compare_means(comparisons = list(c("Control","Cancer"))) + # Add pairwise comparisons p-value
  geom_jitter(width = .2) +
  scale_y_continuous(limits = c(0.1, 0.8) , oob = scales::squish) + # 限定范围方便比较
  theme_bw()
ggsave(filename = paste(inputdir,"/boxplot.pdf",sep = ""))


# t test
t_test <- t.test(as.double(Cancer$V2),as.double(Control$V2))
print(t_test)

# ROC   (好像label得指定为1，0，字符或数值型均可，肿瘤较高水平为1, 否则出现很小的auc）
if(!require(devtools)) install.packages("ROCR")
if(!require(devtools)) install.packages("rms")
library(ROCR,quietly = T)
library(rms,quietly = T)
ROC1 <- prediction(score$DeepCAT_Score, score$labels)   #构建ROC预测模型 
ROC2 <- performance(ROC1,"tpr","fpr")   #计算预测模型的TPR/FPR值
AUC <- performance(ROC1,"auc")   #计算曲线下面积(AUC)值
AUC.value <- unlist(AUC@y.values) #查看AUC值
print(AUC.value)

pdf(paste(inputdir,"/auc.pdf",sep = ""))  # 需要提前指定普通绘图路径和文件名
plot(ROC2, 
     col="red",   #曲线的颜色
     xlab="False positive rate", ylab="True positive rate",   # x轴和y轴的名称
     lty=1,lwd=3,
     main=paste("AUC=",AUC.value))
abline(0, 1, lty=2, lwd=3)   #绘制对角线

