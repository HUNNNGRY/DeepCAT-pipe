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
  scale_y_continuous(limits = c(0.1, 0.8) , oob = scales::squish) + # limit to better compare
  theme_bw()
ggsave(filename = paste(inputdir,"/boxplot.pdf",sep = ""))


# t test
t_test <- t.test(as.double(Cancer$V2),as.double(Control$V2))
print(t_test)

# ROC  (tumor: higher)
if(!require(devtools)) install.packages("ROCR")
if(!require(devtools)) install.packages("rms")
library(ROCR,quietly = T)
library(rms,quietly = T)
ROC1 <- prediction(score$DeepCAT_Score, score$labels)   #get ROC model 
ROC2 <- performance(ROC1,"tpr","fpr")   #calculate model's TPR/FPR
AUC <- performance(ROC1,"auc")   #calculate model's AUC
AUC.value <- unlist(AUC@y.values) #check AUC
print(AUC.value)

pdf(paste(inputdir,"/auc.pdf",sep = ""))  # predefine path to save plot 
plot(ROC2, 
     col="red",   #curve color
     xlab="False positive rate", ylab="True positive rate",   # x/y-axis name
     lty=1,lwd=3,
     main=paste("AUC=",AUC.value))
abline(0, 1, lty=2, lwd=3)   # draw diagonal

