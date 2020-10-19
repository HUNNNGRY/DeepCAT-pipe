#!/usr/bin/Rscript

##sample ID : `ls -hl /readonly/Share2/home/lulab/jinyunfan/exSeek-dev/output/PBMC/bam`
##most important columns/headers: CDR3aa	V	frequency	(Group)
## among them : frequency=frequencyCount..../sum(total TCRB template reads number)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
 stop("usage: Rscript filter.R indir outdir sampleid")
} else if (length(args)==3) {
 print(paste("indir:",args[1],";outdir:",args[2],",sampleid:",args[3]))
 print(paste("start",args[3],"filter at",Sys.time()))
}

# read input
input.name <- paste(args[1],"/",args[3],"_report.tsv",sep = "")
TRUST4 <- read.delim2(file = input.name)
TRUST4 <- TRUST4[,c(4,5,1)]
dim0 <- dim(TRUST4)
print(dim0)

# remove/filter poor/low quality CDR3 sequences/calls/records
## 1.rm record gene locus was not solved
TRUST4$V.nchar <- nchar(as.vector(TRUST4$V))
TRUST4 <- TRUST4[TRUST4$V.nchar > 1,] 
dim1 <- dim(TRUST4)
print(dim1)

## 2.rm record's sequence length was <10 or >24
TRUST4$aa.length <- nchar(as.vector(TRUST4$CDR3aa))
TRUST4 <- TRUST4[TRUST4$aa.length >= 10 & TRUST4$aa.length <= 24, ] 
dim2 <- dim(TRUST4)
print(dim2)

###optional: length distribution plot
###frequency.seq <- table(TRUST4$aa.length)
###frequency.seq
###plot(frequency.seq)

## 3.rm record's CDR3.aa has special symbol: _ * . + X ?        NOTE:perl = F
sequence_without_special <- TRUST4$CDR3aa %in% grep("[^A-Z]", TRUST4$CDR3aa, value = TRUE, ignore.case = T, invert = T)
TRUST4 <- subset(x = TRUST4, subset = sequence_without_special)
dim3 <- dim(TRUST4)
print(dim3)

## 4.rm record's CDR3.aa not start with C, end with F
###filter.head <- grep(pattern="^C", x = TRUST4$CDR3aa, fixed = F,ignore.case = T, value = T)
###filter.head
###length(filter.head)
filter.head <- grep(pattern="^C", x = TRUST4$CDR3aa, fixed = F,ignore.case = T)
TRUST4 <- TRUST4[filter.head, ]

###filter.tail
###filter.tail <- grep(pattern="F$", x = TRUST4$CDR3aa, fixed = F,ignore.case = T, value = T)
###length(filter.tail)
filter.tail <- grep(pattern="F$", x = TRUST4$CDR3aa, fixed = F,ignore.case = T)
TRUST4 <- TRUST4[filter.tail, ]

dim4 <- dim(TRUST4)
print(dim4)

## 5.rm not TRBV
TRUST4$v.class <- substr(TRUST4$V,1,4)
TRUST4 <- TRUST4[TRUST4$v.class == "TRBV", ]
dim5 <- dim(TRUST4)
print(dim5)

## change frequency after filter
TRUST4$frequency <- TRUST4$X.count/sum(TRUST4$X.count)

###optional: plot receptor class
###table(TRUST4$v.class)
###plot(table(TRUST4$v.class))

# export filtered file
output.name <- paste(args[2],"/",args[3],"_report_filter.tsv",sep = "")
log.name <- paste(args[2],"/",args[3],"_filter.log",sep = "")
## save filter output
write.table(x = TRUST4, file = output.name,sep = "\t",row.names = F,quote = F)
## save filter info
TRUST4.log <- t(data.frame(dim0,dim1,dim2,dim3,dim4,dim5))
write.table(x = TRUST4.log, file = log.name, sep = "\t",row.names = T,quote = F)
