library("DESeq2")





directory<-"............."
sampleFiles <- grep("??",list.files(directory),value=TRUE)

time <- factor(c("T12", "T13", "T14", "T1", "T11", "T5", "T7", "T2", "T9", "T3", "T6", "T8", "T4", "T10", "T12", "T13", "T14", "T1", "T11", "T5", "T7", "T2", "T9", "T3", "T6", "T8", "T4", "T10" ))
genotype <- factor(c("GF","GF","GF","GF","GF","GF","GF","GF","GF","GF","GF","GF","GF","GF","VB","VB","VB","VB","VB","VB","VB","VB","VB","VB","VB","VB","VB","VB"))
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles,  genotype=genotype, time=time)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~genotype+time+time:genotype)

ddsHTSeq$genotype<-factor(ddsHTSeq$genotype, levels=c("GF","VB"))
ddsHTSeq$time <- factor(ddsHTSeq$time, levels=c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12", "T13", "T14"))

dds <- DESeq(ddsHTSeq)
res <- results(dds)

ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
ddsHTSeq <- estimateDispersions(ddsHTSeq)
ddsLRT <- nbinomLRT(dds, reduced = formula(~genotype+time)