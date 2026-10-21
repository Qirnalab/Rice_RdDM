args <- commandArgs(T)
name=basename(args[1])
id <- strsplit(name, ".", fixed=TRUE)
name <- id[[1]]

library(DESeq)
sampleTable <- read.table(args[1],row.names=1,header=F,sep="\t")
#new_name <- paste(name,"_wt",sep="") 
new_name <- newCountDataSet(sampleTable, c( "WT", "WT", "WT", "pol4", "pol4", "pol4" ))
new_name <- estimateSizeFactors(new_name)
new_name <- estimateDispersions(new_name)
diff <- nbinomTest( new_name, "WT","pol4")
write.table(diff,file="WT_pol4.diff",quote=F,sep="\t",row.names=F,col.names=T)
