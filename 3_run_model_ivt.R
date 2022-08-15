#!/usr/bin/env Rscript

# usage: run_model_ivt.R ftom.ivt.count.table.txt ftom.ivtcount.table.out.txt 50000 T 123 N 10 0.05 0.1

# ftom.ivt.count.table.txt format (two FTO- replicates vs two IVT replicates):
# require at least 11 columns including the following column names: "pos","motif","type","ftom_G_1","ftom_A_1","ftom_G_2","ftom_A_2","ivt_G_1","ivt_A_1","ivt_G_2","ivt_A_2"
# example:
# pos            motif  type      ftom_G_1  ftom_A_1  ftom_G_2  ftom_A_2  ivt_G_1  ivt_A_1  ivt_G_2  ivt_A_2
# chr1_595009_+  UUAUC  nonDRACH  24        0         24        0         22       1        16       0      

# ftom.ivt.count.table.txt format (one FTO- replicate vs one IVT replicate):
# require at least seven columns including the following column names: "pos","motif","type","ftom_G_count","ftom_A_count","ivt_G_count","ivt_A_count"
# example:
# pos            motif  type      ftom_G_count  ftom_A_count  ivt_G_count  ivt_A_count
# chr1_187451_-  GCAUU  nonDRACH  22            1             8            3          

start <- Sys.time()
cat("Starting Time is ",format(start,format="%Y-%m-%d %H:%M:%S"),"\n\n",sep="")

argv <- commandArgs(TRUE)

# count table: e.g., ftom.ivt.count.table.txt
in_file <- argv[1]
# output results: e.g., ftom.ivtcount.table.out.txt
out_file <- argv[2]
# site number for conversion rate estimation: e.g., 50000
est_num <- as.numeric(argv[3])
# enable linear model or not: logical
linear_fit <- as.logical(argv[4])
# specify seeds for reproducible tests
seed <- argv[5]
# either Y or N. If Y, sum up counts from two replicates in the count table (the two FTO- replicates vs two IVT replicates format).  
rep <- argv[6]
# site coverage cutoff for estimation
read_count <- as.numeric(argv[7])
# pre-set cutoff for the identification of m6A sites
fdr <- as.numeric(argv[8])
accessible_methylation <- as.numeric(argv[9])
accessibility <- 80

if (seed == "NULL") {
	seed <- NULL
} else {
	seed <- as.numeric(seed)
}

library(eTAMseq)

# read count table
in_file_data<-read.delim(in_file,header=T)
# process count table
if (rep == "Y") {
	in_file_data$ftom_G_count<-in_file_data$ftom_G_1+in_file_data$ftom_G_2
	in_file_data$ftom_A_count<-in_file_data$ftom_A_1+in_file_data$ftom_A_2
	in_file_data$ivt_G_count<-in_file_data$ivt_G_1+in_file_data$ivt_G_2
	in_file_data$ivt_A_count<-in_file_data$ivt_A_1+in_file_data$ivt_A_2
	in_file_data$ftom_total_count<-in_file_data$ftom_G_count+in_file_data$ftom_A_count
	in_file_data$ivt_total_count<-in_file_data$ivt_G_count+in_file_data$ivt_A_count
	in_file_data$ftom_G_rate<-round(in_file_data$ftom_G_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ftom_A_rate<-round(in_file_data$ftom_A_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ivt_G_rate<-round(in_file_data$ivt_G_count/in_file_data$ivt_total_count*100,digits=4)
	in_file_data$ivt_A_rate<-round(in_file_data$ivt_A_count/in_file_data$ivt_total_count*100,digits=4)
	in_file_data<-subset(in_file_data,ftom_total_count >= read_count & ivt_total_count >= read_count,select=c("pos","motif","type","ftom_G_count","ftom_A_count","ivt_G_count","ivt_A_count","ftom_total_count","ivt_total_count","ftom_G_rate","ftom_A_rate","ivt_G_rate","ivt_A_rate"))
} else if (rep=="N") {
	in_file_data$ftom_total_count<-in_file_data$ftom_G_count+in_file_data$ftom_A_count
	in_file_data$ivt_total_count<-in_file_data$ivt_G_count+in_file_data$ivt_A_count
	in_file_data$ftom_G_rate<-round(in_file_data$ftom_G_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ftom_A_rate<-round(in_file_data$ftom_A_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ivt_G_rate<-round(in_file_data$ivt_G_count/in_file_data$ivt_total_count*100,digits=4)
	in_file_data$ivt_A_rate<-round(in_file_data$ivt_A_count/in_file_data$ivt_total_count*100,digits=4)
	in_file_data<-subset(in_file_data,ftom_total_count >= read_count & ivt_total_count >= read_count,select=c("pos","motif","type","ftom_G_count","ftom_A_count","ivt_G_count","ivt_A_count","ftom_total_count","ivt_total_count","ftom_G_rate","ftom_A_rate","ivt_G_rate","ivt_A_rate"))
} else {
	cat("rep is either Y or N!!!\n")
	q(status = 1)
}
# model test
in_file_data_ivt<-with(in_file_data,analysis_for_FTO_minus_ivt(ivt_A_rate,ftom_A_rate,ivt_A_count,ivt_total_count,ftom_A_count,ftom_total_count,ivt_G_count,ftom_G_count,est.num=est_num,linear.fit=linear_fit,verbose=T,seed=seed))
colnames(in_file_data_ivt$all.chr.estimates)<-c("IVT_minus_pvalue","IVT_minus_est_beta","IVT_minus_G_rate","fisher_pvalue",colnames(in_file_data_ivt$all.chr.estimates[5:ncol(in_file_data_ivt$all.chr.estimates)]))
in_file_data_ivt$all.chr.estimates$IVT_minus_adjp<-p.adjust(in_file_data_ivt$all.chr.estimates$IVT_minus_pvalue,method="BH")
in_file_data_ivt$all.chr.estimates$fisher_adjp<-p.adjust(in_file_data_ivt$all.chr.estimates$fisher_pvalue,method="BH")
in_file_data_summary<-with(cbind(in_file_data,in_file_data_ivt$all.chr.estimates),data.frame(pos=pos,motif=motif,type=type,FTOm_total_count=ftom_total_count,IVT_total_count=ivt_total_count,IVT_minus_pvalue=signif(IVT_minus_pvalue,4),IVT_minus_adjp=signif(IVT_minus_adjp,4),fisher_pvalue=signif(fisher_pvalue,4),fisher_adjp=signif(fisher_adjp,4),IVT_minus_est_beta=round(IVT_minus_est_beta*100,4),IVT_minus_G_rate=round(IVT_minus_G_rate*100,4),IVT_est_accessibility=round(IVT_EstAccessibility*100,4),IVT_methylation=round(IVT_Methylation*100,4)))
# output the results
write.table(subset(in_file_data_summary,fisher_adjp<fdr & IVT_methylation*IVT_est_accessibility/10000>=accessible_methylation),out_file,sep="\t",quote=F,row.names=F)
write.table(in_file_data_summary,paste0(out_file,".all"),sep="\t",quote=F,row.names=F)

# summarize identified m6A sites
cat("m6A sites (accessible methylation >=",accessible_methylation,"accessibility >=",accessibility,", fisher FDR <",fdr,")\n")
table(subset(in_file_data_summary,fisher_adjp<fdr & IVT_methylation*IVT_est_accessibility/10000>=accessible_methylation & IVT_est_accessibility>=accessibility)$type)
cat("m6A sites (accessible methylation >=",accessible_methylation,", fisher FDR <",fdr,")\n")
table(subset(in_file_data_summary,fisher_adjp<fdr & IVT_methylation*IVT_est_accessibility/10000>=accessible_methylation)$type)

warnings()

end <- Sys.time()
cat("Ending Time is",format(end,format="%Y-%m-%d %H:%M:%S"),"\n")
runtime <- difftime(end,start,unit="mins")
cat("Used Time is",runtime,"mins\n")

