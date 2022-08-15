#!/usr/bin/env Rscript

# usage: run_model_ftop.R ftom.ftop.count.table.txt ftom.ftop.count.table.out.txt 50000 T 10000 123 N 10 6 0.05 0.1

# ftom.ftop.count.table.txt format (two FTO- replicates vs two FTO+ replicates):
# require at least 11 columns including the following column names: "pos","motif","type","ftom_G_1","ftom_A_1","ftom_G_2","ftom_A_2","ftop_G_1","ftop_A_1","ftop_G_2","ftop_A_2"
# example:
# pos            motif  type      ftom_G_1  ftom_A_1  ftom_G_2  ftom_A_2  ftop_G_1  ftop_A_1  ftop_G_2  ftop_A_2
# chr1_14687_-   UUAUC  nonDRACH  1         10        34        16        74        52        85        56      

# ftom.ftop.count.table.txt format (one FTO- replicate vs one FTO+ replicate):
# require at least seven columns including the following column names: "pos","motif","type","ftom_G_count","ftom_A_count","ftop_G_count","ftop_A_count"
# example:
# pos           motif  type      ftom_G_count  ftom_A_count  ftop_G_count  ftop_A_count
# chr1_14580_-  GCACU  nonDRACH  8             3             17            3           

start <- Sys.time()
cat("Starting Time is ",format(start,format="%Y-%m-%d %H:%M:%S"),"\n\n",sep="")

argv <- commandArgs(TRUE)

# count table: e.g., ftom.ftop.count.table.txt
in_file <- argv[1]
# output results: e.g., ftom.ftop.count.table.out.txt
out_file <- argv[2]
# site number for conversion rate estimation: e.g., 50000
est_num <- as.numeric(argv[3])
# enable FTO efficiency estimation or not: logical
est_FTO_efficiency <- as.logical(argv[4])
# site number for accessibility estimation: e.g., 10000
select_num <- as.numeric(argv[5])
# specify seeds for reproducible tests
seed <- argv[6]
# either Y or N. If Y, sum up counts from two replicates in the count table (the two FTO- replicates vs two IVT replicates format).  
rep <- argv[7]
# site coverage cutoff for estimation
read_count <- as.numeric(argv[8])
# threads used to run the script
thread <- as.numeric(argv[9])
# pre-set cutoff for the identification of m6A sites
fdr <- as.numeric(argv[10])
accessible_methylation <- as.numeric(argv[11])
accessibility <- 80
delta <- 30

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
	in_file_data$ftop_G_count<-in_file_data$ftop_G_1+in_file_data$ftop_G_2
	in_file_data$ftop_A_count<-in_file_data$ftop_A_1+in_file_data$ftop_A_2
	in_file_data$ftom_total_count<-in_file_data$ftom_G_count+in_file_data$ftom_A_count
	in_file_data$ftop_total_count<-in_file_data$ftop_G_count+in_file_data$ftop_A_count
	in_file_data$ftom_G_rate<-round(in_file_data$ftom_G_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ftom_A_rate<-round(in_file_data$ftom_A_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ftop_G_rate<-round(in_file_data$ftop_G_count/in_file_data$ftop_total_count*100,digits=4)
	in_file_data$ftop_A_rate<-round(in_file_data$ftop_A_count/in_file_data$ftop_total_count*100,digits=4)
	in_file_data<-subset(in_file_data,ftom_total_count >= read_count & ftop_total_count >= read_count,select=c("pos","motif","type","ftom_G_count","ftom_A_count","ftop_G_count","ftop_A_count","ftom_total_count","ftop_total_count","ftom_G_rate","ftom_A_rate","ftop_G_rate","ftop_A_rate"))
} else if (rep=="N") {
	in_file_data$ftom_total_count<-in_file_data$ftom_G_count+in_file_data$ftom_A_count
	in_file_data$ftop_total_count<-in_file_data$ftop_G_count+in_file_data$ftop_A_count
	in_file_data$ftom_G_rate<-round(in_file_data$ftom_G_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ftom_A_rate<-round(in_file_data$ftom_A_count/in_file_data$ftom_total_count*100,digits=4)
	in_file_data$ftop_G_rate<-round(in_file_data$ftop_G_count/in_file_data$ftop_total_count*100,digits=4)
	in_file_data$ftop_A_rate<-round(in_file_data$ftop_A_count/in_file_data$ftop_total_count*100,digits=4)
	in_file_data<-subset(in_file_data,ftom_total_count >= read_count & ftop_total_count >= read_count,select=c("pos","motif","type","ftom_G_count","ftom_A_count","ftop_G_count","ftop_A_count","ftom_total_count","ftop_total_count","ftom_G_rate","ftom_A_rate","ftop_G_rate","ftop_A_rate"))
} else {
	cat("rep is either Y or N!!!\n")
	q(status = 1)
}
# model test
in_file_data_ftop<-with(in_file_data,analysis_for_FTO_minus_plus_with_shrinkage(ftop_A_count,ftop_A_rate,ftom_A_count,ftom_A_rate,ftom_G_count,ftom_total_count,ftop_G_count,ftop_total_count,est_FTO_efficiency=est_FTO_efficiency,default_FTO_efficiency=0.5,est.num=est_num,verbose=T,selected.num=select_num,seed=seed,thread=thread))
colnames(in_file_data_ftop$all.chr.estimates)<-c("FTO_minus_pvalue","FTO_minus_est_beta","FTO_minus_G_rate","fisher_pvalue",colnames(in_file_data_ftop$all.chr.estimates[5:ncol(in_file_data_ftop$all.chr.estimates)]))
in_file_data_ftop$all.chr.estimates$FTO_minus_adjp<-p.adjust(in_file_data_ftop$all.chr.estimates$FTO_minus_pvalue,method="BH")
in_file_data_ftop$all.chr.estimates$fisher_adjp<-p.adjust(in_file_data_ftop$all.chr.estimates$fisher_pvalue,method="BH")
in_file_data_summary<-with(cbind(in_file_data,in_file_data_ftop$all.chr.estimates),data.frame(pos=pos,motif=motif,type=type,FTOm_A_rate=ftom_A_rate,FTOp_A_rate=ftop_A_rate,FTO_minus_pvalue=signif(FTO_minus_pvalue,4),FTO_minus_adjp=signif(FTO_minus_adjp,4),fisher_pvalue=signif(fisher_pvalue,4),fisher_adjp=signif(fisher_adjp,4),FTO_minus_est_beta=round(FTO_minus_est_beta*100,4),FTO_minus_G_rate=round(FTO_minus_G_rate*100,4),FTO_est_accessibility=round(FTO_EstAccessibility*100,4),FTO_methylation=round(FTO_Methylation*100,4),FTO_final_accessibility=round(final_Accessibility*100,4),FTOm_total_count=ftom_total_count,FTOp_total_count=ftop_total_count))
# output the results
write.table(subset(in_file_data_summary,fisher_adjp<fdr & FTO_methylation*FTO_est_accessibility/10000>=accessible_methylation),out_file,sep="\t",quote=F,row.names=F)
write.table(in_file_data_summary,paste0(out_file,".all"),sep="\t",quote=F,row.names=F)

# summarize identified m6A sites
cat("m6A sites (weighted accessible methylation >=",accessible_methylation,"accessibility >=",accessibility,", fisher FDR <",fdr,")\n")
table(subset(in_file_data_summary,fisher_adjp<fdr & FTO_methylation*FTO_final_accessibility/10000>=accessible_methylation & FTO_final_accessibility>=accessibility)$type)
cat("m6A sites (weighted accessible methylation >=",accessible_methylation,", fisher FDR <",fdr,")\n")
table(subset(in_file_data_summary,fisher_adjp<fdr & FTO_methylation*FTO_final_accessibility/10000>=accessible_methylation)$type)
cat("m6A sites (accessible methylation >=",accessible_methylation,"accessibility >=",accessibility,", fisher FDR <",fdr,")\n")
table(subset(in_file_data_summary,fisher_adjp<fdr & FTO_methylation*FTO_est_accessibility/10000>=accessible_methylation & FTO_est_accessibility>=accessibility)$type)
cat("m6A sites (accessible methylation >=",accessible_methylation,", fisher FDR <",fdr,")\n")
table(subset(in_file_data_summary,fisher_adjp<fdr & FTO_methylation*FTO_est_accessibility/10000>=accessible_methylation)$type)

warnings()

end <- Sys.time()
cat("Ending Time is",format(end,format="%Y-%m-%d %H:%M:%S"),"\n")
runtime <- difftime(end,start,unit="mins")
cat("Used Time is",runtime,"mins\n")

