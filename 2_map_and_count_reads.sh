#!/bin/bash


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

## main setting
# sample fastq file path
fastq_dir=/sample/fastq/dir
# hisat-3n rRNA index name: e.g., rRNA
hisat3n_rep_index_name=rRNA
# hisat-3n rRNA index path
hisat3n_rep_index=/hisat-3n/rRNA/index/path/$hisat3n_rep_index_name
# hisat-3n genome index name: e.g., GRCh38_tran
hisat3n_index_name=GRCh38_tran
# hisat-3n genome index path
hisat3n_index=/hisat-3n/genome/index/path/$hisat3n_index_name
# hisat-3n alignment result output path
hisat3n_out_dir=/hisat-3n/output/dir
# rRNA sequneces in fasta format
rep_fa=/rRNA/sequences/rRNA.fa
# genome sequneces in fasta format
genome_fa=/genome/sequences/genome.fa
# threads to run the bash script
ncpus=24
# minimum percentage of converted As in a read to be considered as a valid read
A2G_percent=0.5
# minimum read coverage in a position to be considered for output
cov=1
# sample name: e.g., hela.polya.wt.ftom.rep1
sample=sample_name
# library type: F for forward, R for reverse
strandness="F"
# here R2 reads were used.
sample_r2_fq_file=$fastq_dir/$sample_name.R2.a3a5trim.umi.fastq.gz


echo -e "\n$sample\n"

if [[ -d $hisat3n_out_dir/$sample ]];then rm -rf $hisat3n_out_dir/$sample;fi
mkdir -p $hisat3n_out_dir/$sample

echo -e "\nhisat3n align -- rRNA mapping\n"
hisat-3n -p $ncpus --time --base-change A,G --no-spliced-alignment --no-softclip --norc --no-unal --rna-strandness $strandness -x $hisat3n_rep_index -U $sample_r2_fq_file --un-gz $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz | \
	samtools sort -@ $ncpus -T $hisat3n_out_dir/$sample -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam -

samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam

samtools stats -r $rep_fa $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam > $sample/$sample.$hisat3n_rep_index_name.align.sorted.stats

# filter reads by the percentage of converted As
samtools view -@ $ncpus -hb -e "([Yf]+[Zf]>0) && ([Yf]/([Yf]+[Zf])>=$A2G_percent)" $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam -
# Count converted As and unconverted As. Here hisat-3n-table is used as an example. Other similar tools can also be used to complete the task. 
#samtools view -@ $ncpus $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam | hisat-3n-table -u -p $ncpus --alignments - --ref $rep_fa --output-name $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.conversion.flt.txt --base-change A,G
# Here, pileup2var is used.
pileup2var -f 524 -s $strandness -g $rep_fa -b $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.align.sorted.flt.bam -o $hisat3n_out_dir/$sample/$sample.$hisat3n_rep_index_name.pileup2var.flt.txt

echo -e "\nhisat3n align -- genome mapping\n"
hisat-3n -p $ncpus --time --base-change A,G --repeat --repeat-limit 1000 --bowtie2-dp 0 --no-unal --rna-strandness $strandness -x $hisat3n_index -U $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz | \
	samtools view -@ $ncpus -Shb - -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam

# filter out multiple-loci mapped reads
samtools view -@ $ncpus -q 60 -hb $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam -
rm $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.raw.bam

samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam

# deduplication
umi_tools dedup --random-seed=123 --method=unique --spliced-is-unique -I $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.bam -S $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam \
	--output-stats=$hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup -L $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.log

samtools index $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam

samtools stats -r $genome_fa $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam > $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.stats

# filter reads by the percentage of converted As
samtools view -@ $ncpus -hb -e "([Yf]+[Zf]>0) && ([Yf]/([Yf]+[Zf])>=$A2G_percent)" $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.bam | samtools sort -T $hisat3n_out_dir/$sample -@ $ncpus -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam -
# Count converted As and unconverted As. Here hisat-3n-table is used as an example. Other similar tools can also be used to complete the task. 
#samtools view -@ $ncpus $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam | hisat-3n-table -p $ncpus --alignments - --ref $genome_fa --output-name $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.conversion.flt.txt --base-change A,G
# Here, pileup2var is used.
pileup2var -t $ncpus -f 524 -a A -c $cov -s $strandness -g $genome_fa -b $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.align.sorted.dedup.flt.bam -o $hisat3n_out_dir/$sample/$sample.$hisat3n_index_name.pileup2var.flt.txt

rm $hisat3n_out_dir/$sample/$sample.rmRep.fastq.gz
wait

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
