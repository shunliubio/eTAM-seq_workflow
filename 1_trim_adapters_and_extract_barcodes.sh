#!/bin/bash


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

# fastq file name format: example.R1.fastq.gz example.R2.fastq.gz
raw_name=example

# sample name: e.g., hela.polya.wt.ftom.rep1
sample=sample_name

# trim R1 5' end adapters and R2 3' end adapters
cutadapt -e 0.1 -n 1 -O 1 -q 6 -m 0:46 -g GACGCTCTTCCGATCT -A AGATCGGAAGAGCGTC \
	-o $sample.R1.a5trim.fastq.gz -p $sample.R2.a3trim.fastq.gz $raw_name.R1.fastq.gz $raw_name.R2.fastq.gz > $sample.a3a5trim.umi.round1.log
# trim R1 3' end adapters and R2 5' end adapters
cutadapt -e 0.1 -n 1 -O 7 -q 6 -m 0:46 -a NNNNNNAGATCGGAAGAGCACA -G TGTGCTCTTCCGATCT \
	-o $sample.R1.a3a5trim.fastq.gz -p $sample.R2.a3a5trim.fastq.gz $sample.R1.a5trim.fastq.gz $sample.R2.a3trim.fastq.gz > $sample.a3a5trim.umi.round2.log

# extract R2 5' end barcodes (6-mer)
umi_tools extract --random-seed=123 --extract-method=regex --bc-pattern=".*" --bc-pattern2="^(?P<umi_1>.{6}).*" -I $sample.R1.a3a5trim.fastq.gz --read2-in=$sample.R2.a3a5trim.fastq.gz \
	--stdout=$sample.R1.a3a5trim.umi.fastq.gz --read2-out=$sample.R2.a3a5trim.umi.fastq.gz -L $sample.a3a5trim.umi.log

rm $sample.R1.a5trim.fastq.gz $sample.R2.a3trim.fastq.gz
rm $sample.R[12].a3a5trim.fastq.gz

echo -e "\nAll done\n"

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
