#!/bin/bash

##Trims raw fastq files using BBTools(bbduk). 
##Runs FastQC on trimmed files. 

#fqlist= list of raw fq files to trim, Note: further fastqc selection is made by calling specific lanes within for loop. 
#BBmap= filepath to bbduk script in BBmap suite
#resource= filepath to Illumina adapter and poly(A) trimming reference files. 
#out_fastqc= filepath to FastQC output files.
#filename= removes fastq.gz from original filename

fqlist=`ls /ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/data/L12345678/lane*_R1_001.fastq.gz`
BBmap=/ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/bbmap/bbduk.sh
resource=/ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/bbmap/resources
out_fastqc=/ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/data/L12345678/fastqc_trim/

mkdir -p ${out_fastqc}

for fq in ${fqlist}
do
if [[ "${fq}" =~ lane(7|8).* ]]
then
	echo true
	echo ${fq} 
	filename="$(basename -s .fastq ${fq})"

	${BBmap} \
	in=${fq} \
	out="${filename}_trim.fastq.gz" \
	ref=${resource}/polyA.fa.gz,${resource}/truseq_rna.fa.gz \
	k=13 ktrim=r useshortkmers=t \
	mink=5 qtrim=r trimq=10 minlength=20

	fastqc -o ${out_fastqc} -t 8 --nogroup ${fq} "${fq}_trim.fastq.gz"

else
	echo False
fi
done 
