#!/bin/bash

##### For loop aligns a list of trimmed fastq files using STAR aligning tool.
##### Creates sam and bam output files for each alignment. Use these files to create count tables. 
###See STAR Manual for STAR parameters. 


##FILEPATHS
#trimfq= list of trimmed fq files to align in specified directory
#STAR= filepath to STAR alignment program directory
#genomeDir= filepath to alignment reference genomes
#outfile=filepath to directory for bam output files and outfile prefix

trimfq=`ls /ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/data/L12345678/trimmed/e*/lane*`
STAR=/local/cluster/STAR/bin/Linux_x86_64/STAR
genomeDir=/ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/Brian/genomedir
outfile=/ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/data/L12345678/${filename}

for fq in ${trimfq}; do \
	echo "processing file: ${fq}"
	filename="$(basename -s .fq ${fq})"
		
	${STAR} --runThreadN 8 \
	--genomeDir ${genomeDir} \
	--readFilesIn ${fq}	\
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverLmax 0.6 \ 
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \ 
	--alignMatesGapMax 1000000 \
	--outSAMattributes NH HI NM MD \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ${outfile} \
	--quantMode GeneCounts
done
