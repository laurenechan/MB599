#!/bin/bash                                                                                                                \
                                                                                                                            

CPU=8
memory=64G

#access STAR software available on the CGRB cluster                                                                         

/local/cluster/STAR/bin/Linux_x86_64/STAR --runThreadN $CPU \

#generate a zebrafish genome index                                                                                          
--runMode genomeGenerate \

#create genome in this directory that was already created with mkdir                                                        
--genomeDir /ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/Brian/genomedir \

#primary zfish assembly retrived by wget from ensembl                                                                       
--genomeFastaFiles /ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/data/zebrafish_genome/Danio_rerio.GRC\
z11.dna.primary_assembly.fa \

#arbitrary limit? This came from an error so idk how to comment this one out                                                
--limitGenomeGenerateRAM 41387679999 \

#annotations were retrieved by wget from ensembl in gtf format                                                              
--sjdbGTFfile /ACTF/Course/mb599_bds_s20/data/share/Group5_Lauren_Chloe_Brian_MK/data/zebrafish_genome/Danio_rerio.GRCz11.1\
00.gtf \

#splice junction database overhang should be the number of bp/read - 1 = 99                                                 
--sjdbOverhang 99


