#!/bin/bash                                                                                                                                      
# dependencies - files must be copied into folder that this processing occurs in

# rename all the files so that only the lane and sample are apparent, add a header for each column                                               
#for i in `ls *out.tab`                                                                                                                          
#do                                                                                                                                              
#    lane=`echo -n $i | cut -f1 -d-`                                                                                                             
#    sample=`echo -n $i | cut -f6 -d- | cut -f1 -d_`                                                                                             
#    echo -e 'transcript_id\tcounts\tcol3\tcol4' > header.tsv                                                                                    
#    cat header.tsv $i > tmpfile; mv tmpfile sample_${sample}_${lane}.tsv                                                                        
#done                                                                                                                                            

# get a list of the transcript ids and put into a new file all-results.tsv                                                                       
#cut -f1 `find . -name "sample*" | head -1` > all-results.tsv                                                                                    

# for each .tsv result file with header extract the first two columns, join them to a new finished results fil                                   
#for i in `find . -name "sample*" | sort -V`                                                                                                     
#do                                                                                                                                              
#    cat $i | cut -f1,2 > $i-counts                                                                                                              
#    join -t $'\t' all-results.tsv $i-counts > /tmp/new-results.tsv                                                                              
#    mv /tmp/new-results.tsv all-results.tsv                                                                                                     
#    rm *-counts                                                                                                                                 
#done                                                                                                                                            

#                                                                                                                                                
#header=`find . -name "sample*" | sort -V | cut -f2 -d/ | sed 's/.tsv//g' `                                                                      
#echo -n -e "transcript_id\t${header}" | tr "\n" "\t" > header                                                                                   
#echo "" >> header                                                                                                                               
#grep -v transcript_id all-results.tsv > /tmp/all-results.tsv                                                                                    
#cat header /tmp/all-results.tsv > all-results.tsv                                                                                               

