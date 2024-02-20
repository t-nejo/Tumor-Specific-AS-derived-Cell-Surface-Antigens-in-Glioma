#!/bin/sh

echo "Written by Takahide Nejo, MD, PhD"
echo "UC San Francisco, Neurological Surgery, Okada Lab"
echo ' '
date
echo ' '
echo 'SET VARIABLES' 
export FASTA1=path_to_file/UP000005640_9606.fasta # downloaded from UniProt (https://www.uniprot.org/proteomes/UP000005640) in November 2020
export FASTA2=path_to_file/aaseq_alt_n17_from_14events.fa
export FASTA3=path_to_output_directory/swissprot_n20600_and_usr_prot_n17.fa
echo ' '
echo 'CONCATENATE'
cat ${FASTA1} ${FASTA2} > ${FASTA3}
echo 'DONE'
date
## end
