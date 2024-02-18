# Tumor-Specific AS-derived Cell-Surface Antigens in Glioma

## Introduction  

This repository contains the code used to analyze the data in Nejo et al. bioRxiv (2023) "Challenges in the discovery of tumor-specific alternative splicing-derived cell-surface antigens in glioma" (doi: [10.1101/2023.10.26.564156](https://www.biorxiv.org/content/10.1101/2023.10.26.564156v2.full)). In this study, we investigated tumor-specific alternative splicing-derived, cell-surface antigens in glioma, through the analyses of transcriptome data of TCGA-GBM/LGG, GTEx, Mayo-GBM-PDX, as well as the clinical samples of the University of California, San Francisco (UCSF) Brain Tumor Center. Moreover, we conducted nanopore full-length transcript sequencing and proteomics analysis of the CPTAC-GBM project. Our investigation illustrated the diverse characteristics of the tumor-specific AS events and the challenges of antigen exploration due to their notable spatiotemporal heterogeneity and elusive nature at the protein levels. 


## Key analysis tools:  

[STAR](https://github.com/alexdobin/STAR): Short-read RNA-seq read alignment and junction detection  
[IRFinder](https://github.com/williamritchie/IRFinder): Intron-retention analysis  
[Guppy basecaller](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview): Nanopore long read sequencing base calling  
[FLAIR](https://github.com/BrooksLabUCSC/flair): Nanopore long read sequencing FASTQ analysis  
[MSGF+](https://github.com/MSGFPlus/msgfplus): Peptide spectra database search  
[MASIC](https://github.com/PNNL-Comp-Mass-Spec/MASIC): Quantification of 11-TMT reporter ions  

## Key data visualization tools:  
[BioRender](https://www.biorender.com/)  
[Affinity Designer](https://affinity.serif.com/en-us/designer/)  
[R package ggplot2](https://ggplot2.tidyverse.org/)  

## Citation

**If you use our data or code for your work, please consider citing the following publication:**  

Nejo T, Wang L, Leung KK, Wang A, Lakshmanachetty S, Gallus M, Kwok DW, Hong C, Chen LH, Carrera DA, Zhang MY, Stevers NO, Maldonado GC, Yamamichi A, Watchmaker P, Naik A, Shai A, Phillips JJ, Chang SM, Wiita AP, Wells JA, Costello JF, Diaz AA, Okada H. Challenges in the discovery of tumor-specific alternative splicing-derived cell-surface antigens in glioma. bioRxiv [Preprint](https://www.biorxiv.org/content/10.1101/2023.10.26.564156v2.full). 2023 Nov 2:2023.10.26.564156. doi: 10.1101/2023.10.26.564156. PMID: [37961484](https://pubmed.ncbi.nlm.nih.gov/37961484/); PMCID: PMC10634890.  
  
  
## Primary contact: 
  
**Takahide Nejo, MD, PhD**  
Postdoctoral Scholar  
University of California, San Francisco, Department of Neurological Surgery  
takahide.nejo@ucsf.edu  
  
  
**Hideho Okada, MD, PhD**  
Professor  
University of California, San Francisco, Department of Neurological Surgery  
hideho.okada@ucsf.edu  
  
  
## Questions about the code:  
  
**Takahide Nejo, MD, PhD**  
takahide.nejo@ucsf.edu

