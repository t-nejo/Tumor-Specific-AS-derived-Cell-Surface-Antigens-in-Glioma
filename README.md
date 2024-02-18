# Tumor-Specific AS-derived Cell-Surface Antigens in Glioma

## Introduction  

This repository contains the code used to analyze the data in Nejo et al. bioRxiv (2023) "Challenges in the discovery of tumor-specific alternative splicing-derived cell-surface antigens in glioma" (doi: [10.1101/2023.10.26.564156](https://www.biorxiv.org/content/10.1101/2023.10.26.564156v2.full)). In this study, we investigated tumor-specific alternative splicing-derived, cell-surface antigens in glioma, through the analyses of transcriptome data of TCGA-GBM/LGG, GTEx, Mayo-GBM-PDX, as well as the clinical samples of the University of California, San Francisco (UCSF) Brain Tumor Center. Moreover, we conducted nanopore full-length transcript sequencing and proteomics analysis of the CPTAC-GBM project. Our investigation illustrated the diverse characteristics of the tumor-specific AS events and the challenges of antigen exploration due to their notable spatiotemporal heterogeneity and elusive nature at the protein levels. 

  
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

  
## Key analysis tools:  

[STAR](https://github.com/alexdobin/STAR): Short-read RNA-seq read alignment and junction detection  
[IRFinder](https://github.com/williamritchie/IRFinder): Intron-retention analysis  
[Guppy basecaller](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview): Nanopore long read sequencing base calling  
[FLAIR](https://github.com/BrooksLabUCSC/flair): Nanopore long read sequencing FASTQ analysis  
[MSGF+](https://github.com/MSGFPlus/msgfplus): Peptide spectra database search  
[MASIC](https://github.com/PNNL-Comp-Mass-Spec/MASIC): Quantification of 11-TMT reporter ions  


## Data availability:  
Nanopore full-length transcript amplicon sequencing dataset generated in this study will be made available through the NCBI Gene Expression Omnibus (GEO) website upon manuscript acceptance. RNA-seq data of the TCGA ([phs000178.v10.p8](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v10.p8)) and the GTEx ([phs000424.v9.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v9.p2)) datasets were obtained through the GDC Data Portal with an approved dbGaP data access. RNA-seq data of Mayo Clinic GBM PDX was obtained through the NCBI GEO (accession, [PRJNA548556](https://www.mayo.edu/research/labs/translational-neuro-oncology/mayo-clinic-brain-tumor-patient-derived-xenograft-national-resource/pdx-characteristics/cbioportal)). The CPTAC-GBM proteomics data was obtained through NIH Proteomic Data Commons (Study ID [PDC000204](https://proteomic.datacommons.cancer.gov/pdc/study/PDC000204)). The RNA-seq data of the spatially mapped dataset have been deposited in the European Genome-phenome Archive (EGA) under accession number [EGAS00001003710](https://ega-archive.org/studies/EGAS00001003710). The data for the longitudinally collected dataset have been deposited in the EGA under accession number [EGAS00001002368](https://ega-archive.org/studies/EGAS00001002368) and dbGaP study accession: [phs002034.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002034.v1.p1).


## Key data visualization tools:  
[BioRender](https://www.biorender.com/)  
[Affinity Designer](https://affinity.serif.com/en-us/designer/)  
[R package ggplot2](https://ggplot2.tidyverse.org/)  

## Citation

**If you use our data or code for your work, please consider citing the following publication:**  

Nejo T, Wang L, Leung KK, Wang A, Lakshmanachetty S, Gallus M, Kwok DW, Hong C, Chen LH, Carrera DA, Zhang MY, Stevers NO, Maldonado GC, Yamamichi A, Watchmaker PB, Naik A, Shai A, Phillips JJ, Chang SM, Wiita AP, Wells JA, Costello JF, Diaz AA, Okada H. Challenges in the discovery of tumor-specific alternative splicing-derived cell-surface antigens in glioma. bioRxiv [Preprint](https://www.biorxiv.org/content/10.1101/2023.10.26.564156v2.full). 2023 Nov 2:2023.10.26.564156. doi: 10.1101/2023.10.26.564156. PMID: [37961484](https://pubmed.ncbi.nlm.nih.gov/37961484/); PMCID: PMC10634890.  
  
  


