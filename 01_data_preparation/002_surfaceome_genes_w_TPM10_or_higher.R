# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco

# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# load packages -----------------------------------------------------------

library(tidyverse) ;
library(tidylog) ;


# dir ---------------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.2 <- "~/path/to/output_folder"


# load data -----------------------------------------------------------------

setwd(dir.1)
in.f <- # "TcgaTargetGtex_rsem_gene_tpm_log2_0001_TCGA_n656.tsv"
tcga.exp.data <- read_tsv(in.f)

## "RNA-seq gene expression abundance data of both the TCGA and the GTEx datasets were downloaded through the UCSC Xena-Toil web portal (dataset ID: TcgaTargetGtex_rsem_gene_tpm; filename: ‘TcgaTargetGtex_rsem_gene_tpm.gz’; version: 2016-09-03) [52]."
## 52. Vivian J, Rao AA, Nothaft FA, Ketchum C, Armstrong J, Novak A, et al. Toil enables reproducible, open source, big biomedical data analyses. Nat Biotechnol. 2017;35:314–6.

setwd(dir.1)
in.f <- # "07032020_tcga_glioma_annotation_n663.tsv"
tcga.annot <- read_tsv(in.f)

## "Clinical, pathological, molecular diagnosis information, and ABSOLUTE-estimated tumor purity data were gathered for the TCGA dataset from the reference [32]." 
## 32. Ceccarelli M, Barthel FP, Malta TM, Sabedot TS, Salama SR, Murray BA, et al. Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Cell. 2016;164:550–63. 


setwd(dir.1)
in.f <- # "pnas.1808790115.sd01_11.7.surfaceome.txt"
surfaceome <- read_tsv(in.f)

## To investigate the genes encoding cell-surface proteins, we selected 2,886 “surfaceome” genes, referring to [28]. 
## 28. Bausch-Fluck D, Goldmann U, Müller S, van Oostrum M, Müller M, Schubert OT, et al. The in silico human surfaceome. Proc Natl Acad Sci U S A. 2018;115:E10988–97.


# check------------------------------------------------------------

tcga.exp.data %>% dim()
# [1] 60498   657

tcga.exp.data[, 1:4] %>% head()
# # A tibble: 6 × 4
# sample             `TCGA-19-1787-01` `TCGA-S9-A7J2-01` `TCGA-E1-A7YI-01`
# <chr>                          <dbl>             <dbl>             <dbl>
# 1 ENSG00000242268.2              -9.97             0.300            -0.452
# 2 ENSG00000259041.1              -9.97            -9.97             -9.97 
# 3 ENSG00000270112.3              -3.82            -3.05             -0.735
# 4 ENSG00000167578.16              5.30             4.89              5.76 
# 5 ENSG00000278814.1              -9.97            -9.97             -9.97 
# 6 ENSG00000078237.5               3.51             2.30              4.23 

tcga.exp.data %>% anyNA() %>% print()
# [1] FALSE

tcga.exp.data %>% colnames() %>% head()
# [1] "sample"          "TCGA-19-1787-01" "TCGA-S9-A7J2-01" "TCGA-E1-A7YI-01" "TCGA-06-5412-01"
# [6] "TCGA-DU-7302-01"


tcga.annot %>% dim()
# [1] 663   8

tcga.annot %>% anyNA() %>% print()
# [1] TRUE

tcga.annot %>% head()
# # A tibble: 6 × 8
# case         study    cls   grade histology         IDH    codel   abs.purity
# <chr>        <chr>    <chr> <chr> <chr>             <chr>  <chr>        <dbl>
# 1 TCGA-CS-4938 TCGA-LGG IDH-A G2    astrocytoma       Mutant non-co…       0.79
# 2 TCGA-CS-4941 TCGA-LGG IDHwt G3    astrocytoma       WT     non-co…       0.61
# 3 TCGA-CS-4942 TCGA-LGG IDH-A G3    astrocytoma       Mutant non-co…       0.76
# 4 TCGA-CS-4943 TCGA-LGG IDH-A G3    astrocytoma       Mutant non-co…       0.83
# 5 TCGA-CS-4944 TCGA-LGG IDH-A G2    astrocytoma       Mutant non-co…       0.74
# 6 TCGA-CS-5390 TCGA-LGG IDH-O G2    oligodendroglioma Mutant codel         0.85


surfaceome %>% dim()
# [1] 2886  39

surfaceome %>% anyNA() %>% print()
# [1] TRUE

surfaceome %>% colnames()
# [1] "UniProt name"                                 
# [2] "UniProt accession"                            
# [3] "UniProt description"                          
# [4] "UniProt gene"                                 
# [5] "Surfaceome Label"                             
# [6] "Surfaceome Label Source"                      
# [7] "Comment"                                      
# [8] "length"                                       
# [9] "TM domains"                                   
# [10] "signalpeptide"                                
# [11] "topology"                                     
# [12] "topology source"                              
# [13] "MachineLearning trainingset"                  
# [14] "SURFY score"                                  
# [15] "MachineLearning FPR class (1=1%, 2=5%, 3=15%)"
# [16] "Ensembl gene"                                 
# [17] "Ensembl protein"                              
# [18] "CD number"                                    
# [19] "Membranome Almen main-class"                  
# [20] "Membranome Almen sub-class"                   
# [21] "nxst motifs"                                  
# [22] "noncyt. nxst count"                           
# [23] "glycomineN sites"                             
# [24] "glycomineO sites"                             
# [25] "glycomineC sites"                             
# [26] "CSPA category"                                
# [27] "CSPA peptide count"                           
# [28] "CSPA peptides"                                
# [29] "CSPA N115 sites"                              
# [30] "CSPA id"                                      
# [31] "UniProt subcellular"                          
# [32] "UniProt keywords"                             
# [33] "UniProt uniref"                               
# [34] "COMPARTMENTS link"                            
# [35] "COMPARTMENTS benchmark pos"                   
# [36] "COMPARTMENTS benchmark neg"                   
# [37] "HPA antibody"                                 
# [38] "DrugBank approved drug IDs"                   
# [39] "GeneID"                               


# edit tcga.annot --------------------------------------------------------------

tcga.annot.edit <- tcga.annot %>% 
  dplyr::filter(!is.na(abs.purity)) %>% 
  dplyr::filter(abs.purity >= 0.6)
  

# check -------------------------------------------------------------------

tcga.annot.edit %>% 
  pull(cls) %>% 
  table()
# IDH-A IDH-O IDHwt 
# 140   123   166 

tcga.annot.edit %>% head()
# # A tibble: 6 × 8
# case         study    cls   grade histology         IDH    codel     abs.purity
# <chr>        <chr>    <chr> <chr> <chr>             <chr>  <chr>          <dbl>
# 1 TCGA-CS-4938 TCGA-LGG IDH-A G2    astrocytoma       Mutant non-codel       0.79
# 2 TCGA-CS-4941 TCGA-LGG IDHwt G3    astrocytoma       WT     non-codel       0.61
# 3 TCGA-CS-4942 TCGA-LGG IDH-A G3    astrocytoma       Mutant non-codel       0.76
# 4 TCGA-CS-4943 TCGA-LGG IDH-A G3    astrocytoma       Mutant non-codel       0.83
# 5 TCGA-CS-4944 TCGA-LGG IDH-A G2    astrocytoma       Mutant non-codel       0.74
# 6 TCGA-CS-5390 TCGA-LGG IDH-O G2    oligodendroglioma Mutant codel           0.85


# edit surfaceome ---------------------------------------------------------

surf.edit <- surfaceome %>% 
  dplyr::select(symbol = `UniProt gene`, ensg = `Ensembl gene`) # %>% 
  # distinct(symbol, .keep_all = T)
  

# edit gene expression data ----------------------------------------------------

# prep
colnames(tcga.exp.data)[-1] <- substr(colnames(tcga.exp.data)[-1], 1, 12) 
colnames(tcga.exp.data)[1] <- "ensg"
tcga.exp.data$ensg <- substr(tcga.exp.data$ensg, 1, 15)

tcga.exp.data.2 <- tcga.exp.data %>% 
  dplyr::select(ensg, colnames(tcga.exp.data)[is.element(colnames(tcga.exp.data), tcga.annot.edit$case)]) 

tcga.exp.data.2 <- surf.edit %>% 
  inner_join(tcga.exp.data.2, by = "ensg")
# inner_join: added no columns
# > rows only in x  (   341)
# > rows only in y  (57,966)
# > matched rows      2,545
# >                 ========
#   > rows total        2,545


# check -------------------------------------------------------------------

tcga.exp.data.2 %>% dim()
# [1] 2545  426

tcga.exp.data.2 %>% anyNA()
# [1] FALSE

tcga.exp.data.2[, 1:4] %>% head()
# # A tibble: 6 × 4
# symbol  ensg            `TCGA-19-1787` `TCGA-S9-A7J2`
# <chr>   <chr>                    <dbl>          <dbl>
# 1 SLC12A8 ENSG00000221955           1.74          0.972
# 2 ESYT3   ENSG00000158220          -2.47         -0.997
# 3 SLC5A10 ENSG00000154025          -3.82         -4.61 
# 4 CLRN2   ENSG00000249581          -9.97         -9.97 
# 5 TMEM30C ENSG00000235156          -9.97         -9.97 
# 6 ANO9    ENSG00000185101          -2.55         -5.01 


# transform to calculate average ------------------------------------------------

tcga.exp.data.3 <- tcga.exp.data.2 %>% 
  gather(barcode, log2tpm, -c(1:2)) %>% 
  mutate(tpm = 2^(log2tpm) - 0.001)


# join with annot ---------------------------------------------------------

tcga.exp.data.4 <- tcga.annot.edit %>% 
  dplyr::select(barcode = case, study, cls) %>% 
  inner_join(tcga.exp.data.3, by = "barcode")


# check -------------------------------------------------------------------

tcga.exp.data.4 %>% anyNA()
# [1] FALSE


# average TPM of each gene ------------------------------------------------

tcga.exp.data.mean <- tcga.exp.data.4 %>% 
  group_by(symbol) %>% 
  summarise(mean.tpm = mean(tpm), sd.tpm = sd(tpm)) %>% 
  arrange(desc(mean.tpm))

tcga.exp.data.mean <- tcga.exp.data.mean %>% 
  mutate(n = 1:nrow(tcga.exp.data.mean))


# check -------------------------------------------------------------------

tcga.exp.data.mean %>% nrow()
# [1] 2532

tcga.exp.data.mean %>% 
  dplyr::filter(mean.tpm >= 10) %>% 
  nrow()
# [1] 737

tcga.exp.data.mean %>% 
  dplyr::filter(mean.tpm < 10) %>% 
  nrow()
# [1] 1795

tcga.exp.data.mean %>% head(n = 10)
# # A tibble: 10 × 4
# symbol  mean.tpm sd.tpm     n
# <chr>      <dbl>  <dbl> <int>
# 1 BCAN       1505.   907.     1
# 2 CD74       1247.  1225.     2
# 3 CD63       1059.   732.     3
# 4 ITM2B       893.   312.     4
# 5 PLP1        869.   900.     5
# 6 ITM2C       802.   470.     6
# 7 PCDHGC3     708.   332.     7
# 8 HLA-B       690.   829.     8
# 9 HLA-A       669.   679.     9
# 10 APP         664.   274.    10

genes.of.interest <- c("PTPRZ1", "BCAN", "DLL3", "NCAM1", "NCAM2", "NRCAM", "TREM2")

tcga.exp.data.mean %>% 
  dplyr::filter(is.element(symbol, genes.of.interest)) %>% 
  print()
# # A tibble: 8 × 6
# symbol mean.tpm sd.tpm lab   lab.2     n
# <chr>     <dbl>  <dbl> <chr> <chr> <int>
# 1 BCAN     1505.   907.  show  show      1
# 2 PTPRZ1    589.   357.  show  show     13
# 3 NCAM1     249.   134.  show  show     29
# 4 EGFR      185.   390.  show  show     41
# 5 NRCAM     169.    86.9 show  show     45
# 6 DLL3      144.   151.  show  show     54
# 7 TREM2      82.1   64.0 show  show    106
# 8 NCAM2      22.6   14.9 show  show    428


# save ----------------------------------------------------

setwd(dir.2)
out.f <- # "average_surfaceome_gene_expression_in_TCGA_glioma.tsv"                                        
write_tsv(tcga.exp.data.mean, out.f)


# si ----------------------------------------------------------------------

Sys.time() %>% print() ; 
# 

sessionInfo() %>% print()
