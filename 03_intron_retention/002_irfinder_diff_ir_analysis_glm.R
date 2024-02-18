# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE)) ;


# load packages -----------------------------------------------------------

library(tidyverse) 
library(tidylog)
library(DESeq2)


# dir --------------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.2 <- "~/path/to/output_folder"


# load functions from IRFinder ------------------------------------------

source("/path_to_irfinder/IRFinder-1.2.3/bin/DESeq2Constructor.R")


# load data -----------------------------------------------------

setwd(dir.1)
in.f <- filename # "tcga_glioma_annotation_abs_purity_060_n429.tsv"
tcga.annot <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "experiment.txt"
experiment <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "filepaths.txt" - manually prepared
path.info <- read_tsv(in.f)


# check -------------------------------------------------------------------

tcga.annot %>% dim()
# [1] 429   8

tcga.annot %>% head()
# # A tibble: 6 × 8
# case         study    cls   grade histology         IDH    codel     abs.purity
# <chr>        <chr>    <chr> <chr> <chr>             <chr>  <chr>          <dbl>
# 1 TCGA-CS-4938 TCGA-LGG IDH-A G2    astrocytoma       Mutant non-codel       0.79
# 2 TCGA-CS-4941 TCGA-LGG IDHwt G3    astrocytoma       WT     non-codel       0.61
# 3 TCGA-CS-4942 TCGA-LGG IDH-A G3    astrocytoma       Mutant non-codel       0.76
# 4 TCGA-CS-4943 TCGA-LGG IDH-A G3    astrocytoma       Mutant non-codel       0.83
# 5 TCGA-CS-4944 TCGA-LGG IDH-A G2    astrocytoma       Mutant non-codel       0.74
# 6 TCGA-CS-5390 TCGA-LGG IDH-O G2    oligodendroglioma Mutant codel           0.85

experiment %>% head()
# # A tibble: 6 × 2
# SampleNames Condition
# <chr>       <chr>    
# 1 SRR1069188  gtex_cns 
# 2 SRR1070986  gtex_cns 
# 3 SRR1075433  gtex_cns 
# 4 SRR1079591  gtex_cns 
# 5 SRR1080172  gtex_cns 
# 6 SRR1082262  gtex_cns 

experiment %>% pull(Condition) %>% unique()
# [1] "gtex_cns"  "tcga_n429"

path.info %>% 
  mutate(X1 = gsub(".*/", "path_to_file/", X1)) %>% 
  head()
# # A tibble: 6 × 1
# X1                                            
# <chr>                                         
# 1 path_to_file/SRR1069188_IRFinder-IR-nondir.txt
# 2 path_to_file/SRR1070986_IRFinder-IR-nondir.txt
# 3 path_to_file/SRR1075433_IRFinder-IR-nondir.txt
# 4 path_to_file/SRR1079591_IRFinder-IR-nondir.txt
# 5 path_to_file/SRR1080172_IRFinder-IR-nondir.txt
# 6 path_to_file/SRR1082262_IRFinder-IR-nondir.txt


# extract path info ------------------------------------------------------------

paths <- path.info %>% 
  pull(X1)



# loop --------------------------------------------------------------------

for (i in 1:3){
  print(i)
  
  
  # test.subject ------------------------------------------------------------
  
  test.name.i <- c("IDHwt_n166", "IDH-A_n140", "IDH-O_n123")[i]
  
  cases.i <- list(
    tcga.annot %>% dplyr::filter(cls == "IDHwt") %>% pull(case), 
    tcga.annot %>% dplyr::filter(cls == "IDH-A") %>% pull(case), 
    tcga.annot %>% dplyr::filter(cls == "IDH-O") %>% pull(case)
  )[[i]]
  
  # extract the cases from paths and experiment -----------------------------
  
  paths.i <- path.info %>% 
    mutate(X2 = gsub(".*\\/", "", X1)) %>% 
    mutate(X2 = gsub("_IRFinder-.*", "", X2)) %>% 
    mutate(X2 = substr(X2, 1, 12)) %>% 
    dplyr::filter(grepl("SRR", X2) | is.element(X2, cases.i)) %>% 
    pull(X1)
  
  experiment.i <- experiment %>% 
    dplyr::filter(Condition == "gtex_cns" | is.element(SampleNames, cases.i)) %>% 
    as.data.frame()
  
  experiment.i$Condition = factor(experiment.i$Condition, levels = c("gtex_cns", "tcga_n429")) # the first condition ("gtex_cns") is set as the baseline in the analysis
  rownames(experiment.i) = NULL  # force removing rownames
  
  
  # check -------------------------------------------------------------------
  
  test.name.i %>% print()
  # [1] "IDHwt_n166"
  
  cases.i %>% length() %>% print()
  # [1] 166
  
  
  # run IRFinder-diff step 1 ---------------------------------------------------------------------
  
  # RUN DESeq2
  
  metaList <- DESeqDataSetFromIRFinder(
    filePaths = paths.i, 
    designMatrix = experiment.i, 
    designFormula = ~1
  )
  
  # extract DESeq2 object with normalization factors ready
  
  dds <- metaList$DESeq2Object                       
  
  
  # check -------------------------------------------------------------------
  
  colData(dds) # Check design of matrix
  
  
  # run IRFinder-diff step 2 --------------------------------------------------------------------
  
  # Build a formula of GLM. 
  design(dds) <- ~ Condition + Condition:IRFinder     
  
  # Estimate parameters and fit to model
  dds <- DESeq(dds)                                  
  
  
  # check -------------------------------------------------------------------
  
  # Check the actual variable name assigned by DESeq2
  resultsNames(dds)                                 
  # [1] "Intercept"                       "Condition_tcga_n429_vs_gtex_cns"
  # [3] "Conditiongtex_cns.IRFinderIR"    "Conditiontcga_n429.IRFinderIR" 
  
  
  # result ------------------------------------------------------------------
  
  # IR-Ratio - GTEx
  
  res.gtex_cns <- results(dds, name = "Conditiongtex_cns.IRFinderIR")
  
  gtex_cns.IR_vs_Splice <- 2^res.gtex_cns$log2FoldChange
  
  IRratio.gtex_cns <- gtex_cns.IR_vs_Splice/(1 + gtex_cns.IR_vs_Splice)
  
  res.gtex_cns.2 <- tibble(
    annot = res.gtex_cns %>% 
      rownames() 
  ) %>% 
    bind_cols(res.gtex_cns %>% as_tibble() ) %>% 
    mutate(IR_vs_Splice = gtex_cns.IR_vs_Splice) %>% 
    mutate(IRratio = IRratio.gtex_cns) %>% 
    dplyr::filter(!is.na(pvalue)) %>% 
    arrange(pvalue)
  
  
  # IR-Ratio - TCGA
  
  res.tcga_n429 <- results(dds, name = "Conditiontcga_n429.IRFinderIR")
  
  tcga_n429.IR_vs_Splice <- 2^res.tcga_n429$log2FoldChange
  
  IRratio.tcga_n429 <- tcga_n429.IR_vs_Splice/(1 + tcga_n429.IR_vs_Splice)
  
  res.tcga_n429.2 = tibble(
    annot = res.tcga_n429 %>% 
      rownames() 
  ) %>% 
    bind_cols(res.tcga_n429 %>% as_tibble() ) %>% 
    mutate(IR_vs_Splice = tcga_n429.IR_vs_Splice) %>% 
    mutate(IRratio = IRratio.tcga_n429) %>% 
    dplyr::filter(!is.na(pvalue)) %>% 
    arrange(pvalue)
  
  
  # diff --------------------------------------------------------------------
  
  res.diff <- results(dds, contrast = list("Conditiontcga_n429.IRFinderIR", "Conditiongtex_cns.IRFinderIR")) ;
  
  IR.change <- IRratio.tcga_n429 - IRratio.gtex_cns ;
  
  annot <- res.diff %>% rownames()
  
  res.diff.2 <- tibble(annot = annot) %>% 
    bind_cols(res.diff %>% as_tibble() ) %>% 
    mutate(IR.change = IR.change) %>% 
    dplyr::filter(!is.na(pvalue)) %>% 
    arrange(pvalue)
  
  
  # save output ---------------------------------------------------------------------
  
  test.name.i.out <- gsub("-", "", test.name.i)
  I <- formatC(i, width = 2, flag = 0)
  
  setwd(dir.2) ;
  out.f <- filename # paste0("irfinder_diff_glm_", I, "_", test.name.i.out, "_01_gtex_cns.tsv")
  write_tsv(res.gtex_cns.2, out.f) 
  
  out.f <- filename # paste0("irfinder_diff_glm_", I, "_", test.name.i.out, "_02_tcga_n429.tsv")
  write_tsv(res.tcga_n429.2, out.f) 
  
  out.f <- filename # paste0("irfinder_diff_glm_", I, "_", test.name.i.out, "_03_diff.tsv")
  write_tsv(res.diff.2, out.f) 
}


# si ----------------------------------------------------------------------

Sys.time()

sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /software/c4/cbi/software/R-4.1.3-gcc8/lib64/R/lib/libRblas.so
# LAPACK: /software/c4/cbi/software/R-4.1.3-gcc8/lib64/R/lib/libRlapack.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
# [1] DESeq2_1.32.0               SummarizedExperiment_1.22.0 Biobase_2.52.0             
# [4] MatrixGenerics_1.4.3        matrixStats_0.62.0          GenomicRanges_1.44.0       
# [7] GenomeInfoDb_1.28.1         IRanges_2.26.0              S4Vectors_0.30.2           
# [10] BiocGenerics_0.38.0         forcats_0.5.1               stringr_1.4.0              
# [13] dplyr_1.0.9                 purrr_0.3.4                 readr_2.1.2                
# [16] tidyr_1.2.0                 tibble_3.1.7                ggplot2_3.3.6              
# [19] tidyverse_1.3.1            
# 
# loaded via a namespace (and not attached):
# [1] bitops_1.0-7           fs_1.5.2               lubridate_1.9.2       
# [4] bit64_4.0.5            RColorBrewer_1.1-3     httr_1.4.3            
# [7] tools_4.1.3            backports_1.2.1        utf8_1.2.2            
# [10] R6_2.5.1               DBI_1.1.2              colorspace_2.0-3      
# [13] withr_2.5.0            tidyselect_1.1.2       bit_4.0.4             
# [16] compiler_4.1.3         cli_3.6.1              rvest_1.0.1           
# [19] xml2_1.3.2             DelayedArray_0.18.0    scales_1.2.0          
# [22] genefilter_1.74.0      XVector_0.32.0         pkgconfig_2.0.3       
# [25] dbplyr_2.1.1           fastmap_1.1.0          rlang_1.1.0           
# [28] readxl_1.3.1           rstudioapi_0.13        RSQLite_2.2.8         
# [31] generics_0.1.2         jsonlite_1.8.0         vroom_1.5.7           
# [34] BiocParallel_1.26.2    RCurl_1.98-1.3         magrittr_2.0.3        
# [37] GenomeInfoDbData_1.2.6 Matrix_1.6-0           Rcpp_1.0.8.3          
# [40] munsell_0.5.0          fansi_1.0.3            lifecycle_1.0.3       
# [43] stringi_1.7.6          zlibbioc_1.38.0        grid_4.1.3            
# [46] blob_1.2.2             crayon_1.5.1           lattice_0.20-44       
# [49] Biostrings_2.60.2      haven_2.4.3            splines_4.1.3         
# [52] annotate_1.70.0        hms_1.1.0              KEGGREST_1.32.0       
# [55] locfit_1.5-9.4         pillar_1.7.0           geneplotter_1.70.0    
# [58] reprex_2.0.1           XML_3.99-0.7           glue_1.6.2            
# [61] modelr_0.1.8           vctrs_0.6.1            png_0.1-7             
# [64] tzdb_0.1.2             cellranger_1.1.0       gtable_0.3.0          
# [67] assertthat_0.2.1       cachem_1.0.6           xtable_1.8-4          
# [70] broom_0.7.9            survival_3.2-12        AnnotationDbi_1.54.1  
# [73] memoise_2.0.0          timechange_0.2.0       ellipsis_0.3.2        


# end ---------------------------------------------------------------------


