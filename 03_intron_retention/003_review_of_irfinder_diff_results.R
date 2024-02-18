# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE)) ;


# load packages -----------------------------------------------------------

library(tidyverse) 
library(tidylog)


# dir --------------------------------------------------------------------

# dir.1 <- "~/path/to/data"
# dir.2 <- "~/path/to/output_folder"


# load data -----------------------------------------------------

setwd(dir.1)
in.f <- filename # "average_surfaceome_gene_expression_in_TCGA_glioma.tsv"     
surfaceome.tcga.exp.data.mean <- read_tsv(in.f)


# load the results of irfinder-diff  --------------------------------------

res.list <- list()
for(i in 1:3){
  print(i)
  test.name.i <- c("IDHwt_n166", "IDH-A_n140", "IDH-O_n123")[i]
  test.name.i.out <- gsub("-", "", test.name.i) 
  
  setwd(dir.1)
  files.i <- list.files()[grepl(test.name.i.out, list.files())]
  
  in.f <- files.i[1]
  res.gtex <- read_tsv(in.f) # e.g., "..._01_gtex_cns.tsv"
  
  in.f <- files.i[2]
  res.tcga <- read_tsv(in.f) # "..._02_tcga_n429.tsv"
  
  in.f <- files.i[3]
  res.diff <- read_tsv(in.f) # "..._03_diff.tsv"
  
  
  # filter: focus on IR events with statistical significance ---------------------
  
  PADJ <- 0.05
  
  res.diff.1 <- res.diff %>% 
    dplyr::filter(padj < PADJ & log2FoldChange > 0) %>% # TCGA > GTEx
    arrange(padj, desc(IR.change))
  
  
  # edit the table -------------------------------------------------------
  
  res.diff.2 <- res.diff.1 %>% 
    mutate(symbol = sapply(str_split(annot, "/"), "[[", 1)) %>% # edit - add new columns
    mutate(ensg = sapply(str_split(annot, "/"), "[[", 2)) %>% 
    mutate(coord = gsub(" ", "", sapply(str_split(annot, "/"), "[[", 4))) %>% 
    mutate(junc.id = paste0( 
      "chr", 
      sapply(str_split(coord, ":"), "[[", 1), 
      ":", 
      sapply(str_split(coord, ":"), "[[", 3),
      ":", 
      sapply(str_split(coord, ":"), "[[", 2)
    )) %>% 
    left_join(res.tcga %>% dplyr::select(annot, IRR.tcga = IRratio), by = "annot") %>% 
    left_join(res.gtex %>% dplyr::select(annot, IRR.gtex = IRratio), by = "annot") %>% 
    dplyr::select(symbol, ensg, junc.id, log2FC = log2FoldChange, pvalue, padj, IR.change, IRR.tcga, IRR.gtex)
  
  res.list[[i]] <- res.diff.2
}


# check -------------------------------------------------------------------

for(i in 1:3){
  I <- formatC(i, width = 2, flag = 0) ;
  test.name.i <- c("IDHwt_n166", "IDH-A_n140", "IDH-O_n123")[i] ;
  test.name.i.out <- gsub("-", "", test.name.i) ; 
  
  print(I) ; print(test.name.i) ; 
  
  res.list[[i]] %>% dim() %>% print() ; 
}
# [1] "01"
# [1] "IDHwt_n166"
# [1] 875   9
# [1] "02"
# [1] "IDH-A_n140"
# [1] 872   9
# [1] "03"
# [1] "IDH-O_n123"
# [1] 841   9


# integrate 4 tables ------------------------------------------------------

for(i in 1:4){
  I <- formatC(i, width = 2, flag = 0)
  test.name.i <- c("IDHwt", "IDH-A", "IDH-O", "GTEx")[i]
  test.name.i.out <- gsub("-", "", test.name.i)
  
  print(I) ; print(test.name.i)
  
  if(i < 4){
    res.i <- res.list[[i]] %>% 
      dplyr::select(symbol, ensg, junc.id, IRR.tcga)
    colnames(res.i)[4] <- test.name.i 
  }else if(i == 4){
    res.i <- bind_rows(res.list[[1]], res.list[[2]], res.list[[3]]) %>% 
      dplyr::select(symbol, ensg, junc.id, IRR.gtex) %>% 
      arrange(junc.id) %>% 
      distinct(symbol, ensg, junc.id, .keep_all = T)
    colnames(res.i)[4] <- test.name.i 
  };
  
  if(i == 1){
    res.merged <- res.i
  }else{
    res.merged <- res.merged %>% 
      full_join(res.i, by = c("symbol" = "symbol", "ensg" = "ensg", "junc.id" = "junc.id"))
  }
}


# check -------------------------------------------------------------------

res.merged %>% dim() %>% print() ; 
# [1] 893   7

res.merged %>% summary() %>% print() ;
# symbol              ensg             junc.id              IDHwt       
# Length:893         Length:893         Length:893         Min.   :0.2720  
# Class :character   Class :character   Class :character   1st Qu.:0.4699  
# Mode  :character   Mode  :character   Mode  :character   Median :0.4819  
# Mean   :0.4764  
# 3rd Qu.:0.4880  
# Max.   :0.5000  
# NA's   :18      
#      IDH-A            IDH-O             GTEx        
#  Min.   :0.4361   Min.   :0.4194   Min.   :0.01577  
#  1st Qu.:0.4929   1st Qu.:0.4959   1st Qu.:0.12761  
#  Median :0.4964   Median :0.5000   Median :0.18446  
#  Mean   :0.4962   Mean   :0.4975   Mean   :0.17702  
#  3rd Qu.:0.5000   3rd Qu.:0.5000   3rd Qu.:0.23002  
#  Max.   :0.5000   Max.   :0.5000   Max.   :0.32252  
#  NA's   :21       NA's   :52       


# narrow down to the expressed, surfaceome genes (mean TPM ≥ 10, n = 737 genes) -----------------------------------------

surfaceome.exp <- surfaceome.exp %>% 
  dplyr::filter(mean.tpm >= 10)

# check
surfaceome.exp %>% nrow()
# [1] 737


# filter ------------------------------------------------------------------

res.merged.surf.exp <- res.merged %>% 
  semi_join(surfaceome.exp, by = "symbol")

# semi_join: added no columns
# > rows only in x  (  0)
# > rows only in y  (711)
# > matched rows      65
# >                 =====
# > rows total        65

res.merged.surf.exp <- res.merged.surf.exp %>% 
  arrange(GTEx) %>% 
  arrange(desc(IDHwt), desc(`IDH-A`), desc(`IDH-O`)) %>% 
  arrange(symbol)


# check -------------------------------------------------------------------

res.merged.surf.exp %>% dim()
# [1] 55  7


# save --------------------------------------------------------------------

setwd(dir.2) ;
out.f <- filename # "irfinder_diff_res_merged_surf_exp_filtered.tsv"
write_tsv(res.merged.surf.exp, out.f)


# examine the table further -----------------------------------------------

res.merged.surf.exp %>% 
  pull(symbol) %>% 
  unique()
# [1] "ABCA2"    "APLP1"    "APP"      "ATP1A1"   "ATP1A2"   "ATP1A3"   "ATP1B1"  
# [8] "ATP1B2"   "BAI2"     "BSG"      "CALY"     "CD63"     "CD74"     "FAIM2"   
# [15] "HLA-B"    "HLA-C"    "HLA-DRB1" "ITM2B"    "ITM2C"    "LRP1"     "PLP1"    
# [22] "SHISA4"   "SLC22A17" "TSPAN7"  

res.merged.surf.exp %>% 
  pull(GTEx) %>% 
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1123  0.1611  0.2021  0.1983  0.2351  0.2841 

res.merged.surf.exp %>% 
  dplyr::filter(GTEx < 0.01) %>% 
  nrow()
# [1] 0


# si ----------------------------------------------------------------------

Sys.time()
# 

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
# [1] tidylog_1.0.2               DESeq2_1.32.0               SummarizedExperiment_1.22.0
# [4] Biobase_2.52.0              MatrixGenerics_1.4.3        matrixStats_0.62.0         
# [7] GenomicRanges_1.44.0        GenomeInfoDb_1.28.1         IRanges_2.26.0             
# [10] S4Vectors_0.30.2            BiocGenerics_0.38.0         forcats_0.5.1              
# [13] stringr_1.4.0               dplyr_1.0.9                 purrr_0.3.4                
# [16] readr_2.1.2                 tidyr_1.2.0                 tibble_3.1.7               
# [19] ggplot2_3.3.6               tidyverse_1.3.1            
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
# [55] locfit_1.5-9.4         pillar_1.7.0           clisymbols_1.2.0      
# [58] geneplotter_1.70.0     reprex_2.0.1           XML_3.99-0.7          
# [61] glue_1.6.2             modelr_0.1.8           vctrs_0.6.1           
# [64] png_0.1-7              tzdb_0.1.2             cellranger_1.1.0      
# [67] gtable_0.3.0           assertthat_0.2.1       cachem_1.0.6          
# [70] xtable_1.8-4           broom_0.7.9            survival_3.2-12       
# [73] AnnotationDbi_1.54.1   memoise_2.0.0          timechange_0.2.0      
# [76] ellipsis_0.3.2        


# end ---------------------------------------------------------------------
