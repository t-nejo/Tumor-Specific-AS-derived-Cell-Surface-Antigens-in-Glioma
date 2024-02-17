# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco

# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# load packages -----------------------------------------------------------

library(tidyverse)
library(tidylog)


# dir ---------------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.2 <- "~/path/to/output_folder"


# load data -----------------------------------------------------

setwd(dir.1)
in.f <- filename # "tcga_glioma_annotation_abs_purity_060_n429.tsv" # 429 TCGA cases to analyze
tcga.annot <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # e.g., "all_the_distinct_junc_ids_min_10_in_TCGA.rds"
junc.ids.total.1 <- readRDS(out.f)


# check -------------------------------------------------------------------

tcga.annot %>% dim() %>% print()
# [1] 429   8

tcga.annot %>% head() %>% print()
# # A tibble: 6 × 8
# case         study    cls   grade histology  IDH   codel abs.purity
# <chr>        <chr>    <chr> <chr> <chr>      <chr> <chr>      <dbl>
# 1 TCGA-CS-4938 TCGA-LGG IDH-A G2    astrocyto… Muta… non-…       0.79
# 2 TCGA-CS-4941 TCGA-LGG IDHwt G3    astrocyto… WT    non-…       0.61
# 3 TCGA-CS-4942 TCGA-LGG IDH-A G3    astrocyto… Muta… non-…       0.76
# 4 TCGA-CS-4943 TCGA-LGG IDH-A G3    astrocyto… Muta… non-…       0.83
# 5 TCGA-CS-4944 TCGA-LGG IDH-A G2    astrocyto… Muta… non-…       0.74
# 6 TCGA-CS-5390 TCGA-LGG IDH-O G2    oligodend… Muta… codel       0.85

junc.ids.total.1 %>% nrow()
# [1] 366734


# load sj.out.tab files of all the TCGA cases -------------------------------------------------------

list.sj.count.tcga <- list() # list object to store junction info identified in each case in TCGA
for(i in 1:nrow(tcga.annot)){
  print(i) ; 
  
  CASE.BARCODE <- tcga.annot %>% 
    dplyr::slice(i) %>% 
    pull(case) ; # e.g., TCGA-CS-4938
  
  setwd(dir.1) # where sj.out.tab files are located. 
  in.files <- list.files()[grep("SJ.out.tab", list.files())] 
  
  in.f <- in.files[grep(CASE.BARCODE, in.files)] 
  sj.i <- read_tsv(in.f, na = c("", "NA"), col_names = F, col_types = cols(X1 = col_character()))
  
  sj.i <- sj.i %>%
    dplyr::rename(chr = X1, int.start = X2, int.end = X3, strand = X4, int.motif = X5, annot = X6, n.uniq.map = X7, n.mult.map = X8, max.spl.overhang = X9) %>%
    mutate(strand = ifelse(strand == 2, "-", ifelse(strand == 1, "+", "undefined"))) %>%
    dplyr::filter(strand !=	"undefined") %>%
    dplyr::filter(chr != "M") %>%
    mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>%
    dplyr::filter(n.uniq.map > 0) 
  
  list.sj.count.tcga[[i]] <- sj.i %>% 
    dplyr::select(junc.id, n.uniq.map) ;
  names(list.sj.count.tcga)[i] = CASE.BARCODE
}


# prep tcga junc.count table ---------------------------------------------------

sj.count.tcga <- junc.ids.total.1 %>% 
  dplyr::select(junc.id)


# make tcga junc.count table -------------------------------------------------------------------

for(i in 1:length(list.sj.count.tcga)){
  print( round(i / length(list.sj.count.tcga), 4)) ;
  
  sj.i <- list.sj.count.tcga[[i]] 
  colnames(sj.i)[2] <- names(list.sj.count.tcga)[i]
  
  sj.count.tcga <- sj.count.tcga %>% 
    left_join(sj.i, by = "junc.id")
}

sj.count.tcga[is.na(sj.count.tcga)] <- 0


# check -------------------------------------------------------------------

sj.count.tcga %>% dim()
# [1] 366734    430

sj.count.tcga %>% anyNA()
# [1] FALSE


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # "junc_count_table_tcga_429_366734.tsv"
write_tsv(sj.count.tcga, out.f)


# si ----------------------------------------------------------------------

Sys.time() ;
# 

sessionInfo() ;
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS:   /software/c4/cbi/software/R-4.1.3-gcc8/lib64/R/lib/libRblas.so
# LAPACK: /software/c4/cbi/software/R-4.1.3-gcc8/lib64/R/lib/libRlapack.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
# [1] ggsci_2.9                         readxl_1.3.1                     
# [3] AnnotationHub_3.0.2               BiocFileCache_2.0.0              
# [5] dbplyr_2.1.1                      Rsamtools_2.8.0                  
# [7] BSgenome.Hsapiens.UCSC.hg19_1.4.3 BSgenome_1.60.0                  
# [9] rtracklayer_1.52.1                Biostrings_2.60.2                
# [11] XVector_0.32.0                    EnsDb.Hsapiens.v75_2.99.0        
# [13] ensembldb_2.16.4                  AnnotationFilter_1.16.0          
# [15] GenomicFeatures_1.44.2            AnnotationDbi_1.54.1             
# [17] Biobase_2.52.0                    GenomicRanges_1.44.0             
# [19] GenomeInfoDb_1.28.1               IRanges_2.26.0                   
# [21] S4Vectors_0.30.2                  BiocGenerics_0.38.0              
# [23] tidylog_1.0.2                     forcats_0.5.1                    
# [25] stringr_1.4.0                     dplyr_1.0.9                      
# [27] purrr_0.3.4                       readr_2.1.2                      
# [29] tidyr_1.2.0                       tibble_3.1.7                     
# [31] ggplot2_3.3.6                     tidyverse_1.3.1                  
# 
# loaded via a namespace (and not attached):
# [1] colorspace_2.0-3              rjson_0.2.20                 
# [3] ellipsis_0.3.2                fs_1.5.2                     
# [5] rstudioapi_0.13               bit64_4.0.5                  
# [7] interactiveDisplayBase_1.30.0 fansi_1.0.3                  
# [9] lubridate_1.9.2               xml2_1.3.2                   
# [11] cachem_1.0.6                  jsonlite_1.8.0               
# [13] broom_0.7.9                   png_0.1-7                    
# [15] shiny_1.7.1                   BiocManager_1.30.16          
# [17] compiler_4.1.3                httr_1.4.3                   
# [19] backports_1.2.1               assertthat_0.2.1             
# [21] Matrix_1.6-0                  fastmap_1.1.0                
# [23] lazyeval_0.2.2                cli_3.6.1                    
# [25] later_1.3.0                   htmltools_0.5.5              
# [27] prettyunits_1.1.1             tools_4.1.3                  
# [29] gtable_0.3.0                  glue_1.6.2                   
# [31] GenomeInfoDbData_1.2.6        rappdirs_0.3.3               
# [33] Rcpp_1.0.8.3                  cellranger_1.1.0             
# [35] vctrs_0.6.1                   rvest_1.0.1                  
# [37] mime_0.12                     timechange_0.2.0             
# [39] lifecycle_1.0.3               restfulr_0.0.13              
# [41] XML_3.99-0.7                  zlibbioc_1.38.0              
# [43] scales_1.2.0                  vroom_1.5.7                  
# [45] promises_1.2.0.1              clisymbols_1.2.0             
# [47] hms_1.1.0                     MatrixGenerics_1.4.3         
# [49] ProtGenerics_1.25.2           SummarizedExperiment_1.22.0  
# [51] yaml_2.3.5                    curl_4.3.2                   
# [53] memoise_2.0.0                 biomaRt_2.48.3               
# [55] stringi_1.7.6                 RSQLite_2.2.8                
# [57] BiocVersion_3.13.1            BiocIO_1.2.0                 
# [59] filelock_1.0.2                BiocParallel_1.26.2          
# [61] rlang_1.1.0                   pkgconfig_2.0.3              
# [63] bitops_1.0-7                  matrixStats_0.62.0           
# [65] lattice_0.20-44               GenomicAlignments_1.28.0     
# [67] bit_4.0.4                     tidyselect_1.1.2             
# [69] magrittr_2.0.3                R6_2.5.1                     
# [71] generics_0.1.2                DelayedArray_0.18.0          
# [73] DBI_1.1.2                     pillar_1.7.0                 
# [75] haven_2.4.3                   withr_2.5.0                  
# [77] KEGGREST_1.32.0               RCurl_1.98-1.3               
# [79] modelr_0.1.8                  crayon_1.5.1                 
# [81] utf8_1.2.2                    tzdb_0.1.2                   
# [83] progress_1.2.2                grid_4.1.3                   
# [85] blob_1.2.2                    reprex_2.0.1                 
# [87] digest_0.6.29                 xtable_1.8-4                 
# [89] httpuv_1.6.5                  munsell_0.5.0       


# end ---------------------------------------------------------------------

