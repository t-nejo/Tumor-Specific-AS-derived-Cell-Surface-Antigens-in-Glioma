# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# load packages -----------------------------------------------------------

library(tidyverse) 
library(tidylog)


# dir --------------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.2 <- "~/path/to/output_folder"


# load data -----------------------------------------------------

setwd(dir.1) ;
in.f <- filename # all_surfaceome_junctions_in_TCGA.rds" # the table containing all the pass-filter junc ids
res.orig <- readRDS(in.f) # or read_tsv(in.f)

setwd(dir.1)
in.f <- filename # e.g., "all_the_distinct_junc_ids_min_10_in_TCGA.rds"
junc.ids.total.1 <- readRDS(in.f) # or read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "junc_count_table_tcga_429_366734.tsv"
sj.count.tcga <- read_tsv(in.f)


# check -------------------------------------------------------------------

res.orig %>% dim()
# [1] 9041    7

res.orig %>% head()
# # A tibble: 6 × 7
# symbol ensg            junc.id                    chr   int.start   int.end strand
# <chr>  <chr>           <chr>                      <chr>     <dbl>     <dbl> <chr> 
#   1 PTPRF  ENSG00000142949 chr1:+:44086903-44087611   1      44086904  44087611 +     
#   2 F3     ENSG00000117525 chr1:-:94997968-94998645   1      94997969  94998645 -     
#   3 S1PR1  ENSG00000170989 chr1:+:101702655-101707244 1     101702656 101707244 +     
#   4 S1PR1  ENSG00000170989 chr1:+:101707328-101707583 1     101707329 101707583 +     
#   5 CELSR2 ENSG00000143126 chr1:+:109816292-109816782 1     109816293 109816782 +     
#   6 CELSR2 ENSG00000143126 chr1:+:109816292-109817037 1     109816293 109817037 +   


# prep --------------------------------------------------------------------

junc.ids.to.test <- res.orig %>% 
  pull(junc.id)

junc.ids.to.test %>% length()
# [1] 9041


# make overlap.table step 1----------------------------------------------

list.junc.ids.to.test <- list()
for(i in 1:length(junc.ids.to.test)){
  print( round(i / length(junc.ids.to.test), 4))
  
  JUNC.ID <- junc.ids.to.test[i]
  
  CHR <- gsub(":.*", "", JUNC.ID)
  STRAND <- gsub(":.*", "", gsub("chr..:", "", gsub("chr.:", "", JUNC.ID)))
  START <- gsub("-.*", "", gsub(".*:", "", JUNC.ID)) %>% as.integer()
  END <- gsub(".*-", "", gsub(".*:", "", JUNC.ID)) %>% as.integer()
  
  JUNC.ID.O <- sj.annot %>% 
    dplyr::filter(chr == gsub("chr", "", CHR) & strand == STRAND & int.start < END & START < int.end) %>% 
    pull(junc.id) ; 

  if(length(JUNC.ID.O) == 0){
    tmp <- tibble(
      junc.id = JUNC.ID, 
      junc.id.overlap = NA
    )
  }else{
    tmp <- tibble(
      junc.id = JUNC.ID, 
      junc.id.overlap = JUNC.ID.O
    )
  }
  list.junc.ids.to.test[[i]] <- tmp
}


# make overlap.table step 1----------------------------------------------

junc.ids.to.check <- NULL
for(i in 1:length(list.junc.ids.to.test)){
  print(i)
  
  junc.ids.i <- list.junc.ids.to.test[[i]] %>% 
    gather("label", "junc.id", 1:2) %>% 
    pull(junc.id)
  junc.ids.to.check <- c(junc.ids.to.check, junc.ids.i) 
}


# check -------------------------------------------------------------------

junc.ids.to.check %>% length()
# [1] 38390

junc.ids.to.check <- junc.ids.to.check[!is.na(junc.ids.to.check)] %>% unique()
junc.ids.to.check %>% length()
# [1] 15354

sj.count.tcga %>% dim()
# [1] 366734    430

sj.count.tcga %>% anyNA()
# [1] FALSE


# retain the necessary lines only -----------------------------------------

sj.count.tcga.2 <- sj.count.tcga %>% 
  semi_join(tibble(junc.id = junc.ids.to.check), by = "junc.id")


# check -------------------------------------------------------------------

sj.count.tcga.2 %>% dim()
# [1] 14028   430


# identify the “most-dominant” junctions among multiple overlaps for each candidate event --------------


overlap.table <- NULL # blank object
for(i in 1:length(junc.ids.to.test)){
  print( round(i/length(junc.ids.to.test), 4))
  
  overlap.table.pre.i <- list.junc.ids.to.test[[i]]
  
  if(nrow(overlap.table.pre.i ) == 1){
    overlap.table <- overlap.table %>% 
      bind_rows(overlap.table.pre.i)
  }else{
    sj.tmp <- overlap.table.pre.i %>% 
      left_join(sj.count.tcga.2, by = c("junc.id.overlap" = "junc.id"))
    
    sj.tmp <- sj.tmp %>%
      gather("barcode", "count", -c(1:2)) %>% 
      group_by(junc.id, junc.id.overlap) %>% 
      summarize(sum = sum(count)) %>% 
      arrange(desc(sum)) %>% 
      head(n = 1) %>% 
      dplyr::select(1:2) %>% 
      ungroup()
    
    overlap.table = overlap.table %>% 
      bind_rows(sj.tmp)
  }
}


# check -------------------------------------------------------------------

overlap.table %>% dim()
# [1] 9041    2


# edit sj.count --------------------------------------------------------------------

sj.count.tcga.3 <- sj.count.tcga.2 %>% # to be used later
  semi_join(overlap.table %>% 
               gather("label", "junc.id", 1:2) %>% 
               dplyr::filter(!is.na(junc.id)) %>% 
               distinct(junc.id), 
             by = "junc.id"
  )


# check -------------------------------------------------------------------

sj.count.tcga.2 %>% dim()
# [1] 14028   430

sj.count.tcga.3 %>% dim()
# [1] 12271   430


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # "junc_count_table_tcga_429_14028.tsv"
write_tsv(sj.count.tcga.2, out.f)

out.f <- filename # "junc_count_table_tcga_429_12271.tsv"
write_tsv(sj.count.tcga.3, out.f)

out.f <- filename # "overlap_table_9041.tsv"
write_tsv(overlap.table, out.f)


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

