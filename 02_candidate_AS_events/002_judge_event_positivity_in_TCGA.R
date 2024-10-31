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

setwd(dir.1)
in.f <- filename # "tcga_429_n9041_junc_ids_sj_count.tsv"
sj.count.tcga <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "overlap_table_9041.tsv"
overlap.table <- read_tsv(in.f)

setwd(dir.1) ; 
in.f <- filename # "tcga_glioma_annotation_abs_purity_060.txt" 
tcga.annot <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "gtex_9166_n9041_psr.tsv"
gtex.psr <- read_tsv(in.f)


# check -------------------------------------------------------------------

overlap.table %>% dim() %>% print()
# [1] 9041    2

overlap.table %>% head() %>% print()
# # A tibble: 6 × 2
# junc.id                    junc.id.overlap           
# <chr>                      <chr>                     
# 1 chr1:+:44086903-44087611   chr1:+:44086903-44087605  
# 2 chr1:-:94997968-94998645   chr1:-:94998036-94998645  
# 3 chr1:+:101702655-101707244 chr1:+:101702655-101704377
# 4 chr1:+:101707328-101707583 NA                        
# 5 chr1:+:109816292-109816782 chr1:+:109816292-109816643
# 6 chr1:+:109816292-109817037 chr1:+:109816292-109816643  

tcga.annot %>% dim()
# [1] 429   8

tcga.annot %>% head()

gtex.psr %>% dim()
# [1] 9041    4


# define the threshold values -------------------------------------------

THR.C <- 10
THR.D <- 20
THR.F <- 0.01


# obtain count, depth, freq, and judge -------------------------------------------

for(i in 1:nrow(overlap.table)){
  print(paste0(i, " (", round(i / nrow(overlap.table), 3), ")" ))
  
  JUNC.ID <- overlap.table %>% 
    dplyr::slice(i) %>% 
    pull(junc.id) ; 
  JUNC.ID.O <- overlap.table %>% 
    dplyr::slice(i) %>% 
    pull(junc.id.overlap) ; 

  if(is.na(JUNC.ID.O)){
    res.i <- sj.count.tcga %>% 
      dplyr::filter(junc.id == JUNC.ID) %>% 
      gather(case, value, -1) %>% 
      spread(junc.id, value) %>% 
      dplyr::rename(count = JUNC.ID) %>% 
      mutate(count = ifelse(is.na(count), 0, count)) %>% 
      mutate(depth = count) %>%
      mutate(freq = round(count / depth, 4)) %>% 
      mutate(freq = ifelse(is.na(freq), 0, freq)) %>% 
      mutate(judge = ifelse(count >= THR.C & depth >= THR.D & freq >= THR.F, 1, 0)) %>% 
      mutate(junc.id = JUNC.ID) %>% 
      dplyr::select(junc.id, case, count, depth, freq, judge)
  }else{
    res.i <- bind_rows(
      sj.count.tcga %>% dplyr::filter(junc.id == JUNC.ID), 
      sj.count.tcga %>% dplyr::filter(junc.id == JUNC.ID.O)
      ) %>% 
      gather(case, value, -1) %>% 
      spread(junc.id, value) %>% 
      dplyr::select(case, count = JUNC.ID, count.o = JUNC.ID.O) %>% 
      mutate(count = ifelse(is.na(count), 0, count)) %>% 
      mutate(count.o = ifelse(is.na(count.o), 0, count.o)) %>% 
      mutate(depth = count + count.o) %>% 
      mutate(freq = round(count / depth, 4)) %>% 
      mutate(freq = ifelse(is.na(freq), 0, freq)) %>% 
      mutate(judge = ifelse(count >= THR.C & depth >= THR.D & freq >= THR.F, 1, 0)) %>% 
      mutate(junc.id = JUNC.ID) %>% 
      dplyr::select(junc.id, case, count, depth, freq, judge)
  }
  
  count.i <- res.i %>% 
    dplyr::select(junc.id, case, count) %>% 
    spread(case, count)
  depth.i <- res.i %>% 
    dplyr::select(junc.id, case, depth) %>% 
    spread(case, depth)
  freq.i <- res.i %>% 
    dplyr::select(junc.id, case, freq) %>% 
    spread(case, freq)
  judge.i <- res.i %>% 
    dplyr::select(junc.id, case, judge) %>% 
    spread(case, judge)
  
  if(i == 1){
    count <- count.i 
    depth <- depth.i 
    freq <- freq.i  
    judge <- judge.i
  }else{
    count <- bind_rows(count, count.i) 
    depth <- bind_rows(depth, depth.i)
    freq <- bind_rows(freq, freq.i) 
    judge <- bind_rows(judge, judge.i)
  }
}


# check -------------------------------------------------------------------

count %>% dim() %>% print()
depth %>% dim() %>% print() 
freq %>% dim() %>% print() 
judge %>% dim() %>% print() 
# [1] 9041  430

count %>% anyNA() %>% print()
depth %>% anyNA() %>% print() 
freq %>% anyNA() %>% print() 
judge %>% anyNA() %>% print() 
# [1] FALSE


# out ---------------------------------------------------------------------

setwd(dir.2)
# count
out.f <- filename # "tcga_429_n9041_count.tsv"
write_tsv(count, out.f)
# depth
out.f <- filename # "tcga_429_n9041_depth.tsv"
write_tsv(depth, out.f)
# freq
out.f <- filename # "tcga_429_n9041_freq.tsv"
write_tsv(freq, out.f)
# judge
out.f <- filename # "tcga_429_n9041_judge.tsv"
write_tsv(judge, out.f)


# prep 3 subgroups --------------------------------------------------------------------

cases <- list(
  tcga.annot$case[tcga.annot$cls == "IDHwt"],
  tcga.annot$case[tcga.annot$cls == "IDH-A"],
  tcga.annot$case[tcga.annot$cls == "IDH-O"]
  )
names(cases) <- c("IDHwt_n166", "IDH-A_n140", "IDH-O_n123")

  
# calculate psr -----------------------------------------------------------

for(i in 1:3){
  print(i) ;
  CASES <- cases[[i]] ;
  
  judge.i <- judge %>% 
    dplyr::select(junc.id, colnames(judge)[is.element(colnames(judge), CASES)])
  
  psr.i <- judge.i %>% 
    mutate(pos = apply(judge.i[, -1], 1, sum)) %>% 
    mutate(psr = round(pos / length(CASES), 4)) %>% 
    dplyr::select(junc.id, psr)

  colnames(psr.i)[2] <- names(cases)[i]

  if(i == 1){
    res.psr <- psr.i
  }else{
    res.psr <- res.psr %>% 
      bind_cols(psr.i %>% dplyr::select(- junc.id) )
  }
} 


# edit -------------------------------------------------------------------

res.psr <- res.psr %>% 
  mutate(max = apply(res.psr[, -1], 1, max))


# check -------------------------------------------------------------------

res.psr %>% dim() %>% print()
# [1] 9041    5


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # "tcga_429_n9041_psr.tsv"
write_tsv(res.psr, out.f)


# combine with psr-gtex ---------------------------------------------------------------------

res.psr.merged <- res.psr %>%
  left_join(gtex.psr, by = "junc.id")

# check
res.psr.merged %>% dim()
# [1] 9041    8

res.psr.merged %>% head()
# # A tibble: 6 × 8
#   junc.id             IDHwt_n166 `IDH-A_n140` `IDH-O_n123`    max GTEx_all_n9166
#   <chr>                    <dbl>        <dbl>        <dbl>  <dbl>          <dbl>
# 1 chr1:+:44086903-44…     0.120        0.0643       0.0163 0.120          0.127 
# 2 chr1:-:94997968-94…     0            0.0071       0      0.0071         0.0109
# 3 chr1:+:101702655-1…     0.217        0.157        0.0163 0.217          0.390 
# 4 chr1:+:101707328-1…     0.0904       0.114        0      0.114          0.0919
# 5 chr1:+:109816292-1…     0.892        0.979        0.935  0.979          0.359 
# 6 chr1:+:109816292-1…     0.151        0.564        0.138  0.564          0.189 
# # … with 2 more variables: GTEx_brain_n1436 <dbl>, GTEx_non_brain_n7730 <dbl>


# filter ---------------------------------------------------------------------

res.psr.merged.2 <- res.psr.merged %>%
  mutate(filter.01 = ifelse(GTEx_all_n9166 < 0.01, "tumor-specific", "not")) %>% 
  mutate(filter.02 = ifelse(GTEx_all_n9166 < 0.01 & max >= 0.1, "shared", "not"))


# check ---------------------------------------------------------------------

res.psr.merged.2 %>% pull(filter.01) %>% table() %>% print()
          #  not tumor-specific 
          # 3094           5947 

res.psr.merged.2 %>% pull(filter.02) %>% table() %>% print()
  #  not shared 
  # 8975     66 


# filter ---------------------------------------------------------------------

res.psr.merged.3 <- res.psr.merged.2 %>%
  dplyr::filter(filter.02 == "shared") %>% 
  dplyr::select(- filter.01, -filter.02)


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # "tcga_429_n_gtex_9166_n9041_psr_merged_w_label.tsv"
write_tsv(res.psr.merged.2, out.f)

out.f <- filename # "tcga_429_n_gtex_9166_n9041_psr_merged_filtered.tsv"
write_tsv(res.psr.merged.3, out.f)


# si ----------------------------------------------------------------------

Sys.time() %>% print()
# 

sessionInfo() %>% print()
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

