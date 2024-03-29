# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco

# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# load packages -----------------------------------------------------------

library(tidyverse) ;


# dir ---------------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.2 <- "~/path/to/output_folder"


# load data -----------------------------------------------------

setwd(dir.1)
in.f <- filename # "07032020_tcga_glioma_annotation_abs_purity_060_n429.tsv" # 429 TCGA cases to analyze
tcga.annot <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "20201128_sjdbList.fromGTF.out.tab_to_analyze.tsv" ; # sj.annot.ref, made from "sjdbList.fromGTF.out.tab" generated by STAR-alinger index
sj.ref <- read_tsv(in.f, na = c("", "NA"), col_names = T, col_types = cols(chr = col_character())) ;

setwd(dir.1)
in.f <- filename # "pnas.1808790115.sd01_11.7.surfaceome.txt"
surfaceome <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "average_surfaceome_gene_expression_in_TCGA_glioma.tsv"     
surfaceome.tcga.exp.data.mean <- read_tsv(in.f)

setwd(dir.1)
in.f <- filename # "gencode.v29lift37.annotation_protein_coding_and_more_n20723.gtf" # GTF file focusing on protein-coding gene info
pc.gtf <- read_tsv(in.f)




# check -------------------------------------------------------------------

tcga.annot %>% dim() %>% print() ;
# [1] 429   8

tcga.annot %>% head() %>% print() ;
# # A tibble: 6 × 8
# case         study    cls   grade histology  IDH   codel abs.purity
# <chr>        <chr>    <chr> <chr> <chr>      <chr> <chr>      <dbl>
# 1 TCGA-CS-4938 TCGA-LGG IDH-A G2    astrocyto… Muta… non-…       0.79
# 2 TCGA-CS-4941 TCGA-LGG IDHwt G3    astrocyto… WT    non-…       0.61
# 3 TCGA-CS-4942 TCGA-LGG IDH-A G3    astrocyto… Muta… non-…       0.76
# 4 TCGA-CS-4943 TCGA-LGG IDH-A G3    astrocyto… Muta… non-…       0.83
# 5 TCGA-CS-4944 TCGA-LGG IDH-A G2    astrocyto… Muta… non-…       0.74
# 6 TCGA-CS-5390 TCGA-LGG IDH-O G2    oligodend… Muta… codel       0.85


# find all the unique junctions with minimum counts ≥ 10 counts each sj file -------------------------------------------------------

sj.list <- list() # list object to store junction info identified in each case in TCGA
junc.ids.total.1 <- NULL # blank object to store all the distinct junction info identified in all TCGA cases

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
  
  sj.list[[i]] <- sj.i %>% 
    dplyr::select(junc.id, n.uniq.map) ;
  names(sj.list)[i] = CASE.BARCODE
  
  junc.ids.total.1 <- junc.ids.total.1 %>%
    bind_rows(
      sj.i %>% 
        dplyr::filter(n.uniq.map >= 10) %>%
        dplyr::select(chr, int.start, int.end, strand, junc.id)
    ) %>%
    distinct(junc.id, .keep_all = T)
}


# check -------------------------------------------------------------------

sj.list %>% length()
# [1] 429

junc.ids.total.1 %>% nrow()
# [1] 366734


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # e.g., "sj.list.rds"
saveRDS(sj.list, out.f)

setwd(dir.2)
out.f <- filename # e.g., "all_the_distinct_junc_ids_min_10_in_TCGA.rds"
saveRDS(junc.ids.total.1, out.f)


# 2. Non-annotated --------------------------------------------------------

# check -------------------------------------------------------------------

sj.ref %>% nrow()
# [1] 293310


sj.ref %>% head(3)
# # A tibble: 3 x 4
#      X1    X2    X3 X4
#   <dbl> <dbl> <dbl> <chr>
# 1     1 12058 12178 +
# 2     1 12228 12594 +
# 3     1 12228 12612 +

# prep --------------------------------------------------------------------

colnames(sj.ref) <- c("chr", "int.start", "int.end", "strand") ;

sj.ref <- sj.ref %>% 
  mutate(junc.id = paste0("chr", chr, ":", strand, ":", (int.start - 1), "-", int.end)) %>% 
  dplyr::select(junc.id, chr, strand, int.start, int.end) ; 


# identify the "non-annotated" junctions ----------------------------------

junc.ids.total.2 <- junc.ids.total.1 %>% 
  anti_join(sj.ref, by = "junc.id")

# anti_join: added no columns
# > rows only in x   168,286
# > rows only in y  ( 94,860)
# > matched rows    (198,448)
# >                 =========
# > rows total       168,286


# check -------------------------------------------------------------------

junc.ids.total.2 %>% nrow()
# [1] 168286


# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # e.g., "all_non_annotated_junctions_in_TCGA.rds"
saveRDS(junc.ids.total.2, out.f)


# 3. surfaceome genes --------------------------------------------------------

surfaceome %>% nrow()
# [1] 2886

surfaceome.tcga.exp.data.mean %>% nrow()
# [1] 2532

surfaceome.tcga.exp.data.mean %>% head(n = 10)
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


# edit (focus on surfaceome genes with mean TPM ≥ 10) --------------------------------------------------------------------

surfaceome.tcga.exp.data.mean <- surfaceome.tcga.exp.data.mean %>% 
  dplyr::filter(mean.tpm >= 10)

pc.gtf.2 <- pc.gtf %>% 
  mutate(ensg = substr(ensg, 1, 15)) %>% 
  dplyr::select(X1, X4, X5, X7, ensg, symbol)
colnames(pc.gtf.2)[1:4] <- c("chr", "start", "end", "strand")


# check -------------------------------------------------------------------

surfaceome %>% nrow()
# [1] 2886

surfaceome.tcga.exp.data.mean %>% nrow()
# [1] 737


# edit (add Ensembl gene ID)--------------------------------------------------------------------

surfaceome.tcga.exp.data.mean.2 <- surfaceome.tcga.exp.data.mean %>% 
  left_join(
    surfaceome %>% dplyr::select(ensg = `Ensembl gene`, symbol = `UniProt gene`), 
    by = "symbol"
  ) %>% 
  distinct(symbol, .keep_all = T) %>% 
  inner_join(
    pc.gtf.2 %>% dplyr::rename(ensg.symbol = symbol), 
    by = "ensg"
  )
# left_join: added one column (ensg)
# > rows only in x       0
# > rows only in y  (2,068)
# > matched rows       818    (includes duplicates)
# >                 =======
# > rows total         818
# distinct: removed 81 rows (10%), 737 rows remaining
# inner_join: added 6 columns (symbol.x, chr, start, end, strand, …)
# > rows only in x  (     2)
# > rows only in y  (19,987)
# > matched rows        736    (includes duplicates)
# >                 ========
# > rows total          736


# check -------------------------------------------------------------------

surfaceome.tcga.exp.data.mean.2 %>% nrow()
# [1] 736

surfaceome.tcga.exp.data.mean.2 %>% 
  pull(ensg) %>% 
  table() %>% 
  as_tibble() %>% 
  arrange(desc(n)) %>% 
  head()
# # A tibble: 6 × 2
# .                   n
# <chr>           <int>
# 1 ENSG00000198223     2
# 2 ENSG00000000003     1
# 3 ENSG00000001461     1
# 4 ENSG00000003056     1
# 5 ENSG00000003436     1
# 6 ENSG00000004399     1

surfaceome.tcga.exp.data.mean.2 %>% 
  dplyr::filter(ensg == "ENSG00000198223")
# # A tibble: 2 × 8
# symbol mean.tpm ensg            chr     start     end strand ensg.symbol
# <chr>     <dbl> <chr>           <chr>   <dbl>   <dbl> <chr>  <chr>      
# 1 CSF2RA     14.3 ENSG00000198223 chrX  1387693 1429274 +      CSF2RA     
# 2 CSF2RA     14.3 ENSG00000198223 chrY  1337693 1379274 +      CSF2RA     


# edit (remove duplication) -----------------------------------------------

surfaceome.tcga.exp.data.mean.2 <- surfaceome.tcga.exp.data.mean.2 %>% 
  distinct(ensg, .keep_all = T)


# 3. junc.ids in surfaceome genes --------------------------------------------------------

junc.ids.total.3 <- NULL
for(i in 1:nrow(junc.ids.total.2)){
  print( round( i / nrow(junc.ids.total.2), 5))

  junc.ids.i <- junc.ids.total.2 %>% 
    dplyr::slice(i)
  
  CHR <- junc.ids.i %>% pull(chr) ; 
  INT.START <- junc.ids.i %>% pull(int.start) ; 
  INT.END <- junc.ids.i %>% pull(int.end) ;
  STRAND <- junc.ids.i %>% pull(strand)
  
  match.i <- surfaceome.exp.2 %>% 
    dplyr::filter(chr == paste0("chr", CHR)) %>% 
    dplyr::filter(strand == STRAND) %>% 
    dplyr::filter(start < INT.START & INT.END < end)
  
  if(nrow(match.i) == 0){
    next
  }else{
    SYMBOL <- match.i %>% pull(ensg.symbol) %>% paste(collapse = ";")
    ENSG <- match.i %>% pull(ensg) %>% paste(collapse = ";")
    
    junc.ids.i <- junc.ids.i %>% 
      mutate(symbol = SYMBOL) %>% 
      mutate(ensg = ENSG) %>% 
      dplyr::select(symbol, ensg, junc.id, chr, int.start, int.end, strand)
    
    junc.ids.total.3 <- junc.ids.total.3 %>% 
      bind_rows(junc.ids.i)
  }
}


# check -------------------------------------------------------------------

junc.ids.total.3 %>% nrow()
# [1] 9041

junc.ids.total.3 %>% 
  print()
# # A tibble: 9,041 × 7
# symbol ensg            junc.id                chr   int.start int.end strand
# <chr>  <chr>           <chr>                  <chr>     <dbl>   <dbl> <chr> 
#   1 PTPRF  ENSG00000142949 chr1:+:44086903-44087… 1      44086904  4.41e7 +     
#   2 F3     ENSG00000117525 chr1:-:94997968-94998… 1      94997969  9.50e7 -     
#   3 S1PR1  ENSG00000170989 chr1:+:101702655-1017… 1     101702656  1.02e8 +     
#   4 S1PR1  ENSG00000170989 chr1:+:101707328-1017… 1     101707329  1.02e8 +     
#   5 CELSR2 ENSG00000143126 chr1:+:109816292-1098… 1     109816293  1.10e8 +     
#   6 CELSR2 ENSG00000143126 chr1:+:109816292-1098… 1     109816293  1.10e8 +     
#   7 BCAN   ENSG00000132692 chr1:+:156617902-1566… 1     156617903  1.57e8 +     
#   8 BCAN   ENSG00000132692 chr1:+:156617902-1566… 1     156617903  1.57e8 +     
#   9 BCAN   ENSG00000132692 chr1:+:156618653-1566… 1     156618654  1.57e8 +     
#   10 BCAN   ENSG00000132692 chr1:+:156618653-1566… 1     156618654  1.57e8 +     
#   # … with 9,031 more rows
  

# save --------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # e.g., "all_surfaceome_junctions_in_TCGA.rds"
saveRDS(junc.ids.total.3, out.f)


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
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] tidylog_1.0.2   readxl_1.3.1    forcats_0.5.1   stringr_1.4.0  
# [5] dplyr_1.0.9     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
# [9] tibble_3.1.7    ggplot2_3.3.6   tidyverse_1.3.1
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.8.3     cellranger_1.1.0 pillar_1.7.0     compiler_4.1.3  
# [5] dbplyr_2.1.1     tools_4.1.3      bit_4.0.4        timechange_0.2.0
# [9] jsonlite_1.8.0   lubridate_1.9.2  lifecycle_1.0.3  gtable_0.3.0    
# [13] pkgconfig_2.0.3  rlang_1.1.0      reprex_2.0.1     DBI_1.1.2       
# [17] cli_3.6.1        rstudioapi_0.13  haven_2.4.3      xml2_1.3.2      
# [21] withr_2.5.0      httr_1.4.3       fs_1.5.2         generics_0.1.2  
# [25] vctrs_0.6.1      hms_1.1.0        bit64_4.0.5      grid_4.1.3      
# [29] tidyselect_1.1.2 glue_1.6.2       R6_2.5.1         fansi_1.0.3     
# [33] vroom_1.5.7      tzdb_0.1.2       modelr_0.1.8     magrittr_2.0.3  
# [37] clisymbols_1.2.0 backports_1.2.1  scales_1.2.0     ellipsis_0.3.2  
# [41] rvest_1.0.1      assertthat_0.2.1 colorspace_2.0-3 utf8_1.2.2      
# [45] stringi_1.7.6    munsell_0.5.0    broom_0.7.9      crayon_1.5.1   


# end ---------------------------------------------------------------------
