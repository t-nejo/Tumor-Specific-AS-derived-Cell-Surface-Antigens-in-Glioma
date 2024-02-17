# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE)) ;


# load packages -----------------------------------------------------------

# https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
library(tidyverse) ;
library(tidylog)
library(ensembldb)
library(EnsDb.Hsapiens.v75)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
library(AnnotationHub)


# dir --------------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.2 <- "~/path/to/output_folder"


# load data -----------------------------------------------------

setwd(dir.1)
in.f <- filename # "20231010_psr_tcga_429_n_gtex_9166_n65_in_surfaceome_exp_gnm_coord_xc_or_not.tsv" # a table containing the junction ids of the candidate events.
res.psr.orig <- read_tsv(in.f)


# prep -------------------------------------------------------------

edb <- EnsDb.Hsapiens.v75
edb

organism(edb) # check
supportedFilters(edb)

bsg <- BSgenome.Hsapiens.UCSC.hg19
unique(genome(bsg)) # check

ah <- AnnotationHub() # retrieve genomic seq
query(ah, c("Homo Sapiens", "TwoBit", "hg19"))
dna = ah[["AH13964"]]
head(seqlevels(dna))


# check -------------------------------------------------------------------

res.psr.orig %>% dim()
# [1] 65 13


# edit (focus on the "extra-cellular" events) -------------------------------------------------------------------

res.psr <- res.psr.orig %>% 
  dplyr::filter(xc.or.not == "xc") %>% 
  dplyr::select(junc.id, ensg.symbol, ensg, tx_id, tm.seq, n.tm.domains = TM.domains) %>% 
  mutate(tx_id = ifelse(ensg.symbol == "PTPRZ1", "ENST00000393386;ENST00000449182", tx_id))
# PTPRZ1: 2 canonical transcripts


# check -------------------------------------------------------------------

res.psr %>% dim()
# [1] 42 4

res.psr %>% head()
# # A tibble: 6 × 6
# junc.id                     ensg.symbol ensg        tx_id tm.seq n.tm.domains
# <chr>                       <chr>       <chr>       <chr> <chr>         <dbl>
# 1 chr3:-:66433823-66434414    LRIG1       ENSG000001… ENST… VGIFT…            1
# 2 chr7:+:31127310-31132302    ADCYAP1R1   ENSG000000… ENST… ALYTV…            7
# 3 chr11:+:113076388-113076776 NCAM1       ENSG000001… ENST… LSTGA…            1
# 4 chr19:+:39993685-39994710   DLL3        ENSG000000… ENST… LLPPA…            1
# 5 chr1:+:156615937-156616594  BCAN        ENSG000001… ENST… NA                0
# 6 chr1:+:156615937-156616766  BCAN        ENSG000001… ENST… NA                0


# analysis ----------------------------------------------------------------

res.2 <- NULL
for(i in 1:nrow(res.psr)){
  print( paste0("i = ", i) );
  
  res.i <- res.psr %>% 
    dplyr::slice(i)
  
  JUNC.ID <- res.i %>% pull(junc.id);
  CHR <- gsub("chr", "", gsub(":.*", "", JUNC.ID))
  STRAND <- gsub(":.*", "", gsub("chr..:", "", gsub("chr.:", "", JUNC.ID)))
  START <- gsub("-.*", "", gsub(".*:", "", JUNC.ID)) %>% as.integer()
  END <- gsub(".*-", "", gsub(".*:", "", JUNC.ID)) %>% as.integer()
  
  ENSG <- res.i %>% pull(ensg);
  ENSTs <- res.i %>% pull(tx_id) %>% str_split(";") %>% unlist()
  
  N.TM <- res.i %>% pull(n.tm.domains)
  TM.SEQ <- res.i %>% pull(tm.seq)
  
  # based on the ENSG, list up all the candidate ENSTs ---------------------------
  
  tbl.tx <- transcripts(
    edb, 
    columns = c("tx_id", "protein_sequence"), 
    filter = list(GeneIdFilter(ENSG), TxBiotypeFilter("protein_coding"))
  ) %>% 
    as_tibble() %>% 
    mutate(aa.length = nchar(protein_sequence)) %>% 
    mutate(seqnames = paste0("chr", seqnames)) %>% 
    dplyr::filter(is.element(tx_id, ENSTs))
  
  if(tbl.tx %>% nrow() != length(ENSTs)){
    print("check"); 
    stop
  }
  
  
  # investigate each candidate tx ----------------------------------------
  
  for(j in 1:nrow(tbl.tx)){
    print( paste0("j = ", j) );
    
    ENST.j <- tbl.tx %>% 
      dplyr::slice(j) %>% 
      pull(tx_id)
    
    
    # get the list of exons within the tx ---------------------------------
    
    exons <- exons(
      edb, 
      columns = c("tx_id", "exon_idx"), 
      filter = list(TxIdFilter(ENST.j))
    ) %>% # filter
      as_tibble() %>% 
      mutate(seqnames = paste0("chr", seqnames)) # edit chr names
    
    exons <- exons %>% 
      mutate(seq = as.vector(getSeq(dna, GRanges(exons) ) ) ) %>% # get the dna seq
      arrange(exon_idx) %>% # re-order
      dplyr::select(- exon_id) # get rid of unnecessary col
    
    
    # identify the starting and ending positions of the junctions within the exons --------
    
    START.EXON.ID <- exons %>% 
      dplyr::filter(start <= START & START <= end) %>% 
      pull(exon_idx)
    
    if(START.EXON.ID %>% length() == 0){
      START.EXON.ID <- exons %>% 
        dplyr::select(start, end, exon_idx) %>% 
        dplyr::filter(end < START) %>% 
        arrange(desc(start)) %>% 
        dplyr::slice(1) %>% 
        pull(exon_idx)
    }
    
    END.EXON.ID <- exons %>% 
      dplyr::select(start, end, exon_idx) %>% 
      dplyr::filter(start <= (END + 1) & (END + 1) <= end) %>% 
      pull(exon_idx)
    
    if(END.EXON.ID %>% length() == 0){
      END.EXON.ID <- exons %>% 
        dplyr::select(start, end, exon_idx) %>% 
        dplyr::filter(END + 1 < start) %>% 
        arrange(start) %>% 
        dplyr::slice(1) %>% 
        pull(exon_idx)
    }
    
    
    # identify the non-affected exons ----------------------------------------------
    
    INTACT.EXONS <- exons %>% 
      dplyr::filter(exon_idx < min(START.EXON.ID, END.EXON.ID) | max(START.EXON.ID, END.EXON.ID) < exon_idx)
    
    
    # edit the affected exons ----------------------------------------------
    
    NEW.START.EXON <- GRanges(
      seqnames = paste0("chr", CHR),
      strand = STRAND, 
      ranges = IRanges(
        exons %>% dplyr::filter(exon_idx == START.EXON.ID) %>% pull(start), 
        START
      )
    ) %>% 
      as_tibble() %>% 
      mutate(tx_id = NA) %>% 
      mutate(exon_idx = START.EXON.ID + 0.5) %>% 
      mutate(seqnames = as.character(seqnames)) %>% 
      mutate(seq = as.vector(
        getSeq(
          dna, 
          GRanges(
            seqnames = paste0("chr", CHR),
            strand = STRAND, 
            ranges = IRanges(
              exons %>% dplyr::filter(exon_idx == START.EXON.ID) %>% pull(start), 
              START
            )
          )
        )
      )
      )
    
    NEW.END.EXON <- GRanges(
      seqnames = paste0("chr", CHR),
      strand = STRAND, 
      ranges = IRanges(
        END + 1, 
        exons %>% dplyr::filter(exon_idx == END.EXON.ID) %>% pull(end) 
      )
    ) %>% 
      as_tibble() %>% 
      mutate(tx_id = NA) %>% 
      mutate(exon_idx = END.EXON.ID + 0.5) %>% 
      mutate(seqnames = as.character(seqnames)) %>% 
      mutate(seq = as.vector(
        getSeq(
          dna, 
          GRanges(
            seqnames = paste0("chr", CHR),
            strand = STRAND, 
            ranges = IRanges(
              END + 1, 
              exons %>% dplyr::filter(exon_idx == END.EXON.ID) %>% pull(end)
            )
          )
        )
      )
      )
    
    
    # build the structure of the exons of the altered tx ---------------------------------
    
    exons.new <- exons %>% 
      dplyr::filter(exon_idx < min(START.EXON.ID, END.EXON.ID) | max(START.EXON.ID, END.EXON.ID) < exon_idx) %>% 
      bind_rows(NEW.START.EXON) %>% 
      bind_rows(NEW.END.EXON) %>% 
      arrange(exon_idx)
    
    
    # get tx.seq ------------------------------------------------------------------
    
    tx.seq.wt <- paste(exons$seq, collapse = "")
    tx.seq.alt <- paste(exons.new$seq, collapse = "")
    
    
    # rm. 5'UTR / 3’UTR -----------------------------------------------------------
    
    # 5'UTR
    utr.5 <- fiveUTRsByTranscript(edb, 
                                  columns = c("tx_id", "exon_idx"), 
                                  filter = list(TxIdFilter(ENST.j))
    )
    if(is.null(utr.5)){
      ln.utr.5 <- 0
    }else{
      ln.utr.5 <- utr.5 %>% as_tibble() %>% 
        pull(width) %>% 
        sum()
    }
    
    # 3'UTR
    utr.3 <- threeUTRsByTranscript(edb, 
                                   columns = c("tx_id", "exon_idx"), 
                                   filter = list(TxIdFilter(ENST.j))
    )
    if(is.null(utr.3)){
      ln.utr.3 <- 0
    }else{
      ln.utr.3 <- utr.3 %>% as_tibble() %>% 
        pull(width) %>% 
        sum()
    }
    
    
    # get CDS seq ------------------------------------------------------------------
    
    cds.seq.wt <- substr(tx.seq.wt, (ln.utr.5 + 1), (nchar(tx.seq.wt) - ln.utr.3))
    cds.seq.alt <- substr(tx.seq.alt, (ln.utr.5 + 1), (nchar(tx.seq.alt) - ln.utr.3))
    
    
    # translation ------------------------------------------------------------------
    
    aa.seq.wt <- as.character( translate( DNAString(cds.seq.wt) ) )
    aa.seq.alt <- as.character( translate( DNAString(cds.seq.alt) ) )
    
    aa.seq.wt <- ifelse(str_sub(aa.seq.wt, -1, -1) == "*", str_sub(aa.seq.wt, 1, -2), aa.seq.wt)
    aa.seq.alt <- ifelse(str_sub(aa.seq.alt, -1, -1) == "*", str_sub(aa.seq.alt, 1, -2), aa.seq.alt)    
    
    
    # check -------------------------------------------------------------------
    
    aa.seq.wt == proteins(edb, 
                          columns = c("protein_sequence"), 
                          filter = TxIdFilter(ENST.j)) %>% 
      as_tibble() %>% 
      pull(protein_sequence) 
    
    
    # check frameshift and PTCs --------------------------------------------------------------------
    
    FS <- ifelse(abs(nchar(cds.seq.wt) - nchar(cds.seq.alt) ) %% 3 == 0, "in-frame", "fs")
    SC <- ifelse(grepl("\\*", aa.seq.alt), "sc", "no.sc")
    
    
    # check tm status ---------------------------------------------------------
    
    if(N.TM == 0){
      TM.STATUS <- ifelse(
        substr(aa.seq.wt, nchar(aa.seq.wt) - 4, nchar(aa.seq.wt)) == substr(aa.seq.alt, nchar(aa.seq.alt) - 4, nchar(aa.seq.alt)), 
        "last5.intact", 
        "last5.lost"
      )
    }else if(N.TM == 1 & !is.na(TM.SEQ)){
      TM.STATUS <-ifelse(grepl(TM.SEQ, aa.seq.alt), "tm.intact", "tm.lost")    
    }else{
      TM.STATUS <- "check"
    }
    
    
    # prep out ----------------------------------------------------------------
    
    res.i2 <- res.i %>% 
      mutate(enst.model = ENST.j) %>% 
      mutate(change = ifelse(nchar(aa.seq.wt) > nchar(aa.seq.alt), "loss", "gain")) %>% 
      mutate(aa.seq.wt = aa.seq.wt) %>% 
      mutate(aa.seq.alt = aa.seq.alt) %>% 
      mutate(ln.aa.wt = nchar(aa.seq.wt)) %>% 
      mutate(ln.aa.alt = nchar(aa.seq.alt)) %>% 
      mutate(ln.aa.diff = ln.aa.alt - ln.aa.wt) %>% 
      mutate(fs = FS) %>% 
      mutate(sc = SC) %>% 
      mutate(tm.status = TM.STATUS)
    
    res.2 <- res.2 %>% 
      bind_rows(res.i2)
  }
}


# check -------------------------------------------------------------------

res.2 %>% nrow()
# [1] 50

res.2 %>% 
  dplyr::select(junc.id, ensg.symbol, tx_id, enst.model, fs, sc, tm.status) %>% 
  dplyr::filter(fs == "in-frame" & sc == "no.sc") %>% 
  nrow()
# [1] 20


# manual check ------------------------------------------------------------

res.2 %>% 
  dplyr::filter(fs == "in-frame" & sc == "no.sc") %>% 
  dplyr::filter(n.tm.domains >= 2)
# # A tibble: 2 × 16
# junc.id         ensg.symbol ensg  tx_id tm.seq n.tm.domains enst.model change
# <chr>           <chr>       <chr> <chr> <chr>         <dbl> <chr>      <chr> 
# 1 chr7:+:3112731… ADCYAP1R1   ENSG… ENST… ALYTV…            7 ENST00000… loss  
# 2 chr5:+:8997202… ADGRV1      ENSG… ENST… AFFTS…            7 ENST00000… loss  
# # … with 8 more variables: aa.seq.wt <chr>, aa.seq.alt <chr>, ln.aa.wt <int>,
# #   ln.aa.alt <int>, ln.aa.diff <int>, fs <chr>, sc <chr>, tm.status <chr>


res.2 %>% 
  dplyr::filter(junc.id == "chr7:+:31127310-31132302") %>% 
  dplyr::select(junc.id, n.tm.domains, tm.seq)
# # A tibble: 1 × 3
# junc.id                  n.tm.domains tm.seq                                 
# <chr>                           <dbl> <chr>                                  
# 1 chr7:+:31127310-31132302            7 ALYTVGYSTSLVTLTTAMVILC;IHMNLFVSFMLRAIS…

# to be excluded

res.2 %>% 
  dplyr::filter(junc.id == "chr5:+:89972026-89975368") %>% 
  dplyr::select(junc.id, n.tm.domains, tm.seq)
# # A tibble: 1 × 3
# junc.id                  n.tm.domains tm.seq                                 
# <chr>                           <dbl> <chr>                                  
# 1 chr5:+:89972026-89975368            7 AFFTSGFICISGLCLAVLSHIF;LLTHMMAASLGTQIL…

# to be excluded


# filtering based on the manual review --------------------------------------------------

res.3 <- res.2 %>% 
  dplyr::filter(fs == "in-frame" & sc == "no.sc") %>% 
  dplyr::filter(junc.id != "chr5:+:89972026-89975368" & junc.id != "chr7:+:31127310-31132302")



# check -------------------------------------------------------------------

res.3 %>% dim()
# [1] 18 24

res.3 %>% distinct(junc.id, .keep_all = T) %>% dim()
# [1] 14 24

res.2 %>% 
  anti_join(res.3, by = "junc.id") %>% 
  distinct(junc.id, ensg.symbol, fs, sc) 

# 26 events - sc
# 2 events - multi-pass

26/42
# [1] 0.6190476


# out ---------------------------------------------------------------------

setwd(dir.2)
out.f <- filename # "res_aaseq_pred_tcga_429_n_gtex_9166_n42.tsv"
write_tsv(res.2, out.f)

setwd(dir.2)
out.f <- filename # "res_aaseq_pred_tcga_429_n_gtex_9166_n14_final_candidates.tsv"
write_tsv(res.3, out.f)


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
# [1] AnnotationHub_3.0.2               BiocFileCache_2.0.0              
# [3] dbplyr_2.1.1                      Rsamtools_2.8.0                  
# [5] BSgenome.Hsapiens.UCSC.hg19_1.4.3 BSgenome_1.60.0                  
# [7] rtracklayer_1.52.1                Biostrings_2.60.2                
# [9] XVector_0.32.0                    EnsDb.Hsapiens.v75_2.99.0        
# [11] ensembldb_2.16.4                  AnnotationFilter_1.16.0          
# [13] GenomicFeatures_1.44.2            AnnotationDbi_1.54.1             
# [15] Biobase_2.52.0                    GenomicRanges_1.44.0             
# [17] GenomeInfoDb_1.28.1               IRanges_2.26.0                   
# [19] S4Vectors_0.30.2                  BiocGenerics_0.38.0              
# [21] tidylog_1.0.2                     forcats_0.5.1                    
# [23] stringr_1.4.0                     dplyr_1.0.9                      
# [25] purrr_0.3.4                       readr_2.1.2                      
# [27] tidyr_1.2.0                       tibble_3.1.7                     
# [29] ggplot2_3.3.6                     tidyverse_1.3.1                  
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
# [85] readxl_1.3.1                  blob_1.2.2                   
# [87] reprex_2.0.1                  digest_0.6.29                
# [89] xtable_1.8-4                  httpuv_1.6.5                 
# [91] munsell_0.5.0                


# end ---------------------------------------------------------------------
