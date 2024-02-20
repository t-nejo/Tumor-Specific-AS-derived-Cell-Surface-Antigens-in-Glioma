# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco

# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# load packages -----------------------------------------------------------

library(tidyverse)
library(tidylog)
library(Spectra)
library(doParallel) 
library(PSMatch) 


# dir --------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.2 <- "~/path/to/output_folder"


# import experiment ----------------------------------------------------------------

setwd(dir.1)
in.f <- filename # "cptac_msgfplus_masic_combined_cell_surface_all.tsv"
res.0.all <- read_tsv(in.f)

in.f <- filename # "cptac_msgfplus_masic_combined_qvalue0.05_passed_curated.tsv"
res.0.curated <- read_tsv(in.f)


# obtain the necessary file info -----------------------------------------------

res.0.all %>% dim()
# [1] 43 35

res.0.all %>% colnames()
#  [1] "symbol_juncid"           "#SpecFile"               "SpecID"                 
#  [4] "ScanNum"                 "FragMethod"              "Precursor"              
#  [7] "IsotopeError"            "PrecursorError(ppm)"     "Charge"                 
# [10] "Peptide"                 "Formula"                 "Protein"                
# [13] "DeNovoScore"             "MSGFScore"               "SpecEValue"             
# [16] "EValue"                  "QValue"                  "PepQValue"              
# [19] "alias"                   "ProteinPeptideFormula"   "ParentIonMZ"            
# [22] "BasePeakIntensity"       "BasePeakMZ"              "ReporterIonIntensityMax"
# [25] "TMT11-126"               "TMT11-127N"              "TMT11-127C"             
# [28] "TMT11-128N"              "TMT11-128C"              "TMT11-129N"             
# [31] "TMT11-129C"              "TMT11-130N"              "TMT11-130C"             
# [34] "TMT11-131N"              "TMT11-131C"     

res.0.curated %>% 
  dplyr::select(symbol, peptide, detected.peptide.signal) %>% print()
# # A tibble: 3 × 3
# symbol   peptide                             detected.peptide.signal         
# <chr>    <chr>                               <chr>                           
# 1 EGFRvIII GNYVVTDHGSCVR                       +229.163GNYVVTDHGSC+57.021VR    
# 2 EGFRvIII KGNYVVTDHGSCVR                      +229.163K+229.163GNYVVTDHGSC+57…
# 3 BCAN     GAIYSIPIMEDGGGGSSTPEDPAEAPRTLLEEEGK +229.163GAIYSIPIMEDGGGGSSTPEDPA…


# edit --------------------------------------------------------------

res.0.curated.b <- res.0.curated %>% 
  dplyr::select(peptide, detected.peptide.signal) 

res.0.all.b <- res.0.all %>% 
  inner_join(res.0.curated.b, by = c("Peptide" = "detected.peptide.signal")) %>% 
  dplyr::select(symbol_juncid, `#SpecFile`, SpecID, ScanNum, EValue, QValue, peptide)


# check -------------------------------------------------------------------

res.0.all.b %>% 
  dplyr::filter(grepl("BCAN", symbol_juncid)) %>% print()
# # A tibble: 3 × 6
# symbol_juncid                   `#SpecFile`      SpecID ScanNum EValue QValue
# <chr>                           <chr>            <chr>    <dbl>  <dbl>  <dbl>
# 1 BCAN_chr1:+:156621481-156622081 09CPTAC_GBM_W_P… contr…   39254  0.256 0.0188
# 2 BCAN_chr1:+:156621481-156622081 10CPTAC_GBM_W_P… contr…   39957 28.1   0.345 
# 3 BCAN_chr1:+:156621481-156622081 10CPTAC_GBM_W_P… contr…   39986  3.33  0.136 

res.0.all.b %>% 
  dplyr::filter(!grepl("BCAN", symbol_juncid)) %>% 
  distinct(peptide, `#SpecFile`, .keep_all = T) %>% 
  dplyr::select(`#SpecFile`, peptide, EValue, QValue) %>% 
  print(n = Inf)
# distinct: removed one row (5%), 21 rows remaining
# # A tibble: 21 × 4
# `#SpecFile`                               peptide          EValue QValue
# <chr>                                     <chr>             <dbl>  <dbl>
# 1 01CPTAC_GBM_W_PNNL_20190123_B1S1_f04.mzML KGNYVVTDHGSCVR 5.26e-11      0
# 2 01CPTAC_GBM_W_PNNL_20190123_B1S1_f05.mzML GNYVVTDHGSCVR  1.56e-11      0
# 3 02CPTAC_GBM_W_PNNL_20190123_B1S2_f04.mzML KGNYVVTDHGSCVR 8.15e-10      0
# 4 02CPTAC_GBM_W_PNNL_20190123_B1S2_f05.mzML GNYVVTDHGSCVR  1.56e-11      0
# 5 03CPTAC_GBM_W_PNNL_20190123_B1S3_f05.mzML GNYVVTDHGSCVR  2.38e- 7      0
# 6 04CPTAC_GBM_W_PNNL_20190123_B1S4_f03.mzML KGNYVVTDHGSCVR 1.37e- 8      0
# 7 04CPTAC_GBM_W_PNNL_20190123_B1S4_f04.mzML GNYVVTDHGSCVR  1.56e-11      0
# 8 04CPTAC_GBM_W_PNNL_20190123_B1S4_f05.mzML GNYVVTDHGSCVR  1.63e- 4      0
# 9 05CPTAC_GBM_W_PNNL_20190306_B2S1_f03.mzML KGNYVVTDHGSCVR 4.40e-10      0
# 10 05CPTAC_GBM_W_PNNL_20190306_B2S1_f04.mzML KGNYVVTDHGSCVR 1.59e- 9      0
# 11 05CPTAC_GBM_W_PNNL_20190306_B2S1_f04.mzML GNYVVTDHGSCVR  1.22e-10      0
# 12 05CPTAC_GBM_W_PNNL_20190306_B2S1_f05.mzML GNYVVTDHGSCVR  3.64e-11      0
# 13 07CPTAC_GBM_W_PNNL_20190306_B2S3_f04.mzML KGNYVVTDHGSCVR 2.88e-10      0
# 14 07CPTAC_GBM_W_PNNL_20190306_B2S3_f05.mzML GNYVVTDHGSCVR  2.01e-10      0
# 15 08CPTAC_GBM_W_PNNL_20190306_B2S4_f04.mzML KGNYVVTDHGSCVR 1.38e- 9      0
# 16 09CPTAC_GBM_W_PNNL_20190501_B3S1_f05.mzML GNYVVTDHGSCVR  1.87e-10      0
# 17 10CPTAC_GBM_W_PNNL_20190501_B3S2_f03.mzML KGNYVVTDHGSCVR 9.03e-11      0
# 18 10CPTAC_GBM_W_PNNL_20190501_B3S2_f04.mzML GNYVVTDHGSCVR  2.21e-11      0
# 19 10CPTAC_GBM_W_PNNL_20190501_B3S2_f05.mzML GNYVVTDHGSCVR  2.44e-10      0
# 20 11CPTAC_GBM_W_PNNL_20190501_B3S3_f04.mzML GNYVVTDHGSCVR  7.42e-10      0
# 21 11CPTAC_GBM_W_PNNL_20190501_B3S3_f05.mzML GNYVVTDHGSCVR  8.26e- 5      0


# filter ------------------------------------------------------------------

res.0.all.e <- res.0.all.b %>% 
  dplyr::filter(!grepl("BCAN", symbol_juncid)) %>% 
  arrange(EValue, QValue) %>% 
  distinct(peptide, .keep_all = T)

res.0.all.e[1, ] %>% 
  pull(2)
# [1] "01CPTAC_GBM_W_PNNL_20190123_B1S1_f05.mzML"

res.0.all.e[1, ] %>% 
  pull(3)
# [1] "controllerType=0 controllerNumber=1 scan=10989"

res.0.all.e[2, ] %>% 
  pull(2)
# [1] "01CPTAC_GBM_W_PNNL_20190123_B1S1_f04.mzML"

res.0.all.e[2, ] %>% 
  pull(3)
# [1] "controllerType=0 controllerNumber=1 scan=11119"

res.0.all.e[1, ] %>% pull(7)
# [1] "GNYVVTDHGSCVR"
res.0.all.e[2, ] %>% pull(7)
# [1] "KGNYVVTDHGSCVR"
res.0.all.b[1, ] %>% pull(7)
# [1] "GAIYSIPIMEDGGGGSSTPEDPAEAPRTLLEEEGK"

SEQUENCES <- c("GNYVVTDHGSCVR", "KGNYVVTDHGSCVR", "GAIYSIPIMEDGGGGSSTPEDPAEAPRTLLEEEGK")


# load the mass-spec files ----------------------------------------------------------


files <- c(
  "01CPTAC_GBM_W_PNNL_20190123_B1S1_f04.mzid", 
  "01CPTAC_GBM_W_PNNL_20190123_B1S1_f05.mzid", 
  "09CPTAC_GBM_W_PNNL_20190501_B3S1_f22.mzid"
)

list.ids <- list()
for(i in 1:length(files)){
  print(i)
  
  setwd(dir.2)
  in.f <- files[i]
  list.ids[[i]] <- readPSMs(in.f)
}


# 3.1.2 Spectra from mzML files -------------------------------------------

# setwd(dir.3)
files <- c(
  "01CPTAC_GBM_Proteome_PNNL_20190123/01CPTAC_GBM_W_PNNL_20190123_B1S1_f04.mzML.gz", 
  "01CPTAC_GBM_Proteome_PNNL_20190123/01CPTAC_GBM_W_PNNL_20190123_B1S1_f05.mzML.gz", 
  "09CPTAC_GBM_Proteome_PNNL_20190501/09CPTAC_GBM_W_PNNL_20190501_B3S1_f22.mzML.gz"
)

list.sp <- list()
for(i in 1:length(files)){
  print(i)
  
  setwd(dir.1) # directory where the mzML.gz files are stored. 
  in.f <- files[i]
  list.sp[[i]] <- Spectra(in.f)
}


# check sp ----------------------------------------------------------------

for(i in 1:3){
  print(i)
  list.sp[[i]] %>% print()
}

# [1] 1
# MSn data (Spectra) with 52240 spectra in a MsBackendMzR backend:
#   msLevel     rtime scanIndex
# <integer> <numeric> <integer>
# 1             1 0.0272359         1
# 2             1 0.3111044         2
# 3             1 0.5897405         3
# 4             1 0.8727232         4
# 5             1 1.1895906         5
# ...         ...       ...       ...
# 52236         2   6599.81     52236
# 52237         2   6599.93     52237
# 52238         2   6600.06     52238
# 52239         2   6600.19     52239
# 52240         2   6600.31     52240
# ... 33 more variables/columns.
# 
# file(s):
#   01CPTAC_GBM_W_PNNL_20190123_B1S1_f04.mzML.gz
# [1] 2
# MSn data (Spectra) with 52528 spectra in a MsBackendMzR backend:
#   msLevel     rtime scanIndex
# <integer> <numeric> <integer>
# 1             1 0.0271159         1
# 2             1 0.3097275         2
# 3             1 0.5878591         3
# 4             1 0.8702158         4
# 5             1 1.1555794         5
# ...         ...       ...       ...
# 52524         1   6599.75     52524
# 52525         2   6599.88     52525
# 52526         2   6600.01     52526
# 52527         2   6600.15     52527
# 52528         2   6600.28     52528
# ... 33 more variables/columns.
# 
# file(s):
#   01CPTAC_GBM_W_PNNL_20190123_B1S1_f05.mzML.gz
# [1] 3
# MSn data (Spectra) with 51157 spectra in a MsBackendMzR backend:
#   msLevel     rtime scanIndex
# <integer> <numeric> <integer>
# 1             1 0.0272253         1
# 2             1 0.3112476         2
# 3             1 0.5923460         3
# 4             1 0.8760788         4
# 5             1 1.1616845         5
# ...         ...       ...       ...
# 51153         1   6599.81     51153
# 51154         2   6599.93     51154
# 51155         2   6600.07     51155
# 51156         2   6600.20     51156
# 51157         2   6600.33     51157
# ... 33 more variables/columns.
# 
# file(s):
#   09CPTAC_GBM_W_PNNL_20190501_B3S1_f22.mzML.gz


# 4.2 Keeping all matches -------------------------------------------------

spectrums <- c(
  "controllerType=0 controllerNumber=1 scan=11119", 
  "controllerType=0 controllerNumber=1 scan=10989", 
  "controllerType=0 controllerNumber=1 scan=39254"
)

list.id2 <- list()
for(i in 1:length(list.ids)){
  print(i)
  id <- list.ids[[i]]
  I <- which(id$spectrumID == spectrums[i])
  
  # reduce
  id2 <- QFeatures::reduceDataFrame(id, id$spectrumID)
  rownames(id2) <- NULL 
  list.id2[[i]] <- id2
}


# check -------------------------------------------------------------------

for(i in 1:3){
  print(i)
  list.ids[[i]] %>% dim() %>% print()
  list.id2[[i]] %>% dim() %>% print()
}
# [1] 1
# [1] 266767     36
# [1] 38401    36
# [1] 2
# [1] 207109     36
# [1] 40017    36
# [1] 3
# [1] 258869     36
# [1] 35697    36



# filter IDs --------------------------------------------------------------

list.id.filtered <- list()
for(i in 1:3){
  print(i)
  id <- list.ids[[i]]
  
  id_filtered <- filterPSMs(id)
  list.id.filtered[[i]] <- id_filtered
}


# check --------------------------------------------------------------------

list.id.filtered[[i]] %>% names()
#  [1] "sequence"                 "spectrumID"              
#  [3] "chargeState"              "rank"                    
#  [5] "passThreshold"            "experimentalMassToCharge"
#  [7] "calculatedMassToCharge"   "peptideRef"              
#  [9] "modNum"                   "isDecoy"                 
# [11] "post"                     "pre"                     
# [13] "start"                    "end"                     
# [15] "DatabaseAccess"           "DBseqLength"             
# [17] "DatabaseSeq"              "DatabaseDescription"     
# [19] "scan.number.s."           "scan.start.time"         
# [21] "acquisitionNum"           "spectrumFile"            
# [23] "idFile"                   "MS.GF.RawScore"          
# [25] "MS.GF.DeNovoScore"        "MS.GF.SpecEValue"        
# [27] "MS.GF.EValue"             "MS.GF.QValue"            
# [29] "MS.GF.PepQValue"          "modPeptideRef"           
# [31] "modName"                  "modMass"                 
# [33] "modLocation"              "subOriginalResidue"      
# [35] "subReplacementResidue"    "subLocation"  


# join spectra data -------------------------------------------------------

list.sp2 <- list()
for(i in 1:3){
  print(i)
  id_filtered <- list.id.filtered[[i]] 
  sp <- list.sp[[i]]
  sp <- joinSpectraData(sp, id_filtered,
                        by.x = "spectrumId",
                        by.y = "spectrumID")
  list.sp2[[i]] <- sp
}


# check -------------------------------------------------------------------

spectraVariables(list.sp2[[2]])

#  [1] "msLevel"                  "rtime"                   
#  [3] "acquisitionNum"           "scanIndex"               
#  [5] "dataStorage"              "dataOrigin"              
#  [7] "centroided"               "smoothed"                
#  [9] "polarity"                 "precScanNum"             
# [11] "precursorMz"              "precursorIntensity"      
# [13] "precursorCharge"          "collisionEnergy"         
# [15] "isolationWindowLowerMz"   "isolationWindowTargetMz" 
# [17] "isolationWindowUpperMz"   "peaksCount"              
# [19] "totIonCurrent"            "basePeakMZ"              
# [21] "basePeakIntensity"        "ionisationEnergy"        
# [23] "lowMZ"                    "highMZ"                  
# [25] "mergedScan"               "mergedResultScanNum"     
# [27] "mergedResultStartScanNum" "mergedResultEndScanNum"  
# [29] "injectionTime"            "filterString"            
# [31] "spectrumId"               "ionMobilityDriftTime"    
# [33] "scanWindowLowerLimit"     "scanWindowUpperLimit"    
# [35] "sequence"                 "chargeState"             
# [37] "rank"                     "passThreshold"           
# [39] "experimentalMassToCharge" "calculatedMassToCharge"  
# [41] "peptideRef"               "modNum"                  
# [43] "isDecoy"                  "post"                    
# [45] "pre"                      "start"                   
# [47] "end"                      "DatabaseAccess"          
# [49] "DBseqLength"              "DatabaseSeq"             
# [51] "DatabaseDescription"      "scan.number.s."          
# [53] "scan.start.time"          "acquisitionNum.y"        
# [55] "spectrumFile"             "idFile"                  
# [57] "MS.GF.RawScore"           "MS.GF.DeNovoScore"       
# [59] "MS.GF.SpecEValue"         "MS.GF.EValue"            
# [61] "MS.GF.QValue"             "MS.GF.PepQValue"         
# [63] "modPeptideRef"            "modName"                 
# [65] "modMass"                  "modLocation"             
# [67] "subOriginalResidue"       "subReplacementResidue"   
# [69] "subLocation"    


all(is.na(filterMsLevel(list.sp2[[3]], 1)$sequence))
# [1] TRUE


# visualize ------------------------------------------------------------------

for(i in 1:3){
  print(i)
  sp <- list.sp2[[i]]
  j <- which(sp$MS.GF.RawScore > 100)[1]

  setwd(dir.2)
  out.f <- filename
  pdf(out.f, w = 8, h = 4.5)
  plotSpectra(sp[j], labels = addFragments,
              labelPos = 3, labelCol = "steelblue", 
              main = SEQUENCES[i])
  dev.off()
  
  
  # intensity-based filtering applied
  
  pdf(out.f, w = 8, h = 4.5)
  filterIntensity(sp[j], 5e3) %>% 
    plotSpectra(labels = addFragments,
                labelPos = 3, labelCol = "steelblue", 
                main = SEQUENCES[i])
  dev.off()
}


# save data ----------------------------------------------------------------

setwd(dir.2)
out.f <- filename
saveRDS(list.sp2, out.f)


# si ----------------------------------------------------------------------

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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] ggplotify_0.1.0      PSMatch_0.7.1        doParallel_1.0.16   
# [4] iterators_1.0.13     foreach_1.5.1        Spectra_1.5.1       
# [7] ProtGenerics_1.25.2  BiocParallel_1.26.2  readxl_1.3.1        
# [10] tidyverse_1.3.1      RColorBrewer_1.1-3   rtracklayer_1.52.1  
# [13] GenomicRanges_1.44.0 GenomeInfoDb_1.28.1  IRanges_2.26.0      
# [16] S4Vectors_0.30.2     BiocGenerics_0.38.0  magrittr_2.0.3      
# [19] ggtranscript_0.99.9  devtools_2.4.2       usethis_2.0.1       
# [22] ggrepel_0.9.1        cowplot_1.1.1        patchwork_1.1.1     
# [25] clustree_0.4.3       ggraph_2.0.5         harmony_0.1.0       
# [28] Rcpp_1.0.8.3         glmGamPoi_1.4.0      sctransform_0.3.5   
# [31] SeuratObject_4.1.3   Seurat_4.3.0.1       viridis_0.6.2       
# [34] viridisLite_0.4.0    paletteer_1.5.0      ggthemes_4.2.4      
# [37] ggsci_2.9            ggpubr_0.4.0         tidylog_1.0.2       
# [40] forcats_0.5.1        stringr_1.4.0        dplyr_1.0.9         
# [43] purrr_0.3.4          readr_2.1.2          tidyr_1.2.0         
# [46] tibble_3.1.7         ggplot2_3.3.6       
# 
# loaded via a namespace (and not attached):
#   [1] scattermore_0.7             ragg_1.2.5                 
# [3] irlba_2.3.5                 DelayedArray_0.18.0        
# [5] data.table_1.14.2           AnnotationFilter_1.16.0    
# [7] RCurl_1.98-1.3              generics_0.1.2             
# [9] callr_3.7.0                 RANN_2.6.1                 
# [11] future_1.25.0               tzdb_0.1.2                 
# [13] spatstat.data_3.0-1         xml2_1.3.2                 
# [15] lubridate_1.9.2             httpuv_1.6.5               
# [17] SummarizedExperiment_1.22.0 assertthat_0.2.1           
# [19] hms_1.1.0                   promises_1.2.0.1           
# [21] fansi_1.0.3                 restfulr_0.0.13            
# [23] dbplyr_2.1.1                igraph_1.3.1               
# [25] DBI_1.1.2                   htmlwidgets_1.6.2          
# [27] spatstat.geom_3.2-4         ellipsis_0.3.2             
# [29] QFeatures_1.5.1             backports_1.2.1            
# [31] deldir_1.0-6                MatrixGenerics_1.4.3       
# [33] vctrs_0.6.1                 Biobase_2.52.0             
# [35] remotes_2.4.2.1             ROCR_1.0-11                
# [37] abind_1.4-5                 cachem_1.0.6               
# [39] withr_2.5.0                 ggforce_0.3.3              
# [41] progressr_0.10.0            GenomicAlignments_1.28.0   
# [43] MultiAssayExperiment_1.18.0 prettyunits_1.1.1          
# [45] goftest_1.2-2               cluster_2.1.2              
# [47] lazyeval_0.2.2              crayon_1.5.1               
# [49] spatstat.explore_3.2-1      labeling_0.4.2             
# [51] pkgconfig_2.0.3             tweenr_1.0.2               
# [53] nlme_3.1-152                pkgload_1.2.1              
# [55] rlang_1.1.0                 globals_0.15.0             
# [57] lifecycle_1.0.3             miniUI_0.1.1.1             
# [59] modelr_0.1.8                cellranger_1.1.0           
# [61] rprojroot_2.0.2             polyclip_1.10-0            
# [63] matrixStats_0.62.0          lmtest_0.9-40              
# [65] Matrix_1.6-0                carData_3.0-4              
# [67] zoo_1.8-10                  reprex_2.0.1               
# [69] ggridges_0.5.3              processx_3.5.2             
# [71] mzR_2.26.1                  png_0.1-7                  
# [73] rjson_0.2.20                clisymbols_1.2.0           
# [75] bitops_1.0-7                KernSmooth_2.23-20         
# [77] Biostrings_2.60.2           parallelly_1.31.1          
# [79] spatstat.random_3.1-5       gridGraphics_0.5-1         
# [81] rstatix_0.7.0               ggsignif_0.6.2             
# [83] scales_1.2.0                memoise_2.0.0              
# [85] plyr_1.8.7                  ica_1.0-2                  
# [87] zlibbioc_1.38.0             compiler_4.1.3             
# [89] BiocIO_1.2.0                clue_0.3-60                
# [91] fitdistrplus_1.1-5          Rsamtools_2.8.0            
# [93] cli_3.6.1                   XVector_0.32.0             
# [95] listenv_0.8.0               pbapply_1.5-0              
# [97] ps_1.6.0                    MASS_7.3-54                
# [99] tidyselect_1.1.2            stringi_1.7.6              
# [101] textshaping_0.3.6           yaml_2.3.5                 
# [103] grid_4.1.3                  tools_4.1.3                
# [105] timechange_0.2.0            future.apply_1.8.1         
# [107] rio_0.5.27                  rstudioapi_0.13            
# [109] MsCoreUtils_1.5.1           foreign_0.8-81             
# [111] gridExtra_2.3               farver_2.1.0               
# [113] Rtsne_0.16                  digest_0.6.29              
# [115] shiny_1.7.1                 car_3.0-11                 
# [117] broom_0.7.9                 later_1.3.0                
# [119] ncdf4_1.17.1                RcppAnnoy_0.0.19           
# [121] httr_1.4.3                  colorspace_2.0-3           
# [123] rvest_1.0.1                 XML_3.99-0.7               
# [125] fs_1.5.2                    tensor_1.5                 
# [127] reticulate_1.25             splines_4.1.3              
# [129] yulab.utils_0.0.4           uwot_0.1.16                
# [131] rematch2_2.1.2              spatstat.utils_3.0-3       
# [133] graphlayouts_0.7.1          sp_2.0-0                   
# [135] systemfonts_1.0.4           plotly_4.10.0              
# [137] sessioninfo_1.1.1           xtable_1.8-4               
# [139] jsonlite_1.8.0              tidygraph_1.2.0            
# [141] testthat_3.0.4              R6_2.5.1                   
# [143] pillar_1.7.0                htmltools_0.5.5            
# [145] mime_0.12                   glue_1.6.2                 
# [147] fastmap_1.1.0               codetools_0.2-18           
# [149] pkgbuild_1.2.0              utf8_1.2.2                 
# [151] lattice_0.20-44             spatstat.sparse_3.0-2      
# [153] curl_4.3.2                  leiden_0.3.9               
# [155] zip_2.2.0                   openxlsx_4.2.4             
# [157] survival_3.2-12             desc_1.3.0                 
# [159] munsell_0.5.0               GenomeInfoDbData_1.2.6     
# [161] haven_2.4.3                 reshape2_1.4.4             
# [163] gtable_0.3.0     


# end ---------------------------------------------------------------------
