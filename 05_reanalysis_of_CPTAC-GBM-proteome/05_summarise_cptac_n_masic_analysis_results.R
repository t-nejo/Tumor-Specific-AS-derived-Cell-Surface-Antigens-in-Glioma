# Takahide Nejo, MD, PhD 
# Okada Lab, Department of Neurological Surgery, UC San Francisco


# rm all ------------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# load packages -----------------------------------------------------------

library(tidyverse)
library(tidylog)
library(readxl)


# dir ---------------------------------------------------------------------

dir.1 <- "~/path/to/data"
dir.1b <- "~/path/to/MASIC_output_data"
dir.2 <- "~/path/to/output_folder"


# load data -----------------------------------------------------------------

setwd(dir.1)
in.f <- "junc_id_of_the_final_14_events.tsv"
res.0 <- read_tsv(in.f)

# annotation
  # downloaded from CPTAC-GBM project webpage
setwd(dir.1)
in.f <- "S048_CPTAC_GBM_Discovery_Cohort_TMT11_CaseID_SampleID_AliquotID_Map_Dec2019_r1.xlsx" ;
annot.orig <- read_xlsx(in.f, skip = 7)


# check/edit annot --------------------------------------------------------------

annot.orig %>% dim() %>% print()
# [1] 121  15
annot.orig %>% pull(`TMT plex`) %>% table()
#  1  2  3  4  5  6  7  8  9 10 11 
# 11 11 11 11 11 11 11 11 11 11 11 


# edit --------------------------------------------------------------------

annot <- annot.orig %>% 
  dplyr::select(alias = Alias, TMT.channel = `TMT channel`, case.id = `Case ID (Participant ID)`, sample.type = `Sample type`) %>% 
  mutate(TMT.channel = paste0("TMT11-", TMT.channel)) %>% 
  mutate(sample.type = ifelse( (is.na(sample.type) | sample.type == "NA") & case.id == "ref", "ref", sample.type)) %>% 
  dplyr::filter(!grepl("Withdrawn", case.id))


# check -------------------------------------------------------------------

annot %>% dim() %>% print() 
# [1] 120   4

annot %>% head()
# # A tibble: 6 × 4
# alias TMT.channel case.id                      sample.type
# <chr> <chr>       <chr>                        <chr>      
# 1 B1S1  TMT11-126   ref                          ref        
# 2 B1S1  TMT11-127N  GTEX-Y8DK-0011-R10A-SM-HAKY1 normal     
# 3 B1S1  TMT11-127C  C3N-03183                    tumor      
# 4 B1S1  TMT11-128N  C3N-01505                    tumor      
# 5 B1S1  TMT11-128C  C3N-03188                    tumor      
# 6 B1S1  TMT11-129N  C3L-02984                    tumor  

annot %>% tail()
# # A tibble: 6 × 4
# alias TMT.channel case.id                      sample.type
# <chr> <chr>       <chr>                        <chr>      
# 1 B3S3  TMT11-129N  C3N-01515                    tumor      
# 2 B3S3  TMT11-129C  C3L-00104                    tumor      
# 3 B3S3  TMT11-130N  C3N-02769                    tumor      
# 4 B3S3  TMT11-130C  C3N-02785                    tumor      
# 5 B3S3  TMT11-131N  GTEX-NPJ7-0011-R10A-SM-HAKXW normal     
# 6 B3S3  TMT11-131C  C3N-02190                    tumor      

annot %>% anyNA()
# [1] FALSE

annot %>% pull(sample.type) %>% table() %>% print()
# normal    ref  tumor 
#     10     11     99 


# load files --------------------------------------------------------------

list.res <- list() 
for(i in 1:11){
  I <- formatC(i, width = 2, flag = 0)
  print(I)
  

  # load MSGF+ output -------------------------------------------------------

  setwd(dir.1) 
  files <- list.files()[grepl(paste0(I, "CPTAC"), list.files())]
  files <- files[grepl("tsv", files)]
  
  setwd(dir.1b)
  dir.i <- list.files()[grepl(paste0(I, "CPTAC"), list.files())] ;
  
  res.i <- NULL
  for(j in 1:24){
    J = formatC(j, width = 2, flag = 0)
    print( paste0(I, " - ", J))
    
    setwd(dir.1)
    in.f <- files[grepl(paste0("_f", J), files)]
    
    res.j <- read_tsv(in.f, na = c("", "NA"), col_names = T) %>% 
      dplyr::filter(!grepl("XXX", Protein)) %>% # remove the all lines for decoys
      arrange(ScanNum) %>% 
      mutate(alias = substr(in.f, 29, 32)) %>% 
      mutate(ProteinPeptideFormula = paste0(Protein, ";", Peptide, ";", Formula)) %>% 
      dplyr::filter(grepl("alt_prot", ProteinPeptideFormula)) %>% 
      dplyr::filter(!grepl("sp", ProteinPeptideFormula)) %>% 
      dplyr::filter(!grepl("tr", ProteinPeptideFormula)) 
    
    res.j <- res.j %>% 
      mutate(junc.name = gsub("_ENST.*", "", ProteinPeptideFormula)) %>% 
      semi_join(res.0, by = c("junc.name" = "name")) %>% 
      dplyr::select(- junc.name)
    
    if(nrow(res.j) > 0){
      setwd( paste0(dir.1b, "/", dir.i))
      in.f <- list.files()[grepl(paste0("_f", J, "_ReporterIons.txt"), list.files())]
      res.masic.j <- read_tsv(in.f, na = c("", "NA"), col_names = T) %>% 
        dplyr::select(ScanNumber, ParentIonMZ:Ion_131.144)
      
      colnames(res.masic.j)[6:16] = paste0("TMT11-", c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C"))
      
      res.j <- res.j %>% 
        left_join(res.masic.j, by = c("ScanNum" = "ScanNumber"))
      
      res.i <- res.i %>% 
        bind_rows(res.j)
    }
  }
  list.res[[i]] <- res.i
}


# check -------------------------------------------------------------------

for(i in 1:11){
  # print(i)
  list.res[[i]] %>% nrow() %>% print()
}
# [1] 10
# [1] 7
# [1] 3
# [1] 7
# [1] 12
# [1] 8
# [1] 5
# [1] 10
# [1] 9
# [1] 16
# [1] 8


# filtering based on Q-values ---------------------------------------------------------------

res.1 <- NULL
for(i in 1:11){
  print(i)
  res.i <- list.res[[i]] %>% 
    dplyr::filter(PepQValue < 0.2)
  
  if( nrow(res.i) > 0){
    res.1 <- res.1 %>% 
      bind_rows(res.i)
  }
}


# check -------------------------------------------------------------------

res.1 %>% dim()
# [1] 43 34

res.1 %>% head()
# # A tibble: 6 × 34
# `#SpecFile`    SpecID ScanNum FragMethod Precursor IsotopeError `PrecursorErro…` Charge
# <chr>          <chr>    <dbl> <chr>          <dbl>        <dbl>            <dbl>  <dbl>
# 1 01CPTAC_GBM_W… contr…   36746 HCD             921.            1           -0.649      5
# 2 01CPTAC_GBM_W… contr…   36757 HCD             921.            0            2.98       5
# 3 01CPTAC_GBM_W… contr…   11119 HCD             513.            0           -0.595      4
# 4 01CPTAC_GBM_W… contr…   10989 HCD             565.            0            0.108      3
# 5 01CPTAC_GBM_W… contr…   15497 HCD             460.            1            8.21       3
# 6 01CPTAC_GBM_W… contr…   50255 HCD             782.            0           11.6        2
# # … with 26 more variables: Peptide <chr>, Formula <chr>, Protein <chr>,
# #   DeNovoScore <dbl>, MSGFScore <dbl>, SpecEValue <dbl>, EValue <dbl>, QValue <dbl>,
# #   PepQValue <dbl>, alias <chr>, ProteinPeptideFormula <chr>, ParentIonMZ <dbl>,
# #   BasePeakIntensity <dbl>, BasePeakMZ <dbl>, ReporterIonIntensityMax <dbl>,
# #   `TMT11-126` <dbl>, `TMT11-127N` <dbl>, `TMT11-127C` <dbl>, `TMT11-128N` <dbl>,
# #   `TMT11-128C` <dbl>, `TMT11-129N` <dbl>, `TMT11-129C` <dbl>, `TMT11-130N` <dbl>,
# #   `TMT11-130C` <dbl>, `TMT11-131N` <dbl>, `TMT11-131C` <dbl>


# edit --------------------------------------------------------------------

res.1b <- res.1 %>% 
  mutate(symbol_juncid = gsub("_alt_prot.*", "", ProteinPeptideFormula)) %>% 
  arrange(symbol_juncid) %>% 
  dplyr::select("symbol_juncid", colnames(res.1))


# check -------------------------------------------------------------------

res.1b %>% dim()
# [1] 43 35

res.1b %>% distinct(symbol_juncid)
# distinct: removed 38 rows (88%), 5 rows remaining
# # A tibble: 5 × 1
# symbol_juncid            
# <chr>                    
# 1 BCAN-03_ENST00000329117  
# 2 EGFRvIII_ENST00000275493 
# 3 NCAM1-01_ENST00000316851 
# 4 NRCAM-01_ENST00000379028 
# 5 PTPRZ1-02_ENST00000449182

res.1b %>% 
  dplyr::select(symbol_juncid, EValue, QValue) %>% 
  arrange(QValue) %>% 
  distinct(symbol_juncid, .keep_all = T)
# # A tibble: 5 × 3
# symbol_juncid               EValue   QValue
# <chr>                        <dbl>    <dbl>
# 1 EGFRvIII_ENST00000275493  5.26e-11 0       
# 2 PTPRZ1-02_ENST00000449182 7.77e-28 0       
# 3 NCAM1-01_ENST00000316851  6.47e- 4 0.000130
# 4 BCAN-03_ENST00000329117   2.56e- 1 0.0188  
# 5 NRCAM-01_ENST00000379028  3.58e+ 0 0.162   

res.1b %>% 
  pull(ProteinPeptideFormula) %>% 
  unique() %>% 
  length() %>% print()
# [1] 11

res.1b %>% 
  pull(ProteinPeptideFormula) %>% 
  unique() 
# [1] "BCAN-03_ENST00000329117_alt_prot(pre=R,post=A);+229.163TLLEEEGK+229.163;H(67) C(39) N(9) O(16) 458.3258"                                                                                                                                              
# [2] "BCAN-03_ENST00000329117_alt_prot(pre=R,post=A);+229.163GAIYSIPIMEDGGGGSSTPEDPAEAPRTLLEEEGK+229.163;H(241) C(152) N(39) O(58) S 458.3258"                                                                                                              
# [3] "EGFRvIII_ENST00000275493_alt_prot(pre=K,post=A);+229.163K+229.163GNYVVTDHGSC+57.021VR;H(106) C(66) N(22) O(22) S 458.3258"                                                                                                                            
# [4] "EGFRvIII_ENST00000275493_alt_prot(pre=K,post=A);+229.163GNYVVTDHGSC+57.021VR;H(94) C(60) N(20) O(21) S 229.1629"                                                                                                                                      
# [5] "NCAM1-01_ENST00000316851_alt_prot(pre=K,post=G);+229.163DVLSNNYLQIR;H(95) C(58) N(17) O(19) 229.1629"                                                                                                                                                 
# [6] "NCAM1-01_ENST00000316851_alt_prot(pre=K,post=L);+229.163TQP+15.995VPGEPSAPK+229.163;H(86) C(53) N(14) O(19) 458.3258"                                                                                                                                 
# [7] "NCAM1-01_ENST00000316851_alt_prot(pre=K,post=L);+229.163TQP+15.995VP+15.995GEPSAPK+229.163;H(86) C(53) N(14) O(20) 458.3258"                                                                                                                          
# [8] "NCAM1-01_ENST00000316851_alt_prot(pre=K,post=L);+229.163TQPVP+15.995GEPSAPK+229.163;H(86) C(53) N(14) O(19) 458.3258"                                                                                                                                 
# [9] "NCAM1-01_ENST00000316851_alt_prot(pre=-,post=F);+229.163MFC+57.021FVFSVSLQVDIVP+15.995SQGEISVGESK+229.163;H(209) C(134) N(31) O(43) S(2) 458.3258"                                                                                                    
# [10] "NRCAM-01_ENST00000379028_alt_prot(pre=R,post=Q);+229.163VFNTP+15.995EGAM+15.995ASR;H(86) C(54) N(16) O(20) S 229.1629"                                                                                                                                
# [11] "PTPRZ1-02_ENST00000449182_alt_prot(pre=K,post=H);PTPRZ1-03_ENST00000449182_alt_prot(pre=K,post=H);PTPRZ1-01_ENST00000449182_alt_prot(pre=K,post=H);+229.163HVADLHASSGFTEEFEEVQSC+57.021TVDLGITADSSNHPDNK+229.163;H(263) C(174) N(49) O(67) S 458.3258"


# check -------------------------------------------------------------------

res.2 <- NULL
for(i in 1:11){
  print(i)
  res.i <- list.res[[i]] %>% 
    dplyr::filter(PepQValue < 0.2)
  
  if( nrow(res.i) > 0){
    res.i.2 <- res.i %>% 
      dplyr::group_by(ProteinPeptideFormula) %>% 
      dplyr::select(ProteinPeptideFormula, 24:34) %>% 
      dplyr::summarise(
        ProteinPeptideFormula = ProteinPeptideFormula,  ######
        `TMT11-126` = sum(`TMT11-126`, na.rm = T), 
        `TMT11-127N` = sum(`TMT11-127N`, na.rm = T), 
        `TMT11-127C` = sum(`TMT11-127C`, na.rm = T), 
        `TMT11-128N` = sum(`TMT11-128N`, na.rm = T), 
        `TMT11-128C` = sum(`TMT11-128C`, na.rm = T), 
        `TMT11-129N` = sum(`TMT11-129N`, na.rm = T), 
        `TMT11-129C` = sum(`TMT11-129C`, na.rm = T), 
        `TMT11-130N` = sum(`TMT11-130N`, na.rm = T), 
        `TMT11-130C` = sum(`TMT11-130C`, na.rm = T), 
        `TMT11-131N` = sum(`TMT11-131N`, na.rm = T), 
        `TMT11-131C` = sum(`TMT11-131C`, na.rm = T)
      ) %>% 
      distinct(ProteinPeptideFormula, .keep_all = T) ;
    
    res.i.2 <- res.i.2 %>% 
      gather("TMT.channel", "intensity", `TMT11-127N`:`TMT11-131C`) %>% 
      mutate(relative.intensity = (intensity + 1000)  / (`TMT11-126` + 1000) ) %>%  # calculate relative intensity to ref
      mutate(alias = substr(res.i$`#SpecFile`[1], 29, 32)) %>% 
      left_join(annot, by = c("alias" = "alias", "TMT.channel" = "TMT.channel")) %>%
      dplyr::select(ProteinPeptideFormula, alias, case.id, relative.intensity) %>% 
      mutate(sample.type = ifelse(grepl("GTEX", case.id), "normal", "tumor"))

    res.2 <- res.2 %>% 
      bind_rows(res.i.2)
  }
}

res.2 <- res.2 %>% 
  arrange(ProteinPeptideFormula)


# check -------------------------------------------------------------------

res.2 %>% dim()
# [1] 320   5

res.2 %>% head()
# # A tibble: 6 × 5
# # Groups:   ProteinPeptideFormula [1]
# ProteinPeptideFormula                        alias case.id relative.intens… sample.type
# <chr>                                        <chr> <chr>              <dbl> <chr>      
# 1 BCAN-03_ENST00000329117_alt_prot(pre=R,post… B3S1  GTEX-W…             1.51 normal     
# 2 BCAN-03_ENST00000329117_alt_prot(pre=R,post… B3S1  C3N-01…             1.41 tumor      
# 3 BCAN-03_ENST00000329117_alt_prot(pre=R,post… B3S1  C3N-00…             1.14 tumor      
# 4 BCAN-03_ENST00000329117_alt_prot(pre=R,post… B3S1  C3L-03…             1.15 tumor      
# 5 BCAN-03_ENST00000329117_alt_prot(pre=R,post… B3S1  C3L-03…             1.69 tumor      
# 6 BCAN-03_ENST00000329117_alt_prot(pre=R,post… B3S1  C3L-02…             1.68 tumor     

res.2 %>% 
  dplyr::select(ProteinPeptideFormula, sample.type) %>% 
  table() %>% 
  as_tibble() %>% 
  spread(sample.type, n) %>% 
  dplyr::select(tumor, normal, ProteinPeptideFormula)
# # A tibble: 11 × 3
# tumor normal ProteinPeptideFormula                                                               
# <int>  <int> <chr>                                                                               
# 1     18      2 BCAN_chr1:+:156621481-156622081_alt_prot(pre=R,post=A);+229.163GAIYSIPIMEDGGGGSSTPE…
# 2      9      1 BCAN_chr1:+:156621481-156622081_alt_prot(pre=R,post=A);+229.163TLLEEEGK+229.163;H(6…
# 3     82      8 EGFRvIII_chr7:+:55087058-55223522_alt_prot(pre=K,post=A);+229.163GNYVVTDHGSC+57.021…
# 4     64      6 EGFRvIII_chr7:+:55087058-55223522_alt_prot(pre=K,post=A);+229.163K+229.163GNYVVTDHG…
# 5      9      1 NCAM1_chr11:+:113076388-113076776_alt_prot(pre=-,post=F);+229.163MFC+57.021FVFSVSLQ…
# 6      9      1 NCAM1_chr11:+:113076388-113076776_alt_prot(pre=K,post=G);+229.163DVLSNNYLQIR;H(95) …
# 7      9      1 NCAM1_chr11:+:113076388-113076776_alt_prot(pre=K,post=L);+229.163TQP+15.995VP+15.99…
# 8     28      2 NCAM1_chr11:+:113076388-113076776_alt_prot(pre=K,post=L);+229.163TQP+15.995VPGEPSAP…
# 9     27      3 NCAM1_chr11:+:113076388-113076776_alt_prot(pre=K,post=L);+229.163TQPVP+15.995GEPSAP…
# 10     9      1 NRCAM_chr7:-:107800937-107820666_alt_prot(pre=R,post=Q);+229.163VFNTP+15.995EGAM+15…
# 11    27      3 PTPRZ1_chr7:+:121568275-121608007_alt_prot(pre=K,post=H);PTPRZ1_chr7:+:121644714-12…
                                                                                                 

# check -------------------------------------------------------------------

res.1b %>% dim()
# [1] 43 35

res.1b %>% colnames()
# [1] "symbol_juncid"           "#SpecFile"              
# [3] "SpecID"                  "ScanNum"                
# [5] "FragMethod"              "Precursor"              
# [7] "IsotopeError"            "PrecursorError(ppm)"    
# [9] "Charge"                  "Peptide"                
# [11] "Formula"                 "Protein"                
# [13] "DeNovoScore"             "MSGFScore"              
# [15] "SpecEValue"              "EValue"                 
# [17] "QValue"                  "PepQValue"              
# [19] "alias"                   "ProteinPeptideFormula"  
# [21] "ParentIonMZ"             "BasePeakIntensity"      
# [23] "BasePeakMZ"              "ReporterIonIntensityMax"
# [25] "TMT11-126"               "TMT11-127N"             
# [27] "TMT11-127C"              "TMT11-128N"             
# [29] "TMT11-128C"              "TMT11-129N"             
# [31] "TMT11-129C"              "TMT11-130N"             
# [33] "TMT11-130C"              "TMT11-131N"             
# [35] "TMT11-131C"    


# edit --------------------------------------------------------------------

res.1c <- res.1b %>% 
  dplyr::filter(QValue < 0.05) %>% 
  mutate(pep = gsub("[0-9]", "", gsub("-", "", gsub("\\+", "", gsub("\\.", "", Peptide))))) %>% 
  mutate(symbol = gsub("_.*", "", symbol_juncid)) %>% 
  mutate(junc.id = gsub(".*_chr", "chr", symbol_juncid)) %>% 
  mutate(batch = substr(`#SpecFile`, 1, 2)) %>% 
  dplyr::select(symbol, junc.id, pep, Peptide, EValue, QValue, batch, `#SpecFile`, ScanNum, alias) 


# check -------------------------------------------------------------------

res.1c %>% 
  arrange(QValue, EValue) %>% 
  distinct(junc.id, pep, .keep_all = T) %>% 
  dplyr::select(junc.id, pep, QValue, EValue) 
# # A tibble: 5 × 4
# junc.id                   pep                                QValue   EValue
# <chr>                     <chr>                               <dbl>    <dbl>
# 1 PTPRZ1-02_ENST00000449182 HVADLHASSGFTEEFEEVQSCTVDLGITADSS… 0       7.77e-28
# 2 EGFRvIII_ENST00000275493  GNYVVTDHGSCVR                     0       1.56e-11
# 3 EGFRvIII_ENST00000275493  KGNYVVTDHGSCVR                    0       5.26e-11
# 4 NCAM1-01_ENST00000316851  TQPVPGEPSAPK                      1.30e-4 6.47e- 4
# 5 BCAN-03_ENST00000329117   GAIYSIPIMEDGGGGSSTPEDPAEAPRTLLEE… 1.88e-2 2.56e- 1

res.1c %>% 
  arrange(QValue, EValue) %>% 
  distinct(junc.id, pep, .keep_all = T) %>% 
  pull(pep)
# [1] "HVADLHASSGFTEEFEEVQSCTVDLGITADSSNHPDNK" -> wt seq
# [2] "GNYVVTDHGSCVR"                         
# [3] "KGNYVVTDHGSCVR"                        
# [4] "TQPVPGEPSAPK"                  -> wt seq        
# [5] "GAIYSIPIMEDGGGGSSTPEDPAEAPRTLLEEEGK"   -> mt only. 

# NOTE
# the following two aa seq were found to be WT-derived. 
# [1] "HVADLHASSGFTEEFEEVQSCTVDLGITADSSNHPDNK" -> wt seq
# [4] "TQPVPGEPSAPK"                  -> wt seq        


res.1d <- res.1c %>% 
  dplyr::filter(!grepl("PTPRZ1", symbol)) %>%
  dplyr::filter(!grepl("NCAM1", symbol)) %>%
  arrange(pep, Peptide, QValue) %>% 
  distinct(pep, .keep_all = T) %>% 
  dplyr::select(symbol, junc.id, pep, EValue, QValue) %>% 
  arrange(QValue)


res.1e <- NULL
for(i in 1:nrow(res.1d)){
  print(i)
  pep.i <- res.1d[i, ] %>% 
    pull(pep)
  
  res.i <- res.1c %>% 
    dplyr::filter(pep == pep.i)
  
  all.peptides.i <- res.i %>% 
    pull(Peptide) %>% 
    unique()
  
  res.2.i <- NULL
  for(j in 1:length(all.peptides.i)){
    print(j)
    pep.j <- all.peptides.i[j]
    res.j <- res.i %>% 
      dplyr::filter(Peptide == pep.j)
    
    res.j <- tibble(
      Peptide = pep.j, 
      n.batch = res.j %>% distinct(batch) %>% nrow(), 
      n.event = res.j %>% nrow()
    )
    res.2.i <- res.2.i %>% 
      bind_rows(res.j)

  }
  res.3.i <- res.i %>% 
    arrange(QValue, EValue) %>% 
    dplyr::select(- batch, - `#SpecFile`, - ScanNum, - alias) %>% 
    distinct(Peptide, .keep_all = T) %>% 
    left_join(res.2.i, by = "Peptide") %>% 
    dplyr::rename(peptide = pep, detected.peptide.signal = Peptide)
  
  res.1e <- res.1e %>% 
    bind_rows(res.3.i)
}

res.1e
# # A tibble: 3 × 8
# symbol   junc.id    peptide detected.peptid…   EValue QValue n.batch n.event
# <chr>    <chr>      <chr>   <chr>               <dbl>  <dbl>   <int>   <int>
# 1 EGFRvIII EGFRvIII_… GNYVVT… +229.163GNYVVTD… 1.56e-11 0            9      14
# 2 EGFRvIII EGFRvIII_… KGNYVV… +229.163K+229.1… 5.26e-11 0            7       8
# 3 BCAN-03  BCAN-03_E… GAIYSI… +229.163GAIYSIP… 2.56e- 1 0.0188       1       1


# out.table ---------------------------------------------------------------------

setwd(dir.2) ;
out.f <- filename # "cptac_analysis_result_summarizing_msgfplus_n_masic.tsv
write_tsv(res.1e, out.f)


# si ----------------------------------------------------------------------

Sys.time() %>% print() ; 
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
# [1] readxl_1.3.1                viridis_0.6.2               viridisLite_0.4.0          
# [4] paletteer_1.5.0             ggthemes_4.2.4              ggsci_2.9                  
# [7] ggpubr_0.4.0                tidylog_1.0.2               DESeq2_1.32.0              
# [10] SummarizedExperiment_1.22.0 Biobase_2.52.0              MatrixGenerics_1.4.3       
# [13] matrixStats_0.62.0          GenomicRanges_1.44.0        GenomeInfoDb_1.28.1        
# [16] IRanges_2.26.0              S4Vectors_0.30.2            BiocGenerics_0.38.0        
# [19] forcats_0.5.1               stringr_1.4.0               dplyr_1.0.9                
# [22] purrr_0.3.4                 readr_2.1.2                 tidyr_1.2.0                
# [25] tibble_3.1.7                ggplot2_3.3.6               tidyverse_1.3.1            
# 
# loaded via a namespace (and not attached):
# [1] colorspace_2.0-3       ggsignif_0.6.2         ellipsis_0.3.2        
# [4] rio_0.5.27             XVector_0.32.0         fs_1.5.2              
# [7] rstudioapi_0.13        bit64_4.0.5            AnnotationDbi_1.54.1  
# [10] fansi_1.0.3            lubridate_1.9.2        xml2_1.3.2            
# [13] splines_4.1.3          cachem_1.0.6           geneplotter_1.70.0    
# [16] jsonlite_1.8.0         broom_0.7.9            annotate_1.70.0       
# [19] dbplyr_2.1.1           png_0.1-7              compiler_4.1.3        
# [22] httr_1.4.3             backports_1.2.1        assertthat_0.2.1      
# [25] Matrix_1.6-0           fastmap_1.1.0          cli_3.6.1             
# [28] tools_4.1.3            gtable_0.3.0           glue_1.6.2            
# [31] GenomeInfoDbData_1.2.6 Rcpp_1.0.8.3           carData_3.0-4         
# [34] cellranger_1.1.0       vctrs_0.6.1            Biostrings_2.60.2     
# [37] openxlsx_4.2.4         rvest_1.0.1            timechange_0.2.0      
# [40] lifecycle_1.0.3        rstatix_0.7.0          XML_3.99-0.7          
# [43] zlibbioc_1.38.0        scales_1.2.0           vroom_1.5.7           
# [46] clisymbols_1.2.0       hms_1.1.0              rematch2_2.1.2        
# [49] RColorBrewer_1.1-3     curl_4.3.2             memoise_2.0.0         
# [52] gridExtra_2.3          stringi_1.7.6          RSQLite_2.2.8         
# [55] genefilter_1.74.0      zip_2.2.0              BiocParallel_1.26.2   
# [58] rlang_1.1.0            pkgconfig_2.0.3        bitops_1.0-7          
# [61] lattice_0.20-44        bit_4.0.4              tidyselect_1.1.2      
# [64] magrittr_2.0.3         R6_2.5.1               generics_0.1.2        
# [67] DelayedArray_0.18.0    DBI_1.1.2              pillar_1.7.0          
# [70] haven_2.4.3            foreign_0.8-81         withr_2.5.0           
# [73] survival_3.2-12        KEGGREST_1.32.0        abind_1.4-5           
# [76] RCurl_1.98-1.3         modelr_0.1.8           crayon_1.5.1          
# [79] car_3.0-11             utf8_1.2.2             tzdb_0.1.2            
# [82] locfit_1.5-9.4         grid_4.1.3             data.table_1.14.2     
# [85] blob_1.2.2             reprex_2.0.1           xtable_1.8-4          
# [88] munsell_0.5.0    


# end ---------------------------------------------------------------------
