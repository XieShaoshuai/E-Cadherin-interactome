#analysis IDIRT and IDIRT swap ratio for E-Cadherin

#load library
library(dplyr)
library(BayesFactor)
library(tidyverse)
library(chemometrics)


#load data
mix_ratio <- function() {
  pep <- read.delim("./Mix data/Mix_all_PeptideGroups.txt")
  info <- read.delim("./Mix data/metadata.csv",sep=",")
  
  
  #data clean
  pep <- pep %>% select(Sequence,Master.Protein.Accessions,Quan.Info,
                        contains("Abundance.F"))
  
  pep<- pep %>% filter(Quan.Info=="") 
  
  #replace NA with 0
  pep <- replace(pep,is.na(pep),0)
  
  #calculate mix ratio------------------
  mix <- pep %>% select(Sequence,Master.Protein.Accessions,
                        contains(info$colname[which(info$experiment==1)]))
  
  sample <- mix %>% select(contains("Sample"))
  control <- mix %>% select(contains("Control"))
  ratio <- sample/(sample+control)
  
  
  ratio$Sequence <- mix$Sequence
  ratio$Master.Protein.Accessions <- mix$Master.Protein.Accessions
  
  ratio <- ratio %>% group_by(Master.Protein.Accessions) %>% 
    summarise(ratio1 = median(Abundance.F1.Heavy.Sample.Hevay_HA,na.rm=TRUE),
              ratio2 = median(Abundance.F2.Heavy.Sample.Hevay_HA,na.rm=TRUE),
              ratio3 = median(Abundance.F3.Heavy.Sample.Hevay_HA,na.rm=TRUE),
              ratio4 = median(Abundance.F4.Heavy.Sample.Hevay_HA,na.rm=TRUE),
              ratio5 = median(Abundance.F5.Heavy.Sample.Hevay_HA,na.rm=TRUE),
              ratio6 = median(Abundance.F6.Heavy.Sample.Hevay_HA,na.rm=TRUE),)
  
  mix_ratio <- c(median(ratio$ratio1,na.rm=TRUE),
                 median(ratio$ratio2,na.rm=TRUE),
                 median(ratio$ratio3,na.rm=TRUE),
                 median(ratio$ratio4,na.rm=TRUE),
                 median(ratio$ratio5,na.rm=TRUE),
                 median(ratio$ratio6,na.rm=TRUE))
  
  
  
  
  
  #calculate mix_wap ratio
  swap <- pep %>% select(Sequence,Master.Protein.Accessions,
                        contains(info$colname[which(info$experiment==2)]))
  
  sample <- swap %>% select(contains("Light.Sample"))
  control <- swap %>% select(contains("Heavy.Sample"))
  ratio <- sample/(sample+control)
  
  
  ratio$Sequence <- swap$Sequence
  ratio$Master.Protein.Accessions <- mix$Master.Protein.Accessions
  
  ratio <- ratio %>% group_by(Master.Protein.Accessions) %>% 
    summarise(ratio1 = median(Abundance.F13.Light.Sample.Heavy_wt,na.rm=TRUE),
              ratio2 = median(Abundance.F14.Light.Sample.Heavy_wt,na.rm=TRUE),
              ratio3 = median(Abundance.F15.Light.Sample.Heavy_wt,na.rm=TRUE),
              ratio4 = median(Abundance.F16.Light.Sample.Heavy_wt,na.rm=TRUE),
              ratio5 = median(Abundance.F17.Light.Sample.Heavy_wt,na.rm=TRUE),
              ratio6 = median(Abundance.F18.Light.Sample.Heavy_wt,na.rm=TRUE),)
  
  swap_ratio <- c(median(ratio$ratio1,na.rm=TRUE),
                 median(ratio$ratio2,na.rm=TRUE),
                 median(ratio$ratio3,na.rm=TRUE),
                 median(ratio$ratio4,na.rm=TRUE),
                 median(ratio$ratio5,na.rm=TRUE),
                 median(ratio$ratio6,na.rm=TRUE))
  
  res <- list(mix=mix_ratio, mix_swap=swap_ratio)
  return(res)
}


