#src to get interactors
#Buffer <- xxx

#load library
#load library
library(tidyverse)
library(dplyr)
library(BayesFactor)
library(ggplot2)

#load scr
source("./src/MIX ratio.R")
source("./src/IDIRT ratio.R")
source("./src/get_interaction.R")


#load data
pep <- read.delim("./IDIRT data/IDIRTall_PeptideGroups.txt")
info <- read.delim("./IDIRT data/IDIRTall_InputFiles.txt")

#check the mix ratio
mix <- mix_ratio()


buffer <- function(buffer){
  res_all <- idirt_merge(buffer)
  res <- Bayesian_test(res_all,buffer,mix)
  
}

Buffer2 <- buffer("buffer2")
Buffer3 <- buffer("buffer3")
Buffer16 <- buffer("buffer16")
Buffer19 <- buffer("buffer19")
Buffer25 <- buffer("buffer25")
Buffer29 <- buffer("buffer29")
