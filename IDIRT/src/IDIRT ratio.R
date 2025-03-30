#analysis IDIRT and IDIRT swap ratio for E-Cadherin

#load library
library(dplyr)
library(BayesFactor)
library(tidyverse)
library(chemometrics)


#load data
pep <- read.delim("./IDIRT data/IDIRTall_PeptideGroups.txt")
info <- read.delim("./IDIRT data/IDIRTall_InputFiles.txt")

#extract buffer
#Buffer <- "buffer2"
#exp="IDIRT"
get_ratio <- function(Buffer,exp) {
  index <- info %>% filter(buffer==Buffer,experiment==exp) %>% 
    select(File.ID)
  df <- pep %>%  
    select(Annotated.Sequence,"Modifications","Master.Protein.Accessions","Quan.Info",
           starts_with("Abundance.F")) %>% 
    select(Annotated.Sequence,"Modifications","Master.Protein.Accessions","Quan.Info", 
           contains(paste0(index$File.ID,".")))
  
  
  #rename df colnames from Heavy/Light to Sample/control
  if(exp=="IDIRT"){
    colnames(df) <- gsub(".Sample", "", colnames(df))
    colnames(df) <- gsub("Heavy", "Sample", colnames(df))
    colnames(df) <- gsub("Light", "Control", colnames(df))
  }else{
    colnames(df) <- gsub("Sample", "", colnames(df))
    colnames(df) <- gsub("Heavy", "Control", colnames(df))
    colnames(df) <- gsub("Light", "Sample", colnames(df))
  }
  
  
  ##------------data clean--------------
  #select unique peptide
  #1. the unique peptide information were stored in Quan.Info column.
  #2. It can be classfied to "","Redundant","Inconsisitent labeled", "Indistinguishable"，"NoQuan Value" and "Not Unique"
  #ref: page533:  https://assets.thermofisher.com/TFS-Assets/CMD/manuals/Man-XCALI-97808-Proteome-Discoverer-User-ManXCALI97808-EN.pdf
  df <- df %>% filter(Quan.Info=="") 
  
  #delect all NA in four replicates
  df <- df[rowSums(is.na(df))<8,]
  
  #groupby the same peptide together
  df <- df %>% select(Annotated.Sequence,Master.Protein.Accessions,contains("Abundance")) %>% 
    group_by(Annotated.Sequence,Master.Protein.Accessions) %>% 
    summarise(across(everything(),sum),.groups = "drop")
  
  #the peptides found in at least two replicates in sample
  df <- df[rowSums(!is.na(df[,grepl("Sample",colnames(df))]))>=2,]
  

  

  
  ##------------data imputation----------
  #step1, impute control smaples, the hypothesis is that the value in control should be 0, so,a small value is be filled
  #imputate using method 3.1
  #Intnew <- Unif(start,end)
  #start= mean(int, 3xsd(int)), end=mean(int,3xsd(int))
  #The standard deviation was calculated using a 25% trim to mitigate the influence of outliers.
  df[,c(3:ncol(df))] <- log2(df[,c(3:ncol(df))])
    
  control <- grepl("Control",colnames(df))
  sample <- grepl("Sample",colnames(df))
  
  #impute control channel
  con_mean <- apply(df[,control], 2, function(x) mean(x,na.rm=TRUE))
  con_sd <- apply(df[,control], 2, function(x) sd(x,na.rm=TRUE))
  
  #impute 3 or 4 NA in control channel
  for(i in c(1:nrow(df))){
    num_na <- sum(is.na(df[i,control]))
    if(num_na>=3){
      random_values <- runif(num_na, min = con_mean-3*con_sd, max = con_mean-2*con_sd)
      random_values[random_values < 0] <- 0
      df[i,colnames(df[i,control])[is.na(df[i,control])]] <- as.list(random_values)
    }
  }
  
  #impute 1 or 2 NA in control channel
  control_max <- apply(df[,control], 1, max, na.rm=TRUE)
  control_min <- apply(df[,control], 1, min, na.rm=TRUE)
  delta <- (control_max-control_min)/((control_max+control_min)/2)
  mu_delta <- mean(delta)
  sd_delta <- sd(delta)
  cor <- data.frame(df[complete.cases(df), ])
  cor <- cor(cor[,control])
    # fill NA
  for(i in c(1:nrow(df))){
    num_na <- sum(is.na(df[i,control]))
    delta_new <- rnorm(num_na, mean = mu_delta, sd = 2**0.5*sd_delta/mean(cor))
    na_col <- colnames(df)[is.na(df[i,])]
    col <- colnames(df[,3:6])[!is.na(df[i,3:6])]
    #df[i,na_col] <- rowMeans(df[i,col]*abs(1+delta_new))
    df[i,colnames(df[i,control])[is.na(df[i,control])]] <-  as.list(rowMeans(df[i,colnames(df[i,control])[!is.na(df[i,control])]])*abs(1+delta_new))
      
  }
  
  
  
  #impute the sample channel, we have only keeped the peptides that found in at least 2 replicates, so, the impute method is based on the replicate value
  #using method 4.1 in journal NAR 2020
  #delta <- (int1-int2)/mean(int1-int2)
  
  #Build distribution of deltas for all non zero proteins
  sample_max <- apply(df[,sample], 1, max, na.rm=TRUE)
  sample_min <- apply(df[,sample], 1, min, na.rm=TRUE)
  delta <- (sample_max-sample_min)/((sample_max+sample_min)/2)
  mu_delta <- mean(delta)
  sd_delta <- sd(delta)
  cor <- data.frame(df[complete.cases(df), ])
  cor <- cor(cor[,sample])
  # fill NA
  for(i in c(1:nrow(df))){
    num_na <- sum(is.na(df[i,sample]))
    delta_new <- rnorm(num_na, mean = mu_delta, sd = 2**0.5*sd_delta/mean(cor))
    #na_col <- colnames(df)[is.na(df[i,])]
    #col <- colnames(df[,3:6])[!is.na(df[i,3:6])]
    #df[i,na_col] <- rowMeans(df[i,col]*abs(1+delta_new))
    df[i,colnames(df[i,sample])[is.na(df[i,sample])]] <-  as.list(rowMeans(df[i,colnames(df[i,sample])[!is.na(df[i,sample])]])*abs(1+delta_new))
    
  }
  
  ##------------get ratio-----------------
  df[,c(3:ncol(df))] <- 2**(df[,c(3:ncol(df))])
  ratio <- df[,sample]/(df[,sample]+df[,control])
  df$ratio1 <- NA
  df$ratio2 <- NA
  df$ratio3 <- NA
  df$ratio4 <- NA
  
  for (j in c(1:ncol(ratio))){
    df[,(ncol(df)-4+j)] <- ratio[,j]
  }
  
  
  #-------------------------get IDIRT ratio for protein level------------------
  df2 <- df %>% group_by(Master.Protein.Accessions) %>% 
    summarise(ratio1=median(ratio1),
              ratio2=median(ratio2),
              ratio3=median(ratio3),
              ratio4=median(ratio4),
              count=n())
  df2$median_ratio <- apply(df2[, c("ratio1", "ratio2", "ratio3", "ratio4")], 1, median)
  if(exp=="IDIRT"){
    colnames(df2) <- c("Protein","IDRIT_r1","IDIRT_r2","IDIRT_r3","IDIRT_r4","IDIRT_count","IDIRT_median")
  }else{
    colnames(df2) <- c("Protein","swap_r1","swap_r2","swap_r3","swap_r4","swap_count","swap_median")
  }
  return(df2)
}


idirt_merge <- function(Buffer){
  res <- get_ratio(Buffer,"IDIRT")
  res_swap <- get_ratio(Buffer,"IDIRT_swap")
  
  res_all <- full_join(res,res_swap,by="Protein")
  
  res_all <- res_all %>% filter(IDIRT_count>3|swap_count>3)  #important!
  
  
  #write.csv(res_all,file=paste0("./result/",Buffer,"ratio_result.csv"))
  return(res_all)
}


