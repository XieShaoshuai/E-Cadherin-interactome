compare <- function(buffer,data){
  #get intensity in wt
  data <- data %>% filter(Number.of.Unique.Peptides>2) %>% 
    select(Accession,Gene, contains("Abundance."))
  
  #select wt samples in buffer
  d <- data %>% select(Accession,Gene,contains(paste0(".",buffer))) %>% 
    select(Accession,Gene, contains("wt")) %>% 
    mutate(avg_tmtc = apply(select(., contains("wt")), 1, mean, na.rm = TRUE)) %>% 
    select(Accession,avg_tmtc) %>% 
    mutate(log2FC=0,adj.p=0,buffer=buffer,group="wt")
  
  #d <- left_join(interactor,d, by="Accession")
  #d <- d %>% select(-Interactor)
  
  #d$avg_tmtc <- d$avg_tmtc*(1988099660/d$avg_tmtc[d$Gene=="CDH1"])
  
  for(m in c("TMTC2.","TMTC3.","TMTC2_3.","TMTC1_4")){
    data2 <- data %>% select(Accession,Gene,contains(paste0(".",buffer)))
    data2 <- data2 %>% select(Accession,Gene,contains(c("wt",m)))
    
    df <- data2 %>% select(Accession,contains(c("wt",m)))
    
    #remove low confidence value
    control <- grepl("Control",colnames(df))
    sample <- grepl("Sample",colnames(df))
    df <- df[(rowSums(is.na(df[,control]))<=2 | rowSums(!is.na(df[,sample]))>=2),]
    
    
    #NA imputation
    df_filled <- df[,-1]
    df_filled <- log2(df_filled)
    
    #for all NA in control or sample
    #table(rowSums(is.na(df_filled[,1:4])))
    #table(rowSums(is.na(df_filled[,5:8])))
    for(i in c(1:nrow(df_filled))){
      mu <- colMeans(df_filled[,1:4],na.rm=TRUE)
      sd <- apply(df_filled[,1:4], 2, function(x) sd(x, na.rm = TRUE))
      if(sum(is.na(df_filled[i,1:4]))==4){
        for(j in c(1:4)){
          #df_filled[i,j] <- rnorm(1, mu[j]-1.8*sd[j], 0.3*sd[j])
          df_filled[i,j] <- max(0,runif(1, min=mu[j]-3*sd[j], max=mu[j]-2*sd[j]))
        }
        
      }
    }
    
    for(i in c(1:nrow(df_filled))){
      mu <- colMeans(df_filled[,5:ncol(df_filled)],na.rm=TRUE)
      sd <- apply(df_filled[,5:ncol(df_filled)], 2, function(x) sd(x, na.rm = TRUE))
      if(sum(is.na(df_filled[i,5:ncol(df_filled)]))==(ncol(df_filled)-4)){
        for(j in c(1:(ncol(df_filled)-4))){
          #df_filled[i,j+4] <- rnorm(1, mu[j]-1.8*sd[j], 0.3*sd[j])
          df_filled[i,j+4] <- max(0,runif(1, min=mu[j]-3*sd[j], max=mu[j]-2*sd[j]))
        }
        
      }
    }
    
    #impute control-----
    #using method 4.1 in journal NAR 2020
    #delta <- (int1-int2)/mean(int1-int2)
    
    #Build distribution of deltas for partical non zero proteins
    control <- grepl("Control",colnames(df_filled))
    control_max <- apply(df_filled[,control], 1, max, na.rm=TRUE)
    control_min <- apply(df_filled[,control], 1, min, na.rm=TRUE)
    
    delta <- (control_max-control_min)/((control_max+control_min)/2)
    mu_delta <- mean(delta)
    sd_delta <- sd(delta)
    
    cor <- data.frame(df_filled[complete.cases(df_filled), ])
    cor <- cor(cor[,control])
    
    
    # fill NA
    for(i in c(1:nrow(df_filled[,1:4]))){
      num_na <- sum(is.na(df_filled[i,control]))
      if(num_na>0){
        delta_new <- rnorm(num_na, mean = mu_delta, sd = 2**0.5*sd_delta/mean(cor))
        na_col <- colnames(df_filled[,1:4])[is.na(df_filled[i,1:4])]
        col <- colnames(df_filled[,1:4])[!is.na(df_filled[i,1:4])]
        df_filled[i,na_col] <- mean(as.numeric(df_filled[i,col]),na.rm=TRUE)*abs(1+delta_new)
        
      }
    }
    #impute sample-----
    #using method 4.1 in journal NAR 2020
    #delta <- (int1-int2)/mean(int1-int2)
    
    #Build distribution of deltas for all non zero proteins
    sample <- grepl("Sample",colnames(df_filled))
    sample_max <- apply(df_filled[,sample], 1, max, na.rm=TRUE)
    sample_min <- apply(df_filled[,sample], 1, min, na.rm=TRUE)
    
    delta <- (sample_max-sample_min)/((sample_max+sample_min)/2)
    mu_delta <- mean(delta)
    sd_delta <- sd(delta)
    
    cor <- data.frame(df_filled[complete.cases(df_filled), ])
    cor <- cor(cor[,sample])
    
    
    # fill NA
    for(i in c(1:nrow(df_filled[,5:ncol(df_filled)]))){
      num_na <- sum(is.na(df_filled[i,sample]))
      if(num_na>0){
        delta_new <- rnorm(num_na, mean = mu_delta, sd = 2**0.5*sd_delta/mean(cor))
        na_col <- colnames(df_filled[,5:ncol(df_filled)])[is.na(df_filled[i,5:ncol(df_filled)])]
        col <- colnames(df_filled[,5:ncol(df_filled)])[!is.na(df_filled[i,5:ncol(df_filled)])]
        df_filled[i,na_col] <- mean(as.numeric(df_filled[i,col]),na.rm=TRUE)*abs(1+delta_new)
        
      }
    }
    
    df_filled$Accession <- df$Accession
    
    
    
    
    df_filled0 <- df_filled
    df_filled <- df_filled0
    
    #normalziatzion using the bait intensity
    df_filled_nor_bait <- df_filled
    CDH1 <- as.matrix(df_filled[df_filled$Accession=="P12830",1:ncol(df_filled)-1])
    df_filled_nor_bait[,-ncol(df_filled)] <- sweep(df_filled[,-ncol(df_filled)], MARGIN = 2, CDH1, `-`)

    #normalization using the median intensity
    df_filled_nor_median <- df_filled
    median <- apply(df_filled_nor_median[,1:(ncol(df_filled)-1)],2,median)
    df_filled_nor_median[,-ncol(df_filled)] <- sweep(df_filled[,-ncol(df_filled)], MARGIN = 2, median, `-`)
    
    #calculate logFC and p.value df_filed----
    for(k in c(1:nrow(df_filled))){
      group1 <- as.numeric(df_filled[k,] %>% select(contains("wt")))
      group2 <- as.numeric(df_filled[k,] %>% select(contains("TMTC")))
      t_test <- t.test(group1,group2,paired = FALSE)
      df_filled$p.value[k] <- t_test$p.value
      df_filled$log2FC[k] <- mean(group2)-mean(group1)
    }
    
    #adjust pvalue
    adjusted_p_values <- p.adjust(df_filled$p.value, method = "BH")
    df_filled$adj.p <- adjusted_p_values
    
    #add gene name information
    df_filled <- left_join(df_filled,gene, by="Accession")
    #add interactor infomration
    df_filled <- left_join(df_filled,interactor,by=c("Gene"="Stable.interaction"))
    
    #calculate logFC and p.value of df_filled_nor_median----
    for(k in c(1:nrow(df_filled_nor_median))){
      group1 <- as.numeric(df_filled_nor_median[k,] %>% select(contains("wt")))
      group2 <- as.numeric(df_filled_nor_median[k,] %>% select(contains("TMTC")))
      t_test <- t.test(group1,group2,paired = FALSE)
      df_filled_nor_median$p.value[k] <- t_test$p.value
      df_filled_nor_median$log2FC[k] <- mean(group2)-mean(group1)
    }
    
    #adjust pvalue
    adjusted_p_values <- p.adjust(df_filled_nor_median$p.value, method = "BH")
    df_filled_nor_median$adj.p <- adjusted_p_values
    
    #add gene name information
    df_filled_nor_median <- left_join(df_filled_nor_median,gene, by="Accession")
    #add interactor infomration
    df_filled_nor_median <- left_join(df_filled_nor_median,interactor,by=c("Gene"="Stable.interaction"))
    
    #calculate logFC and p.value for df_filled_nor_bait------
    for(k in c(1:nrow(df_filled_nor_bait))){
      group1 <- as.numeric(df_filled_nor_bait[k,] %>% select(contains("wt")))
      group2 <- as.numeric(df_filled_nor_bait[k,] %>% select(contains("TMTC")))
      t_test <- t.test(group1,group2,paired = FALSE)
      df_filled_nor_bait$p.value[k] <- t_test$p.value
      df_filled_nor_bait$log2FC[k] <- mean(group2)-mean(group1)
    }
    
    #adjust pvalue
    adjusted_p_values <- p.adjust(df_filled_nor_bait$p.value, method = "BH")
    df_filled_nor_bait$adj.p <- adjusted_p_values
    
    #add gene name information
    df_filled_nor_bait <- left_join(df_filled_nor_bait,gene, by="Accession")
    #add interactor infomration
    df_filled_nor_bait <- left_join(df_filled_nor_bait,interactor,by=c("Gene"="Stable.interaction"))
    
    
    filename1 <- paste0("./output/plot/all/",m,"_Buffer_",buffer,"_all_noNorm.pdf")
    filename2 <- paste0("./output/plot/all/",m,"_Buffer_",buffer,"_all_Norm_Bait.pdf")
    filename3 <- paste0("./output/plot/all/",m,"_Buffer_",buffer,"_all_Norm_Median.pdf")
    p1<- ggplot(df_filled,aes(x=log2FC,y=-log10(adj.p)))+
      geom_point(color="grey35")+
      ggtitle(paste0(m,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj')+
      theme_bw()
    p2<- ggplot(df_filled_nor_bait,aes(x=log2FC,y=-log10(adj.p)))+
      geom_point(color="grey35")+
      ggtitle(paste0(m,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj')+
      theme_bw()
    p3<- ggplot(df_filled_nor_median,aes(x=log2FC,y=-log10(adj.p)))+
      geom_point(color="grey35")+
      ggtitle(paste0(m,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj')+
      theme_bw()
    
    filename <- paste0("./output/plot/all/",m,"_Buffer_",buffer,"_all_BaitNorm.csv")
    write.table(df_filled_nor_bait,file=filename,row.names = FALSE,sep=",")
    data00 <- df_filled_nor_bait
    
    
    ggsave(p1,filename=filename1,height=8,width =8,dpi=600,units="in")
    ggsave(p2,filename=filename2,height=8,width =8,dpi=600,units="in")
    ggsave(p3,filename=filename3,height=8,width =8,dpi=600,units="in")
    
    
    #volcano plot for no normalizatiom----
    plot_data <- df_filled %>% filter(Class=="Stable"|Class=="Bait"|Class=="Dynamic")
    plot_data$color <- "No sig"
    plot_data$color[plot_data$log2FC>1&plot_data$adj.p<0.05] <- "up in TMTC KO"
    plot_data$color[plot_data$log2FC< -1&plot_data$adj.p<0.05] <- "down in TMTC KO"
    
    plot_data$color <- factor(plot_data$color,c("No sig","up in TMTC KO","down in TMTC KO"))
    filename <- paste0("./output/plot/interactors/",m,"_Buffer_",buffer,"_no_Norm.pdf")
    color_mapping <- c("No sig" = "grey", "up in TMTC KO" = "#fc6800", "down in TMTC KO" = "#0068c4")
    p <- ggplot(data=plot_data,aes(x=log2FC,y=-log10(adj.p),color=color))+
      geom_point(aes(alpha=0.8,size=0.4))+
      scale_color_manual(values=color_mapping)+
      geom_text_repel(data=plot_data %>% filter(abs(log2FC) > 1,adj.p <= 0.05),
                      aes(x=log2FC,y=-log10(adj.p),label=Gene))+
      ggtitle(paste0(m,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj)')+
      theme_bw()+
      theme(legend.position="none",axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))
    
    
    ggsave(p,filename=filename,height=8,width =8,dpi=600,units="in")
    
    filename <- paste0("./output/data/",m,"_Buffer_",buffer,".csv")
    write.table(plot_data,file=filename,row.names = FALSE,sep=",")
    
    
    #volcano plot for normalizatiom with bait----
    plot_data <- df_filled_nor_bait %>% filter(Class=="Stable"|Class=="Bait"|Class=="Dynamic")
    plot_data$color <- "No sig"
    plot_data$color[plot_data$log2FC>1&plot_data$adj.p<0.05] <- "up in TMTC KO"
    plot_data$color[plot_data$log2FC< -1&plot_data$adj.p<0.05] <- "down in TMTC KO"
    
    plot_data$color <- factor(plot_data$color,c("No sig","up in TMTC KO","down in TMTC KO"))
    
    filename <- paste0("./output/plot/",m,"_Buffer_",buffer,"_Norm_bait.pdf")
    color_mapping <- c("No sig" = "grey", "up in TMTC KO" = "#fc6800", "down in TMTC KO" = "#0068c4")
    p <- ggplot(data=plot_data,aes(x=log2FC,y=-log10(adj.p),color=color))+
      geom_point(aes(alpha=0.8,size=0.4))+
      scale_color_manual(values=color_mapping)+
      geom_text_repel(data=plot_data %>% filter(abs(log2FC) > 1,adj.p <= 0.05),
                      aes(x=log2FC,y=-log10(adj.p),label=Gene))+
      ggtitle(paste0(m,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj)')+
      theme_bw()+
      theme(legend.position="none",axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))
    
    
    ggsave(p,filename=filename,height=8,width =8,dpi=600,units="in")
    
    filename <- paste0("./output/data/",m,"_Buffer_",buffer,"_Norm_bait.csv")
    write.table(plot_data,file=filename,row.names = FALSE,sep=",")
    
    
    #volcano plot for normalizatiom with median----
    plot_data <- df_filled_nor_median %>% filter(Class=="Stable"|Class=="Bait"|Class=="Dynamic")
    plot_data$color <- "No sig"
    plot_data$color[plot_data$log2FC>1&plot_data$adj.p<0.05] <- "up in TMTC KO"
    plot_data$color[plot_data$log2FC< -1&plot_data$adj.p<0.05] <- "down in TMTC KO"
    
    plot_data$color <- factor(plot_data$color,c("No sig","up in TMTC KO","down in TMTC KO"))
    filename <- paste0("./output/plot/",m,"_Buffer_",buffer,"_Norm_median.pdf")
    color_mapping <- c("No sig" = "grey", "up in TMTC KO" = "#fc6800", "down in TMTC KO" = "#0068c4")
    p <- ggplot(data=plot_data,aes(x=log2FC,y=-log10(adj.p),color=color))+
      geom_point(aes(alpha=0.8,size=0.4))+
      scale_color_manual(values=color_mapping)+
      geom_text_repel(data=plot_data %>% filter(abs(log2FC) > 1,adj.p <= 0.05),
                      aes(x=log2FC,y=-log10(adj.p),label=Gene))+
      ggtitle(paste0(m,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj)')+
      theme_bw()+
      theme(legend.position="none",axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))
    
    
    ggsave(p,filename=filename,height=8,width =8,dpi=600,units="in")
    
    filename <- paste0("./output/data/",m,"_Buffer_",buffer,"Norm_median.csv")
    write.table(plot_data,file=filename,row.names = FALSE,sep=",")
    
    #return(data00)

  }
  
  #return(d)
  
}
