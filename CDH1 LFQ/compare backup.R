compare <- function(buffer,data){
  #get intensity in wt
  data <- data %>% select(Accession,Gene, contains("Abundance."))
  
  #select wt samples in buffer
  d <- data %>% select(Accession,Gene,contains(paste0(".",buffer))) %>% 
    select(Accession,Gene, contains("wt")) %>% 
    mutate(avg_tmtc = apply(select(., contains("wt")), 1, mean, na.rm = TRUE)) %>% 
    select(Accession,avg_tmtc) %>% 
    mutate(log2FC=0,adj.p=0,buffer=buffer,group="wt")
  
  d <- left_join(interactor,d, by="Accession")
  d <- d %>% select(-Interactor)
  
  #d$avg_tmtc <- d$avg_tmtc*(1988099660/d$avg_tmtc[d$Gene=="CDH1"])
  
  for(i in c("TMTC2.","TMTC3.","TMTC2_3.","TMTC1_4")){
    data2 <- data %>% select(Accession,Gene,contains(paste0(".",buffer)))
    data2 <- data2 %>% select(Accession,Gene,contains(c("wt",i)))
    
    df <- data2 %>% select(Accession,contains(c("wt",i)))
    
    #calculate ratio by CDH1
    #CDH1 <- df[df$Accession=="P12830",2:ncol(df)]/100
    #df[,2:ncol(df)] <- as.data.frame(mapply(`/`, df[,2:ncol(df)], CDH1))
    
    #Remove row that NA in all samples
    df <- df[rowSums(is.na(df))< ncol(df)-1,]
    
    #NA imputation
    df_filled <- df[,-1]
    df_filled <- log2(df_filled)
    replace_na <- function(col) {
      mu <- mean(col, na.rm = TRUE)
      sigma <- sd(col, na.rm = TRUE)
      min_val <- mu - sigma
      max_val <- mu + sigma
      ifelse(is.na(col), rnorm(sum(is.na(col)), mu-1.8*sigma, 0.3*sigma), col)
    }
    df_filled <- data.frame(lapply(df_filled, replace_na))
    
    
    df_filled$Accession <- df$Accession
    
    
    #normalziatzion using median 
    
    
    CDH1 <- as.matrix(df_filled[df_filled$Accession=="P12830",-ncol(df_filled)])
    df_filled[,-ncol(df_filled)] <- sweep(df_filled[,-ncol(df_filled)], MARGIN = 2, CDH1, `-`)
    
    #remove NA
    df_filled <- na.omit(df_filled)  
    
    #calculate logFC and p.value
    for(k in c(1:nrow(df_filled))){
      group1 <- as.numeric(df_filled[k,] %>% select(contains("wt")))
      group2 <- as.numeric(df_filled[k,] %>% select(contains("TMTC")))
      if(df_filled$Accession[k]!="P12830"){
        t_test <- t.test(group1,group2,paired = FALSE)
        df_filled$p.value[k] <- t_test$p.value
        df_filled$log2FC[k] <- mean(group2)-mean(group1)
      }else{
        df_filled$p.value[k] <- 1
        df_filled$log2FC[k] <-0
      }
      
    }
    
    #adjust pvalue
    adjusted_p_values <- p.adjust(df_filled$p.value, method = "BH")
    df_filled$adj.p <- adjusted_p_values
    
    #add gene name information
    df_filled <- left_join(df_filled,gene, by="Accession")
    #add interactor infomration
    df_filled <- left_join(df_filled,interactor,by="Accession")
    
    
    filename <- paste0("./output/plot/",i,"_Buffer_",buffer,"_all.pdf")
    p<- ggplot(df_filled,aes(x=log2FC,y=-log10(adj.p)))+
      geom_point()+
      ggtitle(paste0(i,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj')+
      theme_bw()
    
    ggsave(p,filename=filename,height=8,width =8,dpi=600,units="in")
    
    
    #volcano plot
    plot_data <- df_filled %>% filter(Interactor=="Stable"|Interactor=="Bait"|Interactor=="Unstable")
    plot_data$color <- "No sig"
    plot_data$color[plot_data$log2FC>1&plot_data$adj.p<0.05] <- "up in TMTC KO"
    plot_data$color[plot_data$log2FC< -1&plot_data$adj.p<0.05] <- "down in TMTC KO"
    
    plot_data$color <- factor(plot_data$color,c("No sig","up in TMTC KO","down in TMTC KO"))
    filename <- paste0("./output/plot/",i,"_Buffer_",buffer,".pdf")
    
    p <- ggplot(data=plot_data,aes(x=log2FC,y=-log10(adj.p),color=color))+
      geom_point(aes(alpha=0.8,size=0.4))+
      scale_color_manual(values=c( "grey","red","blue"))+
      geom_text_repel(data=plot_data %>% filter(abs(log2FC) > 1,adj.p <= 0.05),
                      aes(x=log2FC,y=-log10(adj.p),label=Gene.x))+
      ggtitle(paste0(i,"  Buffer:",buffer))+
      xlab('log2(FC)')+
      ylab('-log10(p.adj)')+
      theme_bw()+
      theme(legend.position="none",axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12))
    
    
    ggsave(p,filename=filename,height=8,width =8,dpi=600,units="in")
    
    filename <- paste0("./output/data/",i,"_Buffer_",buffer,".csv")
    write.table(plot_data,file=filename,row.names = FALSE,sep=",")
    
    
    df_temp <- plot_data %>% select(Accession,Gene.x,log2FC,adj.p)
    df_temp$buffer <- buffer
    df_temp$group <- i
    
    
    
    data_temp <- data2 %>% mutate(avg_tmtc = apply(select(., contains("TMTC")), 1, mean, na.rm = TRUE))%>% 
      select(Accession,avg_tmtc)
    
    df_temp <- left_join(df_temp,data_temp,by="Accession")
    x <- interactor[,1:2]
    colnames(x) <- c("Accession","Gene.x")
    df_temp <- left_join(x,df_temp,by=c("Accession","Gene.x"))
    colnames(df_temp) <- c("Accession","Gene","log2FC","adj.p","buffer","group","avg_tmtc")
    df_temp$buffer <- buffer
    df_temp$group <- i
    df_temp$avg_tmtc <- df_temp$avg_tmtc*(1988099660/df_temp$avg_tmtc[df_temp$Gene=="CDH1"])
    
    d <- rbind(d,df_temp)
  }
  
  return(d)
  
}