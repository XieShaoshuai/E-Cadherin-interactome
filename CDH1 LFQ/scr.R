library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)

set.seed(12345678)



#load data
data <- read.table("./data/CDH1_LFQ_Proteins.txt",header = TRUE)
interactor <- read.delim("./data/Summary of interactors.csv", sep=",")
gene <- read.delim("./data/Gene name.csv",sep=",")
newcolname <- read.delim("./data/new name.csv",sep=",")
new <- newcolname$new_name

data <- data %>% select(newcolname$Colname)
colnames(data) <- new

data0 <- dplyr::left_join(data, gene ,by="Accession")




source("./compare.R")

buffer2 <- compare("02",data0)
buffer3 <- compare("3",data0)
buffer16 <-compare("16",data0)
buffer19 <-compare("19",data0)
buffer25 <-compare("25",data0)
buffer29 <-compare("29",data0)







buffer_all <- rbind(buffer2,buffer3,buffer16,buffer19,buffer25,buffer29)

#read data

buffer_all[is.na(buffer_all)]<-0



xx <- buffer_all
xx <- buffer_all %>% filter(Gene %in% interactor$Gene)

xx$label <- ""
xx$label[xx$log2FC>1&xx$adj.p<0.05] <- "+"
xx$label[xx$log2FC< -1&xx$adj.p<0.05] <- "-"

xx$avg_tmtc <- log2(xx$avg_tmtc+1)

xx <- na.omit(xx)

xx$group <- factor(xx$group,c("wt","TMTC2.","TMTC3.","TMTC2_3.","TMTC1_4"))
xx$buffer <- factor(xx$buffer,c("02","3","16","19","25","29"))


xx$log2FC[xx$avg_tmtc==0] <- NA


ggplot(data=xx %>% filter(xx$group!="wt")) + 
  geom_tile(aes(x = group, y =Gene , fill = log2FC),color="grey50") +
  geom_text(aes(x = group, y = Gene, label = label),color="black",
            size = 4)+
  # Colors 
  #scale_fill_gradient(low = "#56B4E9", high = "#D55E00", na.value = "grey50") +
  # scale_fill_gradientn(colors = hcl.colors(20, "Reds", rev = TRUE), 
  #                      na.value = "grey50") +
  scale_fill_gradient2(limits = c(-8,8),
                      high = "red",mid = "white",low = "blue",
                      midpoint = 0,na.value="grey80")  +
  facet_grid(cols=vars(buffer)) +
  # Separate by group0
  # Labels
  labs(x = "", y = "") +
  # Theme
  theme_classic() +
  theme(
    strip.text.y = element_text(angle = 0, size =10),
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

write.csv(buffer_all,file="fold change_20240203.csv")
