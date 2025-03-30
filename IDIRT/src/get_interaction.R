library(ggplot2)
library(dplyr)
library(BayesFactor)
library(ggrepel)


#Bayesian test
Bayesian_test <- function(res_all,buffer,mix) {
  mix_ratio <- mix[["mix"]]
  swap_ratio <- mix[["mix_swap"]]
  
  #sd(mix_ratio)
  #median(mix_ratio)
  
  df <- res_all %>% select(Protein,contains("_r"))
  
  #replace NA with a random value based mix distribution
  replace_na_with_normal <- function(x, mu, sd) {
    is_na <- is.na(x)
    x[is_na] <- rnorm(sum(is_na), mean = mu, sd = sd)
    return(x)
  }
  
  df[, 2:5] <- apply(df[, 2:5], 2, function(x) replace_na_with_normal(x, mean(mix_ratio), sd(mix_ratio)))
  df[, 6:9] <- apply(df[, 6:9], 2, function(x) replace_na_with_normal(x, mean(swap_ratio), sd(swap_ratio)))
  
  
  df$idirt_BF <- 0
  df$swap_BF <-0
  #get bayesian factor
  for (i in c(1:nrow(df))){
    ratio <- na.omit(as.numeric(df[i,2:5]))
    dd <- data.frame(
      value = c(ratio, mix_ratio),
      group = factor(c(rep(1, length(ratio)), rep(2, length(mix_ratio)))))
    bf_result <- ttestBF(formula = value ~ group, data = dd)
    df$idirt_BF[i] <-bf_result@bayesFactor$bf
    
    ratio <- na.omit(as.numeric(df[i,6:9]))
    dd <- data.frame(
      value = c(ratio, swap_ratio),
      group = factor(c(rep(1, length(ratio)), rep(2, length(swap_ratio)))))
    bf_result <- ttestBF(formula = value ~ group, data = dd)
    df$swap_BF[i] <-bf_result@bayesFactor$bf
  }
  
  
  df$IDIRT_median <- apply(df[,c("IDRIT_r1","IDIRT_r2","IDIRT_r3","IDIRT_r4")],1,mean)
  df$swap_median <- apply(df[,c("swap_r1","swap_r2","swap_r3","swap_r4")],1,mean)
  
  
  
  a2 <- 1-(median(mix_ratio)+sd(mix_ratio))    #important!
  b2 <- 1-(median(swap_ratio)+sd(swap_ratio))
  
  df$group <- (df$IDIRT_median-1)**2/a2**2+ (df$swap_median-1)**2/b2**2
  
  
  
  
  
  #check whether the protein is specific interactor
  #filter 1: BF>3 in idirt or idirt_swap
  #filter 2: ratio> mix_ratio-sd, swap_ratio> swap_ratio
  filter1 <- median(mix_ratio)-sd(mix_ratio)
  filter2 <- median(swap_ratio)-sd(swap_ratio)
  #bollen <- (df$idirt_BF>3|df$swap_BF>3)&(df$IDIRT_median>filter1&df$swap_median>filter2)
  
  #df$interactor <- bollen
  
  
  
  #convert protein accession to gene name
  gene <- read.delim("./IDIRT data/Gene name.csv",sep=",")
  
  df <- left_join(df,gene,by=c("Protein"="Accession"))
  
  # 创建一个新列color，根据条件对点进行着色
  df$color <- ifelse(df$idirt_BF > 3 & df$swap_BF > 3 & df$group<1, "Interactor in both",
                     ifelse(df$idirt_BF > 3 & df$swap_BF < 3 & df$group<1, "Interactor in IDIRT",
                            ifelse(df$idirt_BF < 3 & df$swap_BF > 3 & df$group<1, "Interactor in Swap", 
                                   "Non interactor")))
  
  
  
  
  t <- seq(0, 2 * pi, by = 0.01)
  h <- 1
  k <- 1
  
  
  
  x2 <- h + a2 * cos(t)
  y2 <- k + b2* sin(t)
  
  # generat data
  ellipse_data2 <- data.frame(x2, y2)
  
  
  my_colors <- c("Interactor in both" = "#00D300", 
                 "Interactor in IDIRT" = "#89b5f1", 
                 "Interactor in Swap" = "#FF5C53", 
                 "Non interactor" = "grey70")
  # 绘制散点图
  p <- ggplot() +
    geom_point(data=df, aes(x=IDIRT_median,y=swap_median, color = color),size = 4,alpha=0.75) +
    scale_color_manual(values =my_colors) +
    geom_path(data=ellipse_data2, aes(x = x2, y = y2),linetype = "dashed",color="grey20",size=1)+
    geom_text_repel(data = df %>% filter(color!="Non interactor"), 
                    aes(x=IDIRT_median,y=swap_median,label = Gene,color="grey20"))+
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
    coord_fixed()+
    theme_bw()
  
  p <- p+geom_path(data=ellipse_data2, aes(x = x2, y = y2),linetype = "dashed",color="grey20",size=1)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(p,file=paste0("./result/",buffer,".pdf"),
         height =12,
         units = c("in"),
         dpi = 900,)
  write.csv(df,file=paste0("./result/",buffer,"_IDRIT_result.csv"))
  return(df)
}
