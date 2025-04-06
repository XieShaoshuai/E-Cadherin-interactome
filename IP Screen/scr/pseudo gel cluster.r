# Install and load required packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(stats)
library(cluster)
library(tibble)
library(tidyr)
# 在文件开头添加这两个包
library(ggdendro)
library(cowplot)

#load data
data <- read.delim("./Data/IP screening/CHD1 IP screening_Proteins.txt", stringsAsFactors = FALSE) %>% 
  filter(
    Master == "IsMasterProtein",
    Protein.FDR.Confidence.Combined == "High",
    Number.of.Unique.Peptides > 1,
    Contaminant == "False"
  ) %>%
  dplyr::select(
    Accession, MW.in.kDa, starts_with("Abundance.F")
  )


#select external gene name and abundance columns
df <- data %>%
  dplyr::select(Accession, MW.in.kDa, starts_with("Abundance.F"))

# plot the persu gel
#the speed a protein moved in gel is linear to its ln(MW)
#the regress linear was calculated manunually by hand
#distance = -17.6 x ln(MW)+98.954

#convert MW to distance
df$MW.in.kDa <- -17.6 * log(df$MW.in.kDa) + 98.954

#replace NA with 0, and transform the abundance to log2(abundance+1) scale
df <- df %>%
  replace(is.na(.), 0) %>%
  mutate(across(starts_with("Abundance.F"), ~ log2(. + 1)))


#check the values
#median(df$Abundance.F1.Sample,na.rm = TRUE)



colnames(df) <- c("Gene_Name", "MW.in.kDa", c(1:32))

#combine all values to one columns for plotting
df2 <- df  %>% 
  gather(key = "group", value = "value", -Gene_Name, -MW.in.kDa)

#cutoff the intensity
cutoff <- quantile(df2$value,0.9)
max <- round(max(df2$value))
df2$value[df2$value < cutoff] <- 0


df2$group <- as.numeric(df2$group)
df2$group <- df2$group * 2

# 定义特殊基因的颜色映射
#special_color_mapping <- c(
#  "P12830" = "red",
#  "A0A2R8YCH5" = "blue"
#)

# 为特殊基因设置颜色
df2$color <- NA

# 为其他基因根据 value 值设置颜色
df2$color <- ifelse(is.na(df2$color),
                    scales::gradient_n_pal(c("white", "black"))(scales::rescale(df2$value, from = c(cutoff, max))), 
                    "white")

df2$color <- ifelse(df2$Gene_Name == "P12830",
                    scales::gradient_n_pal(c("#ff7676", "#fd0000c6"))(scales::rescale(df2$value, from = c(cutoff, max))), 
                    df2$color)

df2$color <- ifelse(df2$Gene_Name == "A0A2R8YCH5",
                    scales::gradient_n_pal(c("#BA8DFFFF", "#6600FFFF"))(scales::rescale(df2$value, from = c(cutoff, max))), 
                    df2$color)

# 对group进行聚类分析
# 将数据转换为宽格式用于聚类
df_wide <- df2 %>% 
  dplyr::select(Gene_Name, group, value) %>% 
  pivot_wider(names_from = group, values_from = value) %>% 
  column_to_rownames("Gene_Name")

# 计算距离矩阵并进行层次聚类
dist_matrix <- dist(t(df_wide), method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")

# 获取聚类顺序
cluster_order <- hc$order

# 重新映射group编号但保持x轴位置均匀分布
original_groups <- unique(df2$group)
new_group_order <- original_groups[cluster_order]
df2$group <- factor(df2$group, levels = new_group_order)

# 添加新的x_position列用于均匀分布
df2$x_position <- as.numeric(factor(df2$group, levels = new_group_order))

# 修改绘图部分
p_main <- ggplot(df2, aes(x = x_position, y = MW.in.kDa, color = I(color))) +
  geom_segment(aes(x = x_position - 0.3, xend = x_position + 0.3, 
                   y = MW.in.kDa, yend = MW.in.kDa), size = 1) +
  labs(title = "",
       x = "Extrsctant",
       y = "MW (kDa)") +
  scale_y_continuous(
    breaks = c(-17.6 * log(15) + 98.954, -17.6 * log(25) + 98.954, -17.6 * log(50) + 98.954, 
               -17.6 * log(100) + 98.954, -17.6 * log(125) + 98.954, -17.6 * log(150) + 98.954),
    labels = c("15", "25", "50", "100", "125", "150"),
    limits = c(-17.6 * log(150) + 98.954, -17.6 * log(10) + 98.954)
  ) +
  scale_x_continuous(
    breaks = seq_along(unique(df2$group)),
    labels = as.numeric(levels(df2$group)) / 2
  ) +
  theme_bw()+
   theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

# 创建水平方向的树状图
dendro_data <- dendro_data(hc)
p_dendro <- ggplot(segment(dendro_data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_continuous(limits = c(0.5, max(df2$x_position)+0.5)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# 调整主图的边距
p_main <- p_main + 
  theme(plot.margin = margin(0, 0, 0, 0))

# 组合两个图（垂直排列）
combined_plot <- plot_grid(p_dendro, p_main, 
                          align = "v", 
                          axis = "l",
                          rel_heights = c(0.2, 0.8),
                          ncol = 1)

# 输出组合图
plot(combined_plot)
