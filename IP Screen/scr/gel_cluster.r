# 加载必要的库
library(imager)
library(ggplot2)
library(cluster)
library(tidyr)

#读取文件
file <- c(1:32)

profiles <- as.data.frame(matrix(NA, nrow = 236, ncol = 32)) # 创建一个空的数据框来存储所有的profile
images <- list()
for(i in file){
    image_path <- paste0("./lanes/lanes rescale size/", i, ".jpg")
    gel_image <- load.image(image_path)


    #check the imgae size, make sure the size for all images are the same
    image_height <- dim(gel_image)[2]
    image_width <- dim(gel_image)[1]
    cat(image_width, image_height, "\n")
    images[[i]] <- gel_image

    #plot(gel_image, axes = FALSE, main = paste("Gel Image", i))

    # 假设gel_image是矩阵形式，如果不是需要先转换
    if(!is.matrix(gel_image)) {
       gel_matrix <- matrix(gel_image, 
                     nrow = image_height,  # 使用之前定义的图像高度
                     ncol = image_width,   # 使用之前定义的图像宽度
                     byrow = TRUE)
    } else {
    gel_matrix <- gel_image
    }

# 计算每行5-25列的平均值
profile <- data.frame(rowMeans(gel_matrix[, 5:18]))
colnames(profile) <- i  # 1表示按行计算
profiles[,i] <- profile  # 将计算得到的profile添加到profiles数据框中
}

# 层次聚类分析
dist_matrix <- dist(t(profiles), method = "euclidean")  # 计算样本间距离
hc <- hclust(dist_matrix, method = "ward.D2")  # 使用ward方法进行层次聚类

# 绘制聚类图并保存为PDF
pdf("cluster_dendrogram.pdf", width = 10, height = 6)
plot(hc, main = "Gel Image Clustering Dendrogram", 
     xlab = "Sample", sub = "", hang = -1)
dev.off()

