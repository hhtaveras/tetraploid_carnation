library(adegenet)
library(poppr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)


# read in geno data

geno<-read.csv("input/d001.csv", sep = ",", dec = ".", header = T)

# format geno and filter to get genetic distances between tetraploids 
geno0<-geno[,4:ncol(geno)]
rownames(geno0)<-geno$rn
maf <- rowMeans(geno0, na.rm=TRUE) / 4
maf <- ifelse(maf > 0.5, 1-maf, maf)
maf_threshold <- 0.1
keep_snps <- maf >= maf_threshold

geno1 <- geno0[keep_snps,]
geno1<-t(geno1)

# genid format

genid1<-df2genind(geno1, ploidy = 4, ncode = 1)

# calculate distances 
r.dist<- rogers.dist(genid1)

##################################
# Principal coordinates analysis 

d3 <- as.matrix(r.dist)
n <- ncol(d3)
I <- diag(n)
K <- I-1/n
A <- d3 * d3 * -0.5
B <- K %*% A %*% K
E <- eigen(B)
norm <- t(matrix(sqrt(abs(E$values)),ncol=n,nrow=n))
Y <- norm * E$vectors
Y3 <- Y[,1:3]
Y3 <- data.frame(Y3)
Y3$RGS <- rownames(d3)
colnames(Y3) <- c("PCo1","PCo2","PCo3", "car")
# Explained variance
var.exp <- round (100* (E$values / sum (E$values) ),2)


# PCoA plot

theme_set(theme_bw())

plot1 <- ggplot(Y3, aes(x = (PCo1*-1), y = PCo2)) + 
  geom_point(size=2.4, col = "navyblue", shape = 23, bg = "cornflowerblue") +
  xlab(label = paste0("PCo 1 (", var.exp[1], " %)")) +
  ylab(label = paste0("PCo 2 (", var.exp[2], " %)")) +
  coord_fixed() + # square plot
  
  #ggtitle("a") +
  
  theme(panel.border = element_blank(), # remove black border around the figure
        aspect.ratio=0.9,
        text=element_text(size=10,  family="sans"),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=12))

plot1

################################################################################

diag(d3)<-NA

plot2<- Heatmap(d3,
                col = c("orangered", "yellow2"),
                na_col = "orangered3",
                cluster_rows = T,
                cluster_columns = T,
                show_row_dend = T,
                show_column_dend = F,
                row_dend_width = unit(25, "mm"),
                height = unit(115, "mm"),
                name = "MRD",
                heatmap_legend_param = list(border = "black",
                                            at = c(0, 0.029, 0.058, 0.088, 0.12),
                                            labels = c("0.00", "0.02", "0.05", "0.09", "0.12"),
                                            
                                            labels_gp = gpar(col = "black", fontsize = 12)
                                            #ticks_gp = grid::gpar(col = "black", lwd = 2) 
                ),
                row_labels = rep("", nrow(d3)),
                border = F,
                
                column_labels = rep("", ncol(d3)))


plot2
