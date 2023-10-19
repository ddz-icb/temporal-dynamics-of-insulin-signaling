#install.packages('heatmaply',dependencies = TRUE)
library(heatmaply)
library(RColorBrewer)
library(readxl)
library(dplyr)
#library(scales)

cluster_data_import <- read.xlsx2("C:\\add_your_path_to_example_file\\merged_cluster_data_5.xlsx")
rownames(cluster_data_import) <- cluster_data_import$Gene
cluster_data_import<-cluster_data_import[cluster_data_import$assignment!="complex B specific",] #no specific snrp-complex in figure
#reorder rows by sidebar column
cluster_data_import <- cluster_data_import %>% arrange(assignment)
# rename columns
cluster_data_import <-cluster_data_import %>% rename("1" = "Abundance Ratio: (1) / (0)", "2.5" = "Abundance Ratio: (2.5) / (0)",
                                                     "5" ="Abundance Ratio: (5) / (0)","15"= "Abundance Ratio: (15) / (0)",
                                                     "30"= "Abundance Ratio: (30) / (0)","60"= "Abundance Ratio: (60) / (0)")

# cap outliers
col_val_raw <- c("1","2.5","5","15","30","60")
for (temp_col in col_val_raw) {
  cluster_data_import[temp_col][cluster_data_import[temp_col] > 2.5 ] <- 2.5
}
data <- as.matrix(cluster_data_import %>% select(4:9))
data <- log2(data)
# normalize, center must be 1 log-transformde-> log2(1)=0

maxs <- apply(data, 1, max) # rep(2.5, 30) #
mins <- apply(data, 1, min) #rep(-0.5,30 ) #
my_center <- rep(0,length(data)/6)
#abundance <- t(apply(data, 1, rescale, to=c(-1,1)))
abundance <-t(scale(t(data), center=my_center, scale=(maxs-mins)/2)) # transpose for row wise scaling);center=(maxs+mins)/2

# remove infinite values (no standard deviation, all values in one row are the same)
abundance[which(!is.finite(abundance))] <-  2.61357001

rownames(abundance) <- cluster_data_import$Gene
my_group <- as.numeric(factor(cluster_data_import$assignment,levels = unique(cluster_data_import$assignment)))
#colSide <- brewer.pal(8, "Set3")[my_group] #,colorblindFriendly=TRUE
my_color_pal <- c("#7fc97f","#beaed4","#fdc086","#ffff99","#6a3d9a","#b15928","#bf5b17","#666666")[my_group] #assingn color to sidebar

# with heatmaply
heatmaply(abundance,
          #seriate = "none",
          heatmap_layers = theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")),
          fontsize_row = 11,
          fontsize_col = 16,
          column_text_angle = 0,
          row_text_angle = 0,
          Rowv = NA,
          Colv = NA,
          xlab = "time[min]",
          ylab = "Gene",
          colors = rev(colorRampPalette(brewer.pal(9, "RdBu"))(15)), # reverse color order
          RowSideColors = my_color_pal,
          #row_side_colors = cluster_data_import$assignment,
          #row_side_palette = my_color_pal,
          hide_colorbar = TRUE,
          #subplot_widths = c(0.5, 0.5)
)

