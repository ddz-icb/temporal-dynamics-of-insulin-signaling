library(ggplot2)


sample_path = "C:\\your_path_to_file\\"
de <- read.csv(paste(c(sample_path, "volcano_example.csv"), sep = "",collapse=""),sep=",",header = TRUE)
de <- data.frame(de)

# set criteria
# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FC > 0.1 & de$adj_pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FC < -0.1 & de$adj_pvalue < 0.05] <- "DOWN"

#same for point size and alpha values
de$criteria <- 0.5
de$criteria[de$log2FC_betrag > 0.1 & de$adj_pvalue < 0.05]<- 0.7

#gene labels
library(ggrepel)
#gene_info <- read.csv(paste(c(sample_path, "all_DTU_genes_anno.tsv"), sep = "",collapse=""),sep="\t",header = TRUE)
#merged_df <- merge(de, gene_info, by= "gene_id", all.x = TRUE)
merged_df <- de
merged_df$delabel <- NA
merged_df$delabel[merged_df$diffexpressed != "NO"] <- merged_df$hgnc_symbol[merged_df$hgnc_symbol != "NO"]
p <- ggplot(data=merged_df, aes(x=log2FC, y=-log10(adj_pvalue), col=diffexpressed, label=delabel)) + 
  geom_point(alpha = 1/10,aes(size=criteria)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text_repel(box.padding = 0.3,max.overlaps =17)+
  theme(legend.position="none")+
  xlab("log2Foldchange") + ylab("-log10(adjusted_pvalue)")
p2 <- p + geom_vline(xintercept=c(-0.1,0.1), col="grey") +
  geom_hline(yintercept=-log10(0.05), col="grey")+
  xlim(-3, 4)
# set colors
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)+ scale_size_continuous(range = c(2, 3))
p3
png(paste(c(sample_path, "Volcano_2_sided_transcripts.png"), sep = "",collapse=""), 
    width=3600,
    height=3600,
    res=600)
p3
dev.off()
