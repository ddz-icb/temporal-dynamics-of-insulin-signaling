library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
sample_path = "C:\\path_to_your_folder\\"
selected_metadata <- read.csv(paste(c(sample_path, "data_for_graphics_proportions.csv"), sep = "",collapse=""),sep=",",header = TRUE)
long_data <- pivot_longer(data = selected_metadata, 
                     cols = basal_1:insulin_3, 
                     names_to = "condition",
                     values_to = "counts")
#merge_groups
long_data$condition <-gsub('_[[:digit:]]+', '', long_data$condition)
long_data$transcript_id <- factor(long_data$transcript_id)

test1 <- long_data %>% filter(hgnc_symbol == 'RGL4')
p1 <- ggboxplot(test1, x = "transcript_id", y = "counts", color = "black", fill = "condition", 
          palette =c("white", "darkgrey"), 
          width = 0.70,
          xlab = FALSE, ylab = FALSE,
          add.params = list(size = 0.9),
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.15)
ggpar(p1, main = "",legend = 'none',ylim = c(0,0.6))

test1 <- long_data %>% filter(hgnc_symbol == 'SMARCE1')
p1 <- ggboxplot(test1, x = "transcript_id", y = "counts", color = "black", fill = "condition", 
          palette =c("white", "darkgrey"), 
          width = 0.70,
          xlab = FALSE, ylab = FALSE,
          add.params = list(size = 0.9),
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.15)
ggpar(p1, main = "",legend = 'none',ylim = c(0,0.6),xlab = FALSE,ylab = FALSE)

test1 <- long_data %>% filter(hgnc_symbol == 'GMEB1')
p1 <- ggboxplot(test1, x = "transcript_id", y = "counts", color = "black", fill = "condition", 
                palette =c("white", "darkgrey"), 
                width = 0.70,
                xlab = FALSE, ylab = FALSE,
                add.params = list(size = 0.9),
                bxp.errorbar = TRUE, bxp.errorbar.width = 0.15)
ggpar(p1, main = "",legend = 'none',ylim = c(0,0.6))
