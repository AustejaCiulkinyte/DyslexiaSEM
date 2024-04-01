##----------------------------##
# Adapted from from Else Eising's GitHUB, based on script from Kate Doust
#-----------------------------##

# load libraries
library(data.table)
library(dplyr)
library(tidyverse)
library(corrplot)
library(plyr)
library(grid)
library(viridis)

#####################################################################
##################         Data clean-up           ##################
#####################################################################
# Read in data from CTG-VL
# Two analyses were run: batch LDSC of F5 against all publicly available GWAS, and batch LDSC of F5 against
# the 10 GWAS used in this study

data_public <- fread("F5_sumstat_ldsc.csv", header=T, na.strings="\"NA\"")
data_mine <- fread("F5_sumstat_ldsc-2.csv", header=T, na.strings = "\"NA\"")

# Clean up data
data_public <- data_public %>% filter(!if_all(c(rg, se, z, p), is.na))
myPhenotypes <- c("Factor 5", "ADHD_mine", "AN_mine", "ANX_mine",
                  "AUT_mine", "BIP_mine", "DYX_mine", "MDD_mine",
                  "OCD_mine", "SCZ_mine", "TS_mine") # must be in the same order as in the raw dataset
data_mine$phenotype <- myPhenotypes

# Join the two datasets (LDSC with public sumstats and LDSC with sumstats provided manually)
data <- full_join(data_public, data_mine)
fwrite(data, "CTGVL_complete.tsv", col.names=T, row.names=F, quote=F, sep = "\t") # Export all correlations from CTG-VL

# Filter to Bonferroni-significant correlations
p_Bonf <- 0.05/1468 # in total, 1457 comparisons were made using public data and 11 using my own data - including F5 against itself
data_Bonf <- subset(data, data$p < p_Bonf)

fwrite(data_Bonf, "CTGVL_significant.tsv", col.names=T, row.names=F, quote=F, sep = "\t") # Export the table with significant correlations only
# The output was annotated with category data (e.g. "Psychiatric", "Wellbeing" phenotypes etc.)

#####################################################################
##################            Plotting             ##################
#####################################################################

ldsc_sig <- fread("CTGVL_significant_annotated.tsv")

# lock in factor level order to preserve order
ldsc_sig$phenotype <- factor(ldsc_sig$phenotype, levels = ldsc_sig$phenotype)
# Category annotations were added manually
ldsc_sig$category2 <- factor(ldsc_sig$category2, levels = rev(unique(ldsc_sig$category2)))

# Plot trait correlation with F5 and confidence intervals for each trait
ldsc_plot <- ggplot(ldsc_sig, aes(y = phenotype, color = category2)) +
  
  # Add a point range 
  geom_pointrange(size = 1, fatten=2, aes(x = rg, xmin = rg - se, xmax = rg + se), show.legend = FALSE) +
  
  scale_color_manual(palette = viridis) +
  
  scale_y_discrete(labels=ldsc_sig$renamed) +
  
  scale_alpha_discrete(range=c(0.25, 0.8)) +
  
  # fix scaling, limits and breaks   
  scale_x_continuous(expand = c(0, 0), limits = c(-1.1,1.1), breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
  
  # Add a line at zero to divide positive from negative correlations
  geom_vline(xintercept = 0) +
  
  # Change theme to classic (white background and no gridlines)
  theme_classic()+
  
  theme(# Increase thickness of x axis
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    # Remove legend
    legend.position = "none",
    # Remove ticks from y axis
    axis.ticks.y = element_blank()) +
  # specify axis labels
  labs(x="Genetic correlation (rg)", y="") +

# add headings. Adjust y-axis values as needed
  annotate(geom = "text", x = -1, y = 66.5, label = "Psychiatric", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 60, label = "Cognitive", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 53, label = "Education", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 48, label = "Occupation", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 45, label = "Physical health", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 21, label = "Lifestyle", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 13, label = "Wellbeing", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +


  # Add horizontal lines to divide categories. Adjust y-axis values as needed

  geom_segment(aes(x = -1, xend = 1, y = 13.5, yend = 13.5), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 21.5, yend = 21.5), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 45.5, yend = 45.5), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 48.5, yend = 48.5), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 53.5, yend = 53.5), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 60.5, yend = 60.5), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 67, yend = 67), colour = "grey70")


# View plot
ldsc_plot

# save figure as pdf
ggsave("Figure3.pdf", ldsc_plot, bg="white", height = 12, width = 9, scale = 1.2)


