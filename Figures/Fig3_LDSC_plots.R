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
# Read in data
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

p_Bonf <- 0.05/1468 # in total, 1457 comparisons were made using public data and 11 using my own data
data_Bonf <- subset(data, data$p < p_Bonf)

fwrite(data_Bonf, "CTGVL_significant.tsv", col.names=T, row.names=F, quote=F, sep = "\t") # Export the table with significant correlations only


#####################################################################
##################            Plotting             ##################
#####################################################################

ldsc_sig <- fread("CTGVL_significant.csv")

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

# print pic
ggsave("ldsc_v3.svg", ldsc_plot, bg="white", height = 12, width = 9, scale = 1.2)

























##----------------------------##
# From Else Eising's GitHUB

# bonferroni corrected is calculated
ldsc_sig <- filter(ldsc_raw, p <= (0.05/nrow(ldsc_raw)))


heat_map <- ggplot(data = sig, aes(p1, phenotype, fill = rg)) +
  geom_tile (color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1.1, 1.1), space = "Lab", name = "Genetic correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1)) +
  coord_fixed()


#----------------------------------------------------------------------------------------------------------
# script for the plotting of LDhub results
# By Else Eising, based on script from Kate Doust
# The input file is formatted in Excel to have the appropriate columns and to select the traits to plot
#----------------------------------------------------------------------------------------------------------


# data loaded into ldsc_raw

# add column that specifies whether MTAG or CP has significant correlation
ldsc_raw$Significance_MTAG <- "significant"
ldsc_raw$Significance_MTAG[which(ldsc_raw$p > (0.05/(1404)))] <- "not significant"

#----------------------------------------------------------------------------------------------------------
# Plot data
#----------------------------------------------------------------------------------------------------------

ldsc_raw$Short_name <- factor(ldsc_raw$phenotype, levels=unique(ldsc_raw$phenotype)[order(unique(ldsc_raw$phenotype), decreasing=TRUE)])
ldsc_raw$Short_name2 <- sub("blank[1-9]", "", ldsc_raw$Short_name)

# Plot correlation with dyslexia and confidence intervals for each trait
ldsc_plot <- ggplot(ldsc_raw, aes(y = Short_name)) + 
  
  # Add a point range for MTAG and CP results
  geom_pointrange(size = 1.5, fatten=2, aes(x = rg, xmin = rg - se, xmax = rg + se, alpha=Significance_MTAG), color = viridis(n=3)[2], show.legend = FALSE) +
  
  scale_y_discrete(labels=ldsc_raw$Short_name2[order(ldsc_raw$phenotype, decreasing=TRUE)]) +
  
  scale_alpha_discrete(range=c(0.25, 0.8)) +
  # Add a line at zero to divide positive from negative correlations   
  scale_x_continuous(expand = c(0, 0), limits = c(-1.1,1.1), breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1))+
  
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
  labs(x="Genetic correlation (rg)", y="") 
  # use other column for y-axis labels

# This works - think that the order object is a number which is then deleted (sub)
# Original code has an $Order column which I think is for the group by levels
# LDhub$Short_name <- factor(LDhub$Y_axis_labels, levels=unique(LDhub$Y_axis_labels)[order(unique(LDhub$Order), decreasing=TRUE)])

ldsc_plot  
  

# add headings
  annotate(geom = "text", x = -1, y = 86.5, label = "Cognitive", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 80.5, label = "Education", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 69.5, label = "Eyesight", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 65.5, label = "Chronotype", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 63.5, label = "Lifestyle", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 49.5, label = "Wellbeing", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 38.5, label = "Psychiatric", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 34.5, label = "Pain", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 26.5, label = "Physical health and exercise", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  annotate(geom = "text", x = -1, y = 10.5, label = "SES", hjust = 0, vjust = 1, size = 3, colour="grey31", fontface =2) +
  
  # Add horizontal lines to divide categories
  geom_segment(aes(x = -1, xend = 1, y = 11, yend = 11), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 27, yend = 27), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 35, yend = 35), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 39, yend = 39), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 50, yend = 50), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 64, yend = 64), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 66, yend = 66), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 70, yend = 70), colour = "grey70")+
  geom_segment(aes(x = -1, xend = 1, y = 81, yend = 81), colour = "grey70")+
  geom_segment(aes(x = 0, xend = 0, y = 0.5, yend = 86.5), colour = "grey31")


# View plot
ldsc_plot

# store plot
pdf(paste("Plot_MTAG_and_CP_LDhub_results_subset_viridis.pdf", sep=""), width=9, height=12)   
ldhub_plot
dev.off()

