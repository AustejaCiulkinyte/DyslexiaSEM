# load packages
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(scales)

## load files for plotting

# neff values for each dataset
neff <- read.csv("Neff.csv", header=TRUE)
  #    Name       Neff
  #1  ADHD 103135.500
  #2    AN  46321.904
  #3   ANX 248239.000
  #4   AUT  41457.380
  #5   BIP 101574.300
  #6   DYX 113887.000
  #7   MDD 430423.800
  #8   OCD   7281.307
  #9   SCZ  58746.070
  #10   TS  12140.028

# full correlation matrix from Supplementary Table 2
corr <- read.csv("correlations_full.csv", header=TRUE, row.names=1)

# a dataframe containing names of traits in all possible pairs, but otherwise empty (45 pairs for 10 traits)
df <- read.csv("pairs.csv", header=TRUE)

# convert variables
corr <- as.matrix(corr)
t(corr)
corr <- as.vector(t(corr))

# remove correlation values that are equal to 1 (only those where a trait is compared
# to itself)
corr[corr==1] <- NA
corr <- na.omit(corr)

# add pairwise correlation values to the empty dataframe with all possible pairs
df$corr <- corr

# get neff for the first member of each pair
Neff1 <- merge(df, neff, by.x="d1", by.y="Name", sort=FALSE)$Neff
df$Neff1 <- Neff1

# get neff for the second member of each pair
df <- merge(df, neff, by.x="d2", by.y="Name", sort=FALSE)

# clean up the dataframe
df <- df[,c("d1", "d2", "corr", "Neff1", "Neff")]
df <- df[order(df$d1),]
names(df)[names(df) == 'Neff'] <- 'Neff2'

# calculate pairwise Neff
df$sqrtNeff <- sqrt(df$Neff1 * df$Neff2)

# perform a regression of pairwise Neff on genetic correlation
mod1 = lm(corr~sqrtNeff, data=df)

# get statistics
modsum = summary(mod1)
modsum

# define outliers, i.e., which points should be labelled
# here, points where rg is >0.3 or <0.1, or where pairwise effective Neff is
# >100,000 will be labelled
df_outliers <- df %>% filter((corr > 0.3 | corr < 0.1 | sqrtNeff > 100000))

# save the plot as .pdf
pdf(file="rg_vs_neff.pdf", width=5, height=5)

df%>%
  ggplot(aes(x=sqrtNeff, y=corr)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, col="black") +
  geom_label_repel(data=df_outliers,
                   aes(x=sqrtNeff, y=corr, label=paste(d1, d2, sep="/")), 
                   size=5, box.padding = 0.15, label.padding=0.15, max.overlaps=5) +
  theme(plot.title=element_text(size=72,face="bold"),
                axis.text=element_text(size=72),
                axis.title=element_text(size=72)) +
  theme_classic()+
  xlab(expression(Effective~sample~size~group("(",sqrt(Neff[1]*~"\u00D7"~Neff[2]), ")"))) +
  ylab("Genetic correlation") +
  scale_x_continuous(labels=comma)

dev.off()
