# load libraries
library(PolarMorphism)
library(whitening)
library(data.table)

# These libraries are not necessary for PolarMorphism, but are used to make plots later
library(readr)
library(tidyverse)
library(BiocManager)
library(qvalue)

# define folder where sumstats are kept
wdir <- "./"
traits <- c("ADHD", "DLX")

ADHD <- fread(paste0(wdir, "ADHD", "_clean.tsv"))
# change column names
colnames(ADHD) <- c("snpid", "chr", "bp", "freq", "a1","a2","beta","se", "n","pval") # note that PolarMorphism does not need or use the "n" column

DLX <- fread(paste0(wdir, "DLX", "_clean.tsv"))
# change column names
colnames(DLX) <- c("snpid", "chr", "bp", "freq", "a1","a2","beta","se", "pval","pass") 
# note that PolarMorphism does not need or use the "n" column


# We need to choose one of the GWAS as reference, to make sure all GWAS's 
# have the same reference and alternative allele for each SNP
# We will make DLX the reference, and 'flip' the alleles of ADHD so they align
ADHD <- AlleleFlip(sumstats = ADHD, snps = DLX %>% select(snpid, a1, a2), 
                   snpid = "snpid", only.a2 = F)

# because the function Alleleflip not only flips the alleles, but also adds 
# a z-score column, we have to manually do that for DLX
DLX$z <- DLX$beta/DLX$se

# run PolarMorphism
DLX.ADHD <- ConvertToPolar(dfnames = traits, snpid = "snpid", whiten = T, LDcorrect = F)

# p-value & q-value for r
DLX.ADHD$r.pval <- PvalueForR(r = as.numeric(DLX.ADHD$r), p = 2)
DLX.ADHD$r.qval <- qvalue(p = DLX.ADHD$r.pval)$qvalues

# filter on r q-value
PolarMorphies <- DLX.ADHD[DLX.ADHD$r.qval < 0.05,]

# get p- and q-values for theta
PolarMorphies$theta.pval <- PvalueForAngle(angle.trans = as.numeric(PolarMorphies$angle), 
                                           r = as.numeric(PolarMorphies$r))
PolarMorphies$theta.qval <- qvalue(p = PolarMorphies$theta.pval)$qvalues

# filter on theta q-value
PolarMorphies.q <- PolarMorphies[PolarMorphies$theta.qval < 0.05,]

#calculate lambda values for QQ plot
DLX.ADHD$chisq <- qchisq(DLX.ADHD$r.pval,1,lower.tail=FALSE)
case <- 19099+51800 # sum of cases in ADHD and DLX GWASes
control <- 34194+1087070 # sum of controls in ADHD and DLX GWASes
lambda <- median(DLX.ADHD$chisq) / qchisq(0.5,1)
lambda1000 <- 1 + (lambda - 1) * (1 / case + 1 / control) * 500

# print lambda values to console
lambda
lambda1000

# export a table where all SNPs have an r and theta value calculated, as well as
# p- and q- values for both
write.table(DLX.ADHD, "polarmorph.tsv", row.names=F, quote=F, sep="\t")

# export SNPs filtered to r.qval <0.05 and theta.qval <0.05
write.table(PolarMorphies.q, "polarmorph_filtered_theta_r.tsv", row.names=F, quote=F, sep="\t")

## The filtered output was compressed and uploaded to FUMA for gene mapping
