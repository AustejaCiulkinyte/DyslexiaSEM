# This code was run as multiple independent scripts, but I combined them into one R file on Github for conciseness

#################################################
## Step 1: split sumstats into chunks
## Otherwise, time and memory consumption when running userGWAS becomes unreasonable

## Load required packages
library(lspcheminf)

## Load LDSC output and munged summary statistics
load("sumstats.RData") # output from Analysis_1_sumstats.R

## Split sumstats into 5 chunks
sumstats_chunks <- chunk_df(sumstats, 5, seed = NULL)

SumstatsChunk1 <- sumstats_chunks[[1]]
save(SumstatsChunk1, file = "SumstatsChunk1.RData")

SumstatsChunk2 <- sumstats_chunks[[2]]
save(SumstatsChunk2, file = "SumstatsChunk2.RData")

SumstatsChunk3 <- sumstats_chunks[[3]]
save(SumstatsChunk3, file = "SumstatsChunk3.RData")

SumstatsChunk4 <- sumstats_chunks[[4]]
save(SumstatsChunk4, file = "SumstatsChunk4.RData")

SumstatsChunk5 <- sumstats_chunks[[5]]
save(SumstatsChunk5, file = "SumstatsChunk5.RData")

#################################################
## Step 2: run userGWAS
## This code should be run separately for each of the 5 sumstats chunks. Only one example for chunk 1 is shown
## 

## Load required packages
library(GenomicSEM)

## Load LDSC output and chunk of summary statistics
load("SumstatsChunk1.RData")
sumstats <- SumstatsChunk1

load("LDSCoutput.RData")

## Define model syntax
FiveModel <- 'F1 =~ OCD + AN + TS
F2 =~ BIP + SCZ
F3 =~ ANX + MDD
F4 =~ ADHD + ASD
F5 =~ ADHD + DLX
F1 ~~ F3
F1 ~~ F4
F1 ~~ F5
F2 ~~ F3
F2 ~~ F4
F2 ~~ F5
F3 ~~ F4
F3 ~~ F5
F4 ~~ F5
F5 ~ SNP # this line tells userGWAS to compute correlations between F5 and each SNP
'

# run userGWAS. The "sub" variable tells userGWAS we are only interested in the factor-SNP correlations
CorrelatedFactors <- userGWAS(covstruc = LDSCoutput, SNPs = sumstats, estimation = "DWLS", 
                                model = FiveModel, sub=c("F5~SNP"), toler = 1e-30, cores = 7) 

save(CorrelatedFactors, file = "CorrelatedFactors1.RData")

#################################################
## Step 3: combine 5 userGWAS outputs into one file

library(dplyr)
library(data.table)

# load in each chunk of sumstats
load("CorrelatedFactors1.RData")
CorrFact1 <- CorrelatedFactors[1]
rm(CorrelatedFactors)

load("CorrelatedFactors2.RData")
CorrFact2 <- CorrelatedFactors[1]
rm(CorrelatedFactors)

load("CorrelatedFactors3.RData")
CorrFact3 <- CorrelatedFactors[1]
rm(CorrelatedFactors)

load("CorrelatedFactors4.RData")
CorrFact4 <- CorrelatedFactors[1]
rm(CorrelatedFactors)

load("CorrelatedFactors5.RData")
CorrFact5 <- CorrelatedFactors[1]
rm(CorrelatedFactors)

# convert to dataframes
CorrFact1df <- CorrFact1[[1]]
CorrFact2df <- CorrFact2[[1]]
CorrFact3df <- CorrFact3[[1]]
CorrFact4df <- CorrFact4[[1]]
CorrFact5df <- CorrFact5[[1]]

# define vector of dataframes
df_list <- list(CorrFact1df, CorrFact2df, CorrFact3df, CorrFact4df, CorrFact5df)

# bind all dataframes together
CorrFactFull <- do.call(rbind, df_list)

save(CorrFactFull, file="F5sumstats.RData")
write.table(CorrFactFull, file="F5sumstats.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

# Check how many SNPs have a warning
noWarn<-which(CorrFactFull$warning=="0")
length(CorrFactFull$warning)-length(noWarn)
# 121207 SNPs did not converge

# Subset to exclude SNPs that did not converge
CorrFactClean <- subset(CorrFactFull, warning=="0")
save(CorrFactClean, file="F5sumstatsclean.RData")
write.table(CorrFactClean, file="F5sumstatsclean.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

##Calculate Implied Sample Size for Factor 5
#restrict to MAF of 40% and 10%
CorrFactforN<-subset(CorrFactClean, CorrFactClean$MAF <= .4 & CorrFactClean$MAF >= .1)

N_hat_F5<-mean(1/((2*CorrFactforN$MAF*(1-CorrFactforN$MAF))*CorrFactforN$SE^2))
# N_hat_F5 = 254601

## Clean F5 sumstats were compressed and uploaded to CTG-VL for batch LDSC.
