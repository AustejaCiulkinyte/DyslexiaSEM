library(data.table)
library(GenomicSEM)

#create vector of clean summary statistics files
files <- c("ADHD_raw.tsv",
           "AN_raw.tsv",
           "ANX_raw.txt",
           "ASD_raw.tsv",
           "BIP_raw.tsv",
           "DLX_raw.tsv",
           "MDD_raw.tsv",
           "OCD_raw.tsv",
           "SCZ_raw.tsv",
           "TS_raw.tsv"
           )

#define the reference file being used to align alleles across summary stats
#here we are using 1000 genomes
ref <- "reference.1000G.maf.0.005.txt"

# name the traits. Note that these names were updated later - in the paper, 
# we replaced "ASD" with "AUT" and "DLX" with "DYX".
trait.names <- c("ADHD", "AN", "ANX", "ASD", "BIP",
                 "DLX", "MDD", "OCD", "SCZ", "TS")

# a  vector indicating whether the standard errors in 
# each set of summary statistics is on the logit scale
# if SE of logistic beta then TRUE
se.logit <- c(T,T,F,T,T,
              T,T,T,T,T)

# a logical vector indicating whether the GWAS was
# for a continuous trait and used OLS (or a LMM)
OLS <- c(F,F,F,F,F,
         F,F,F,F,F)

# a logical vector indicating whether the GWAS is 
# a binary outcome with only Z-statistics or 
# was analyzed using a linear probability model 
# i.e. a dichotomous trait using OLS (or a LMM)
linprob <- c(F,F,T,F,F,
             F,F,F,F,F)

# A vector of total sample sizes for continuous traits and 
# the sum of effective sample sizes for binary traits
# All SNP-specific sums of effective sample sizes have
# been calculated previously or already present, except DLX
N <- c(NA, NA, NA, NA, NA,
       1138870, NA ,NA ,NA, NA)

# vector of column names of betas for continuous traits that 
# are known to have been standardized prior to running the GWAS
# No continuous traits being analysed
betas <- c(NA, NA, NA, NA, NA,
           NA, NA ,NA ,NA, NA)

#run sumstats
sumstats<-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=OLS,linprob=linprob,N=N,betas=betas,
         info.filter = 0.9,maf.filter=0.01, keep.indel=FALSE,parallel=FALSE,
        cores=NULL,ambig=FALSE,direct.filter=FALSE)

save(sumstats, file="sumstats.RData")

##############################################

#extract formatted lists of SNPs, betas and SEs
#read in raw files and add back Neff column
#this is needed for munge and PolarMorphism later. 
#the munge function does not perform all the alignments that sumstats does, 
#therefore we run sumstats first

library(dplyr)

# name output files
output_files <- c("ADHD_clean.tsv",
                  "AN_clean.tsv",
                  "ANX_clean.tsv",
                  "ASD_clean.tsv",
                  "BIP_clean.tsv",
                  "DLX_clean.tsv",
                  "MDD_clean.tsv",
                  "OCD_clean.tsv",
                  "SCZ_clean.tsv",
                  "TS_clean.tsv")

save_files <- function(){
  for(i in 1:length(files)){
    beta_col <- paste("beta.",trait.names[i], sep="")
    se_col <- paste("se.",trait.names[i], sep="")
    df <- subset(sumstats, select=c("SNP", "CHR", "BP", "MAF", "A1", "A2", beta_col, se_col))

    df_neff <- fread(files[i])
    Neff_col_name <- names(df_neff)[grep(glob2rx("Neff*"), 
                                         names(df_neff), 
                                         ignore.case = TRUE)]
    p_col_name <- names(df_neff)[grep(glob2rx("P*"),
                                         names(df_neff),
                                         ignore.case = TRUE)]

    df_neff <- subset(df_neff, select=c("SNP", Neff_col_name, p_col_name))
    
    df <- inner_join(df, df_neff, by="SNP", copy=TRUE)
    names(df)[names(df)==beta_col] <- "BETA"
    names(df)[names(df)==se_col] <- "SE"
    
    write.table(df, file=output_files[i], sep = "\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  }
  return()
}

save_files()
