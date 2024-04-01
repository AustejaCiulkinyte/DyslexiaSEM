# load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(VennDiagram)

# read files
ADHD.loci <- fread("ADHD_signif_SNPs.csv") # Loci associated with ADHD - from ADHD GWAS paper
DLX.loci <- fread("DLX_signif_SNPs.csv") # Loci associated with DLX - from DLX GWAS paper
double.loci <- fread("GenomicRiskLoci.txt") # loci associated with both ADHD and DLX - from FUMA

### clean up files
DLX.loci$CHR <- as.integer(str_extract(DLX.loci$Cytoband, "[0-9]+")) # get chromosome number
DLX.loci <- separate(data = DLX.loci, col = 'Start-stop positions of credible set', 
                     into = c("start", "end"), sep = ":") # get start-end positions of loci
DLX.loci$start <- as.integer(DLX.loci$start)
DLX.loci$end <- as.integer(DLX.loci$end)
colnames(ADHD.loci)[2:4] <- c("CHR", "BP", "SNP") # rename columns
colnames(DLX.loci)[2] <- "BP" # rename column
colnames(double.loci)[3:5] <- c("SNP", "CHR", "BP") # rename columns

# define function to find number of loci in df1 that are within start-end positions of df2
countInRange = function(df1, df2){
  count=0
  for(i in 1:nrow(df2)){
     for(j in 1:nrow(df1)){
        if(inrange(df1$BP[j], lower=df2$start[i], upper=df2$end[i]) && df1$CHR[j] == df2$CHR[i]){
          count=count+1
        }
     }
  }
  return(count)
}

# find number of loci in ADHD that are within start-end positions of double
countInRange(ADHD.loci, double.loci)
# 4

# find number of loci in DLX that are within start-end positions of double
countInRange(DLX.loci, double.loci)
# 15

# find number of loci in ADHD that are within start-end positions of DLX
countInRange(ADHD.loci, DLX.loci)
# 1


# a function to return names of SNPs from df1 that are within start-end positions of
# loci in df2

returnInRange = function(df1, df2){
  
  matched <- data.frame(CHR=character(), df1_SNP=character(), df2_Risk_locus_start=character(), 
                        df2_Risk_locus_end=character(), stringsAsFactors = FALSE)
  
  for(i in 1:nrow(df2)){
    for(j in 1:nrow(df1)){
      if(inrange(df1$BP[j], lower=df2$start[i], upper=df2$end[i]) && df1$CHR[j] == df2$CHR[i]){
        matched[nrow(matched)+1,] <- c(df1$CHR[j], df1$SNP[j], df2$start[i], df2$end[i])
      }
    }
  }
  return(matched)
}


ADHD.double.loci <- returnInRange(ADHD.loci, double.loci)
#   CHR    df1_SNP df2_Risk_locus_start df2_Risk_locus_end
# 1   3  rs7613360             48719638           50399695
# 2  10 rs11596214            106451870          106832251
# 3  11  rs2582895             28591168           28710118
# 4  20  rs6082363             21154234           21549630

DLX.double.loci <- returnInRange(DLX.loci, double.loci)

  # CHR    df1_SNP df2_Risk_locus_start df2_Risk_locus_end
# 1   2 rs72916919            198146381          198954774
# 2   3  rs2624839             48719638           50399695
# 3   3 rs10511073             85396592           85958954
# 4   3 rs13082684            135619585          136752590
# 5   7  rs3735260             68914449           69896122
# 6   9  rs3122702             15529852           16053793

ADHD.DLX.loci <- returnInRange(ADHD.loci, DLX.loci)
#   CHR   df1_SNP df2_Risk_locus_start df2_Risk_locus_end
# 1   3 rs7613360             49735746           50209053




###############################################################################
## The following code formats lists of loci (from papers and derived by PolarMorphism)
## such that they can be read by the draw.venn command later.
## If this code is run, you MUST reread the original loci files back into R
## before running any further analysis.

## In short, the SNPIDs are changed such that they are all aligned to SNPIDs as
## reported in the ADHD GWAS. This ensures that draw.venn can correctly detect 
## how many loci are shared between each list (ADHD, dyslexia and pleotropic risk
## loci). However, these SNPIDs will no longer reflect those which have been originally
## assigned to risk loci identified in dyslexia GWAS or the pleiotropic analysis.

double.loci$match <- runif(49)

for(i in 1:length(double.loci$SNP)){
  for(j in 1:length(DLX.double.loci$df1_SNP)){

    if(double.loci$start[i] == DLX.double.loci$df2_Risk_locus_start[j] &
     double.loci$end[i] %in% DLX.double.loci$df2_Risk_locus_end[j]){
      double.loci$match[i] <- DLX.double.loci$df1_SNP[j]
  }
  }
}

for(i in 1:length(double.loci$SNP)){
  for(j in 1:length(ADHD.double.loci$df1_SNP)){
    
    if(double.loci$start[i] == ADHD.double.loci$df2_Risk_locus_start[j] &
       double.loci$end[i] %in% ADHD.double.loci$df2_Risk_locus_end[j]){
      double.loci$match[i] <- ADHD.double.loci$df1_SNP[j]
    }
  }
}

for(i in 1:length(DLX.loci$SNP)){
  for(j in 1:length(ADHD.DLX.loci$df1_SNP)){
    
    if(DLX.loci$start[i] == ADHD.DLX.loci$df2_Risk_locus_start[j] &
       DLX.loci$end[i] %in% ADHD.DLX.loci$df2_Risk_locus_end[j]){
      DLX.loci$SNP[i] <- ADHD.DLX.loci$df1_SNP[j]
    }
  }
}

###############################################################################


# draw a Venn diagram and save as .pdf
draw.venn(ADHD.loci$SNP, DLX.loci$SNP, double.loci$match,
          title="", subtitle="", 
          xtitle="ADHD risk loci", xt_f="sans", xt_fb=1, xt_s=3,
          ytitle="DYX risk loci", yt_f="sans", yt_fb=1, yt_s=3,
          ztitle="Pleiotropic risk loci", zt_f="sans", zt_fb=1, zt_s=3,
          nr_f="sans", nr_fb=1, nr_s=3,
          x_c="yellow", y_c="cyan", z_c="magenta",
          output="pdf", filename="loci_Venn.pdf")

