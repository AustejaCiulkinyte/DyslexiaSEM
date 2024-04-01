library(GenomicSEM)

#create vector of the summary statistics files
files<-c("ADHD_clean.tsv",
         "AN_clean.tsv",
         "ANX_clean.tsv",
         "ASD_clean.tsv",
         "BIP_clean.tsv",
         "DLX_clean.tsv",
         "MDD_clean.tsv",
         "OCD_clean.tsv",
         "SCZ_clean.tsv",
         "TS_clean.tsv"
)

#define the reference file being used to align alleles across summary stats
#here we are using 1000 genomes
ref<-"reference.1000G.maf.0.005.txt"

#name the traits 
trait.names<-c("ADHD","AN","ANX", "ASD", "BIP", "DLX", "MDD", "OCD", "SCZ", "TS")

#list the sample sizes. All SNP-specific sums of effective sample sizes have
#been calculated previously or already present, except DLX (and ANX because transfer
#of N after sumstats has not worked there)
N=c(NA, NA, 248239, NA, NA,
    1138870, NA ,NA ,NA, NA)

#definte the imputation quality filter
info.filter=0.9

#define the MAF filter
maf.filter=0.05

#run munge
munge(files=files,hm3=ref,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)
