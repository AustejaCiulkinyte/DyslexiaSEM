# load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(BioVenn)
library(BiocManager)


# read files
ADHD.genes <- fread("ADHD_signif_genes.csv") # Genes associated with ADHD - from ADHD GWAS paper
DLX.genes <- fread("DLX_signif_genes.csv") # Genes associated with DLX - from DLX GWAS paper
double.genes <- fread("genes.txt") # genes associated with both ADHD and DLX - from FUMA (Supp. Table 5)

# clean up files
#colnames(ADHD.genes)[5] <- "symbol" # rename column
colnames(DLX.genes)[9] <- "symbol" # rename column

# function to find number of genes in df1 that are present in df2 and return their names
commonGenes = function(df1, df2){
  
  common <- inner_join(df1, df2, by="symbol")
  count <- nrow(common)
  geneNames <- common$symbol
  
  return(list(count, geneNames))
  
}

# find genes that are associated with both ADHD and DLX
commonGenes(ADHD.genes, DLX.genes)


# find genes that are associated with both ADHD and the ADHD/DLX subset
commonGenes(ADHD.genes, double.genes)


# find genes that are associated with both DLX and the ADHD/DLX subset
commonGenes(DLX.genes, double.genes)


# check if any of the genes shared between ADHD/double and DLX/double are the same
ADHD_shared_genes <- commonGenes(ADHD.genes, double.genes)[2][[1]]
DLX_shared_genes <- commonGenes(DLX.genes, double.genes)[2][[1]]
double_shared_genes <- intersect(ADHD_shared_genes, DLX_shared_genes)
print(double_shared_genes)


# draw a Venn diagram and save as .pdf
draw.venn(ADHD.genes$symbol, DLX.genes$symbol, double.genes$symbol,
          title="", subtitle="", 
          xtitle="ADHD mapped genes", xt_f="sans", xt_fb=1, xt_s=3,
          ytitle="DLX mapped genes", yt_f="sans", yt_fb=1, yt_s=3,
          ztitle="Pleiotropic mapped genes", zt_f="sans", zt_fb=1, zt_s=3,
          nr_f="sans", nr_fb=1, nr_s=3,
          x_c="yellow", y_c="cyan", z_c="magenta",
          output="pdf", filename="gene_Venn.pdf")
