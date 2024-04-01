## 1. Data availability
This repository relates to Ciulkinyte et al. (2023) "Genetic neurodevelopmental clustering and dyslexia", preprint available at https://doi.org/10.1101/2023.10.04.23296530.

All data are publicly available, with the exception of anxiety and dyslexia GWAS. Links to public repositories are below:

- ADHD: https://figshare.com/articles/dataset/adhd2019/14671965
- Anorexia nervosa: https://figshare.com/articles/dataset/an2019/14671980
- Autism: https://figshare.com/articles/dataset/asd2019/14671989
- Bipolar disorder: https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594
- MDD: https://datashare.ed.ac.uk/handle/10283/3203
- OCD: https://figshare.com/articles/dataset/ocd2018/14672103
- Schizophrenia: https://figshare.com/articles/dataset/scz2022/19426775
- Tourette syndrome: https://figshare.com/articles/dataset/ts2019/14672232

Notes on which file to download where multiple are available:
- Bipolar disorder: pgc-bip2021-all.vcf.tsv.gz
- MDD: “Genome-wide summary statistics from a meta-analysis of PGC and UK Biobank (347.4Mb)” (PGC_UKB_depression_genome-wide.txt)
- Schizophrenia: PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz

Anxiety GWAS was obtained through communication with authors of Romero et al. (2022) https://doi.org/10.1038/s41588-022-01245-2. The anxiety summary statistics included in this study is a meta-analysis of two previously published studies: Purves et al. (2020) https://doi.org/10.1038/s41380-019-0559-1 and Levey et al. (2020) https://doi.org/10.1176/appi.ajp.2019.19030256. We have obtained permission from the authors of the meta-analysis and of individual analyses to use these summary statistics in our study.

Dyslexia GWAS summary statistics are available through 23andMe website (https://research.23andme.com/dataset-access/) to qualified researchers under an agreement with 23andMe that protects the privacy of the 23andMe participants.

## 2. Data analysis

### 2.1. Effective sample size calculation
To quote from https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics, effective sample size (**Neff**) should be obtained from, in order of priority:
1.	Available in-file
2.	Calculated from cohort-level info
3.	Calculated from in-file MAF
4.	Calculated from reference MAF

In this study, we use the following sources to calculate Neff:

| Trait | Neff source |
| --- | --- | 
| ADHD	| Calculated from cohort-level info | 
| AN	| Available in-file |
| ANX	| Available in-file |
| ASD	| Calculated from reference MAF |
| BIP	| Available in-file |
| DLX	| Not a meta-analysis, will use total sample size and prevalence |
| MDD	| Calculated from in-file MAF |
| OCD	| Calculated from cohort-level info |
| SCZ	| Available in-file |
| TS	| Calculated from cohort-level info |

Please refer to https://github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-Sample-Size-and-Preparing-GWAS-Summary-Statistics for a tutorial on how to calculate Neff.

### 2.2. Quality control
Before running gSEM, we filter to only use SNPs shared across datasets where MAF > 0.05 and INFO > 0.09 using _sumstats_ and _munge_ functions from the GenomicSEM R package. Please refer to https://github.com/GenomicSEM/GenomicSEM/wiki/ for an in-depth explanation on how to run these and other gSEM functions. In this repository, we provide the following code relevant to our study:

 ` Analysis_1_sumstats.R `

 ` Analysis_2_munge.R `

In addition, we provide the log files from these two functions:

 ` ADHD_AN_ANX_ASD_BIP_DLX_MDD_OCD_SCZ_TS_sumstats.txt ` 
 
 ` ADHD_AN_ANX_ASD_BIP_DLX_MDD_OCD_SCZ_TS_munge.txt ` 

### 2.3. LDSC
Following QC, we ran ld-score regression using the _ldsc_ function in the GenomicSEM R package. In this repository, we provide the code:

` Analysis_3_ldsc.R ` 

the log file

` ADHD.sumstats.gz_AN.sumstats.gz_ANX.sumstats.gz_ASD.sumstats.gz_BIP.sumstats.gz_DLX.sumstats.gz_MDD._ldsc.txt `

and the output

` LDSCoutput.RData `

relevant to this study.

Using the log file, we did the following:
1.	Checked that for Heritability Results for each individual trait, Intercept is close to 1; Ratio is close to 0; h2 Z is >4
2.	Made a matrix containing Genetic Correlation Results (observed scale h2, rg values, standard errors and g_cov P-values) for each pair of traits (Supplementary Table 2)

### 2.4. Exploratory factor analysis

See ` Analysis_4_EFA.R ` .

A complete output is provided as ` EFA_output.xlsx `.

### 2.5. GenomicSEM

See ` Analysis_5_gSEM.R ` .

### 2.6. userGWAS

This function was used to extract summary statistics describing Factor 5 (Attention and learning difficulties latent factor).

See ` Analysis_6_userGWAS.R ` .

### 2.7. PolarMorphism

See ` Analysis_6_polarMorphism.R ` .

## 3. Figures

We provide code used to generate each figure, where relevant, in a separate directory.

## 4. Supplementary tables

We provide supplementary tables as an .xlsx file for ease of use.
