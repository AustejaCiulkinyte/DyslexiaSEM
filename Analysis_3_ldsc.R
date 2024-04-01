library(GenomicSEM)

#vector of munged summary statisitcs
traits<-c("ADHD.sumstats.gz", 
          "AN.sumstats.gz", 
          "ANX.sumstats.gz",
          "ASD.sumstats.gz",
          "BIP.sumstats.gz",
          "DLX.sumstats.gz",
          "MDD.sumstats.gz",
          "OCD.sumstats.gz",
          "SCZ.sumstats.gz",
          "TS.sumstats.gz"
)

#enter sample prevalence of .5 to reflect that all traits except DLX
#were munged using the sum of effective sample size
sample.prev<-c(.5,.5,.5,.5,.5,
               .045,.5,.5,.5,.5)

#vector of population prevalences
population.prev<-c(0.05,  #ADHD
                   0.01,  #AN
                   0.1,   #ASD
                   0.012, #ASD
                   0.02,  #BIP
                   0.1,   #DLX
                   0.302, #MDD
                   0.025, #OCD
                   0.01,  #SCZ
                   0.008  #TS
)

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("ADHD","AN","ANX","ASD", "BIP", "DLX", "MDD", "OCD", "SCZ", "TS")

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names,stand=TRUE)

#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="LDSCoutput.RData")
