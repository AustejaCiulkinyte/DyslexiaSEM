library(GenomicSEM)
library(data.table)


#vector of munged summary statisitcs
traits<-c("sumstatsF5.txt")

#enter sample prevalence of .5 to reflect that all traits except DLX
#were munged using the sum of effective sample size
sample.prev<-NA

#vector of population prevalences
population.prev<-NA

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("F5")

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names,stand=TRUE)

#optional command to save the output as a .RData file for later use
save(LDSCoutput,file="LDSCoutputF5.RData")


