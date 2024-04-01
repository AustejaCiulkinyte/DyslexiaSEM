load("LDSCoutput.RData")

#smooth the S matrix for EFA using the nearPD function in the Matrix package. 
require(Matrix)

Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)

#run EFA with promax rotation and using the factanal function in the stats package
require(stats)

#change number of factors incrementally. We provide EFA output ranging from 1 to 6 factors in a separate file
EFA<-factanal(covmat = Ssmooth, factors = 6, rotation = "promax")

#output results to console
EFA
