# load pacakges
library(GenomicSEM)
library(lavaan)
library(semPlot)
library(semptools)

#get variance/covariance matrix
load("LDSCoutput.RData")

#change abbreviations from ASD & DLX to AUT & DYX
attributes(LDSCoutput[["S"]])$dimnames[[2]][4] <- "AUT"
attributes(LDSCoutput[["S_Stand"]])$dimnames[[2]][4] <- "AUT"
attributes(LDSCoutput[["S"]])$dimnames[[2]][6] <- "DYX"
attributes(LDSCoutput[["S_Stand"]])$dimnames[[2]][6] <- "DYX"


#defines the function to convert GenomicSEM output into an output readable by 
#semPlot(). NOT directly compatible between common factor & user-defined models.

semPlotModel_GSEM=function(gsem.object=GWISoutput, 
                           unstd.measure.label="Unstand_Est",
                           unstd.SE.label="Unstand_SE",
                           std.measure.label="STD_Genotype",
                           std.SE.label="STD_Genotype_SE"){ 
  
  object=gsem.object$results
  object$free=0
  numb=1:length(which(object$op!="~~"))
  object$free[which(object$op!="~~")]=numb
  varNames <- lavaanNames(object, type = "ov")
  factNames <- lavaanNames(object, type = "lv")
  factNames <- factNames[!factNames %in% varNames]
  n <- length(varNames)
  k <- length(factNames)
  # if (is.null(object$label)) 
  #   object$label <- rep("", nrow(object))
  semModel <- new("semPlotModel")
  
  for (i in 1:length(gsem.object$results[,std.measure.label])){
    if(gsem.object$results[,std.SE.label][i]==""){
      object$label[i] <- round(as.numeric(object[,std.measure.label][i]),2)
    }
    else
      object$label[i] <- paste(as.character(round(as.numeric(object[,std.measure.label][i]),2)),
                               paste("(",as.character(round(as.numeric(object[,std.SE.label][i]),2)),")", sep=""),
                               sep="\n")
  }
  
  object$est <- object[,unstd.measure.label]
  object$std <- object[,std.measure.label]
  if (is.null(object$group)) 
    object$group <- ""
  semModel@Pars <- data.frame(label = object$label, lhs = ifelse(object$op == 
                                                                   "~" | object$op == "~1", object$rhs, object$lhs), edge = "--", 
                              rhs = ifelse(object$op == "~" | object$op == "~1", object$lhs, 
                                           object$rhs), est = object$est, std = object$std, group = object$group, 
                              fixed = object$free==0, par = object$free, stringsAsFactors = FALSE)
  semModel@Pars$edge[object$op == "~~"] <- "<->"
  semModel@Pars$edge[object$op == "~*~"] <- "<->"
  semModel@Pars$edge[object$op == "~"] <- "~>"
  semModel@Pars$edge[object$op == "=~"] <- "->"
  semModel@Pars$edge[object$op == "~1"] <- "int"
  semModel@Pars$edge[grepl("\\|", object$op)] <- "|"
  semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge), 
                                       -(3:4)]
  semModel@Pars <- semModel@Pars[!object$op %in% c(":=", "<", 
                                                   ">", "==", "|", "<", ">"), ]
  semModel@Vars <- data.frame(name = c(varNames, factNames), 
                              manifest = c(varNames, factNames) %in% varNames, exogenous = NA, 
                              stringsAsFactors = FALSE)
  semModel@ObsCovs <- list()
  semModel@ImpCovs <- list()
  semModel@Computed <- FALSE
  semModel@Original <- list(object)
  return(semModel)
  
}


# Syntax defining a three correlated factor model

TFModel <- 'F1~~F2
F2~~F3
F1 =~ OCD + AN + TS
F2 =~ BIP + SCZ + MDD + ANX
F3 =~ ADHD + AUT + DYX + MDD + ANX + TS
'
# Create a structural model
ThreeFactor<-usermodel(LDSCoutput, estimation = "DWLS", model = TFModel, std.lv=TRUE)

#  get model fit statistics
ThreeFactor$modelfit
#       chisq df      p_chisq      AIC       CFI      SRMR
# df 329.3954 29 1.178691e-52 381.3954 0.9333533 0.0770417

# remove an aberrant correlation line
results<-ThreeFactor$results
noF1F3<-subset(results, !(results$lhs=="F1" & results$op=="~~" & results$rhs=="F3"))
ThreeFactor$results<-noF1F3

# convert gSEM output to a format readable by the plotting function
formatted_model <- semPlotModel_GSEM(ThreeFactor)

#run semPaths to draw a diagram and save it as .pdf
semPaths(formatted_model,
         layout="tree2",
         whatLabels="label",
         style="ram",
         edge.label.cex=0.9,
         sizeMan=8,
         sizeLat=10,
         nCharNodes=4,
         nCharEdges=0,
         label.norm="OOOO",
         fixedStyle = c("black",1),
         freeStyle = c("black",1),
         width=15,
         height=8,
         mar=c(8,1,8,1),
         edge.label.position= c(rep(0.5, 14), 0.3, rep(0.5, 15)),
         filetype="pdf",
         filename="path_diagram_3F_initial.pdf"
)


# Model suggested by EFA: syntax
EFAModel <- 'F1 =~ OCD + AN + TS
F2 =~ BIP + SCZ
F3 =~ ANX + MDD
F4 =~ ADHD + AUT + DYX
F1 ~~ F3
F1 ~~ F4
F2 ~~ F3
F2 ~~ F4
F3 ~~ F4
'

# create a structural model
CFAofEFA<-usermodel(LDSCoutput, estimation = "DWLS", model = EFAModel, std.lv=TRUE)

# get model fit statistics
CFAofEFA$modelfit
#       chisq df      p_chisq      AIC       CFI       SRMR
# df 172.5952 29 2.326762e-22 224.5952 0.9681415 0.06195346

#remove aberrant correlation lines
results<-CFAofEFA$results
results$formula <- paste(results$lhs,results$op, results$rhs)
unwanted <- c("F1 ~~ F2")
results_clean<-subset(results, !(results$formula %in% unwanted))
CFAofEFA$results<-results_clean

#convert gSEM output to a format readable by semPaths
formatted_model <- semPlotModel_GSEM(CFAofEFA)

g<-10

#save the model to a variable
semPaths_plot<-semPaths(formatted_model,
         layout="tree2",
         whatLabels="label",
         style="ram",
         curvature = 5,
         edge.label.cex=0.8,
         reorder=TRUE,
         sizeMan=8,
         sizeLat=10,
         nCharNodes=4,
         nCharEdges=0,
         label.norm="OOOO",
         optimizeLatRes = TRUE,
         optimPoints = 0.4,
         fixedStyle = c("black",1),
         freeStyle = c("black",1),
         width=45,
         height=25,
         levels=c(1,2,4,5),
         mar=c(10,2,15,1)
         #edge.label.position=c(rep(0.5,9), 0.3, rep(0.5,26))
         
)

# save the model as a .pdf. Saving the model as a variable first allows me to
# apply the rotate_resid() function from the semptools package. This function
# rotates factor residual arrows  to prevent overlapping with factor correlation
# arrows

pdf(file="CFAofEFA.pdf", width=5, height = 3)

plot(rotate_resid(semPaths_plot, rotate_resid_list = c(F1=225,F2=225,F3=225,F4=225)))

dev.off()

# Five factor model: syntax
FiveModel <- 'F1 =~ OCD + AN + TS
F2 =~ BIP + SCZ
F3 =~ ANX + MDD
F4 =~ ADHD + AUT
F5 =~ ADHD + DYX
F1 ~~ F3
F1 ~~ F4
F1 ~~ F5
F2 ~~ F3
F2 ~~ F4
F2 ~~ F5
F3 ~~ F4
F3 ~~ F5
F4 ~~ F5
'

#get structural model
FiveFactor<-usermodel(LDSCoutput, estimation = "DWLS", model = FiveModel, std.lv=TRUE)

#get model fit statistics
FiveFactor$modelfit
#       chisq df      p_chisq      AIC       CFI       SRMR
# df 96.56931 24 1.148296e-10 158.5693 0.9838995 0.04777765

#remove aberrant correlation lines
results<-FiveFactor$results
results$formula <- paste(results$lhs,results$op, results$rhs)
unwanted <- c("F1 ~~ F2")
results_clean<-subset(results, !(results$formula %in% unwanted))
FiveFactor$results<-results_clean

#convert gsem output to format readable by semPaths
formatted_model <- semPlotModel_GSEM(FiveFactor)

# save model as a variable
semPaths_plot<-semPaths(formatted_model,
         layout="tree2",
         whatLabels="label",
         style="ram",
         curvature=5,
         edge.label.cex=0.8,
         reorder = TRUE,
         sizeMan=8,
         sizeLat=10,
         nCharNodes=4,
         nCharEdges=0,
         label.norm="OOOO",
         optimizeLatRes = TRUE,
         optimPoints = 0.4,
         fixedStyle = c("black",1),
         freeStyle = c("black",1),
         width=45,
         height=25,
         levels=c(1,2,6,7),
         mar=c(10,2,15,1),
         edge.label.position=0.5
         #filetype="png",
         #filename="path_diagram_5F.png"
)

#save model as .pdf
pdf(file="5fact.pdf", width=5, height = 3)

plot(rotate_resid(semPaths_plot, rotate_resid_list = c(F1=225,F2=225,F3=225,F4=225,F5=225)))

dev.off()
