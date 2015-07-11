## Start Date: May 21, 2014
## Package{drc}

#################################################################################
## GENERIC EC50 4PL BLOCK, modified for each fitting (last edit: 2014.06.09)
#################################################################################
## Set work directory, path depends on the location of the PC (work or home)
setwd("C:/Users/amedvedev/Google Drive/Current_Work/DB_work/2014-08-04-LOR/data")
## Load dataset: substitute "xxxxxxx.csv" with specific file name such as "pparg.csv". 
## File structure: 1st column is "Conc", 2nd column is "Fold"
dr1 <- read.csv("lor1.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
dr2 <- as.data.frame(dr1)
## Loading Package (drc)
library(drc)
################################
## Fitting a four-parameter logistic model (LL.4), if problem try LL.3 or LL.5.
pxrT.m1 <- drm(Fold ~ Conc, data = dr2, fct = LL.4())

## Displaying a summary of the model fit
summary(pxrT.m1)
## Plotting the fitted curve together with the original data
plot(pxrT.m1, main = "pxrT: 4PL Model", xlab="Conc (uM)", ylab="Induction (fold)") 
## New plot
plot(pxrT.m1, broken=TRUE, main = "PXR by T0901317: 4PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.pxrT.m1 <- ED(pxrT.m1, c(10, 50, 90), interval = "delta", level = 0.95)
## Save EC50 as .csv
write.csv(ed.pxrT.m1, file = "ec50-pxrT-m1.csv")

###############################
## Fitting a five-parameter logistic model (LL.5)
gr.m2 <- drm(Fold ~ Conc, data = dr2, fct = LL.5())

## Displaying a summary of the model fit
summary(gr.m2)
## Plotting the fitted curve together with the original data
plot(gr.m2, main = "GR: 5PL Model", xlab="Conc (uM)", ylab="Induction (fold)") 
## New plot
plot(gr.m2, broken=TRUE, main = "GR: 5PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.gr.m2 <- ED(gr.m2, c(10, 50, 90), interval = "delta", level = 0.95)
## Save EC50 as .csv
write.csv(ed.gr.m2, file = "ec50-gr-m2.csv")




## Set working directory
setwd("C:/Users/Alex/Google Drive/Current_Work/DB_work/2014-05-16-EC50/data")

## Load datasets
dr1 <- read.csv("pparg.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
dr2 <- as.data.frame(dr1)
## Load modified dataset where value 1 was subtracted from fold-induction values
dr3 <- read.csv("pparg1m.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
dr4 <- as.data.frame(dr3)

## Loading Package (drc)
library(drc)

################################
## Fitting a four-parameter logistic model (LL.4)
pparg.m1 <- drm(Fold ~ Conc, data = dr2, fct = LL.4())

## Displaying a summary of the model fit
summary(pparg.m1)

## Capturing (and appending) model output with time stamp in the .txt file 
##from http://stackoverflow.com/questions/1734896/r-capturing-elements-of-r-output-into-text-files
out<-date()
cat(out,file="pparg-out.txt",sep="\n",append=TRUE)
out<-capture.output(summary(pparg.m1))
cat(out,file="pparg-out.txt",sep="\n",append=TRUE)

## Plotting the fitted curve together with the original data
plot(pparg.m1, broken=TRUE, main = "PPARg: 4PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.pparg.m1 <- ED(pparg.m1, c(10, 50, 90), interval = "delta", level = 0.95)
zz <- ED(pparg.m1, c(50), interval = "delta", level = 0.95)
x <- zz[1]
abline(v = x )
## Save EC50 as .csv
write.csv(ed.pparg.m1, file = "ec50-pparg-m1.csv")
write.csv(qq, file = "model-pparg-m1.csv")
###############################
## Using modified dataset (subtracting value 1 from fold-induction values, so lowest values will be close to 0 not to 1)
## Fitting a four-parameter logistic model (LL.4)
pparg1m.m1 <- drm(Fold ~ Conc, data = dr4, fct = LL.4())

## Displaying a summary of the model fit
summary(pparg1m.m1)

## Plotting the fitted curve together with the original data
plot(pparg1m.m1, broken=TRUE, main = "PPARg: 4PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.pparg1m.m1 <- ED(pparg1m.m1, c(10, 50, 90), interval = "delta", level = 0.95)
## Save EC50 as .csv
write.csv(ed.pparg1m.m1, file = "ec50-pparg1m-m1.csv")

##############################
## Fitting a five-parameter logistic model (LL.5)
## Using original (unmodified) dataset
pparg.m2 <- drm(Fold ~ Conc, data = dr2, fct = LL.5())

## Displaying a summary of the model fit
summary(pparg.m2)

## Plotting the fitted curve together with the original data
plot(pparg.m2, broken=TRUE, main = "PPARg: 5PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.pparg.m2 <- ED(pparg.m2, c(10, 50, 90), interval = "delta", level = 0.95)
## Save EC50 as .csv
write.csv(ed.pparg.m2, file = "ec50-pparg-m2.csv")

##############################
## Fitting a three-parameter logistic model (LL.3)
## Using original (unmodified) dataset
pparg.m3 <- drm(Fold ~ Conc, data = dr2, fct = LL.3())

## Displaying a summary of the model fit
summary(pparg.m3)

## Plotting the fitted curve together with the original data
plot(pparg.m3, broken=TRUE, main = "PPARg: 3PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.pparg.m3 <- ED(pparg.m3, c(10, 50, 90), interval = "delta", level = 0.95)

## Save EC50 as .csv
write.csv(ed.pparg.m3, file = "ec50-pparg-m3.csv")

##############################

## All models Summary plot
plot(pparg.m1, broken=TRUE, main = "PPARg: 4PL, 5PL and 3PL Models", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)
minor.tick()
plot(pparg.m2, add=TRUE, broken=TRUE, lty=3, lwd=2)
plot(pparg.m3, add=TRUE, broken=TRUE, lty=3, lwd=2)
plot(, add=TRUE, broken=TRUE, lty=3, lwd=2)
arrows(0.1, 15, 0.05, 15, 0.15, lwd=2)
text(0.12,15, "EC50=0.042uM", pos=4, cex=1.2)
abline(v = 0.42)

##################
## NEW TEST BLOCK
##################
## Fitting a four-parameter logistic with coefficient names defined
pparg.m4 <- drm(Fold ~ Conc, data = dr2, fct =  LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
## Displaying a summary of the model fit
summary(pparg.m4)
## New plot
plot(pparg.m4, broken=TRUE, xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)
plot(pparg.m4, add=TRUE, broken=TRUE, lty=3, lwd=2)
##################
## NEW TEST BLOCK END
##################

#################################################################################
## GENERIC EC50 4PL BLOCK, modified for each fitting (last edit: 2014.06.09)
#################################################################################
## Set work directory, path depends on the location of the PC (work or home)
setwd("C:/Users/amedvedev/Google Drive/Current_Work/DB_work/2014-05-16-EC50/data")
## Load dataset: substitute "xxxxxxx.csv" with specific file name such as "pparg.csv". 
## File structure: 1st column is "Conc", 2nd column is "Fold"
dr1 <- read.csv("pxrT.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA",""))
dr2 <- as.data.frame(dr1)
## Loading Package (drc)
library(drc)
################################
## Fitting a four-parameter logistic model (LL.4), if problem try LL.3 or LL.5.
pxrT.m1 <- drm(Fold ~ Conc, data = dr2, fct = LL.4())

## Displaying a summary of the model fit
summary(pxrT.m1)
## Plotting the fitted curve together with the original data
plot(pxrT.m1, main = "pxrT: 4PL Model", xlab="Conc (uM)", ylab="Induction (fold)") 
## New plot
plot(pxrT.m1, broken=TRUE, main = "PXR by T0901317: 4PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.pxrT.m1 <- ED(pxrT.m1, c(10, 50, 90), interval = "delta", level = 0.95)
## Save EC50 as .csv
write.csv(ed.pxrT.m1, file = "ec50-pxrT-m1.csv")
###############################
## Fitting a five-parameter logistic model (LL.5)
gr.m2 <- drm(Fold ~ Conc, data = dr2, fct = LL.5())

## Displaying a summary of the model fit
summary(gr.m2)
## Plotting the fitted curve together with the original data
plot(gr.m2, main = "GR: 5PL Model", xlab="Conc (uM)", ylab="Induction (fold)") 
## New plot
plot(gr.m2, broken=TRUE, main = "GR: 5PL Model", xlab="Conc (uM)", ylab="Induction (fold)", lwd=2, 
     cex=1.2, cex.axis=1.2, cex.lab=1.2)

## ED estimates effective doses (ECp/EDp/ICp) for given reponse levels.
ed.gr.m2 <- ED(gr.m2, c(10, 50, 90), interval = "delta", level = 0.95)
## Save EC50 as .csv
write.csv(ed.gr.m2, file = "ec50-gr-m2.csv")


##################
## GENERIC EC50 4PL BLOCK END
##################




################################################
## Fitting a four-parameter Weibull model (type 1)
ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
plot(ryegrass.m2)

## Fitting a four-parameter log-logistic model
## with user-defined parameter names
ryegrass.m3 <- drm(rootl ~ conc, data = ryegrass, 
fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(ryegrass.m3)

## Comparing log-logistic and Weibull models
## (Figure 2 in Ritz (2009))
ryegrass.m0 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4())

plot(ryegrass.m0, broken=TRUE, xlab="Dose (mM)", ylab="Root length (cm)", lwd=2, 
cex=1.2, cex.axis=1.2, cex.lab=1.2)
plot(ryegrass.m1, add=TRUE, broken=TRUE, lty=2, lwd=2)
plot(ryegrass.m2, add=TRUE, broken=TRUE, lty=3, lwd=2)

arrows(3, 7.5, 1.4, 7.5, 0.15, lwd=2)
text(3,7.5, "Weibull-2", pos=4, cex=1.2)

arrows(2.5, 0.9, 5.7, 0.9, 0.15, lwd=2)
text(3,0.9, "Weibull-1", pos=2, cex=1.2)
