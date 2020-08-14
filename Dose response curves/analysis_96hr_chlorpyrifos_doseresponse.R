#RUN THIS FILE THIRD

# Upated for Daniel Frank 2015, updated by Paige Mundy, 2019

##############################################################
###############################################################

Db = read.csv('DM_96h.csv')

##############################################################
###############################################################

Db <- Db %>% 
 dplyr::select(-X)

View(Db)
#install.packages('stats')
library(bbmle)
library(lme4)
library(stats)
library(multcomp)

# Manipulations on columns to prepare the data
# Constrain some of the columns to act as factors, not numerics
Db$Trt = Db$Treatment
Db$Trt = log10(Db$Trt + 0.05) # add to avoid taking log(0) and to make the values evenly spaced on a log scale

Db$Trt_f = as.factor(Db$Treatment)


Cycles = colnames(Db)
Cycles = Cycles[3:7]
View(Db)

# Dose-response curve fitting

#--------------------------------------------------------------
# Do the auto-curve-fitting for all the Cycles
# This will loop over each Cycle, fit all possible curves,
# Then export the p-values for each fit to the Excel file listed below.
# It will also put all a figure showing the best curve into a directory
# named "dose_response_figures" (need to make that folder)
Curve.pvals = matrix(NA,length(Cycles),12)
for (g in 1:length(Cycles)){
  X = fit.curves(Db,Cycles[g])
  Curve.pvals[g, ] = X$Pval 
}

Curve.pvals = data.frame(Curve.pvals,row.names=Cycles)
names(Curve.pvals) = c('Linear','Quad','Sigmoid','Unimodal1','Unimodal2')

write.csv(Curve.pvals,'Curve_DM_96h.csv')
#--------------------------------------------------------------

 