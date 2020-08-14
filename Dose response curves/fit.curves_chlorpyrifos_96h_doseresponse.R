#To use this code there are 3 files. FIRST, Open up fit.curves.r, select-all and run that code. Do the same thing for dose_response_functions.r. Once you have run both of those files (which just contain all of the helper functions) you can go through the analysis file to execute the analyses for each Cycle.

# Function to automatically fit all of the different D-R curves for a given Cycle

# Modified 19 Aug 2015 to have different starting param values for some, to improve fitting.

#Modified 2019 for behavior data

install.packages('bbmle')
library(bbmle)
library(lme4)
library(stats)

#libraries for running both ANOVA (linear mixed effects model) and dose-response curves

fit.curves = function(Db,Cycle){
  # Read in the correct data
  Col = colnames(Db)==Cycle
  OKrows = !is.na(Db[,Col]) #& Db[,Plate]=='1'
  X = Db$Trt[OKrows]
  Y = Db[OKrows,Col]
  
  
  # Intercept-only
  mBase = mle2(Int,start=c(g = mean(Y),s=1),data=list(Y=Y))
  
  # Linear
  mLin = mle2(Linear,start=c(b0 = mean(Y),b1=1,s=1),data=list(X=X,Y=Y))
  
  # Sigmoidal
  if (any(Cycle == c('Dark1', 'Light1', 'Dark2', 'Light2', 'Dark3'))){
    A = -1; B = -10;}else{
      A = 1; B = 10}
  mSigm = mle2(regDR,start=c(g=max(Y),h=min(Y),a=A,b=B,s=1),data=list(X=X,Y=Y))
  
  # Unimodal (complicated)
  G = 0 # Difference between maximum and minimum values (i.e., y-range of data
  H = -1 # Minimum value
  # parameters for the increasing portion:
  A1 = 2 # Intercept (larger positive numbers move this to the left)
  B1 = 2 # Rate of increase (must be > 1)
  # parameters for the decreasing portion
  A2 = 2 # Intercept (larger positive numbers move this to the right)
  B2 = 2 # Rate of decrease (must be < 1)
  mUni1 = mle2(uniDR,start=c(g=G,h=H,a1=A1,b1=B1,a2=A2,b2=B2,s=1),data=list(X=X,Y=Y), 
               control=list(maxit=5000, reltol=1e-4)) #new HERE ########
  
  # Unimodal (reduced)
  if (any(Cycle == c('Dark1', 'Light1', 'Dark2', 'Light2', 'Dark3'))){
    G = 2 # 0 # Difference between maximum and minimum values (i.e., y-range of data
    H = 0 #-1 # Minimum value
    # parameters for the increasing portion:
    B1 = 10 #2 # Rate of increase (must be > 1)
    # parameters for the decreasing portion
    A2 = 2 # Intercept (larger positive numbers move this to the right)
    B2 = -2 # Rate of decrease (must be < 1)
  }else{
    if (any(Cycle== 'Dark1', 'Light1', 'Dark2', 'Light2', 'Dark3')){
      G = 2 # 0 # Difference between maximum and minimum values (i.e., y-range of data
      H = 0 #-1 # Minimum value
      # parameters for the increasing portion:
      B1 = 2 #2 # Rate of increase (must be > 1)
      # parameters for the decreasing portion
      A2 = -1 # Intercept (larger positive numbers move this to the right)
      B2 = -1 # Rate of decrease (must be < 1)     
    }else{
      G = 2 # 0 # Difference between maximum and minimum values (i.e., y-range of data
      H = 0 #-1 # Minimum value
      # parameters for the increasing portion:
      B1 = 2 #2 # Rate of increase (must be > 1)
      # parameters for the decreasing portion
      A2 = 0 # Intercept (larger positive numbers move this to the right)
      B2 = -1 # Rate of decrease (must be < 1)   
    }}
  mUni2 = mle2(uniDR2,start=c(g=G,h=H,b1=B1,a2=A2,b2=B2,s=1),data=list(X=X,Y=Y), 
               control=list(maxit=5000, reltol=1e-4)) #HERE
  
  # Unimodal (reduced #2)
  #G = 0 # Difference between maximum and minimum values (i.e., y-range of data
  #H = -1 # Minimum value
  # parameters for the increasing portion:
  # B1 = 2 # Rate of increase (must be > 1)
  # parameters for the decreasing portion
  #  B2 = -4 # Rate of decrease (must be < 1)
  # mUni3 = mle2(uniDR3,start=c(g=G,h=H,b1=B1,b2=B2,s=1),data=list(X=X,Y=Y))
  
  # Quadratic
  mQuad = mle2(Quadratic,start=c(b0 = mean(Y),b1=1,b2=mean(X),s=1),data=list(X=X,Y=Y), 
               control=list(maxit=5000, reltol=1e-4)) #HERE
  
  AIC.table = AICctab(mBase,mLin,mSigm,mUni1,mUni2,mQuad,nobs=length(X))
  
  ANOVA = list()
  ANOVA[[1]] = anova(mBase,mLin)
  ANOVA[[2]] = anova(mBase,mQuad)
  ANOVA[[3]] = anova(mBase,mSigm)
  ANOVA[[4]] = anova(mBase,mUni2)
  ANOVA[[5]] = anova(mBase,mUni1)
  # ANOVA.6 = anova(mBase,mUni3)
  
  Pval = rep(1,3) # number of replicates
  Dev = rep(1,3) # number of replicates
  for (a in 1:3){
    Pval[a]=ANOVA[[a]][10] # pvalue
    Dev[a]=ANOVA[[a]][4]} # deviance (-2 * log likelihood)
  
  # Pick the best one and make a plot
  #Best= match(min(Pval),Pval)
  Best = match(min(Pval),Pval)
  
  YL = paste(Cycle)
  plot(X,Y,xaxp=c(-3,3,12),
       ylab=YL,
       yaxt='n', xaxt='n',
       ylim=c(0,1))
  axis(2, seq(0,1.5,0.5))
  X_temp = seq(-5,5,length.out=1000) # dummy x-values for plotting the curve
  
  B = coef(mLin)
  P = B[1] + B[2]*X_temp
  # lines(X_temp,P,col='black')
  
  B = coef(mSigm)
  P = (B[2]-B[1])/(1+exp(B[3] + B[4]*X_temp))+B[1]
  #  lines(X_temp,P,col='blue')
  
  B = coef(mUni2)
  P = B[1]+B[2]*(1 + B[3]*X_temp)/(1+exp(-(B[4] + B[5]*X_temp)))
  # lines(X_temp,P,col='green')
  
  B = coef(mQuad)
  P = B[1]+B[2]*(X_temp-B[3])^2
  #  lines(X_temp,P,col='red')
  
  B = coef(mUni1)
  P = B[1]+B[2]*(1 + exp(-(B[3]+B[4]*X_temp)))/(1+exp(-(B[5] + B[6]*X_temp)))
  #  lines(X_temp,P,col='black')
  
  ##B = coef(mUni3)
  #P = B[1]+B[2]*(1 + B[3]*X_temp)/(1+B[4]*X_temp)
  #lines(X_temp,P,col='cyan')
  
  if (Pval[Best]<0.05){Lty = 1}else{Lty=2}
  
  if (Best==0){ # Baseline
    lines(c(X_temp[1],X_temp[1000]),c(mean(Y),mean(Y)))
    File.name = paste("dose_response_figures/",Cycle,"_baseline_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==1){ # Linear
    B = coef(mLin)
    P = B[1]+B[2]*X_temp
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Cycle,"_linear_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==3){ # Sigmoidal
    B = coef(mSigm)
    P = (B[2]-B[1])/(1+exp(B[3] + B[4]*X_temp))+B[1]
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Cycle,"_sigmoidal_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==5){ # Uni1
    B = coef(mUni1)
    P = B[1]+B[2]*(1 + exp(-(B[3]+B[4]*X_temp)))/(1+exp(-(B[5] + B[6]*X_temp)))
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Cycle,"_unimodal1_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==4){ # Uni2
    B = coef(mUni2)
    P = B[1]+B[2]*(1 + B[3]*X_temp)/(1+exp(-(B[4] + B[5]*X_temp)))
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Cycle,"_unimodal2_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  if (Best==2){ # Quadratic
    B = coef(mQuad)
    P = B[1]+B[2]*(X_temp-B[3])^2
    lines(X_temp,P,lty=Lty)
    File.name = paste("dose_response_figures/",Cycle,"_quadratic_dose_response.pdf",sep="")
    dev.print(device=pdf,file=File.name,useDingbats=FALSE)
  }
  
  return(list("baseline"=mBase,"sigmoidal"=mSigm,"linear"=mLin,
              "unimodal1"=mUni1,"unimodal2"=mUni2,"quadratic"=mQuad,
              "AIC"=AIC.table,"ANOVA"=ANOVA,'FN'=File.name,'B'=Best,
              "Pval"=Pval))
}
parnames(fit.curves) = c('Db','Cycle')

