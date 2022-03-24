# ----------------------------------
# Simulation studies for "Regression tree credibility model" by L Diao and C Weng
# NAAJ 2019
# Coding by Liqun Diao
# This file contains the simulation code for one of many simulation settings
# Homogeneous - 5 observations per subject - lognormal distribution
# ----------------------------------

# loading in rpart package - classification and regression trees
library(rpart)
set.seed(50000)
file.nam="homo_5_1.dat"
cat(" ", "\n", append=T, file=file.nam)

nSim = 1000 # number of simulations
nu = 300 # number of subjects
nfold = 5 # number of fold of cross validation
dist = 1 # distribution 1 = lognormal;2 = exponential;3 = pareto
n = 5 #number of observations per subject
 

# --- user-written loss function for rpart package -----


# Initialization function 

init <- function(y,offset,parms, wt)
{
  
  sfun <- function(yval, dev, wt, ylevel, digits) {
  
  paste(" mean=", format(signif(yval, digits)), ", MSE=", format(signif(dev/wt, digits)), sep = "")
  
  }
  
  environment(sfun) <- .GlobalEnv
  
  list(y = y, parms = parms, numresp = 1, numy = 1, summary = sfun)  
}


# ________________________________________________________
# 
# ---------- 4 CREDIBILITY LOSS FUNCTIONS ---------------
#
# ________________________________________________________


# -----------------------------------------------------
#       credibility loss one 
# -----------------------------------------------------

cred1_e = function(y, wt, parms) {

nu = parms[length(parms)] # number of subjects in the data
n  = parms[length(parms)-1] # number of observation of each subject
wm = parms[1:nu] # mean of each subjects in the data

num = length(y) # number of subjects in this node

rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
wmean = y # mean of each subject in this node

# calculating tau and sigma
sigma = mean(wvar)
tau = var(wmean) - sigma/n
tau = max(tau,0)

alpha=n*tau/(n*tau+sigma)
# calculating loss function
rss = num*sigma*tau/(n*tau+sigma)

list(label=alpha , deviance=rss)
}
  
  cred1_s = function(y, wt, x, parms, continuous) {
  
  nu = parms[length(parms)] # number of subjects in the data
  n  = parms[length(parms)-1] # number of observation of each subject
  wm = parms[1:nu] # mean of each subjects in the data
  
  num = length(y) # number of subjects in parent node
  
  
  rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
  wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
  wmean = y # mean of each subject in this node
  
  
  # -- parameters and loss in the parent node	--
  
  # calculating tau and sigma
  sigma = mean(wvar)
  tau = var(wmean)-sigma/n
  tau = max(tau,0)
  # calculating loss function
  SSA = num*sigma*tau/(n*tau+sigma)
  
  if (continuous) {
  # Continuous X variable
  # calculating node sizes of left and right child nodes
  wtL = cumsum(rep(1,num))[-num]
  wtR = num-wtL
  
  # -- parameters and loss in the child nodes	--
  # calculating tau and sigma
  sigmaL = cumsum(wvar)[-num]/wtL
  sigmaR = (sum(wvar)-cumsum(wvar)[-num])/wtR
  
  tauL  = cumsum(wmean^2)[-num]/(wtL-1)-(cumsum(wmean)[-num])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(wmean^2)-cumsum(wmean^2)[-num])/(wtR-1)-(sum(wmean)-cumsum(wmean)[-num])^2/(wtR-1)/wtR-sigmaR/n
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  # calculating loss function
  SSL = wtL*sigmaL*tauL/(n*tauL+sigmaL)
  SSR = wtR*sigmaR*tauR/(n*tauR+sigmaR)
  
  goodness1 = SSA-SSL-SSR
  
  list(goodness= goodness1, direction=sign(cumsum(wmean)[-num]/wtL-mean(wmean)))
  
  }else{
  
  
  # Categorical X variable
  ux <- sort(unique(x))
  wtsum <- tapply(wt, x, sum)
  yvarsum <- tapply(wvar, x, sum)
  sigmac <- yvarsum/wtsum
  tauc = tapply(wmean, x, var) -sigmac/n
  tauc = pmax(tauc,0)
  tauc = ifelse(is.na(tauc)==T,0,tauc)
  
  alphac = n/(n+sigmac/tauc)
  
  ymeansum <- tapply(wmean, x, sum)
  
  
  # We can order the categories by their alpha
  # then use the same code as for a non-categorical
  ord <- order(alphac)
  nn <- length(ord)
  
  # calculating node sizes of left and right child nodes
  wtL = cumsum(wtsum[ord])[-nn]
  wtR = num - wtL
  
  sigmaL = cumsum(yvarsum[ord])[-nn]/wtL
  sigmaR = (sum(yvarsum)-cumsum(yvarsum[ord])[-nn])/wtR
  
  tauL  = cumsum(ymeansum[ord]^2)[-nn]/(wtL-1)-(cumsum(ymeansum[ord])[-nn])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(ymeansum^2)-cumsum(ymeansum[ord]^2)[-nn])/(wtR-1)-(sum(ymeansum)-cumsum(ymeansum[ord])[-nn])^2/(wtR-1)/wtR-sigmaR/n
  
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  # calculating loss function
  SSL = wtL*sigmaL*tauL/(n*tauL+sigmaL)
  SSR = wtR*sigmaR*tauR/(n*tauR+sigmaR)
  
  goodness2 = SSA-SSL-SSR
  
  list(goodness= goodness2, direction=ux[ord])
  
  }
  
  
  }
  
  
  
  
  # -----------------------------------------------------
  #       credibility loss two 
  # -----------------------------------------------------
  
  cred2_e = function(y, wt, parms) {
  
  nu = parms[length(parms)] # number of subjects in the data
  n  = parms[length(parms)-1] # number of observation of each subject
  wm = parms[1:nu] # mean of each subjects in the data
  
  num = length(y) # number of subjects in this node
  
  rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
  wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
  wmean = y # mean of each subject in this node
  
  # calculating tau and sigma
  sigma = mean(wvar)
  tau = var(wmean) - sigma/n
  tau = max(tau,0)
  
  # calculating loss function
  rss = num*sigma*((n+1)*tau+sigma)/(n*tau+sigma)
  
  alpha=n*tau/(n*tau+sigma)
  
  list(label= alpha, deviance=rss)
  }
  
  cred2_s = function(y, wt, x, parms, continuous) {
  
  nu = parms[length(parms)] # number of subjects in the data
  n  = parms[length(parms)-1] # number of observation of each subject
  wm = parms[1:nu] # mean of each subjects in the data
  
  num = length(y) # number of subjects in parent node
  
  rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
  wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
  wmean = y # mean of each subject in this node
  
  # -- parameters and loss in the parent node	--
  
  # calculating tau and sigma
  sigma = mean(wvar)
  tau = var(wmean)-sigma/n
  tau = max(tau,0)
  # calculating loss function
  SSA = num*sigma*((n+1)*tau+sigma)/(n*tau+sigma)
  
  
  if (continuous) {
  
  # calculating node sizes of left and right child nodes
  wtL = cumsum(rep(1,num))[-num]
  wtR = num-wtL
  
  # -- parameters and loss in the child nodes	--
  
  # calculating tau and sigma
  sigmaL = cumsum(wvar)[-num]/wtL
  sigmaR = (sum(wvar)-cumsum(wvar)[-num])/wtR
  
  tauL  = cumsum(wmean^2)[-num]/(wtL-1)-(cumsum(wmean)[-num])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(wmean^2)-cumsum(wmean^2)[-num])/(wtR-1)-(sum(wmean)-cumsum(wmean)[-num])^2/(wtR-1)/wtR-sigmaR/n
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  # calculating loss function
  SSL = wtL*sigmaL*((n+1)*tauL+sigmaL)/(n*tauL+sigmaL)
  SSR = wtR*sigmaR*((n+1)*tauR+sigmaR)/(n*tauR+sigmaR)
  
  goodness1 = SSA-SSL-SSR
  
  list(goodness= goodness1, direction=sign(cumsum(wmean)[-num]/wtL-mean(wmean)))
  }else{
  
  
  # Categorical X variable
  ux <- sort(unique(x))
  wtsum <- tapply(wt, x, sum)
  yvarsum <- tapply(wvar, x, sum)
  sigmac <- yvarsum/wtsum
  tauc = tapply(wmean, x, var) -sigmac/n
  tauc = pmax(tauc,0)
  tauc = ifelse(is.na(tauc)==T,0,tauc)
  
  alphac = n/(n+sigmac/tauc)
  
  ymeansum <- tapply(wmean, x, sum)
  
  
  # We can order the categories by their alpha
  # then use the same code as for a non-categorical
  ord <- order(alphac)
  nn <- length(ord)
  
  # calculating node sizes of left and right child nodes
  wtL = cumsum(wtsum[ord])[-nn]
  wtR = num - wtL
  
  sigmaL = cumsum(yvarsum[ord])[-nn]/wtL
  sigmaR = (sum(yvarsum)-cumsum(yvarsum[ord])[-nn])/wtR
  
  tauL  = cumsum(ymeansum[ord]^2)[-nn]/(wtL-1)-(cumsum(ymeansum[ord])[-nn])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(ymeansum^2)-cumsum(ymeansum[ord]^2)[-nn])/(wtR-1)-(sum(ymeansum)-cumsum(ymeansum[ord])[-nn])^2/(wtR-1)/wtR-sigmaR/n
  
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  # calculating loss function
  SSL = wtL*sigmaL*((n+1)*tauL+sigmaL)/(n*tauL+sigmaL)
  SSR = wtR*sigmaR*((n+1)*tauR+sigmaR)/(n*tauR+sigmaR)
  
  goodness2 = SSA-SSL-SSR
  
  list(goodness= goodness2, direction=ux[ord])
  
  }
  
  }
  
  
  
  # -----------------------------------------------------
  #       credibility loss three 
  # -----------------------------------------------------
  
  cred3_e = function(y, wt, parms) {
  
  nu = parms[length(parms)] # number of subjects in the data
  n  = parms[length(parms)-1] # number of observation of each subject
  wm = parms[1:nu] # mean of each subjects in the data
  
  num = length(y) # number of subjects in this node
  
  rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
  wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
  wmean = y # mean of each subject in this node
  
  # calculating tau and sigma
  sigma = mean(wvar)
  tau = var(wmean) - sigma/n
  tau = max(tau,0)
  alpha=n*tau/(n*tau+sigma)
  
  # calculating loss function
  rss = num*(1-alpha)*tau+(1-alpha)^2*(n*tau+sigma)/n
  
  list(label= alpha, deviance=rss)
  }
  
  cred3_s = function(y, wt, x, parms, continuous) {
  
  nu = parms[length(parms)] # number of subjects in the data
  n  = parms[length(parms)-1] # number of observation of each subject
  wm = parms[1:nu] # mean of each subjects in the data
  
  num = length(y) # number of subjects in parent node
  
  
  rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
  wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
  wmean = y # mean of each subject in this node
  
  # -- parameters and loss in the parent node	--
  
  # calculating tau and sigma
  sigma = mean(wvar)
  tau = var(wmean)-sigma/n
  tau = max(tau,0)
  alpha=n*tau/(n*tau+sigma)
  
  # calculating loss function
  SSA = num*(1-alpha)*tau+(1-alpha)^2*(n*tau+sigma)/n
  
  
  if (continuous) {
  # calculating node sizes of left and right child nodes
  wtL = cumsum(rep(1,num))[-num]
  wtR = num-wtL
  
  # -- parameters and loss in the child nodes	--
  
  # calculating tau and sigma
  sigmaL = cumsum(wvar)[-num]/wtL
  sigmaR = (sum(wvar)-cumsum(wvar)[-num])/wtR
  
  tauL  = cumsum(wmean^2)[-num]/(wtL-1)-(cumsum(wmean)[-num])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(wmean^2)-cumsum(wmean^2)[-num])/(wtR-1)-(sum(wmean)-cumsum(wmean)[-num])^2/(wtR-1)/wtR-sigmaR/n
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  alphaL=n*tauL/(n*tauL+sigmaL)
  alphaR=n*tauR/(n*tauR+sigmaR)
  
  # calculating loss function
  SSL = wtL*(1-alphaL)*tauL+(1-alphaL)^2*(n*tauL+sigmaL)/n
  SSR = wtR*(1-alphaR)*tauR+(1-alphaR)^2*(n*tauR+sigmaR)/n
  
  goodness1 = SSA-SSL-SSR
  
  list(goodness= goodness1, direction=sign(cumsum(wmean)[-num]/wtL-mean(wmean)))
  }else{
  
  
  # Categorical X variable
  ux <- sort(unique(x))
  wtsum <- tapply(wt, x, sum)
  yvarsum <- tapply(wvar, x, sum)
  sigmac <- yvarsum/wtsum
  tauc = tapply(wmean, x, var) -sigmac/n
  tauc = pmax(tauc,0)
  tauc = ifelse(is.na(tauc)==T,0,tauc)
  
  alphac = n/(n+sigmac/tauc)
  
  ymeansum <- tapply(wmean, x, sum)
  
  
  # We can order the categories by their alpha
  # then use the same code as for a non-categorical
  ord <- order(alphac)
  nn <- length(ord)
  
  # calculating node sizes of left and right child nodes
  wtL = cumsum(wtsum[ord])[-nn]
  wtR = num - wtL
  
  sigmaL = cumsum(yvarsum[ord])[-nn]/wtL
  sigmaR = (sum(yvarsum)-cumsum(yvarsum[ord])[-nn])/wtR
  
  tauL  = cumsum(ymeansum[ord]^2)[-nn]/(wtL-1)-(cumsum(ymeansum[ord])[-nn])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(ymeansum^2)-cumsum(ymeansum[ord]^2)[-nn])/(wtR-1)-(sum(ymeansum)-cumsum(ymeansum[ord])[-nn])^2/(wtR-1)/wtR-sigmaR/n
  
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  alphaL=n*tauL/(n*tauL+sigmaL)
  alphaR=n*tauR/(n*tauR+sigmaR)
  
  # calculating loss function
  SSL = wtL*(1-alphaL)*tauL+(1-alphaL)^2*(n*tauL+sigmaL)/n
  SSR = wtR*(1-alphaR)*tauR+(1-alphaR)^2*(n*tauR+sigmaR)/n
  
  goodness2 = SSA-SSL-SSR
  
  list(goodness= goodness2, direction=ux[ord])
  
  }
  
  }
  
  
  # -----------------------------------------------------
  #       credibility loss four 
  # -----------------------------------------------------
  
  cred4_e = function(y, wt, parms) {
  
  nu = parms[length(parms)] # number of subjects in the data
  n  = parms[length(parms)-1] # number of observation of each subject
  wm = parms[1:nu] # mean of each subjects in the data
  
  num = length(y) # number of subjects in this node
  
  rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
  wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
  wmean = y # mean of each subject in this node
  
  # calculating tau and sigma
  sigma = mean(wvar)
  tau = var(wmean) - sigma/n
  tau = max(tau,0)
  alpha=n*tau/(n*tau+sigma)
  
  # calculating loss function
  rss = num*(1-alpha)*tau+(1-alpha)^2*(n*tau+sigma)/n+sigma*num
  
  list(label= alpha, deviance=rss)
  }
  
  cred4_s = function(y, wt, x, parms, continuous) {
  
  nu = parms[length(parms)] # number of subjects in the data
  n  = parms[length(parms)-1] # number of observation of each subject
  wm = parms[1:nu] # mean of each subjects in the data
  
  num = length(y) # number of subjects in parent node
  # calculating node sizes of left and right child nodes
  wtL = cumsum(rep(1,num))[-num]
  wtR = num-wtL
  
  rsett = mapply(function(t){(1:nu)[wm==t]},t=y) # which individuals fall into this node
  wvar  = parms[(nu+1):(2*nu)][rsett] #variance of each subject in this node
  wmean = y # mean of each subject in this node
  
  if (continuous) {
  
  # -- parameters and loss in the parent node	--
  
  # calculating tau and sigma
  sigma = mean(wvar)
  tau = var(wmean)-sigma/n
  tau = max(tau,0)
  alpha=n*tau/(n*tau+sigma)
  
  # calculating loss function
  SSA = num*(1-alpha)*tau+(1-alpha)^2*(n*tau+sigma)/n+sigma*num
  
  # -- parameters and loss in the child nodes	--
  
  # calculating tau and sigma
  sigmaL = cumsum(wvar)[-num]/wtL
  sigmaR = (sum(wvar)-cumsum(wvar)[-num])/wtR
  
  tauL  = cumsum(wmean^2)[-num]/(wtL-1)-(cumsum(wmean)[-num])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(wmean^2)-cumsum(wmean^2)[-num])/(wtR-1)-(sum(wmean)-cumsum(wmean)[-num])^2/(wtR-1)/wtR-sigmaR/n
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  alphaL=n*tauL/(n*tauL+sigmaL)
  alphaR=n*tauR/(n*tauR+sigmaR)
  
  
  # calculating loss function
  SSL = wtL*(1-alphaL)*tauL+(1-alphaL)^2*(n*tauL+sigmaL)/n+wtL*sigmaL
  SSR = wtR*(1-alphaR)*tauR+(1-alphaR)^2*(n*tauR+sigmaR)/n+wtR*sigmaR
  
  
  goodness1 = SSA-SSL-SSR
  
  list(goodness= goodness1, direction=sign(cumsum(wmean)[-num]/wtL-mean(wmean)))
  }else{
  
  
  # Categorical X variable
  ux <- sort(unique(x))
  wtsum <- tapply(wt, x, sum)
  yvarsum <- tapply(wvar, x, sum)
  sigmac <- yvarsum/wtsum
  tauc = tapply(wmean, x, var) -sigmac/n
  tauc = pmax(tauc,0)
  tauc = ifelse(is.na(tauc)==T,0,tauc)
  
  alphac = n/(n+sigmac/tauc)
  
  ymeansum <- tapply(wmean, x, sum)
  
  
  # We can order the categories by their alpha
  # then use the same code as for a non-categorical
  ord <- order(alphac)
  nn <- length(ord)
  
  # calculating node sizes of left and right child nodes
  wtL = cumsum(wtsum[ord])[-nn]
  wtR = num - wtL
  
  sigmaL = cumsum(yvarsum[ord])[-nn]/wtL
  sigmaR = (sum(yvarsum)-cumsum(yvarsum[ord])[-nn])/wtR
  
  tauL  = cumsum(ymeansum[ord]^2)[-nn]/(wtL-1)-(cumsum(ymeansum[ord])[-nn])^2/wtL/(wtL-1)-sigmaL/n
  tauR  = (sum(ymeansum^2)-cumsum(ymeansum[ord]^2)[-nn])/(wtR-1)-(sum(ymeansum)-cumsum(ymeansum[ord])[-nn])^2/(wtR-1)/wtR-sigmaR/n
  
  tauL[is.na(tauL)==T]=0
  tauR[is.na(tauR)==T]=0
  tauL[tauL<0]=0
  tauR[tauR<0]=0
  
  alphaL=n*tauL/(n*tauL+sigmaL)
  alphaR=n*tauR/(n*tauR+sigmaR)
  
  # calculating loss function
  SSL = wtL*(1-alphaL)*tauL+(1-alphaL)^2*(n*tauL+sigmaL)/n+wtL*sigmaL
  SSR = wtR*(1-alphaR)*tauR+(1-alphaR)^2*(n*tauR+sigmaR)/n+wtR*sigmaR
  
  goodness2 = SSA-SSL-SSR
  
  list(goodness= goodness2, direction=ux[ord])
  
  }
  
  }


# ------------------------------------------
# longitudinal cross validation - since interests will be in predicting next outcome of each subject with given observations of the subject 
# ------------------------------------------  
  longi_cv = function(fitt,mlist,dat,yy,xx,n,nfold,fold)
  {
  
  nn = length(fitt$cp[,1]) # number of cp values
  nu = nrow(dat)
  cp = (c(1000,fitt$cp[,1])/2 + c(fitt$cp[,1],0)/2) # the average of adjacent cp values
  nsize = n/nfold
  # matrix - row represents each subject - column represents each cp value
  xfit=matrix(0,nu,nn)
  
  for (i in 1:nn)
  {
  cpi = cp[i] # the ith cp value
  
  for (j in 1:nfold)
  {
  
  wmeanj = rowMeans(yy[,fold!=j]) # individual mean without the jth observation
  wvarj  = rowSums((yy[,fold!=j]-wmeanj)^2)/(n-nsize-1) # individual var without the jth observation denominator (n-1)
  
  datj = data.frame(y=wmeanj,xx=xx) #data set without the jth obser
  
  parmsj=c(wmeanj,wvarj,n-nsize,nu) # parameters to be put in rpart
  
  fittj = rpart(y~.,method=mlist,parms=parmsj,maxsurrogate=0,cp=cpi,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=datj)
  
  predj = predict(fittj,data=datj)
  
  nsetj=unique(fittj$where) # terminal nodes in tree fittj
  
  for (k in nsetj)
  {
  wmeanjk = wmeanj[fittj$where==k]
  alphajk = predj[fittj$where==k][1]
  pred_jk = alphajk * wmeanjk + (1-alphajk)* mean(wmeanjk)
  
  xfit[fittj$where==k,i] = xfit[fittj$where==k,i] + rowSums(as.matrix((yy[fittj$where==k,fold==j] - pred_jk)^2))
  }
  }
  }
  
  xerror = colMeans(xfit)
  
  overall.min = which.min(xerror)
  #first.min = which(xerror-c(xerror[-1],max(xerror)+1)<0)[1]
  
  ncp1 = (1:nn)[overall.min] #the ncp th cp value selected
  #ncp2 = (1:nn)[first.min] #the ncp th cp value selected
  
  CP_selected1 = mean(c(cp[ncp1],cp[ncp1+1]))
  #CP_selected2 = mean(c(cp[ncp2],cp[ncp2+1]))
  
  return (CP_selected1)
  
  }
  
  


  # naive sets
  generate_naive_sets=function(x1,x2,x3,x4,x5,nu){
  set2=ifelse(x2<=50,1,0)
  set4=ifelse(x4<=50,1,0)
  
  set123=numeric(nu)
  set123[(x1<=50)&(x2<=50)&(x3<=50)]=1
  set123[(x1<=50)&(x2<=50)&(x3>50)]=2
  set123[(x1<=50)&(x2>50)&(x3<=50)]=3
  set123[(x1<=50)&(x2>50)&(x3>50)]=4
  set123[(x1>50)&(x2<=50)&(x3<=50)]=5
  set123[(x1>50)&(x2<=50)&(x3>50)]=6
  set123[(x1>50)&(x2>50)&(x3<=50)]=7
  set123[(x1>50)&(x2>50)&(x3>50)]=8
  
  set124=numeric(nu)
  set124[(x1<=50)&(x2<=50)&(x4<=50)]=1
  set124[(x1<=50)&(x2<=50)&(x4>50)]=2
  set124[(x1<=50)&(x2>50)&(x4<=50)]=3
  set124[(x1<=50)&(x2>50)&(x4>50)]=4
  set124[(x1>50)&(x2<=50)&(x4<=50)]=5
  set124[(x1>50)&(x2<=50)&(x4>50)]=6
  set124[(x1>50)&(x2>50)&(x4<=50)]=7
  set124[(x1>50)&(x2>50)&(x4>50)]=8
  
  set423=numeric(nu)
  set423[(x4<=50)&(x2<=50)&(x3<=50)]=1
  set423[(x4<=50)&(x2<=50)&(x3>50)]=2
  set423[(x4<=50)&(x2>50)&(x3<=50)]=3
  set423[(x4<=50)&(x2>50)&(x3>50)]=4
  set423[(x4>50)&(x2<=50)&(x3<=50)]=5
  set423[(x4>50)&(x2<=50)&(x3>50)]=6
  set423[(x4>50)&(x2>50)&(x3<=50)]=7
  set423[(x4>50)&(x2>50)&(x3>50)]=8
  
  
  set134=numeric(nu)
  set134[(x1<=50)&(x3<=50)&(x4<=50)]=1
  set134[(x1<=50)&(x3<=50)&(x4>50)]=2
  set134[(x1<=50)&(x3>50)&(x4<=50)]=3
  set134[(x1<=50)&(x3>50)&(x4>50)]=4
  set134[(x1>50)&(x3<=50)&(x4<=50)]=5
  set134[(x1>50)&(x3<=50)&(x4>50)]=6
  set134[(x1>50)&(x3>50)&(x4<=50)]=7
  set134[(x1>50)&(x3>50)&(x4>50)]=8
  
  
  set1234=numeric(nu)
  set1234[(x3<=50)&(x1<=50)&(x2<=50)&(x4<=50)]=1
  set1234[(x3<=50)&(x1<=50)&(x2<=50)&(x4>50)]=2
  set1234[(x3<=50)&(x1<=50)&(x2>50)&(x4<=50)]=3
  set1234[(x3<=50)&(x1<=50)&(x2>50)&(x4>50)]=4
  set1234[(x3<=50)&(x1>50)&(x2<=50)&(x4<=50)]=5
  set1234[(x3<=50)&(x1>50)&(x2<=50)&(x4>50)]=6
  set1234[(x3<=50)&(x1>50)&(x2>50)&(x4<=50)]=7
  set1234[(x3<=50)&(x1>50)&(x2>50)&(x4>50)]=8
  set1234[(x3>50)&(x1<=50)&(x2<=50)&(x4<=50)]=9
  set1234[(x3>50)&(x1<=50)&(x2<=50)&(x4>50)]=10
  set1234[(x3>50)&(x1<=50)&(x2>50)&(x4<=50)]=11
  set1234[(x3>50)&(x1<=50)&(x2>50)&(x4>50)]=12
  set1234[(x3>50)&(x1>50)&(x2<=50)&(x4<=50)]=13
  set1234[(x3>50)&(x1>50)&(x2<=50)&(x4>50)]=14
  set1234[(x3>50)&(x1>50)&(x2>50)&(x4<=50)]=15
  set1234[(x3>50)&(x1>50)&(x2>50)&(x4>50)]=16
  return(list(set2=set2,set4=set4,set123=set123,set124=set124,set423=set423,set134=set134,set1234=set1234))
  }
  
  # base simulation scheme
  data_generation1=function(n=n,nu=nu,dist=dist,nfold=nfold)	
  {
  
  
  xx = mapply(function(i){sample(1:100, size=nu, replace = TRUE) },i=1:10)
  
  x1 = xx[,1]
  x2 = xx[,2]
  x3 = xx[,3]
  x4 = xx[,4]
  
  xbeta = 0.01*(x1+2*x2-x3+2*sqrt(x1*x3)-sqrt(x2*x4))
  rbeta = abs(2*x1-x2+sqrt(x1*x2))/102

  # - fix variance -
if (dist==1){
  
  lmd = 1.6487 # mean parameter in the exp dist
  mu = exp(xbeta)+lmd
  yy=mapply(function(i){exp(xbeta)+rexp(n=nu,rate=1/lmd)},i=1:(n+1))	
  
}

  if (dist==2){
  
  delta = 1 # parameter in the lognormal distribution
  mu = exp(xbeta)+exp(delta^2/2)
  yy=mapply(function(i){exp(xbeta)+exp(delta*rnorm(nu,0,1))},i=1:(n+1))	
  
  }
  
  
  if (dist==3){
  
  theta = 3.2974 # parameter in the pareto dist
  mu = exp(xbeta)+theta/2
  yy=mapply(function(i){exp(xbeta)+((1-runif(nu,0,1))^(-1/3)-1)*theta},i=1:(n+1))
  
  }
  
 
  
  wmean = rowMeans(yy[,1:n]) # individual mean
  wvar  = rowSums((yy[,1:n]-wmean)^2)/(n-1) # individual variance
  
  parms = c(wmean,wvar,n,nu) # a stack of parameters to be put in rpart
  
  dat = data.frame(y=wmean,xx=xx)
  
  fold=sample(rep(1:nfold, length = n), n,  replace = FALSE)
  
  return(list(dat=dat,yy=yy[,1:n],xx = xx, fold=fold,parms=parms,mu=mu,wmean=wmean,wvar=wvar,ynext= yy[,n+1]))
  }
  
  
  # random variance
  data_generation2=function(n=n,nu=nu,dist=dist,nfold=nfold)	
  {
  
  
  xx = mapply(function(i){sample(1:100, size=nu, replace = TRUE) },i=1:10)
  
  x1 = xx[,1]
  x2 = xx[,2]
  x3 = xx[,3]
  x4 = xx[,4]
  
  xbeta = 0.01*(x1+2*x2-x3+2*sqrt(x1*x3)-sqrt(x2*x4))
  rbeta = abs(2*x1-x2+sqrt(x1*x2))/102
 
  
  # - random variance -
  

  if (dist==1){
  
lmd = exp(rbeta/2) # mean parameter in the exp dist
mu = exp(xbeta)+lmd
yy=mapply(function(i){exp(xbeta)+rexp(n=nu,rate=1/lmd)},i=1:(n+1))	

  }


  if (dist==2){
  
  mu = exp(xbeta)+exp(rbeta/2)
  yy=mapply(function(i){exp(xbeta)+exp(sqrt(rbeta)*rnorm(nu,0,1))},i=1:(n+1))
  
  }
  

  
  
  if (dist==3){
  
  theta = 2*exp(rbeta/2) # parameter in the pareto dist
  mu = exp(xbeta)+theta/2
  yy=mapply(function(i){exp(xbeta)+((1-runif(nu,0,1))^(-1/3)-1)*theta},i=1:(n+1))
  
  }
  
  
  wmean = rowMeans(yy[,1:n]) # individual mean
  wvar  = rowSums((yy[,1:n]-wmean)^2)/(n-1) # individual variance
  
  parms = c(wmean,wvar,n,nu) # a stack of parameters to be put in rpart
  
  dat = data.frame(y=wmean,xx=xx)
  
  fold=sample(rep(1:nfold, length = n), n,  replace = FALSE)
  
  return(list(dat=dat,yy=yy[,1:n],xx = xx, fold=fold,parms=parms,mu=mu,wmean=wmean,wvar=wvar,ynext= ynext))
  }


  # re1
  data_generation3=function(n=n,nu=nu,dist=dist,nfold=nfold)	
  {
  
  
  xx = mapply(function(i){sample(1:100, size=nu, replace = TRUE) },i=1:10)
  
  x1 = xx[,1]
  x2 = xx[,2]
  x3 = xx[,3]
  x4 = xx[,4]
  
  xbeta = 0.01*(x1+2*x2-x3+2*sqrt(x1*x3)-sqrt(x2*x4))
  rbeta = abs(2*x1-x2+sqrt(x1*x2))/102
  random_effect=runif(n=nu,min=0.9,max=1.1)
  
  
  
  
  if (dist==1){
  
lmd = exp(rbeta/2) # mean parameter in the exp dist
mu = (exp(xbeta)+lmd)*random_effect
yy=mapply(function(i){exp(xbeta)+rexp(n=nu,rate=1/lmd)},i=1:(n+1))	
yy=yy*random_effect

}

  if (dist==2){
  
  mu = (exp(xbeta)+exp(rbeta/2))*random_effect
  yy=mapply(function(i){exp(xbeta)+exp(sqrt(rbeta)*rnorm(nu,0,1))},i=1:(n+1))	
  yy=yy*random_effect
  
  }
  

  
  if (dist==3){
  
  theta = 2*exp(rbeta/2) # parameter in the pareto dist
  mu = (exp(xbeta)+theta/2)*random_effect
  yy=mapply(function(i){exp(xbeta)+((1-runif(nu,0,1))^(-1/3)-1)*theta},i=1:(n+1))
  yy=yy*random_effect
  
  }
  
  wmean = rowMeans(yy[,1:n]) # individual mean
  wvar  = rowSums((yy[,1:n]-wmean)^2)/(n-1) # individual variance
  
  parms = c(wmean,wvar,n,nu) # a stack of parameters to be put in rpart
  
  dat = data.frame(y=wmean,xx=xx)
  
  fold=sample(rep(1:nfold, length = n), n,  replace = FALSE)
  
  return(list(dat=dat,yy=yy[,1:n],xx = xx, fold=fold,parms=parms,mu=mu,wmean=wmean,wvar=wvar,ynext= yy[,n+1]))
  }
  
  
  # re2
  data_generation4=function(n=n,nu=nu,dist=dist,nfold=nfold)	
  {
  
  
  xx = mapply(function(i){sample(1:100, size=nu, replace = TRUE) },i=1:10)
  
  x1 = xx[,1]
  x2 = xx[,2]
  x3 = xx[,3]
  x4 = xx[,4]
  
  random_effect1=runif(n=nu,min=0.9,max=1.1)
  random_effect2=runif(n=nu,min=0.9,max=1.1)
  xbeta = 0.01*(x1+2*x2-x3+2*sqrt(x1*x3)-sqrt(x2*x4))*random_effect1
  rbeta = abs(2*x1-x2+sqrt(x1*x2))/102*random_effect2
  
  
  
  if (dist==2){
  
  mu = (exp(xbeta)+exp(rbeta/2))
  yy=mapply(function(i){exp(xbeta)+exp(sqrt(rbeta)*rnorm(nu,0,1))},i=1:(n+1))	
  
  }
  
  if (dist==1){
  
  lmd = exp(rbeta/2) # mean parameter in the exp dist
  mu = (exp(xbeta)+lmd)
  yy=mapply(function(i){exp(xbeta)+rexp(n=nu,rate=1/lmd)},i=1:(n+1))	
  
  }
  
  
  if (dist==3){
  
  theta = 2*exp(rbeta/2) # parameter in the pareto dist
  mu = (exp(xbeta)+theta/2)
  yy=mapply(function(i){exp(xbeta)+((1-runif(nu,0,1))^(-1/3)-1)*theta},i=1:(n+1))
  
  }
  
  
  wmean = rowMeans(yy[,1:n]) # individual mean
  wvar  = rowSums((yy[,1:n]-wmean)^2)/(n-1) # individual variance
  
  parms = c(wmean,wvar,n,nu) # a stack of parameters to be put in rpart
  
  dat = data.frame(y=wmean,xx=xx)
  
  fold=sample(rep(1:nfold, length = n), n,  replace = FALSE)
  
  return(list(dat=dat,yy=yy[,1:n],xx = xx, fold=fold,parms=parms,mu=mu,wmean=wmean,wvar=wvar,ynext= yy[,n+1]))
  }
  # calculating prediction error
  PE=function(sset,mu,wmean,wvar,nu,ynext)
  {
  nset = unique(sset) # terminal nodes in tree fittj
  
  pe1 = numeric(nu)
  pe2 = numeric(nu)
  
  for (k in nset)
  {
  wvark  = wvar[sset==k]
  wmeank = wmean[sset==k]
  
  sigmak = mean(wvark)
  tauk = var(wmeank)-sigmak/n
  tauk = max(0,tauk)
  alphak = n*tauk/(n*tauk+sigmak)
  
  pred_k = alphak * wmeank + (1-alphak)* mean(wmeank)
  
  pe1[sset==k] = pe1[sset==k] + rowSums(as.matrix(mu[sset==k] - pred_k)^2)
  pe2[sset==k] = pe2[sset==k] + rowSums(as.matrix(ynext[sset==k] - pred_k)^2)
  
  }
  
  return(list(pe1 = mean(pe1), pe2 = mean(pe2)))
  }
  
  
  test = function(fit,mu,wmean,wvar,nu,ynext)
  {
  
  tem = PE(fit$where,mu,wmean,wvar,nu,ynext)
  pe1 = tem$pe1
  pe2 = tem$pe2
  
  size = length(unique(fit$where))
  
  noise=sum(unique(fit$frame$var)!="xx.1"&unique(fit$frame$var)!="xx.2"&unique(fit$frame$var)!="xx.3"&unique(fit$frame$var)!="xx.4")-1
  
  inflv=sum(c(unique(fit$frame$var)=="xx.1",unique(fit$frame$var)=="xx.2", unique(fit$frame$var)=="xx.3",unique(fit$frame$var)=="xx.4"))
  
  return(c(pe1,pe2,size,noise,inflv))
  
  }
  
  result = function(data_list){
  
  dat = data_list$dat
  yy = data_list$yy
  xx = data_list$xx
  fold = data_list$fold
  parms = data_list$parms
  mu = data_list$mu
  wmean = data_list$wmean
  wvar = data_list$wvar
  ynext = data_list$ynext
  
  # ---------  Tree Building -----------------------
  
  # L2 Tree
  fitt1 = rpart(y~.,method="anova",maxsurrogate=0,cp=0.0,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  CP1 = longi_cv(fitt=fitt1,mlist="anova",dat=dat,yy=yy,xx=xx,n=n,nfold=nfold,fold=fold)
  fit1 = rpart(y~.,method="anova",maxsurrogate=0,cp=CP1,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  
  
  # Credibility Trees - Tree one
  fitt2 = rpart(y~.,method=clist1,parms=parms,maxsurrogate=0,cp=0.0,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  CP2 = longi_cv(fitt=fitt2,mlist=clist1,dat=dat,yy=yy,xx=xx,n=n,nfold=nfold,fold=fold)
  fit2 = rpart(y~.,method=clist1,parms=parms,maxsurrogate=0,cp=CP2,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  
  
  # Credibility Trees - Tree two
  fitt3 = rpart(y~.,method=clist2,parms=parms,maxsurrogate=0,cp=0.0,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  CP3 = longi_cv(fitt=fitt3,mlist=clist2,dat=dat,yy=yy,xx=xx,n=n,nfold=nfold,fold=fold)
  fit3 = rpart(y~.,method=clist2,parms=parms,maxsurrogate=0,cp=CP3,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  
  
  # Credibility Trees - Tree three
  fitt4 = rpart(y~.,method=clist3,parms=parms,maxsurrogate=0,cp=0.0,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  CP4 = longi_cv(fitt=fitt4,mlist=clist3,dat=dat,yy=yy,xx=xx,n=n,nfold=nfold,fold=fold)
  fit4 = rpart(y~.,method=clist3,parms=parms,maxsurrogate=0,cp=CP4,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  
  
  # Credibility Trees - Tree four
  fitt5 = rpart(y~.,method=clist4,parms=parms,maxsurrogate=0,cp=0.0,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  CP5 = longi_cv(fitt=fitt5,mlist=clist4,dat=dat,yy=yy,xx=xx,n=n,nfold=nfold,fold=fold)
  fit5 = rpart(y~.,method=clist4,parms=parms,maxsurrogate=0,cp=CP5,usesurrogate=0,xval=0,minbucket=20,minsplit=7,data=dat)
  
  list_set=generate_naive_sets(dat$xx.1,dat$xx.2,dat$xx.3,dat$xx.4,dat$xx.5,nu)
  
  # no subgroups - all together
  res0=PE(rep(0,nu),mu,wmean,wvar,nu,ynext)
  resn2=PE(list_set$set2,mu,wmean,wvar,nu,ynext)
  resn4=PE(list_set$set4,mu,wmean,wvar,nu,ynext)
  
  # naive credibility error without using tree
  res123=PE(list_set$set123,mu,wmean,wvar,nu,ynext)
  res124=PE(list_set$set124,mu,wmean,wvar,nu,ynext)
  res423=PE(list_set$set423,mu,wmean,wvar,nu,ynext)
  res134=PE(list_set$set134,mu,wmean,wvar,nu,ynext)
  res1234=PE(list_set$set1234,mu,wmean,wvar,nu,ynext)
  
  temp_pe1 = c(res0$pe1,resn2$pe1,resn4$pe1,res123$pe1,res124$pe1,res423$pe1,res134$pe1,res1234$pe1)
  temp_pe2 = c(res0$pe2,resn2$pe2,resn4$pe2,res123$pe2,res124$pe2,res423$pe2,res134$pe2,res1234$pe2)
  temp_size=c(0,2,2,8,8,8,8,16)
  temp_noise=c(999,0,0,0,0,0,0,0)
  temp_inflv=c(999,1,1,3,3,3,3,4)
  
  res1=test(fit1,mu,wmean,wvar,nu,ynext)
  res2=test(fit2,mu,wmean,wvar,nu,ynext)
  res3=test(fit3,mu,wmean,wvar,nu,ynext)
  res4=test(fit4,mu,wmean,wvar,nu,ynext)
  res5=test(fit5,mu,wmean,wvar,nu,ynext)
  temp=c(res1,res2,res3,res4,res5)
  
  return(as.vector(cbind(rbind(temp_pe1,temp_pe2,temp_size,temp_noise,temp_inflv),matrix(temp,5,))))
  
  }
  
  # user-written functions
  
  clist1 = list(init = init, split = cred1_s, eval = cred1_e) # credibility loss one
  clist2 = list(init = init, split = cred2_s, eval = cred2_e) # credibility loss two
  clist3 = list(init = init, split = cred3_s, eval = cred3_e) # credibility loss three
  clist4 = list(init = init, split = cred4_s, eval = cred4_e) # credibility loss four
  
  

for (runtime in 1:nSim){


# Data Simulation

data_list1=data_generation1(n=n,nu=nu,dist=dist,nfold=nfold)
data_list2=data_generation2(n=n,nu=nu,dist=dist,nfold=nfold)
data_list3=data_generation3(n=n,nu=nu,dist=dist,nfold=nfold)
data_list4=data_generation4(n=n,nu=nu,dist=dist,nfold=nfold)

re1 = result(data_list1)
re2 = result(data_list2)
re3 = result(data_list3)
re4 = result(data_list4)

temp = c(re1,re2,re3,re4)

cat(temp, sep=" ", "\n", append=T, file=file.nam)

print(runtime)}

