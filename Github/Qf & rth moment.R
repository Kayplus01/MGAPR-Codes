##----The MGAPR  ----------------------------------------------------## 
rm(list=ls())
library(stats4)
library(bbmle)
library(stats)
library(numDeriv)
library(rootSolve)

### MGAPR CDF Graph ###   
CDF <- function (α, δ, λ, z){
  A <- (1+δ)*(1-exp(-(z^2)/2*(λ^2)))
  B <- (δ*α^(-exp(-(z^2)/2*(λ^2)))) + 1 -exp(-(z^2)/2*(λ^2))
  cdf<- A/B
  return(cdf)}



### MGAPR PDF Graph ###   
PDF <- function (α, δ,λ, z){
  A <- exp(-(z^2)/(2*(λ^2)))
  B <- z*δ*(1+δ)*A*(1-log(α)*(1-A))*α^(-A)
  C <- (λ^2)*(δ*α^(-A) + 1- A )^2
  pdf<- B/C
  return(pdf)}


#-----------------------------------------------------------------#
MGAPR_quantile=function(parameter,p){
  α = parameter[1]
  δ = parameter[2]
  λ = parameter[3]
  f=function(z){
    CDF(α, δ, λ, z) - p
  }
  z=min(uniroot.all(f,lower=0,upper=1000000000000,tol=0.0001))
  return(z)
}


MGAPR_QuantileTable=function(parameter_matrix){
  p=seq(0.1,0.9,0.1)
  size=dim(parameter_matrix)[1]
  Table_Quantile=matrix(NA,nrow=length(p),ncol=size)
  row.names(Table_Quantile)=p
  colnames(Table_Quantile)=apply(parameter_matrix,1,function(z){paste0('(',paste0(z,collapse=','),")")})
  Table_Quantile
  for(iter in 1:size){
    parameter=parameter_matrix[iter,]
    for(i in 1:length(p)){
      Table_Quantile[i,iter] = MGAPR_quantile(parameter,p[i])
    }
  }
  return(Table_Quantile)
}

#--------------------------------------------------------------------------
## table of quantile
parameter_matrix=as.matrix(rbind(
  par1 = c(0.2, 7.2, 0.9),
  par2 = c(0.9, 2.4, 0.6),
  par3 = c(0.5, 0.09, 0.2),
  par4 = c(0.7, 1.5, 0.6),
  par5 = c(0.6, 1.5, 1.2)
))

MGAPR_QuantileTable(parameter_matrix)
parameter_matrix[,c(0.2, 7.2, 0.9)]

plot(n,(matriz_par[,1]),type="l", col="green", lty=1, lwd=2, xlab="Sample Size",
     ylab="Parameter Estimate", ylim=c(0,1.5))
lines(n, (matriz_par[,2]), col="blue", lty=5, lwd=2, type="l")
lines(n, (matriz_par[,3]), col="red", lty=8, lwd=2, type="l")
title("")
legend('topright', c(expression(
  hat(α),
  hat(δ),
  hat(λ),
)),lty=c(1,5,8), cex=1, lwd=2, col=c("green", "blue", "red"),
ncol = 2)
#------------------------------------------------------------------------#


MGAPR_moments = function(parameter, r){
  α = parameter[1]
  δ = parameter[2]
  λ = parameter[3]
  r = r
  
  f=function(z){
    (z^r)*(PDF(α, δ, λ, z))
  }
  Moments = integrate(f,lower=0,upper=1,subdivisions = 10000000)
  return(Moments$value)
}


MGAPR_momentsTable = function(parameter_matrix,r){
  list = c(paste0('Moment', 1 : r),'SD','CV','CS','CK')
  Table_Moments = matrix(NA,nrow=length(list),ncol=dim(parameter_matrix)[1])
  row.names(Table_Moments) = list
  colnames(Table_Moments) = apply(parameter_matrix, 1, function(z){paste0('(',paste0(z,collapse=','),")")})
  Table_Moments
  for(iter in 1:dim(parameter_matrix)[1]){
    parameter=parameter_matrix[iter,]
    for(i in 1: r){
      Table_Moments[i,iter] = MGAPR_moments(parameter,i)
    }
    Table_Moments['SD',iter] = sqrt(Table_Moments['Moment2',iter]-Table_Moments['Moment1',iter]^2)
    Table_Moments['CV',iter] = sqrt(Table_Moments['Moment2',iter]/Table_Moments['Moment1',iter]^2-1)
    Table_Moments['CS',iter] = (Table_Moments['Moment3',iter]-3*Table_Moments['Moment1',iter]*
                                  Table_Moments['Moment2',iter]+2*Table_Moments['Moment1',iter]^3)/
      ((Table_Moments['Moment2',iter]-Table_Moments['Moment1',iter]^2)^(3/2))
    Table_Moments['CK',iter] = (Table_Moments['Moment4',iter]-4*Table_Moments['Moment1',iter]*
                                  Table_Moments['Moment3',iter]+6*Table_Moments['Moment1',iter]^2*Table_Moments['Moment2',iter]-
                                  3*Table_Moments['Moment1',iter]^4)/((Table_Moments['Moment2',iter]-Table_Moments['Moment1',iter]^2)^2)
  }
  return(Table_Moments)
}


#--------------------------------------------------------------------------#
## Moments Table
parameter_matrix=as.matrix(rbind(
  par1 = c(0.2, 7.2, 0.9),
  par2 = c(0.9, 2.4, 0.6),
  par3 = c(0.5, 0.09, 0.2),
  par4 = c(0.7, 1.5, 0.6),
  par5 = c(0.6, 1.5, 1.2)
))
r = 5
MGAPR_momentsTable(parameter_matrix, r)
#----------------------------------------END------------------------------#




