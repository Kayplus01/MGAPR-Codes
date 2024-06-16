#required packages
library(AdequacyModel)
library(gsl)
library(ggplot2)
library(numDeriv)
library(bbmle)
library(fBasics)
library(zoo)
library(lmtest)
library(Newdistns)
library(EnvStats)
library(lmomco)
par(mfrow=c(1,1))

### Data_Set ####
#### information on the remission times, measured in months, for a group of 128 
Remission_time <- c(3.88, 5.32, 7.39, 10.34, 14.83, 34.26, 0.90, 2.69, 4.18, 5.34, 7.59, 
                 10.66, 15.96, 36.66, 1.05, 2.69, 4.23, 5.41, 7.62, 10.75, 16.62, 43.01, 1.19, 2.75, 4.26, 
                 5.41, 7.63, 17.12, 46.12, 1.26, 2.83, 4.33, 5.49, 7.66, 11.25, 17.14, 79.05, 1.35, 2.87, 
                 5.62, 7.87, 11.64, 17.36, 1.40, 3.02, 4.34, 5.71, 7.93, 0.08, 2.09, 3.48, 4.87, 6.94, 8.66, 
                 13.11, 23.63, 0.20, 2.23, 3.5, 4.98, 6.97, 9.02, 13.29, 0.40, 2.26, 3.57, 5.06, 7.09, 9.22, 
                 13.80, 25.74, 0.50, 2.46, 3.64, 5.09, 7.26, 9.47, 14.24, 25.82, 0.51, 2.54, 3.70, 5.17, 
                 7.28, 9.74, 14.76, 26.31, 0.81, 2.62, 3.82, 5.32, 7.32, 10.06, 14.77, 32.15, 2.64, 11.79, 
                 18.10, 1.46, 4.40, 5.85, 8.26, 11.98, 19.13, 1.76, 3.25, 4.50, 6.25, 8.37, 12.02, 2.02, 
                 3.31, 4.51, 6.54, 8.53, 12.03, 20.28, 2.02, 3.36, 6.76, 12.07, 21.73, 2.00, 3.36, 6.93, 
                 8.65, 12.63, 22.69)


#The data gives 100 observations on breaking stress of carbon fibres (in Gba):  
carbon_fibres <- c(3.70, 2.74, 2.73, 2.50, 3.6, 3.11, 3.27, 2.87, 1.47, 3.11, 4.42, 2.41, 3.19, 3.22, 1.69, 3.28, 3.09, 1.87, 
        3.15, 4.9, 3.75, 2.43, 2.95, 2.97, 3.39, 2.96, 2.53, 2.67, 2.93, 3.22, 3.39, 2.81, 4.20, 3.33, 2.55, 3.31, 
        3.31, 2.85, 2.56, 3.56, 3.15, 2.35, 2.55, 2.59, 2.38, 2.81, 2.77, 2.17, 2.83, 1.92, 1.41, 3.68, 2.97, 1.36, 
        0.98, 2.76, 4.91, 3.68, 1.84, 1.59, 3.19, 1.57, 0.81, 5.56, 1.73, 1.59, 2.00, 1.22, 1.12, 1.71, 2.17, 1.17, 
        5.08, 2.48, 1.18, 3.51, 2.17, 1.69, 1.25, 4.38, 1.84, 0.39, 3.68, 2.48, 0.85, 1.61, 2.79, 4.70, 2.03, 1.80, 
        1.57, 1.08, 2.03, 1.61, 2.12, 1.89, 2.88 ,2.82, 2.05, 3.65)

       
x=Remission_time
z=x

###TTT plots###
TTT(z, col = "blue", lwd = 2, grid = TRUE, lty = 2)

### Histogram ###
hist(z,col="blue",main="")

### Boxplot ###
boxplot(z, col = "blue")


### Descriptive statistics ###
summaryFull(x)

### Create the function for the mode ###
getmode <- function(x){
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}
mode <- getmode(x)
print(mode)


##############################################################################

### MGAPR CDF Graph ### 
### CDF ###
x=z
CDF<-function(par,x){
  α <-par[1]
  δ <-par[2]
  λ <-par[3]
  K <- α*(1+δ)*(1-exp(-(x^2)/(2*(λ^2))))
  L <- (δ*(α^(1-exp(-(x^2))/(2*(λ^2))))) + α*(1 -exp(-(x^2)/(2*(λ^2))))
  c   <- K/L
  return(c)}
 


### PDF of the MGAPR Distribution ###
PDF<-function(par,x){ 
  α <-par[1]
  δ <-par[2]
  λ<-par[3]
  A <- exp(-(x^2)/(2*(λ^2)))
  B <- x*δ*(1+δ)*A*(1-(log(α)*(1-A)))*(α^(2-A))
  C <- (λ*((δ*(α^(1-A))) + (α*(1- A))))^2
  p <- B/C
  return(p)}


####################################################################################
### Other models ### 
### GAPR cdf Distribution ### 
cdf_GAPR <-function(pars,x){
  η <-pars[1]
  λ <-pars[2]
  G <- η*(1- exp(-(x^2)/(2*(λ^2))))
  H <- η^(1- exp(-(x^2)/(2*(λ^2))))
  d <- G/H
  return(d)}


### GAPR pdf Distribution ###
pdf_GAPR<-function(pars,x){ 
  η <-pars[1]
  λ <-pars[2]
  I <- 1-exp(-(x^2)/(2*(λ^2)))
  J <- ((η*x)*(exp(-(x^2)/(2*(λ^2)))))/((λ^2)*(η^I))
  K <- 1 - (I*log(η))
  e <- J *K
  return(e)}


### APGIR cdf Distribution ### 
cdf_APGIR <-function(pars,x){
  α <- pars[1]
  β <- pars[2]
  λ <- pars[3]
  L <- (1- (exp(-(λ*x)^(-2))))^β
  M <- (α^(1-L)) - 1
  f <- M/(α-1)
  return(f)}


### APGIR pdf Distribution ###
pdf_APGIR <-function(pars,x){ 
  α <- pars[1]
  β <- pars[2]
  λ <- pars[3]
  N <- (log(α)*2*β)/((α-1)*(λ^2)*(x^3))*(exp(-(λ*x)^(-2)))
  O <- (1- (exp(-(λ*x)^(-2))))^(β-1)
  P <- (α^(1-((1- (exp(-(λ*x)^(-2))))^β)))
  g <- N*O*P
  return(g)}


### EIR cdf ###
cdf_EIR <- function(parss,x){
  α <-parss[1]
  θ <-parss[2]
  A <- (1-exp(-(θ/x)^2))^α
  i<- 1- A
  return(i)}

### EIR pdf ###
pdf_EIR <- function(parss,x){
  α <-parss[1]
  θ <-parss[2]
  A <- ((2*α*θ^2)/(x^3))*exp(-(θ/x)^2)
  B <- (1-exp(-(θ/x)^2))^(α-1)
  j <- A*B
  return(j)}


### ERayleigh cdf ###
cdf_ERay <- function(pars,x){
  μ <-pars[1]
  λ <-pars[2]
  A <- (1-exp(-(x^2)/(2*(λ^2))))^μ}

### ERayleigh pdf ###
pdf_ERay <- function(pars,x){
  μ <-pars[1]
  λ<-pars[2]
  A<- ((x*μ)/(λ^2))*(exp(-(x^2)/(2*(λ^2))))
  B <- (1-exp(-(x^2)/(2*(λ^2))))^(μ-1)
  t <- A*B
  return(t)}


### Rayleigh cdf ###
cdf_Ray <- function(pars,x){
  λ<-pars[1]
  A <- 1-exp(-(x^2)/(2*(λ^2)))}

### Rayleigh pdf ###
pdf_Ray <- function(pars,x){
  λ<-pars[1]
  A<- (x/(λ^2))*exp(-(x^2)/(2*(λ^2)))}

#--------------------------END-----------------------------------------------#


### Goodness of fit statistic ###
#### MGAPR Distribution ####
x=sort(z)
set.seed(123)
result = goodness.fit(pdf=PDF, cdf=CDF, 
                      starts = c(α=α, δ=δ, λ=λ), 
                      data = x, method = "B",domain = c(0,Inf),mle = NULL)

result


#--------------------------END-----------------------------------------------#

### Goodness of fit statistic of Other models###
#### APGIR ###
result1 = goodness.fit(pdf=pdf_APGIR, cdf=cdf_APGIR, 
                      starts = c(α=α, β=β, λ=λ), 
                      data = x, method = "S",domain = c(0,Inf),mle = NULL)

result1


#### GAPR ####
result2 = goodness.fit(pdf=pdf_GAPR, cdf=cdf_GAPR, 
                       starts = c(η=η, λ=λ), 
                       data = x, method = "S",domain = c(0,Inf),mle = NULL)

result2


### EIR ###  
result3 = goodness.fit(pdf=pdf_EIR, cdf=cdf_EIR, starts = c(α=α, θ=θ),
              data = x, method = "S",domain = c(0,Inf),mle = NULL)
result3



### ERayleigh ###
result4 = goodness.fit(pdf=pdf_ERay,cdf=cdf_ERay, starts = c(μ=μ  ,λ=λ),
                       data = x, method = "S",domain = c(0,Inf),mle = NULL)
result4



### Rayleigh ###
result5 = goodness.fit(pdf=pdf_Ray,cdf=cdf_Ray, starts = c(λ=λ),
                       data = x, method = "S",domain = c(0,Inf),mle = NULL)
result5

#--------------------------END-----------------------------------------------#


### PLOTS of fitted densities ###
x <- sort(x)
x
hist(x,col="gray", freq=FALSE,  
     ylim=c(0,.11), 
     ylab="pdf",xlab="z", main="")
lines(x, PDF(x, par = result$mle), col = "blue", lwd = 2, lty = 1)
legend("topright", c("MGAPR"), col=c("blue"), lwd=2, ncol = 1, cex = 1)


hist(x,col="white", freq=FALSE,   
     ylim=c(0,.11), ylab="Density",xlab="z", main="")



### Add fit lines####
lines(x, PDF(x, par = result$mle), col = "blue", lwd = 2, lty = 1)
lines(x, pdf_APGIR(x, par = result1$mle), col = "cyan", lwd = 1, lty = 2)
lines(x, pdf_GAPR(x, par = result2$mle), col = "red", lwd = 1, lty = 2)
lines(x, pdf_EIR(x, par = result3$mle), col = "green", lwd = 1, lty = 2)
lines(x, pdf_ERay(x, par = result4$mle), col = "purple", lwd = 1, lty = 2)
lines(x, pdf_Ray(x, par = result5$mle), col = "brown", lwd = 1, lty = 2)
legend("topright", c("MGAPR", "APGIR", "GAPR", "EIR", "ER", "R"), 
       col=c("blue", "cyan", "red", "green", "purple", "brown"), lty=c(1,2,2,2,2,2), 
       lwd=2, ncol = 1, cex = 1)

