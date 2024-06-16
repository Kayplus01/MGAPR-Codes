#R codes for Monte Carlo Simulation for MGAPR Distribution ########
####Defining the MGAPR quantile and negative log-likelihood #####
##### Specifying MGAPR Quantile Function ########

rm(list=ls())
library(stats4)
library(bbmle)
library(stats)
library(numDeriv)
library(Matrix)
library(rootSolve)
library(gsl)
set.seed(123)


### MGAPR CDF Graph ###   
CDF <- function(α, δ, λ, z){
  A <- (1+δ)*(1-exp(-(z^2)/2*(λ^2)))
  B <- (δ*α^(-exp(-(z^2)/2*(λ^2)))) + 1 -exp(-(z^2)/2*(λ^2))
  cdf<- A/B
  return(cdf)}


### MGAPR PDF Graph ###   
PDF <- function(α, δ, λ, z){
  C <- exp(-(z^2)/(2*(λ^2)))
  D <- z*δ*(1+δ)*C*(1-log(α)*(1-C))*α^(-C)
  E <- (λ*(δ*α^(-C) + 1- C))^2
  pdf<- D/E
  return(pdf)}


### MGAPR Qf Graph ###
Qf <- function(α, δ, λ, lambda){
  f <- function(z){
    CDF(α, δ, λ, z)-lambda}
  z = min(uniroot.all(f,lower=0,upper=1000000000000,tol=0.0001))
  return(z)}


#### Negative log-likelihood of MGAPR ####
NLL_pdf <- function(par){
  C <- exp(-(z^2)/(2*(par[3]^2)))
  D <- z*par[2]*(1+par[2])*C*(1-log(par[1])*(1-C))*par[1]^(-C)
  E <- (par[3]*(par[2]*par[1]^(-C) + 1- C))^2
  pdf <- D/E
  NLL <- -sum(log(pdf))
  return(NLL)}


### Generating a Random Sample of Size n from the MGAPR ####
RS_MGAPR = function(par,n){
  α = par[1]
  δ = par[2]
  λ = par[3]  
  z <- c(rep(0,n))
  lambda <- 0
  lambda <- runif(n, min=0, max=1)
  
  for (i in 1:n){
    z[i] = Qf(α, δ, λ, lambda[i])}
  return(z)}

### Monte Carlo Simulation of MGAPR ###
par <- c(α=δ, δ=δ, λ=λ)   
n_replicas = 2000 
matriz_par <- matrix(0, 10, 3)
matriz_bias <- matrix(0, 10, 3)
matriz_MSE <- matrix(0, 10, 3)
matriz_RMSE <- matrix(0, 10, 3)
matriz_std <- matrix(0, 10, 3)

colnames(matriz_par) <- c("α", "δ", "λ")
colnames(matriz_bias) <- c("α", "δ", "λ")
colnames(matriz_MSE) <- c("α", "δ", "λ")
colnames(matriz_RMSE) <- c("α", "δ", "λ")
colnames(matriz_std) <- c("α", "δ", "λ")
cont = 1
n = 50
while(n <= 500){
  par_mean <- c(0,0,0)
  std_mean <- c(0,0,0)
  bias <- c(0,0,0)
  MSE <- c(0,0,0)
  replica = 1
  while(replica <= n_replicas){
    print(paste("n=",n,",replica=",replica))
    z <- RS_MGAPR(par,n)
    Data <- z
    α = α
    δ = δ
    λ = λ
    
    #### Optimization and Generating the Simulation results ####
    result <- nlminb(c(α, δ, λ), NLL_pdf, lower=0, upper=Inf)
    if (class(result) !="try-error" && result$convergence==0){
      par_mean <- par_mean+result$par
      bias = bias+(result$par-par)
      MSE = MSE+(result$par-par)^2
      replica = replica+1}}
  
  par_mean = par_mean/n_replicas
  bias = bias/n_replicas
  MSE = MSE/n_replicas
  RMSE = sqrt(MSE)
  
  matriz_par[cont,] =par_mean
  matriz_std[cont,] =std_mean
  matriz_bias[cont,] =bias
  matriz_MSE[cont,] =MSE
  matriz_RMSE[cont,] =RMSE
  
  print("mean= ")
  print(par_mean)
  print("bias= ")
  print(bias)
  print("MSE= ")
  print(MSE)
  print("RMSE= ")
  print(RMSE)
  n =n+50
  cont = cont+1}

print(matriz_par)
print(matriz_bias)
print(matriz_MSE)
print(matriz_RMSE)

n=seq(50, 500, by=50)
########### END ############