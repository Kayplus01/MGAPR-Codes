### Required Packages ###
library(gsl)
library(ggplot2)
library(gridExtra)
library(rootSolve)

### MGAPR CDF Graph ###   
CDF <- function (α, δ, λ, z){
  A <- α*(1+δ)*(1-exp(-(z^2)/2*(λ^2)))
  B <- (δ*α^(1-exp(-(z^2)/2*(λ^2)))) + (α*(1-exp(-(z^2)/2*(λ^2))))
  cdf<- A/B
  return(cdf)}


#-----------------------------------------------------------------#
MGAPR_quantile=function(α, δ, λ, p){
  f=function(z){
    CDF(α, δ, λ, z) - p
  }
  z=min(uniroot.all(f,lower=0,upper=1000000000000,tol=0.0001))
  return(z)
}


### Skewness of the MGAPR Distribution ###
Skw <- function(α, δ)
{
  Q25 <- MGAPR_quantile(α, δ, 2.7, 0.25)
  Q50 <- MGAPR_quantile(α, δ, 2.7, 0.50)
  Q75 <- MGAPR_quantile(α, δ, 2.7, 0.75)
  valuesk <- (Q75 - 2*Q50 + Q25)/(Q75 - Q25)
  return(valuesk)
}

Skw(0.9,0.03)
 

α <- seq(0.03,0.72, by=0.05)
δ <- seq(0.36,1.05, by=0.05)

Q25 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q50 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q75 <- matrix(NA, ncol = length(α), nrow = length(δ))
Values_skewness <- matrix(NA, ncol = length(α), nrow = length(δ))
for(i in 1:length(α))
{
  for(j in 1:length(δ))
  {
    
    Q25[i,j] <- MGAPR_quantile(α[i], δ[j], 2.7, 0.25)
    Q50[i,j] <- MGAPR_quantile(α[i], δ[j], 2.7, 0.50)
    Q75[i,j] <- MGAPR_quantile(α[i], δ[j], 2.7, 0.75)
    Values_skewness[i,j] <- (Q75[i,j] - 2*Q50[i,j] + Q25[i,j])/(Q75[i,j] - Q25[i,j])
  }
}

par(mfrow=c(1,1))
### Mesh Plot of Skewness ###
persp(α, δ, Values_skewness, theta = 45, phi = 35, r = 4, axes = T, scale = T, 
      box = TRUE, nticks=4, ticktype="detailed", zlim = c(0.01,0.16),
      col="blue", xlab="α", ylab="δ", zlab="Skew", main ="")



### Kurtosis of the MGAPR Distribution ###
α <- seq(0.03,0.72, by=0.05)
δ <- seq(0.36,1.05, by=0.05)

Q125 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q25 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q375 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q50 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q625 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q75 <- matrix(NA, ncol = length(α), nrow = length(δ))
Q875 <- matrix(NA, ncol = length(α), nrow = length(δ))
Values_kurtosis <- matrix(NA, ncol = length(α), nrow = length(δ))
for(i in 1:length(α))
{
  for(j in 1:length(δ))
  {
    Q125[i,j] <- MGAPR_quantile(α[i], δ[j], 2.7, 0.125)
    Q25[i,j]  <- MGAPR_quantile(α[i], δ[j], 2.7, 0.25)
    Q375[i,j] <- MGAPR_quantile(α[i], δ[j], 2.7, 0.375)
    Q50[i,j]  <- MGAPR_quantile(α[i], δ[j], 2.7, 0.50)
    Q625[i,j] <- MGAPR_quantile(α[i], δ[j], 2.7, 0.625)
    Q75[i,j]  <- MGAPR_quantile(α[i], δ[j], 2.7, 0.75)
    Q875[i,j] <- MGAPR_quantile(α[i], δ[j], 2.7, 0.875)
    Values_kurtosis[i,j] <- (Q875[i,j] - Q625[i,j] + Q375[i,j] 
                             - Q125[i,j])/(Q75[i,j] - Q25[i,j])
  }
}


### Mesh Plot of Kurtosis ###
persp(α, δ, Values_kurtosis, theta = 45, phi = 35, r = 5, axes = T, scale = T, 
      box = TRUE, nticks=3,ticktype="detailed", zlim = c(1.2, 1.3),
      col="cyan", xlab="α", ylab="δ", zlab="Kurt", main ="")

