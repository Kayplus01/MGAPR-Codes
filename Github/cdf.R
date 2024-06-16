### R Packages needed
###R code to combine multiple plots###
par(mfrow=c(3,3))
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork) 
library(hrbrthemes)


### MGAPR CDF Graph ###   
CDF <- function (α, δ, λ, z){
  A <- (1+δ)*(1-exp(-(z^2)/2*(λ^2)))
  B <- (δ*α^(-exp(-(z^2)/2*(λ^2)))) + 1 -exp(-(z^2)/2*(λ^2))
  cdf<- A/B
  return(cdf)}
