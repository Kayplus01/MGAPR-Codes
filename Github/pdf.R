### Required Packages###
###R code to combine multiple plots###
par(mfrow=c(3,3))
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork) 
library(hrbrthemes)


### MGAPR PDF Graph ###   
PDF <- function (α, δ,λ, z){
  A <- exp(-(z^2)/(2*(λ^2)))
  B <- z*δ*(1+δ)*A*(1-log(α)*(1-A))*α^(-A)
  C <- (λ*(δ*α^(-A) + 1- A ))^2
  pdf<- B/C
  return(pdf)}