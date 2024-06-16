###R code to combine multiple plots###
par(mfrow=c(3,3))
library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork) 
library(hrbrthemes)


### MGAPR HF Graph ###   
HF <- function (α, δ, λ, z){
  A <- 1-exp(-(z^2)/(2*(λ^2)))
  B <- z*(1+δ)*(exp(-(z^2)/2*(λ^2)))*(1-log(α)*A)*(α^(-exp(-(z^2)/2*(λ^2))))
  C <- (λ^2)*((δ*(α^(-exp(-(z^2)/(2*(λ^2)))))) + A)*((α^(-exp(-(z^2)/(2*(λ^2)))))-A)
  hf<- B/C
  return(hf)}
