#script used to use GenSA to solve optimization problem (Thomson's problem)
#Fang Yu, 03/07/2016
#
#clear the memory
rm(list=ls()) 
#set the working directory (folder)
setwd("C:/Users/fang/Documents/FangRepo/GlobalOptimizationAlgorithms/src")

#load the library GenSA
library(GenSA)

#
n.particles <- 12
lower.T <- rep(0, 2 * n.particles)
upper.T <- c(rep(pi, n.particles), rep(2 * pi, n.particles))

#set a random seed so that the results can be reproduced later on
set.seed(1234)

#set the number of function calls as 0
fn.call <<- 0

# source the Thomson's function 
source("Thomson.R")

#Perform simulation using General Simulated Annealing
fn.call.DEoptim = 10000

out.GenSA <- GenSA(par = NULL, lower = lower.T, upper = upper.T,
                     fn = Thomson.fn, control = list(max.call = fn.call.DEoptim))

#Check the results
out.GenSA[c("value", "counts")]

#Minimum Value
minValue <- out.GenSA["value"]
#The corresponding Theta and Phi which generated the minValue
ThetaPhiFinal <- unlist(out.GenSA["par"])
ThetaFinal <- ThetaPhiFinal[1:n.particles]
PhiFinal <- ThetaPhiFinal[(1:n.particles)+n.particles]


