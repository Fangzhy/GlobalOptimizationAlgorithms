# script used to solve optimization problems (Thomson's problem and 2D Rastrigin function)
# using methods of Generalized Simulated Annealing algorithm (GenSA) and Deferential Evolution algorithm (DEoptm).
# Fang Yu | 04/11/2016

###
# load library DEoptim
library("DEoptim")
library("GenSA")
library("ggplot2")
options(digits = 10)
set.seed(1234)

rm(list=ls()) # delete all objects
setwd("C:/Users/fang/Documents/FangRepo/GlobalOptimizationAlgorithms")



### Rastrigin's function
Rastrigin <- function(x) {
  fn.call <<- fn.call + 1
  sum(x^2 - 10 * cos(2 * pi * x)) + 10 * length(x)}

### Thomson’s function
#Function of Thomson's problem
Thomson.fn <- function(ThetaPhi) {
  fn.call <<- fn.call + 1
  ThetaPhi <- matrix(ThetaPhi, ncol = 2)
  xyz <- t(apply(ThetaPhi, 1, function(thetaphi) {
    c(sin(thetaphi[1]) * cos(thetaphi[2]),
      sin(thetaphi[1]) * sin(thetaphi[2]), cos(thetaphi[1]))}))
  #print(xyz)
  n <- nrow(ThetaPhi)
  tmp <- matrix(NA, nrow = n, ncol = n)
  index <- cbind(as.vector(row(tmp)), as.vector(col(tmp)))
  index <- index[index[, 1] < index[, 2], , drop=F]
  rdist <- apply(index, 1, function(idx) {
    tmp <- 1/sqrt(sum((xyz[idx[1], ] - xyz[idx[2], ])^2))})
  res <- sum(rdist)
  return(res)
}


### Visualization of  2D Rastrigin function
# Adapated from https://jamesmccaffrey.wordpress.com/2015/10/07/graphing-rastrigins-function-in-3d-color-gradient-using-r/

x0 <- seq(-5.12, 5.12, length=100)
x1 <- seq(-5.12, 5.12, length=100)

Rastrigin2D <- function(x0, x1) { 20 + (x0^2 - 10 *
                                          cos(2 * 3.14 *x0)) + (x1^2 - 10 *
                                                                  cos(2 * 3.14 *x1)) }
z <- outer(x0, x1, Rastrigin2D)
jet.colors <- colorRampPalette(c("midnightblue",
                                 "blue", "cyan", "green", "yellow", "orange", "red", "darkred"))
nbcol <- 64
color <- jet.colors(nbcol)
nrz <- nrow(z)
ncz <- ncol(z)
zfacet <- z[-1,-1] + z[-1,-ncz] +
  z[-nrz,-1] + z[-nrz,-ncz]
facetcol <- cut(zfacet, nbcol)
persp(x0, x1, z, col=color[facetcol],
      phi=15, theta=-35, ticktype="detailed",
      d=10, r=1, shade=0.1, expand=0.7,zlab = "f(x0,x1)")

### Use DEoptm to solve 2D Rastrigin function
# Adapted from Xiang Yang (2013)
dimension <- 2
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
sink("tmp.txt")
fn.call <<- 0
out.DEoptim <- DEoptim(fn = Rastrigin, lower = lower, upper = upper,
                       control = list(storepopfrom = 1))
sink(NULL)

#Check output parameters and results
out.DEoptim$optim[c(1, 2)]
fn.call


### Use GenSA to solve 2D Rastrigin function
# Adapted from Xiang Yang (2013)
expected.val <- 0
absTol <- 1e-13
fn.call <- 0
out.GenSA <- GenSA(par = NULL, lower = lower, upper = upper, fn = Rastrigin,
                   control = list(threshold.stop = expected.val + absTol))
#Check output
out.GenSA[c("value", "par", "counts")]

cat("GenSA call functions", fn.call, "times.\n")

### Comparison between two methods for solving 2D Rastrigin function
# This script is used to compare GenSA and DEoptim for solving 2D Rastrigin function

dimension <- 2
NumRun <- 100
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
#
expected.val <- 0
absTol <- c(1e-05, 1e-7, 1e-9, 1e-11)
NumSuccessRun <- matrix(0,nrow=2,ncol=length(absTol))
NumFnCallSA <- matrix(NaN,nrow=NumRun,ncol=length(absTol))
NumFnCallDE <- matrix(NaN,nrow=NumRun,ncol=length(absTol))
sink("tmp.txt") #direct intermidiate output into tmp.txt
for (i in 1:length(absTol)){
  for (j in 1:NumRun){
    #Initialize fn.call as 0 for GenSA
    fn.call <<- 0
    out.GenSA <- GenSA(par = NULL, lower = lower, upper = upper, fn = Rastrigin,
                       control = list(threshold.stop = expected.val + absTol[i]))
    if (out.GenSA["value"] <= (expected.val + absTol[i])) NumSuccessRun[1,i] = NumSuccessRun[1,i] + 1  
    NumFnCallSA[j,i] <- fn.call
    #Initialize fn.call as 0 for DEoptim
    fn.call <<- 0
    out.DEoptim <- DEoptim(fn = Rastrigin, lower = lower, upper = upper,
                           control = list(VTR = expected.val + absTol[i]))
    if (out.DEoptim$optim["bestval"] <= (expected.val + absTol[i])) NumSuccessRun[2,i] = NumSuccessRun[2,i] + 1
    NumFnCallDE[j,i] <- fn.call
  }
}
sink(NULL) # show output in console
NumFnCallAveSA <- round(colMeans(NumFnCallSA))
NumFnCallSdSA <- round(apply(NumFnCallSA,2,sd))
NumFnCallAveDE <- round(colMeans(NumFnCallDE))
NumFnCallSdDE <- round(apply(NumFnCallDE,2,sd))
#Test significance in Number of function calls between GenSA and DEoptim
pvalue <- array(NaN,dim=c(1,length(absTol)))
for (i in 1:length(absTol)){
  ttestRes <- t.test(NumFnCallSA[,i],NumFnCallDE[,i])
  pvalue[i] <- ttestRes["p.value"]
}
#Visualize the results
#successfull run
DataSuccessRun <- data.frame(Methods = rep(c("GenSA","DEoptim"),each=length(absTol)),
                             absTol = as.factor(rep(absTol,2)), NumSuccessRun = as.vector(t(NumSuccessRun)))
ggplot(data=DataSuccessRun, aes(x=absTol, y=NumSuccessRun, fill=Methods)) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_text(aes(label=NumSuccessRun), vjust=1.6, color="white",
            position = position_dodge(0.9), size=5.0) +
  xlab("Absolute Tolerance") +
  ylab("Successfull Runs %")
#theme_minimal()
#Average Number of function calls
DataFnCall <- data.frame(Methods = rep(c("GenSA","DEoptim"),each=length(absTol)),
                         absTol = as.factor(rep(absTol,2)), FnCallAve = c(NumFnCallAveSA,NumFnCallAveDE),
                         FnCallSd = c(NumFnCallSdSA,NumFnCallSdDE))
limits <- aes(ymax = FnCallAve + FnCallSd, ymin=FnCallAve - FnCallSd)
ggplot(DataFnCall, aes(fill=Methods, y=FnCallAve, x=absTol)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(limits, position='dodge', width=0.25) +
  xlab("Absolute Tolerance") +
  ylab("Average Number of Function Calls") + 
  geom_text(aes(y=FnCallAve, label=c("***","***","***","***","\t","\t","\t","\t")), vjust=-1.6, color="black",
            position = position_dodge(0.9), size=5.0)

### Comparison between two methods for solving Thomson’s problem
n.particles <- 12
lower.T <- rep(0, 2 * n.particles)
upper.T <- c(rep(pi, n.particles), rep(2 * pi, n.particles))
#set a random seed so that the results can be reproduced later on
set.seed(1234)
####Solve Thompson problem using GenSA
#Minimum Value
minValue <- out.GenSA["value"]
#The corresponding Theta and Phi which generated the minValue
out.GenSA <- GenSA(par = NULL, lower = lower.T, upper = upper.T, fn = Thomson.fn)
ThetaPhiFinal <- unlist(out.GenSA["par"])
ThetaFinal <- ThetaPhiFinal[1:n.particles]
PhiFinal <- ThetaPhiFinal[(1:n.particles)+n.particles]
data.frame(theta=as.numeric(ThetaFinal),phi= as.numeric(PhiFinal))
out.GenSA["value"]
####Compare with DEoptim (take long time to run)
itermax <- seq(5,250,by=5)
FminGenSA <- rep(NaN,length(itermax))
FminDEoptim <- rep(NaN,length(itermax))
FnCall <-rep(NaN,length(itermax))
sink("tmp.txt")
for (i in (1:length(itermax))){
  #set the number of function calls as 0 for DEoptim
  fn.call <<- 0
  out.DEoptim <- DEoptim(fn = Thomson.fn, lower = lower.T, upper = upper.T, 
                         control = list(itermax=itermax[i]))
  FminDEoptim[i] <- as.numeric(out.DEoptim$optim["bestval"])
  FnCall[i] <- fn.call
  #set the number of function calls as 0 for GenSA
  fn.call <<- 0
  out.GenSA <- GenSA(par = NULL, lower = lower.T, upper = upper.T,
                     fn = Thomson.fn,control = list(max.call=FnCall[i]))
  FminGenSA[i] <- as.numeric(out.GenSA["value"])
}
# save.image("ThomsonCompareData")
sink(NULL)
#Visualize the results
df <- data.frame(x=rep(FnCall,2), y=c(FminGenSA, FminDEoptim), 
                 Method=c(rep("GenSA", length(itermax)), rep("DEoptim", length(itermax))))
ggplot(df, aes(x=x, y=y, color=Method)) + geom_line(size=1.5) +
  xlab("Number of function calls") +
  ylab("Minimum Value of Objective Function") 


