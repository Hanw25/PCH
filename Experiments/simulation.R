################################################################################
# This document includes codes for conducting simulation studies
################################################################################
rm(list=ls())
gc()
library(Matrix)
library(MASS)
library(NbClust)
library(cvxclustr)
source('main-function.R')
source('sim-function.R')
################################################################################
para <- list(lambda1 = seq(0.1,1.6,length=7), lambda2 = seq(0.1,0.6,length=5))
lambda.x <- 10*seq(1,15,length=10)
lambda.z <- 10^seq(11,20,length=10)
################################################################################
n <- 120                      
q <- 6                      
p <- 30
rho <- 0.3                     
va <- 1                        
s <- 0.8                      
tau <- c(0.04,0.08)            
mu1 <- 1.2; mu2 <- 1.6          
stru1 <- 'diag'; stru2 <- 'AR'; stru3 <- 'banded'
sim1 <- 'sim1'; sim2 <- 'sim2'; sim3 <- 'sim3'
# Example
# Gaussian clusters case: Simulation1, mu1, diagonal covariance structure
wholedata <- Generatedata(n, q, p, mu = mu1, stru = stru1, sim = sim1, tau, rho, va, s)
init <- initFun(wholedata)
res.p0 <- tuninglambda(wholedata, init$theta.init0, para, prior = '0', line = F, merge = T)
res.p1 <- tuninglambda(wholedata, init$theta.init1, para, prior = '1', line = F, merge = T)
res.p2 <- tuninglambda(wholedata, init$theta.init2, para, prior = '2', line = F, merge = T)
res.l1 <- hiercvx(wholedata, lambda.x, lambda.z, penal.type = '1')
res.l2 <- hiercvx(wholedata, lambda.x, lambda.z, penal.type = '2')
################################################################################

################################################################################
mu <- 1                      
radius <- 4
va1 <- 0.1^2
va2 <- 1
center1 <- c(-2,0.2); center2 <- c(-3,0.2)
sim4 <- 'sim4'; sim5 <- 'sim5'
# Example
# Two-half moon clusters case: Simulation4, center1
moondata <-  Generatemoondata(n, q, p, mu, sim = sim4, center = center1, radius, tau, va1, va2, s)
init <- initFun(moondata)
res.p0 <- tuninglambda(moondata, init$theta.init0, para, prior = '0', line = F, merge = T)
res.p1 <- tuninglambda(moondata, init$theta.init1, para, prior = '1', line = F, merge = T)
res.p2 <- tuninglambda(moondata, init$theta.init2, para, prior = '2', line = F, merge = T)
res.l1 <- hiercvx(moondata, lambda.x, lambda.z, penal.type = '1')
res.l2 <- hiercvx(moondata, lambda.x, lambda.z, penal.type = '2')
################################################################################

################################################################################
# Case violating hierarchy
# Example
wholedata <- nohierGeneratedata(n, q, p, mu = mu1, stru = stru1, sim = "nohier", tau, rho, va, s)
init <- initFun(wholedata)
res.p0 <- tuninglambda(wholedata, init$theta.init0, para, prior = '0', line = F, merge = T)
res.p1 <- tuninglambda(wholedata, init$theta.init1, para, prior = '1', line = F, merge = T)
res.p2 <- tuninglambda(wholedata, init$theta.init2, para, prior = '2', line = F, merge = T)
res.l1 <- hiercvx(wholedata, lambda.x, lambda.z, penal.type = '1')
res.l2 <- hiercvx(wholedata, lambda.x, lambda.z, penal.type = '2')
################################################################################

################################################################################
# Case with mis-specified prior
# Example
wholedata <- Generatedata(n, q, p, mu = mu1, stru = stru1, sim = sim3, tau, rho, va, s)
misdata <- misGeneratedata(wholedata, c(10,10))
misinit <- initFun(misdata)
misres.p0 <- tuninglambda(misdata, misinit$theta.init0, para, prior = '0', line = F, merge = T)
misres.p1 <- tuninglambda(misdata, misinit$theta.init1, para, prior = '1', line = F, merge = T)
misres.p2 <- tuninglambda(misdata, misinit$theta.init2, para, prior = '2', line = F, merge = T)
################################################################################







