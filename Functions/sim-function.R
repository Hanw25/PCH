#####################################################################################
## All functions for simulation studies.
#####################################################################################

###########################################################################################################################
## Part I #################################################################################################################
##################################### Functions for generating simulation data ############################################
Generategroup <- function(n, K1, K2){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate the hierarchical structure
  ## -----------------------------------------------------------------------------------------------------------------
  group.main <- vector(mode = 'list', length = K1)
  group.hier <- vector(mode = 'list', length = K2)
  div <- K2/K1; group.hier0 <- group.hier
  ind <- sample(K2, n, replace = T)
  min <- vector(mode = 'integer', length = K2)
  for (k2 in 1:K2) {
    group.hier0[[k2]] <- which(ind == k2)
    min[k2] <- group.hier0[[k2]][1]
  }
  rank <- rank(min)
  for (k2 in 1:K2) {
    group.hier[[k2]] <- group.hier0[[which(rank == k2)]]
  }
  for (k1 in 1:K1) {
    vec <- NULL
    for (i in 1:div) {
      vec <- c(vec, group.hier[[(k1-1)*div+i]])
    }
    group.main[[k1]] <- sort(vec)
  }
  return(list(group.main = group.main, group.hier = group.hier))
}

Generateprior <- function(group.hier, tau){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate the prior information
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(group.hier))
  pairind <- pairindFun(group.hier)
  if(tau > 0){
    pairind.p <- pairind[, sort(sample(ncol(pairind), floor(ncol(pairind)*tau), replace = F))]
  }else{
    pairind.p <- matrix(c(0,0), 2, 1)
  }
  group.p <- groupFun(pairmatFun(n, pairind.p), merge = T)
  return(group.p)
}

Generatedata <- function(n, q, p, mu, stru, sim, tau, rho, va, s){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate wholedata
  ## -----------------------------------------------------------------------------------------------------------------
  sigx <- diag(va, q, q)
  sigz <- matrix(0, nrow = p, ncol = p)
  if(stru == 'diag'){
    diag(sigz) <- va
  }
  if(stru == 'AR'){
    for (i in 1:p) {
      for (j in i:p) {
        sigz[i,j] <- va*(rho^abs(i-j))
        sigz[j,i] <- sigz[i,j]
      }
    }
    diag(sigz) <- va
  }
  if(stru == 'banded'){
    for (i in 1:(p-1)) {
      sigz[i,i+1] <- va*rho
    }
    sigz <- sigz+t(sigz)
    diag(sigz) <- va
  }
  
  if(sim == 'sim1'){
    K1 <- 2; K2 <- 4
    one <- rep(1,q/3)
    xi.true.mat <- s*mu*cbind(c(-one,one,one),c(one,one,-one))
    one <- rep(1,p/2)
    alpha.true.mat <- mu*cbind(c(one,one),c(-one,one),
                               c(one,-one),c(-one,-one))
  }
  if(sim == 'sim2'){
    K1 <- 2; K2 <- 6
    one <- rep(1,q/3)
    xi.true.mat <- s*mu*cbind(c(-one,one,one),c(one,one,-one))
    one <- rep(1,p/3)
    alpha.true.mat <- mu*cbind(c(-one,one,one),c(one,-one,one),c(one,one,-one),
                               c(one,-one,-one),c(-one,one,-one),c(-one,-one,one))
  }
  if(sim == 'sim3'){
    K1 <- 3; K2 <- 6
    one <- rep(1,q/3)
    xi.true.mat <- s*mu*cbind(c(-one,one,one),c(one,-one,one),c(one,one,-one))
    one <- rep(1,p/3)
    alpha.true.mat <- mu*cbind(c(-one,one,one),c(one,-one,-one),
                               c(one,-one,one),c(-one,one,-one),
                               c(one,one,-one),c(-one,-one,one))
  }
  
  group <- Generategroup(n,K1,K2)
  group.main <- group$group.main
  group.hier <- group$group.hier
  group.p0 <- Generateprior(group.hier, 0)
  group.p1 <- Generateprior(group.hier, tau[1])
  group.p2 <- Generateprior(group.hier, tau[2])
  
  data.x <- matrix(0, nrow = q, ncol = n)
  beta.true.mat <- matrix(0, nrow = q, ncol = n)
  for (k1 in 1:K1) {
    data.x[,group.main[[k1]]] <- t(MASS::mvrnorm(length(group.main[[k1]]), xi.true.mat[,k1], sigx))
    beta.true.mat[,group.main[[k1]]] <- matrix(rep(xi.true.mat[,k1], length(group.main[[k1]])), nrow = q)
  }
  
  data.z <- matrix(0, nrow = p, ncol = n)
  gamma.true.mat <- matrix(0, nrow = p, ncol = n)
  for (k2 in 1:K2) {
    data.z[,group.hier[[k2]]] <- t(MASS::mvrnorm(length(group.hier[[k2]]), alpha.true.mat[,k2], sigz))
    gamma.true.mat[,group.hier[[k2]]] <- matrix(rep(alpha.true.mat[,k2], length(group.hier[[k2]])), nrow = p)
  }
  
  return(list(data.x = data.x, data.z = data.z, group.main = group.main, group.hier = group.hier, 
              group.p0 = group.p0, group.p1 = group.p1, group.p2 = group.p2,
              xi.true.mat = xi.true.mat, alpha.true.mat = alpha.true.mat, 
              beta.true.mat = beta.true.mat, gamma.true.mat = gamma.true.mat))
}

Generatemoondata <- function(n, q, p, mu, sim, center, radius, tau, va1, va2, s){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate half-moon data
  ## -----------------------------------------------------------------------------------------------------------------
  K1 <- 2
  one <- rep(1,(q-2)/2)
  xi2.true.mat <- s*mu*cbind(c(-one,one),c(one,-one))
  
  if(sim == 'sim4'){
    K2 <- 4
    one <- rep(1,p/2)
    alpha.true.mat <- mu*cbind(c(one,one),c(-one,one),
                               c(one,-one),c(-one,-one))
  }
  if(sim == 'sim5'){
    K2 <- 6
    one <- rep(1,p/3)
    alpha.true.mat <- mu*cbind(c(-one,one,one),c(one,-one,one),c(one,one,-one),
                               c(one,-one,-one),c(-one,one,-one),c(-one,-one,one))
  }
  
  group <- Generategroup(n,K1,K2)
  group.main <- group$group.main
  group.hier <- group$group.hier
  group.p0 <- Generateprior(group.hier, 0)
  group.p1 <- Generateprior(group.hier, tau[1])
  group.p2 <- Generateprior(group.hier, tau[2])
  
  data.x2 <- matrix(0, nrow = q-2, ncol = n)
  beta2.true.mat <- matrix(0, nrow = q-2, ncol = n)
  for (k1 in 1:K1) {
    data.x2[,group.main[[k1]]] <- t(MASS::mvrnorm(length(group.main[[k1]]), xi2.true.mat[,k1], diag(q-2)*va2))
    beta2.true.mat[,group.main[[k1]]] <- matrix(rep(xi2.true.mat[,k1], length(group.main[[k1]])), nrow = q-2)
  }
  
  data.z <- matrix(0, nrow = p, ncol = n)
  gamma.true.mat <- matrix(0, nrow = p, ncol = n)
  for (k2 in 1:K2) {
    data.z[,group.hier[[k2]]] <- t(MASS::mvrnorm(length(group.hier[[k2]]), alpha.true.mat[,k2], diag(p)*va2))
    gamma.true.mat[,group.hier[[k2]]] <- matrix(rep(alpha.true.mat[,k2], length(group.hier[[k2]])), nrow = p)
  }
  
  center1 <- center
  center2 <- -center
  n1 <- length(group.main[[1]])
  n2 <- length(group.main[[2]])
  pi1 <- runif(n1, 0, pi)
  pi2 <- runif(n2, pi, 2*pi)
  beta1.true.mat <- matrix(0, nrow = 2, ncol = n)
  beta1.true.mat[,group.main[[1]]] <- rbind(cos(pi1)*radius+center1[1], sin(pi1)*radius+center1[2])
  beta1.true.mat[,group.main[[2]]] <- rbind(cos(pi2)*radius+center2[1], sin(pi2)*radius+center2[2])
  xi1.true.mat <- matrix(0, nrow = 2, ncol = 2)
  xi1.true.mat[,1] <- apply(beta1.true.mat[,group.main[[1]]], 1, mean)
  xi1.true.mat[,2] <- apply(beta1.true.mat[,group.main[[2]]], 1, mean)
  data.x1 <- beta1.true.mat+t(MASS::mvrnorm(n,c(0,0),diag(2)*va1))
  xi.true.mat <- rbind(xi1.true.mat,xi2.true.mat)
  beta.true.mat <- rbind(beta1.true.mat,beta2.true.mat)
  data.x <- rbind(data.x1,data.x2)
  
  return(list(data.x = data.x, data.z = data.z, group.main = group.main, group.hier = group.hier, 
              group.p0 = group.p0, group.p1 = group.p1, group.p2 = group.p2,
              xi.true.mat = xi.true.mat, alpha.true.mat = alpha.true.mat, 
              beta.true.mat = beta.true.mat, gamma.true.mat = gamma.true.mat))
}


###########################################################################################################################
## Part II ################################################################################################################
######################################## Functions for evaluating performances  ###########################################
grperFun <- function(gr.gf, group.true){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate estimated group structure performance (RI, TPR, FPR)
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(gr.gf))
  cap.gf.matrix <- pairmatFun(n, pairindFun(gr.gf))
  cap.true.matrix <- pairmatFun(n, pairindFun(group.true))
  gr.tp <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.gf.matrix) == 1)))
  gr.fp <- length(intersect(which(as.vector(cap.true.matrix) == 0), which(as.vector(cap.gf.matrix) == 1)))
  gr.fn <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.gf.matrix) == 0)))
  gr.tn <- n*(n-1)/2-gr.tp-gr.fp-gr.fn
  gr.tpr <- gr.tp/(gr.tp+gr.fn)
  gr.fpr <- gr.fp/(gr.fp+gr.tn)
  gr.ri <- (gr.tp+gr.tn)/(n*(n-1)/2)
  gr.df <- data.frame('gr.tpr' = gr.tpr, 'gr.fpr' = gr.fpr, 'gr.ri' = gr.ri)
  return(gr.df)
}

mseFun <- function(beta.gf.mat, beta.true.mat, gamma.gf.mat, gamma.true.mat){
  mse.beta <- sqrt(sum((beta.true.mat-beta.gf.mat)^2)/length(beta.true.mat))
  mse.gamma <- sqrt(sum((gamma.true.mat-gamma.gf.mat)^2)/length(gamma.true.mat))
  mse.df <- data.frame('mse.beta' = mse.beta, 'mse.gamma' = mse.gamma)
  return(mse.df)
}

perFun <- function(result, wholedata){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate performances (BIC, RI, TPR, FPR, MSE)
  ## -----------------------------------------------------------------------------------------------------------------
  q <- nrow(wholedata$data.x)
  p <- nrow(wholedata$data.z)
  n <- ncol(wholedata$data.z)

  mse <- mseFun(result$beta.gf.mat, wholedata$beta.true.mat, result$gamma.gf.mat, wholedata$gamma.true.mat)
  mse.beta <- mse$mse.beta
  mse.gamma <- mse$mse.gamma
  
  gr.main <- result$gr.main
  gr.main.num <- gr.main$gr.num
  gr.main.gf <- gr.main$gr.gf
  gr.main.per <- grperFun(gr.main.gf, wholedata$group.main)
  gr.main.tpr <- gr.main.per$gr.tpr
  gr.main.fpr <- gr.main.per$gr.fpr
  gr.main.ri <- gr.main.per$gr.ri
  
  gr.hier <- result$gr.hier
  gr.hier.num <- gr.hier$gr.num
  gr.hier.gf <- gr.hier$gr.gf
  gr.hier.per <- grperFun(gr.hier.gf, wholedata$group.hier)
  gr.hier.tpr <- gr.hier.per$gr.tpr
  gr.hier.fpr <- gr.hier.per$gr.fpr
  gr.hier.ri <- gr.hier.per$gr.ri
  
  BIC <- result$BIC
  
  per.df <- data.frame('BIC' = BIC,
                       'gr.main.num' = gr.main.num, 'gr.main.ri' = gr.main.ri, 'gr.main.tpr' = gr.main.tpr, 'gr.main.fpr' = gr.main.fpr, 
                       'gr.hier.num' = gr.hier.num, 'gr.hier.ri' = gr.hier.ri, 'gr.hier.tpr' = gr.hier.tpr, 'gr.hier.fpr' = gr.hier.fpr, 
                       'mse.beta' = mse.beta, 'mse.gamma' = mse.gamma)
  
  return(list(per.df = per.df, gr.main.gf = gr.main.gf, gr.hier.gf = gr.hier.gf))
}

RIFun <- function(gr.gf, group.true){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate estimated group structure performance (RI, ARI)
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(gr.gf))
  cap.gf.matrix <- pairmatFun(n, pairindFun(gr.gf))
  cap.true.matrix <- pairmatFun(n, pairindFun(group.true))
  gr.tp <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.gf.matrix) == 1)))
  gr.fp <- length(intersect(which(as.vector(cap.true.matrix) == 0), which(as.vector(cap.gf.matrix) == 1)))
  gr.fn <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.gf.matrix) == 0)))
  gr.tn <- n*(n-1)/2-gr.tp-gr.fp-gr.fn
  gr.ri <- (gr.tp+gr.tn)/(n*(n-1)/2)
  gr.ari <- (n*(n-1)/2*(gr.tp+gr.tn)-(gr.tp+gr.fp)*(gr.tp+gr.fn)-(gr.tn+gr.fp)*(gr.tn+gr.fn))/((n*(n-1)/2)^2-(gr.tp+gr.fp)*(gr.tp+gr.fn)-(gr.tn+gr.fp)*(gr.tn+gr.fn))
  return(list(gr.ri = gr.ri, gr.ari = gr.ari))
}

newperFun <- function(result, wholedata){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate performances (BIC, RI, ARI, MSE)
  ## -----------------------------------------------------------------------------------------------------------------
  q <- nrow(wholedata$data.x)
  p <- nrow(wholedata$data.z)
  n <- ncol(wholedata$data.z)
  mse <- mseFun(result$beta.gf.mat, wholedata$beta.true.mat, result$gamma.gf.mat, wholedata$gamma.true.mat)
  mse.beta <- mse$mse.beta
  mse.gamma <- mse$mse.gamma
  
  gr.main <- result$gr.main
  gr.main.num <- gr.main$gr.num
  gr.main.per <- RIFun(gr.main$gr.gf, wholedata$group.main)
  gr.main.ri <- gr.main.per$gr.ri
  gr.main.ari <- gr.main.per$gr.ari
  
  gr.hier <- result$gr.hier
  gr.hier.num <- gr.hier$gr.num
  gr.hier.per <- RIFun(gr.hier$gr.gf, wholedata$group.hier)
  gr.hier.ri <- gr.hier.per$gr.ri
  gr.hier.ari <- gr.hier.per$gr.ari
  
  BIC <- result$BIC
  
  # performance summary 
  per.df <- data.frame('BIC' = BIC,
                       'gr.main.num' = gr.main.num, 'gr.main.ri' = gr.main.ri, 'gr.main.ari' = gr.main.ari, 
                       'gr.hier.num' = gr.hier.num, 'gr.hier.ri' = gr.hier.ri, 'gr.hier.ari' = gr.hier.ari, 
                       'mse.beta' = mse.beta, 'mse.gamma' = mse.gamma)
  return(per.df)
}

tuninglambda <- function(wholedata, theta.init, para, prior,
                         iter.max = 50, epsi = 0.005, gam.mcp = 3, penal.para = 1, c = 1, merge = F, line = F){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Tuning procedure (perFun or newperFun)
  ## -----------------------------------------------------------------------------------------------------------------
  lambda1 <- para$lambda1; L1 <- length(lambda1)
  lambda2 <- para$lambda2; L2 <- length(lambda2)
  
  if(!line){
    admmres <- vector(mode = "list", length = L1*L2)
    per.df <- data.frame()
    # tuning lambda2 and lambda1 by grid searching
    for (i in 1:L2) {
      for (j in 1:L1) {
        cat('-----------', i, '-th lambda2', lambda2[i] , 'and', j, '-th lambda1', lambda1[j] ,'--------------\n')
        admmres[[(i-1)*L1+j]] <- ADMM(wholedata = wholedata, theta.init = theta.init, lambda1 = lambda1[j], lambda2 = lambda2[i], prior = prior,
                                      iter.max = iter.max, epsi = epsi, gam.mcp = gam.mcp, penal.para = penal.para, c = c, merge = merge)
        peri <- data.frame('lambda1' = lambda1[j], 'lambda2' = lambda2[i], (perFun(admmres[[(i-1)*L1+j]], wholedata))$per.df)
        per.df <- rbind(per.df, peri)
      }
    }
  }else{
    admmres <- vector(mode = "list", length = L1+L2)
    per.df <- data.frame()
    # tuning lambda2 and lambda1 by line searching
    for (i in 1:L2) {
      cat('-----------', i, '-th lambda2', lambda2[i] , 'and', 'fixed lambda1', 0, '--------------\n')
      admmres[[i]] <- ADMM(wholedata = wholedata, theta.init = theta.init, lambda1 = 0, lambda2 = lambda2[i], prior = prior,
                           iter.max = iter.max, epsi = epsi, gam.mcp = gam.mcp, penal.para = penal.para, c = c, merge = merge)
      peri <- data.frame('lambda1' = 0, 'lambda2' = lambda2[i], (perFun(admmres[[i]], wholedata))$per.df)
      per.df <- rbind(per.df, peri)
    }
    lambda2.opt <- per.df[min(which(per.df$BIC == min(per.df$BIC))),]$lambda2
    for (j in 1:L1) {
      cat('-----------', 'fixed optimal lambda2', lambda2.opt, 'and', j, '-th lambda1', lambda1[j], '--------------\n')
      admmres[[L2+j]] <- ADMM(wholedata = wholedata, theta.init = theta.init, lambda1 = lambda1[j], lambda2 = lambda2.opt, prior = prior,
                              iter.max = iter.max, epsi = epsi, gam.mcp = gam.mcp, penal.para = penal.para, c = c, merge = merge)
      peri <- data.frame('lambda1' = lambda1[j], 'lambda2' = lambda2.opt, (perFun(admmres[[L2+j]], wholedata))$per.df)
      per.df <- rbind(per.df, peri)
    }
  }
  
  # summary
  res.tune <- admmres[[min(which(per.df$BIC == min(per.df$BIC)))]]
  per.df.tune <- per.df[min(which(per.df$BIC == min(per.df$BIC))),]
  
  return(list(admmres = admmres, per.df = per.df, res.tune = res.tune, per.df.tune = per.df.tune))
}


###########################################################################################################################
## Part III ###############################################################################################################
######################################## Functions for additional exploration #############################################
nohiergroup <- function(n, K1 = 2, K2 = 6){
  nohier <- Generategroup(n, K1, K2)
  group.main <- nohier$group.main
  group.hier0 <- nohier$group.hier
  group.hier <- vector(mode = 'list', length = 5)
  group.hier[[1]] <- group.hier0[[1]]
  group.hier[[2]] <- group.hier0[[2]]
  group.hier[[3]] <- sort(c(group.hier0[[3]],group.hier0[[4]]))
  group.hier[[4]] <- group.hier0[[5]]
  group.hier[[5]] <- group.hier0[[6]]
  return(list(group.main = group.main, group.hier = group.hier, group.hier0 = group.hier0))
}

nohierGeneratedata <- function(n, q, p, mu, stru, sim = "nohier", tau, rho, va, s){
  sigx <- diag(va, q, q)
  sigz <- matrix(0, nrow = p, ncol = p)
  if(stru == 'diag'){
    diag(sigz) <- va
  }
  if(stru == 'AR'){
    for (i in 1:p) {
      for (j in i:p) {
        sigz[i,j] <- va*(rho^abs(i-j))
        sigz[j,i] <- sigz[i,j]
      }
    }
    diag(sigz) <- va
  }
  if(stru == 'banded'){
    for (i in 1:(p-1)) {
      sigz[i,i+1] <- va*rho
    }
    sigz <- sigz+t(sigz)
    diag(sigz) <- va
  }
  
  if(sim == 'nohier'){
    K1 <- 2; K2 <- 6
    one <- rep(1,q/3)
    xi.true.mat <- s*mu*cbind(c(-one,one,one),c(one,one,-one))
    one <- rep(1,p/2)
    alpha.true.mat <- mu*cbind(c(one,one),c(-one,one),
                               rep(0,p),
                               c(one,-one),c(-one,-one))
  }
  
  group <- nohiergroup(n,K1,K2)
  group.main <- group$group.main
  group.hier <- group$group.hier
  group.hier0 <- group$group.hier0
  group.p0 <- Generateprior(group.hier0, 0)
  group.p1 <- Generateprior(group.hier0, tau[1])
  group.p2 <- Generateprior(group.hier0, tau[2])
  
  data.x <- matrix(0, nrow = q, ncol = n)
  beta.true.mat <- matrix(0, nrow = q, ncol = n)
  for (k1 in 1:K1) {
    data.x[,group.main[[k1]]] <- t(MASS::mvrnorm(length(group.main[[k1]]), xi.true.mat[,k1], sigx))
    beta.true.mat[,group.main[[k1]]] <- matrix(rep(xi.true.mat[,k1], length(group.main[[k1]])), nrow = q)
  }
  
  data.z <- matrix(0, nrow = p, ncol = n)
  gamma.true.mat <- matrix(0, nrow = p, ncol = n)
  for (k2 in 1:5) {
    data.z[,group.hier[[k2]]] <- t(MASS::mvrnorm(length(group.hier[[k2]]), alpha.true.mat[,k2], sigz))
    gamma.true.mat[,group.hier[[k2]]] <- matrix(rep(alpha.true.mat[,k2], length(group.hier[[k2]])), nrow = p)
  }
  
  return(list(data.x = data.x, data.z = data.z, group.main = group.main, group.hier = group.hier, 
              group.hier0 = group.hier0, group.p0 = group.p0, group.p1 = group.p1, group.p2 = group.p2,
              xi.true.mat = xi.true.mat, alpha.true.mat = alpha.true.mat, 
              beta.true.mat = beta.true.mat, gamma.true.mat = gamma.true.mat))
}

misgroup <- function(group, m){
  K <- length(group); n <- length(unlist(group))
  num <- c()
  for (k in 1:K) {
    num[k] <- length(group[[k]])
  }
  gr.min <- group[order(num)[1:m]]
  gr.max <- group[order(num)[(K-m+1):K]]
  gr <- vector(mode = "list", length = m)
  for (i in 1:m) {
    gr[[i]] <- sort(c(gr.min[[i]],gr.max[[(m-i+1)]]))
  }
  group <- c(group[order(num)[(m+1):(K-m)]],gr)
  group <- groupFun(pairmatFun(n,pairindFun(group)))
  return(group)
}

misGeneratedata <- function(wholedata, m){
  group.p1 <- misgroup(wholedata$group.p1, m[1])
  group.p2 <- misgroup(wholedata$group.p2, m[2])
  wholedata$group.p1 <- group.p1
  wholedata$group.p2 <- group.p2
  return(wholedata)
}


###########################################################################################################################
## Part IV ################################################################################################################
##################################### Functions for calculating the initial values ########################################
groupkmeansFun <- function(res){
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Generate the group structure of k-means result
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  group0 <- vector(mode = 'list',length = max(res$cluster))
  group <- vector(mode = 'list',length = max(res$cluster))
  min <- vector(mode = 'integer', length = max(res$cluster))
  for (k in 1:length(group0)) {
    group0[[k]] <- which(res$cluster == k)
    min[k] <- group0[[k]][1]
  }
  rank <- rank(min)
  for (k in 1:length(group)) {
    group[[k]] <- group0[[which(rank == k)]]
  }
  return(list(gr.gf = group, rank = rank))
}

initFun <- function(wholedata, kmin = 2, kmax = 4, divmin = 2, divmax = 3, c = 1){
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Generate the initial value for proposed method via k-means
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  cova.x <- wholedata$data.x; q <- nrow(cova.x)
  cova.z <- wholedata$data.z; p <- nrow(cova.z)
  cova <- rbind(cova.x,cova.z); n <- ncol(cova)
  
  K1 <- NbClust(t(cova.x), distance = "euclidean", min.nc = kmin, max.nc = kmax, method = "kmeans", index = 'ch', alphaBeale=0.1)
  K1 <- max(K1$Best.partition)
  res.x <- kmeans(x = t(cova.x), centers = K1, nstart = 25)
  gr.main.gf <- groupkmeansFun(res.x)$gr.gf
  gr.main.num <- length(gr.main.gf)
  gr.main <- list(gr.num = gr.main.num, gr.gf = gr.main.gf)
  rank.x <- groupkmeansFun(res.x)$rank
  xi.gf.mat <- as.matrix(t(res.x$centers[rank.x,]))
  colnames(xi.gf.mat) <- NULL
  beta.gf.mat <- matrix(0, nrow = q, ncol = n)
  for (k1 in 1:gr.main.num) {
    beta.gf.mat[,gr.main.gf[[k1]]] <- matrix(rep(xi.gf.mat[,k1], length(gr.main.gf[[k1]])), nrow = q)
  }
  
  gamma.gf.mat <- matrix(0, nrow = p, ncol = n)
  gr.hier.gf <- NULL
  for (k in 1:gr.main.num){
    indk <- gr.main.gf[[k]]
    if(length(indk) > 1){
      cova.zk <- cova.z[,indk]
      K2 <- NbClust(t(cova.zk), distance = "euclidean", min.nc = divmin, max.nc = divmax, method = "kmeans", index = 'ch', alphaBeale=0.1)
      K2 <- max(K2$Best.partition)
      res.zk <- kmeans(x = t(cova.zk), centers = K2, nstart = 25)
      gr.gf.zk <- groupkmeansFun(res.zk)$gr.gf
      rank.zk <- groupkmeansFun(res.zk)$rank
      alpha.gf.mat.zk <- as.matrix(t(res.zk$centers[rank.zk,]))
      colnames(alpha.gf.mat.zk) <- NULL
      for (m in 1:length(gr.gf.zk)) {
        gr.hier.gf <- c(gr.hier.gf, list(indk[gr.gf.zk[[m]]]))
        gamma.gf.mat[,indk[gr.gf.zk[[m]]]] <- matrix(rep(alpha.gf.mat.zk[,m], length(gr.gf.zk[[m]])), nrow = p)
      }
    }else{
      gamma.gf.mat[,indk] <- cova.z[,indk]
      gr.hier.gf <- c(gr.hier.gf, list(indk))
    }
  }
  gr.hier.num <- length(gr.hier.gf)
  gr.hier <- list(gr.num = gr.hier.num, gr.gf = gr.hier.gf)
  
  alpha.gf.mat <- matrix(0, nrow = p, ncol = gr.hier.num)
  for (k in 1:gr.hier.num) {
    alpha.gf.mat[,k] <- gamma.gf.mat[,gr.hier.gf[[k]][1]]
  }
  
  BIC <- log(sum((cova-rbind(beta.gf.mat,gamma.gf.mat))^2/n))+c*log(log(n))*log(n)*(gr.main.num+gr.hier.num)/n
  
  kmeansres <- list(beta.gf.mat = beta.gf.mat, gamma.gf.mat = gamma.gf.mat, xi.gf.mat = xi.gf.mat, alpha.gf.mat = alpha.gf.mat,
                    gr.main = gr.main, gr.hier = gr.hier, BIC = BIC)
  
  theta.init.mat <- rbind(beta.gf.mat, gamma.gf.mat)
  
  group.p0 <- wholedata$group.p0
  theta.init0 <- as.matrix(kronecker(LMatrixFun(group.p0), bdiag(diag(q+p)))%*%matrix(theta.init.mat, ncol = 1))
  
  group.p1 <- wholedata$group.p1
  theta.init1.mat <- matrix(0, nrow = q+p, ncol = length(group.p1))
  for (k in 1:length(group.p1)) {
    theta.init1.mat[,k] <- apply(matrix(theta.init.mat[,group.p1[[k]]], nrow = q+p), 1, mean)
  }
  theta.init1 <- as.matrix(kronecker(LMatrixFun(group.p1), bdiag(diag(q+p)))%*%matrix(theta.init1.mat, ncol = 1))
  
  group.p2 <- wholedata$group.p2
  theta.init2.mat <- matrix(0, nrow = q+p, ncol = length(group.p2))
  for (k in 1:length(group.p2)) {
    theta.init2.mat[,k] <- apply(matrix(theta.init.mat[,group.p2[[k]]], nrow = q+p), 1, mean)
  }
  theta.init2 <- as.matrix(kronecker(LMatrixFun(group.p2), bdiag(diag(q+p)))%*%matrix(theta.init2.mat, ncol = 1))
  
  return(list(kmeansres = kmeansres, theta.init0 = theta.init0, theta.init1 = theta.init1, theta.init2 = theta.init2))
}


###########################################################################################################################
## Part V #################################################################################################################
######################################### Functions for applying alternatives #############################################
groupcvxFun <- function(clust){
  group <- vector(mode = 'list', length = length(clust$size))
  for (k in 1:length(group)) {
    group[[k]] <- which(clust$cluster == k)
  }
  return(list(gr.num = length(group), gr.gf = group))
}

hiercvx <- function(wholedata, lambda.x, lambda.z, c = 1, penal.type = '1'){
  cova.x <- wholedata$data.x; q <- nrow(cova.x)
  cova.z <- wholedata$data.z; p <- nrow(cova.z)
  cova <- rbind(cova.x,cova.z); n <- ncol(cova)
  
  w.x <- kernel_weights(X = cova.x, phi = 0.5)
  w.x <- knn_weights(w = w.x, k = 5, n = n)
  res.main <- cvxclust(X = cova.x, w = w.x, gamma = lambda.x, method = 'admm', type = penal.type)
  
  res.list <- vector(mode = 'list',length = length(lambda.x))
  
  for (i1 in 1:length(lambda.x)) {
    beta.gf.mat <- res.main$U[[i1]]
    
    gr.main <- groupcvxFun(find_clusters(create_adjacency(V = res.main$V[[i1]], w = w.x, n = n, method = 'admm')))
    gr.main.num <- gr.main$gr.num
    gr.main.gf <- gr.main$gr.gf
    
    gr.hier.gf.list <- vector(mode = 'list',length = length(lambda.z))
    gamma.gf.mat.list <- vector(mode = 'list',length = length(lambda.z))
    gamma.gf.mat <- matrix(0,p,n)
    gr.hier.gf <- NULL
    BIC.seq <- c()
    for (i2 in 1:length(lambda.z)) {
      for (k in 1:gr.main.num){
        indk <- gr.main.gf[[k]]
        if(length(indk) > 1){
          cova.zk <- cova.z[,indk]
          w <- kernel_weights(X = cova.zk, phi = 0.5)
          if(length(indk) > 5){
            w <- knn_weights(w = w, k = 5, n = length(indk))
          }
          res <- cvxclust(X = cova.zk, w = w, gamma = lambda.z[i2], method = 'admm', type = penal.type)
          gamma.gf.mat[,indk] <- res$U[[1]]
          A <- create_adjacency(V = res$V[[1]], w = w, n = length(indk), method = 'admm')
          gr.gf <- groupcvxFun(find_clusters(A))$gr.gf
          for (m in 1:length(gr.gf)) {
            gr.hier.gf <- c(gr.hier.gf, list(indk[gr.gf[[m]]]))
          }
        }else{
          gamma.gf.mat[,indk] <- cova.z[,indk]
          gr.hier.gf <- c(gr.hier.gf, list(indk))
        }
      }
      residual <- log(sum((cova-rbind(beta.gf.mat,gamma.gf.mat))^2/n))
      BIC.seq <- c(BIC.seq, residual+c*log(log(n))*log(n)*(gr.main.num+length(gr.hier.gf))/n)
      gamma.gf.mat.list[[i2]] <- gamma.gf.mat
      gr.hier.gf.list[[i2]] <- groupFun(pairmatFun(n,pairindFun(gr.hier.gf)))
    }
    res.list[[i1]] <- list(beta.gf.mat = beta.gf.mat, gr.main.gf = gr.main.gf, gamma.gf.mat.list = gamma.gf.mat.list, 
                           gr.hier.gf.list = gr.hier.gf.list, BIC.seq = BIC.seq)
  }
  
  BIC.mat <- matrix(0,length(lambda.x),length(lambda.z))
  for (i1 in 1:length(lambda.x)){
    BIC.mat[i1,] <- res.list[[i1]]$BIC.seq
  }
  ind <- min(which(BIC.mat == min(BIC.mat)))
  indz <- floor((ind-1)/length(lambda.x))+1
  indx <- ind-(indz-1)*length(lambda.x)
  
  beta.gf.mat <- res.list[[indx]]$beta.gf.mat
  gr.main.gf <- res.list[[indx]]$gr.main.gf
  gr.main <- list(gr.num = length(gr.main.gf), gr.gf = gr.main.gf)
  
  gamma.gf.mat <- (res.list[[indx]]$gamma.gf.mat.list)[[indz]]
  gr.hier.gf <- (res.list[[indx]]$gr.hier.gf.list)[[indz]]
  gr.hier <- list(gr.num = length(gr.hier.gf), gr.gf = gr.hier.gf)
  
  xi.gf.list <- vector(mode = 'list', length = gr.main$gr.num)
  for (k in 1:gr.main$gr.num) {
    xi.gf.list[[k]] <- matrix(apply(matrix(beta.gf.mat[,gr.main$gr.gf[[k]]], nrow = q), 1, mean), ncol = 1)
  }
  xi.gf <- do.call(rbind,xi.gf.list)
  xi.gf.mat <- matrix(xi.gf, nrow = q)
  
  alpha.gf.list <- vector(mode = 'list', length = gr.hier$gr.num)
  for (k in 1:gr.hier$gr.num) {
    alpha.gf.list[[k]] <- matrix(apply(matrix(gamma.gf.mat[,gr.hier$gr.gf[[k]]], nrow = p), 1, mean), ncol = 1)
  }
  alpha.gf <- do.call(rbind,alpha.gf.list)
  alpha.gf.mat <- matrix(alpha.gf, nrow = p)
  
  return(list(beta.gf.mat = beta.gf.mat, gamma.gf.mat = gamma.gf.mat, xi.gf.mat = xi.gf.mat, alpha.gf.mat = alpha.gf.mat,
              gr.main = gr.main, gr.hier = gr.hier, BIC.mat = BIC.mat, indx = indx, indz = indz, BIC = BIC.mat[indx,indz], res.list = res.list))
}

































