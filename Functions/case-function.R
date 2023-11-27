#####################################################################################
# Functions for real data analysis.
#####################################################################################
casetuninglambda <- function(wholedata, theta.init, para, prior,
                             iter.max = 50, epsi = 0.005, gam.mcp = 3, penal.para = 1, c = 1, merge = F, line = F){
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
        peri <- data.frame('lambda1' = lambda1[j], 'lambda2' = lambda2[i], 'BIC' = admmres[[(i-1)*L1+j]]$BIC, 
                           'gr.main.num' = admmres[[(i-1)*L1+j]]$gr.main$gr.num, 'gr.hier.num' = admmres[[(i-1)*L1+j]]$gr.hier$gr.num)
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
      peri <- data.frame('lambda1' = 0, 'lambda2' = lambda2[i], 'BIC' = admmres[[i]]$BIC, 
                         'gr.main.num' = admmres[[i]]$gr.main$gr.num, 'gr.hier.num' = admmres[[i]]$gr.hier$gr.num)
      per.df <- rbind(per.df, peri)
    }
    lambda2.opt <- per.df[min(which(per.df$BIC == min(per.df$BIC))),]$lambda2
    for (j in 1:L1) {
      cat('-----------', 'fixed optimal lambda2', lambda2.opt, 'and', j, '-th lambda1', lambda1[j], '--------------\n')
      admmres[[L2+j]] <- ADMM(wholedata = wholedata, theta.init = theta.init, lambda1 = lambda1[j], lambda2 = lambda2.opt, prior = prior,
                              iter.max = iter.max, epsi = epsi, gam.mcp = gam.mcp, penal.para = penal.para, c = c, merge = merge)
      peri <- data.frame('lambda1' = lambda1[j], 'lambda2' = lambda2.opt, 'BIC' = admmres[[L2+j]]$BIC, 
                         'gr.main.num' = admmres[[L2+j]]$gr.main$gr.num, 'gr.hier.num' = admmres[[L2+j]]$gr.hier$gr.num)
      per.df <- rbind(per.df, peri)
    }
  }

  res.tune <- admmres[[min(which(per.df$BIC == min(per.df$BIC)))]]
  per.df.tune <- per.df[min(which(per.df$BIC == min(per.df$BIC))),]
  
  return(list(admmres = admmres, per.df = per.df, res.tune = res.tune, per.df.tune = per.df.tune))
}

caseinitFun <- function(wholedata, kmin = 2, kmax = 3, divmin = 2, divmax = 3, c = 1){
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
  
  return(list(kmeansres = kmeansres, theta.init0 = theta.init0, theta.init1 = theta.init1))
}

grouprm <- function(x, group){
  n <- length(unlist(group))
  new <- setdiff(1:n,x)
  group.new <- NULL
  for (i in 1:length(group)) {
    a <- setdiff(group[[i]],x)
    if(length(a) == 0){
      group.new <- group.new
    }
    if(length(a) > 0){
      group.new <- c(group.new,list(match(a, new)))
    }
  }
  group.new <- groupFun(pairmatFun(length(new),pairindFun(group.new)))
  return(group.new)
}

datarm <- function(a, wholedata){
  allID <- colnames(wholedata$data.x)
  data.x <- wholedata$data.x[,-a]
  data.z <- wholedata$data.z[,-a]
  group.p0 <- grouprm(a, wholedata$group.p0)
  group.p1 <- grouprm(a, wholedata$group.p1)
  wholedatarm <- list(data.x = data.x, data.z = data.z, group.p0 = group.p0, group.p1 = group.p1, rmID = allID[a], rmind = a)
  return(wholedatarm)
}

grriFun <- function(gr.gf, group.true){
  n <- length(unlist(gr.gf))
  cap.gf.matrix <- pairmatFun(n, pairindFun(gr.gf))
  cap.true.matrix <- pairmatFun(n, pairindFun(group.true))
  gr.tp <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.gf.matrix) == 1)))
  gr.fp <- length(intersect(which(as.vector(cap.true.matrix) == 0), which(as.vector(cap.gf.matrix) == 1)))
  gr.fn <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.gf.matrix) == 0)))
  gr.tn <- n*(n-1)/2-gr.tp-gr.fp-gr.fn
  gr.ri <- (gr.tp+gr.tn)/(n*(n-1)/2)
  return(gr.ri)
}








