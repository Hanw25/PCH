#####################################################################################
## Fundamental supporting functions and ADMM function of proposed methods.
#####################################################################################

################################################################################################################################
## Part I ######################################################################################################################
######################################## Some fundamental supporting functions #################################################
grindFun <- function(i, group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Find group index of the i-th subject 
  ## -----------------------------------------------------------------------------------------------------------------
  for (k in 1:length(group)) {
    if(i %in% group[[k]]){
      grind <- k
    }
  }
  return(grind)
}

LMatrixFun <- function(group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Sparse matrix L in computation
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(group)); K <- length(group)
  L <- Matrix(0, nrow = n, ncol = K, sparse = T)
  for (i in 1:n) {
    L[i,grindFun(i, group)] <- 1
  }
  return(L)
}

DMatrixFun <- function(n){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Sparse matrix D in computation
  ## -----------------------------------------------------------------------------------------------------------------
  indx <- combn(n,2)
  D <- Matrix(0, nrow = ncol(indx), ncol = n, sparse = T)
  for (i in 1:ncol(indx)) {
    D[i,indx[1,i]] <- 1
    D[i,indx[2,i]] <- (-1)
  }
  return(D)
}

mcpd <- function(x, lambda, a = 3){
  ## -----------------------------------------------------------------------------------------------------------------
  ## The derivative of MCP function
  ## -----------------------------------------------------------------------------------------------------------------
  if(lambda != 0){
    rho <- lambda*(1 > abs(x)/(lambda*a))*(1-abs(x)/(lambda*a))
  }else{
    rho <- 0
  }
  return(rho)
}

pairindFun <- function(group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate the pair indices of subjects by group structure
  ## -----------------------------------------------------------------------------------------------------------------
  P <- vector(mode = "list")
  if(length(group) == length(unlist(group))){
    P <- matrix(c(0,0), 2, 1)
  }else{
    for (k in 1:length(group)) {
      kgroup <- group[[k]]; m <- length(kgroup); P[[k]] <- vector(mode = "list")
      if(m > 1){
        for (i in 1:(m-1)) {
          P[[k]][[i]] <- rbind(rep(kgroup[i],(m-i)), kgroup[(i+1):m])
        }
        P[[k]] <- do.call(cbind, P[[k]])
      }else{
        P[[k]] <- NULL
      }
    }
    P <- do.call(cbind, P)
  }
  return(P)
}

pairmatFun <- function(n, pairind){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate index upper triangular matrix by the output of pairindFun
  ## -----------------------------------------------------------------------------------------------------------------
  pairmat <- matrix(0, nrow = n, ncol = n)
  if(pairind[1,1] != 0){
    for (i in 1:ncol(pairind)) {
      pairmat[pairind[1,i],pairind[2,i]] <- 1
    }
  }
  return(pairmat)
}

groupFun <- function(pairmat, merge = F){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate a group structure (list) by the output of pairmatFun
  ## -----------------------------------------------------------------------------------------------------------------
  if(merge){
    pairmat2 <- pairmat+t(pairmat); diag(pairmat2) <- 1; n <- nrow(pairmat)
    group <- vector(mode = 'list',length = n)
    for (i in 1:n) {
      group[[i]] <- which(pairmat2[i,] == 1)
    }
    for (i in 1:n) {
      for (j in setdiff(1:n,i)) {
        if(i %in% group[[j]]){
          group[[i]] <- sort(unique(c(group[[i]],group[[j]])))
        }
      }
    }
    group <- unique(group)
    pairmat <- pairmatFun(n, pairindFun(group))
  }
  pairmat2 <- pairmat+t(pairmat); diag(pairmat2) <- 1
  group.unique.matrix <- unique(pairmat2); gr.num <- nrow(group.unique.matrix)
  group <- vector(mode = "list", length = gr.num)
  for (k in 1:gr.num) {
    for (i in 1:nrow(pairmat2)) {
      if(sum(abs(pairmat2[i,]-group.unique.matrix[k,])) == 0){
        group[[k]] <- c(group[[k]],i)
      }
    }
  }
  return(group)
}

get.gr <- function(n, eta.gf, merge = F){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate group structure by estimation of ADMM
  ## -----------------------------------------------------------------------------------------------------------------
  sumdif <- n*(n-1)/2; p <- length(eta.gf)/sumdif; indx <- combn(n,2)
  diff.gf <- apply(matrix(abs(eta.gf), nrow = p, ncol = sumdif), 2, sum)
  cap.mat <- matrix(0,n,n)
  if(length(which(diff.gf == 0)) > 0){
    indx.gf <- as.matrix(indx[,which(diff.gf == 0)])
    for(j in 1:ncol(indx.gf)){
      cap.mat[indx.gf[1,j],indx.gf[2,j]] <- 1
    }
  }
  gr.gf <- groupFun(cap.mat, merge = merge)
  gr.num <- length(gr.gf)
  return(list(gr.num = gr.num, gr.gf = gr.gf))
}


###########################################################################################################################
## Part II ################################################################################################################
########################################### Functions for main algorithms #################################################
ADMM <- function(wholedata, theta.init, lambda1, lambda2, prior,
                 iter.max = 50, epsi = 0.005, gam.mcp = 3, penal.para = 1, c = 1, merge = F){
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The ADMM implementation of proposed method
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ wholedata: The input data analyzed.
  ## @ theta.init: The Initial values of X and Z which associated with the prior.
  ## @ lambda1: The tuning parameter controlling the number of refined subgroups.
  ## @ lambda2: The tuning parameter controlling the number of rough subgroups.
  ## @ prior: The prior setting. '0', '1', '2'.
  ## @ iter.max: int >= 2, maximum number of cycles of the ADMM algorithm, the default setting is 50.
  ## @ epsi: A float value, algorithm termination threshold, the default setting is 0.005.
  ## @ gam.mcp: The regularization parameter in MCP, the default setting is 3.
  ## @ penal.para: The penalty parameter in ADMM algorithm, the default setting is 1.
  ## @ c: The parameter in mBIC, the default setting is 1.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  cova.x <- wholedata$data.x; q <- nrow(cova.x)
  cova.z <- wholedata$data.z; p <- nrow(cova.z)
  cova <- rbind(cova.x,cova.z); n <- ncol(cova)
  
  # Represent prior
  if(prior == '0'){
    group.p <- wholedata$group.p0
  }
  if(prior == '1'){
    group.p <- wholedata$group.p1
  }
  if(prior == '2'){
    group.p <- wholedata$group.p2
  }
  l.matrix <- LMatrixFun(group.p)
  j.matrix <- kronecker(l.matrix, bdiag(diag(q+p)))
  d.matrix <- DMatrixFun(n)
  h.matrix <- kronecker(d.matrix, bdiag(diag(q+p)))
  one.matrix <- bdiag(rep(list(rep(1,q+p)),(n*(n-1)/2)))
  one.matrix.main <- bdiag(rep(list(rep(1,q)),(n*(n-1)/2)))
  one.matrix.hier <- bdiag(rep(list(rep(1,p)),(n*(n-1)/2)))
  
  main.n <- c(1:q)
  hier.n <- setdiff(1:(q+p), main.n)
  a.main.n <- sort(rep((0:(n-1))*(q+p),q))+main.n
  a.hier.n <- sort(rep((0:(n-1))*(q+p),p))+hier.n
  d.main.n <- sort(rep((0:(n*(n-1)/2-1))*(q+p),q))+main.n
  d.hier.n <- sort(rep((0:(n*(n-1)/2-1))*(q+p),p))+hier.n
  
  # initialization  
  zeta.init <- h.matrix%*%theta.init
  w.init <- Matrix(0, nrow = (q+p)*n*(n-1)/2, ncol = 1, sparse = T) 
  
  iter <- 1
  theta.est.list <- vector(mode = "list", length = iter.max); theta.est.list[[iter]] <- theta.init 
  zeta.est.list <- vector(mode = "list", length = iter.max); zeta.est.list[[iter]] <- zeta.init 
  w.est.list <- vector(mode = "list", length = iter.max); w.est.list[[iter]] <- w.init 
  rm(theta.init,zeta.init,w.init)
  gc()
  
  # iteration
  while(iter <= iter.max){
    iter <- iter+1
    if(iter > iter.max){
      cat("The step iteration number exceed the maximum values!\n")
      cat("The eps is now", eps, "\n")
      break
    }
    
    # update theta
    theta.est.list[[iter]] <- (cova+matrix(penal.para*zeta.est.list[[iter-1]]-w.est.list[[iter-1]], nrow = q+p)%*%d.matrix)%*%l.matrix%*%solve(t(l.matrix)%*%l.matrix+penal.para*t(l.matrix)%*%t(d.matrix)%*%d.matrix%*%l.matrix)
    theta.est.list[[iter]] <- j.matrix%*%matrix(theta.est.list[[iter]], ncol = 1)
    
    # hierarchical groupwise thresholding operator
    zeta.i <- h.matrix%*%theta.est.list[[iter]]
    w.i <- w.est.list[[iter-1]]
    zeta.tem1 <- zeta.i+w.i/penal.para
    zeta.tem1.main <- as.matrix(zeta.tem1[d.main.n,], sparse = T)
    zeta.tem1.hier <- as.matrix(zeta.tem1[d.hier.n,], sparse = T)
    zeta.i1 <- zeta.i
    zeta.tem2norm.main <- sqrt(t(one.matrix.main)%*%(zeta.tem1.main^2))
    zeta.tem2norm.hier <- sqrt(t(one.matrix.hier)%*%(zeta.tem1.hier^2))
    zeta.tem2norm <- sqrt(zeta.tem2norm.main^2 + zeta.tem2norm.hier^2)
    

    num.n.n <- which(zeta.tem2norm > gam.mcp*lambda1 & zeta.tem2norm.main > gam.mcp*lambda2)

    oz.coe <- apply(1-lambda1/penal.para/zeta.tem2norm,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
    oo <- oz.coe*zeta.tem2norm.main
    num.s.n <- which(zeta.tem2norm <= gam.mcp*lambda1 & oo > gam.mcp*lambda2)

    oo.coe <- apply(1-lambda2/penal.para/zeta.tem2norm.main,1,function(x) max(x,0)) / (1-1/(gam.mcp*penal.para))
    ee <- sqrt(zeta.tem2norm.hier^2+(oo.coe*zeta.tem2norm.main)^2)
    num.n.s <- which(ee > gam.mcp*lambda1 & zeta.tem2norm.main <= gam.mcp*lambda2)

    num.s.s <- setdiff(1:(n*(n-1)/2), Reduce(union,list(num.n.n,num.s.n,num.n.s)))
    
    if(length(num.n.n) < 2){
      num.n.n2 <- which(rowSums(as.matrix(one.matrix[,num.n.n]))!=0)
    }else{
      num.n.n2 <- which(rowSums(one.matrix[,num.n.n])!=0)
    }
    if(length(num.s.n) < 2){
      num.s.n2 <- which(rowSums(as.matrix(one.matrix[,num.s.n]))!=0)
    }else{
      num.s.n2 <- which(rowSums(one.matrix[,num.s.n])!=0)
    }
    if(length(num.n.s) < 2){
      num.n.s2 <- which(rowSums(as.matrix(one.matrix[,num.n.s]))!=0)
    }else{
      num.n.s2 <- which(rowSums(one.matrix[,num.n.s])!=0)
    }
    if(length(num.s.s) < 2){
      num.s.s2 <- which(rowSums(as.matrix(one.matrix[,num.s.s]))!=0)
    }else{
      num.s.s2 <- which(rowSums(one.matrix[,num.s.s])!=0)
    }
    
    if(length(num.n.n2) > 0){
      zeta.i1[num.n.n2,] <- zeta.i[num.n.n2,]
    }

    num <- num.s.n; num2 <- num.s.n2
    if(length(num) > 0){
      zeta.tem3 <- as.matrix(oz.coe[num])
      zeta.tem4 <- as.vector(apply(zeta.tem3, 1, function(x) rep(x,q+p)))
      zeta.i1[num2,] <- zeta.tem4*zeta.i[num2,]
    }

    num <- num.n.s; num2 <- num.n.s2
    if(length(num) > 0){
      num.main <- num2[sort(rep((1:length(num)-1)*(q+p),q)) + main.n]
      num.hier <- num2[sort(rep((1:length(num)-1)*(q+p),p)) + hier.n]
      zeta.tem3 <- as.matrix(oo.coe[num])
      zeta.tem4 <- as.vector(apply(zeta.tem3, 1, function(x) rep(x,q)))
      zeta.i1[num.main,] <- zeta.tem4*zeta.i[num.main,]
      zeta.i1[num.hier,] <- zeta.i[num.hier,]
    }

    num <- num.s.s; num2 <- num.s.s2
    if(length(num) > 0){
      num.main <- num2[sort(rep((1:length(num)-1)*(q+p),q)) + main.n]
      num.hier <- num2[sort(rep((1:length(num)-1)*(q+p),p)) + hier.n]
      zeta.i0norm <- sqrt(t(one.matrix)%*%(zeta.est.list[[iter-1]]^2))
      zeta.i0norm.main <- sqrt(t(one.matrix.main)%*%(zeta.est.list[[iter-1]][d.main.n,]^2))
      mcp.z <- mcpd(zeta.i0norm[num],lambda1,gam.mcp)/penal.para/(zeta.i0norm[num]+10^(-7))
      mcp.o <- mcpd(zeta.i0norm.main[num],lambda2,gam.mcp)/penal.para/(zeta.i0norm.main[num]+10^(-7))
      zeta.tem3.main <- as.matrix(mcp.z+mcp.o)
      zeta.tem3.hier <- as.matrix(mcp.z)
      zeta.tem4.main <- as.vector(apply(zeta.tem3.main, 1, function(x) rep(x,q)))
      zeta.tem4.hier <- as.vector(apply(zeta.tem3.hier, 1, function(x) rep(x,p)))
      zeta.i1[num.main,] <- zeta.i[num.main,]/(1+zeta.tem4.main)
      zeta.i1[num.hier,] <- zeta.i[num.hier,]/(1+zeta.tem4.hier)
      zeta.i1[abs(zeta.i1) < 10^(-9)] <- 0
    }
    # update zeta
    zeta.est.list[[iter]] <- zeta.i1
    
    # update w
    w.est.list[[iter]] <- w.est.list[[iter-1]]+penal.para*(h.matrix%*%theta.est.list[[iter]]-zeta.est.list[[iter]])
    
    eps <- sqrt(sum((h.matrix%*%theta.est.list[[iter]]-zeta.est.list[[iter]])^2)/((q+p)*n*(n-1)/2))
    if(eps < epsi){
      cat("The iterations finish in", iter, "-th step.\n")
      break
    }
    
    theta.est.list[iter-1] <- list(NULL)
    zeta.est.list[iter-1] <- list(NULL)
    w.est.list[iter-1] <- list(NULL)
  }
  
  ##################
  if(iter > iter.max){
    theta.gf <- as.matrix(theta.est.list[[iter-1]])
    zeta.gf <- as.matrix(zeta.est.list[[iter-1]])
  }else{
    theta.gf <- as.matrix(theta.est.list[[iter]])
    zeta.gf <- as.matrix(zeta.est.list[[iter]])
  }
  
  theta.gf.mat <- matrix(theta.gf, nrow = q+p)
  beta.gf.mat <- theta.gf.mat[1:q,]
  gamma.gf.mat <- theta.gf.mat[(q+1):(q+p),]
  
  zeta.gf[which(abs(zeta.gf[,1])<10^(-5)),1] <- rep(0,length(which(abs(zeta.gf[,1])<10^(-5))))
  omega.gf <- as.matrix(zeta.gf[d.main.n,])
  eta.gf <- as.matrix(zeta.gf[d.hier.n,])
  
  rm(theta.est.list, zeta.est.list, w.est.list)
  gc()
  
  # get rough/refined subgroup structure
  gr.main <- get.gr(n, omega.gf, merge = merge)
  gr.hier <- get.gr(n, eta.gf, merge = merge)
  
  # calculate mBIC
  residual <- log(sum((matrix(cova,ncol=1)-theta.gf)^2/n))
  BIC <- residual+c*log(log(n))*log(n)*(gr.main$gr.num+gr.hier$gr.num)/n
  
  # calculate centroids of rough subgroup structure
  xi.gf.list <- vector(mode = 'list', length = gr.main$gr.num)
  for (k in 1:gr.main$gr.num) {
    xi.gf.list[[k]] <- matrix(apply(matrix(beta.gf.mat[,gr.main$gr.gf[[k]]], nrow = q), 1, mean), ncol = 1)
  }
  xi.gf <- do.call(rbind,xi.gf.list)
  xi.gf.mat <- matrix(xi.gf, nrow = q)
  
  # calculate centroids of refined subgroup structure
  alpha.gf.list <- vector(mode = 'list', length = gr.hier$gr.num)
  for (k in 1:gr.hier$gr.num) {
    alpha.gf.list[[k]] <- matrix(apply(matrix(gamma.gf.mat[,gr.hier$gr.gf[[k]]], nrow = p), 1, mean), ncol = 1)
  }
  alpha.gf <- do.call(rbind,alpha.gf.list)
  alpha.gf.mat <- matrix(alpha.gf, nrow = p)
  
  return(list(beta.gf.mat = beta.gf.mat, gamma.gf.mat = gamma.gf.mat, xi.gf.mat = xi.gf.mat, alpha.gf.mat = alpha.gf.mat,
              omega.gf = omega.gf, eta.gf = eta.gf, gr.main = gr.main, gr.hier = gr.hier, BIC = BIC))
}




