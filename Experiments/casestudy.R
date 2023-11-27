################################################################################
rm(list=ls())
gc()
library(Matrix)
library(MASS)
library(NbClust)
library(cvxclustr)
source("main-function.R")
source("case-function.R")
source("sim-function.R")

################################# Part I #######################################
# Apply the proposed method and the alternatives
load("LUAD.data.RData")

data.x <- LUAD.data$data.x.scale
data.z <- LUAD.data$data.z
group.p <- LUAD.data$group.p
ID <- row.names(data.x)
group.p1 <- group.p
for (i in 1:length(group.p)) {
  group.p1[[i]] <- match(group.p[[i]],ID)
}
group.p1 <- groupFun(pairmatFun(length(ID),pairindFun(group.p1)))
group.p0 <- vector(mode = "list",length = length(ID))
for (i in 1:length(group.p0)) {
  group.p0[[i]] <- i
}

wholedata <- list(data.x = t(data.x), data.z = t(data.z), group.p0 = group.p0, group.p1 = group.p1)

init <- caseinitFun(wholedata)
para <- list(lambda1 = seq(0.4,0.7,length=7), lambda2 = seq(0.1,0.5,length=5))
result <- casetuninglambda(wholedata, theta.init = init$theta.init1, para, prior = "1", merge = T)

lambda.x <- 10*seq(0.8,1.2,length=5)
lambda.z <- 10^seq(8,12,length=5)
other.result1 <- hiercvx(wholedata, lambda.x, lambda.z, penal.type = "1")
other.result2 <- hiercvx(wholedata, lambda.x, lambda.z, penal.type = "2")
################################################################################

################################# Part II ######################################
# Use ANOVA to check the significance of our clustering on some clinical features
# load("wholedata.RData")
# load("result.RData")
result <- result$res.tune
gc()

ID <- colnames(wholedata$data.x)
group.main <- result$gr.main$gr.gf
group.hier <- result$gr.hier$gr.gf
group.main.label <- rep(0,355)
group.hier.label <- rep(0,355)
for (i in 1:2) {
  group.main.label[group.main[[i]]] <- rep(i,length(group.main[[i]]))
}
for (i in 1:4) {
  group.hier.label[group.hier[[i]]] <- rep(i,length(group.hier[[i]]))
}

# all clinical data
alldata <- read.csv("data_bcr_clinical_data_patient_match.csv")
na.split2 <- strsplit(alldata[,1], split = '-')
for (j in 1:length(na.split2)) {
  alldata[j,1] <- paste(na.split2[[j]][1],na.split2[[j]][2],na.split2[[j]][3], sep = ".")
}
alldata <- alldata[match(ID, alldata[,1]),]
rownames(alldata) <- NULL
alldata <- cbind(data.frame('group.main.label' = group.main.label, 'group.hier.label' = group.hier.label), alldata)

for (i in 1:nrow(alldata)) {
  for (j in 1:ncol(alldata)) {
    if(alldata[i,j] == '[Not Available]' || alldata[i,j] == '[Not Applicable]'){
      alldata[i,j] <- NA
    }
  }
}

aov.data1 <- na.omit(alldata[,c('group.main.label','group.hier.label','PATIENT_ID','DFS_MONTHS')])
x1.main <- aov(DFS_MONTHS~group.main.label, data = aov.data1)
x1.hier <- aov(DFS_MONTHS~group.hier.label, data = aov.data1)
summary(x1.main) # 0.0558
summary(x1.hier) # 0.00765
DFS.MONTHS <- list(aov.data1 = aov.data1, x1.main = x1.main, x1.hier = x1.hier)

aov.data2 <- na.omit(alldata[,c('group.main.label','group.hier.label','PATIENT_ID','OS_MONTHS')])
x2.main <- aov(OS_MONTHS~group.main.label, data = aov.data2)
x2.hier <- aov(OS_MONTHS~group.hier.label, data = aov.data2)
summary(x2.main) # 0.0169
summary(x2.hier) # 0.00437
OS.MONTHS <- list(aov.data2 = aov.data2, x2.main = x2.main, x2.hier = x2.hier)
################################################################################

################################ Part III ######################################
# Clustering stability
# load("wholedata.RData")
group.p1 <- wholedata$group.p1
noprior.ind <- NULL
for (i in 1:length(group.p1)) {
  if(length(group.p1[[i]]) == 1){
    noprior.ind <- c(noprior.ind, group.p1[[i]])
  }
}
noprior.ind <- sort(noprior.ind)
m <- 100
para <- list(lambda1 = seq(0.4,0.7,length=7), lambda2 = seq(0.1,0.5,length=5))
lambda.x <- 10*seq(0.8,1.2,length=5)
lambda.z <- 10^seq(8,12,length=5)
################################################################################
## Example: Remove 5 data points
a <- 5
wholedatarm <- vector(mode = "list", length = m)
for (i in 1:m) {
  wholedatarm[[i]] <- datarm(sort(sample(noprior.ind,a)),wholedata)
}
initrm <- vector(mode = "list", length = m)
for (i in 1:m) {
  initrm[[i]] <- caseinitFun(wholedatarm[[i]])
}

resultrm <- vector(mode = "list", length = m)
for (i in 1:m) {
  cat('-----------', i, '-th data' ,'--------------\n')
  resultrm[[i]] <- casetuninglambda(wholedatarm[[i]], theta.init = initrm[[i]]$theta.init1, para, prior = "1", merge = T)[-1]
  gc()
}

other.result1rm <- vector(mode = "list", length = m)
for (i in 1:m) {
  cat('-----------', i, '-th data' ,'--------------\n')
  other.result1rm[[i]] <- hiercvx(wholedatarm[[i]], lambda.x, lambda.z, penal.type = "1")[-11]
  gc()
}

other.result2rm <- vector(mode = "list", length = m)
for (i in 1:m) {
  cat('-----------', i, '-th data' ,'--------------\n')
  other.result2rm[[i]] <- hiercvx(wholedatarm[[i]], lambda.x, lambda.z, penal.type = "2")[-11]
  gc()
}

# Calculate results (rm5)
result <- result$res.tune
gc()
group.main <- result$gr.main$gr.gf
group.hier <- result$gr.hier$gr.gf

similarity.main <- NULL
similarity.hier <- NULL
for (i in 1:m) {
  rmind <- wholedatarm[[i]]$rmind
  similarity.main <- c(similarity.main, grriFun(resultrm[[i]]$res.tune$gr.main$gr.gf, grouprm(rmind, group.main)))
  similarity.hier <- c(similarity.hier, grriFun(resultrm[[i]]$res.tune$gr.hier$gr.gf, grouprm(rmind, group.hier)))
}
similarityrm <- data.frame("similarity.main" = similarity.main, "similarity.hier" = similarity.hier)

similarity.main <- NULL
similarity.hier <- NULL
for (i in 1:m) {
  rmind <- wholedatarm[[i]]$rmind
  similarity.main <- c(similarity.main, grriFun(other.result1rm[[i]]$gr.main$gr.gf, grouprm(rmind, group.main)))
  similarity.hier <- c(similarity.hier, grriFun(other.result1rm[[i]]$gr.hier$gr.gf, grouprm(rmind, group.hier)))
}
other.similarity1rm <- data.frame("similarity.main" = similarity.main, "similarity.hier" = similarity.hier)

similarity.main <- NULL
similarity.hier <- NULL
for (i in 1:m) {
  rmind <- wholedatarm[[i]]$rmind
  similarity.main <- c(similarity.main, grriFun(other.result2rm[[i]]$gr.main$gr.gf, grouprm(rmind, group.main)))
  similarity.hier <- c(similarity.hier, grriFun(other.result2rm[[i]]$gr.hier$gr.gf, grouprm(rmind, group.hier)))
}
other.similarity2rm <- data.frame("similarity.main" = similarity.main, "similarity.hier" = similarity.hier)
################################################################################












