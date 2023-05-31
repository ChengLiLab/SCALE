library(Matrix)
library(MASS)

findVariableGenes <- function(counts, groups, genelist=rownames(counts),
                              method=c("KW","likelihood","glm"),
                              thresh=null){
  method = match.arg(method)
  if (method == "KW"){
    # find differentially distributed genes by Kruskal-Wallis test
    if (is.null(thresh)){
      print("Please specify a p-value threshold!")
      return(-1)
    }
    return(.KWGenes(counts, groups, genelist=genelist, p.thresh = thresh))
    
  }
  
  if (method == "likelihood"){
    # find genes that could best separate age groups
    if (is.null(thresh)){
      print("Please specify an acuracy threshold!")
      return(-1)
    }
    return(.likelihoodGenes(counts, groups, genelist=genelist, acu.thresh = thresh))
  }
  
  if (method == "glm"){
    # find differentially distributed genes by generalized linear model
    if (is.null(thresh)){
      print("Please specify a coeficient threshold!")
      return(-1)
    }
    return(.glmGenes(counts, groups, genelist=genelist, coef.thresh = thresh))
    
  }
  
}

# find differentially distributed genes by Kruskal-Wallis test
.KWGenes <- function(counts, groups, genelist, p.thresh){
  # Kruskal-Wallis test
  pvalue = c()
  genelist = intersect(rownames(counts), genelist)
  pb = txtProgressBar(1, length(genelist), style = 3)
  
  for (i in 1:length(genelist)){
    setTxtProgressBar(pb, i)
    kw = kruskal.test(counts[genelist[i],], groups)
    pvalue = c(pvalue, kw$p.value)
  }
  close(pb)
  
  return(data.frame(
    genes = genelist[pvalue<=p.thresh],
    pvalue = pvalue[pvalue<=p.thresh]))
}

# find genes that could best separate age groups
.likelihoodGenes <- function(counts, groups, genelist, acu.thresh){
  acuracy = c()
  genelist = intersect(rownames(counts), genelist)
  pb = txtProgressBar(1, length(genelist), style = 3)
  
  for (i in 1:length(genelist)) {
    setTxtProgressBar(pb, i)
    loglh = singleGeneLikelihood(counts[genelist[i],], groups)
    pred = apply(loglh, 1, which.max)
    acu.i = sum(diag(table(pred, groups)))/dim(counts)[2]
    acuracy = c(acuracy, acu.i)
  }
  close(pb)
  
  return(data.frame(
    genes = genelist[acuracy>=acu.thresh],
    acuracy = acuracy[acuracy>=acu.thresh]))
  
}

# find differentially distributed genes by generalized linear model
.glmGenes <- function(counts, groups, genelist, coef.thresh){
  require(glmnet)
  require(caTools)
  
  # empty
 
}


# calculate the likelihood of each age level
singleGeneLikelihood <- function(counts, groups){
  ageLevels = levels(factor(groups))
  likelihood = matrix(0, length(counts), length(ageLevels))
  colnames(likelihood) = ageLevels
  
  for (i in 1:length(ageLevels)){
    d = density(as.numeric(counts[groups==ageLevels[i]]))
    likelihood[,i] = approx(d$x, d$y, xout=counts)[2][[1]]
  }
  return(log10(likelihood))
}

# test the accuracy
test.ac.list <- function(matrix, glist, age, return.lh=F){
  level <- levels(age)
  lh <- matrix(0, dim(matrix)[2], length(level))
  for (i in 1:length(glist)){
    lh.i <- SingleGeneLh(matrix, as.character(glist[i]), age)
    rm <- rowSums(is.na(lh.i))!=0
    lh.i[rm,] <- rep(0, length(level))
    lh <- lh + lh.i
  }  
  
  if (return.lh){
    return(as.data.frame(lh))
  }
  
  prediction <- apply(lh, 1, which.max)
  age.copy <- age
  levels(age.copy) <- c(1,2,3)
  s <- sum(age.copy==prediction)
  return(s/length(age))
  #return(lh)
}

# boostraping 
bootstrap <- function(matrix, age, itertimes,numgenes, glist=rownames(matrix)){
  accu <- c()
  for (i in 1:itertimes){
    accu <- c(accu, test.ac.list(matrix, sample(glist,numgenes), age))  
    if (i%%100==0){print(i)}
  }
  return(accu)
}

# pairwise correlation
cor_plot <- function(matrix, age, g1, g2, plot=T){
  ix <- which(rownames(matrix)==g1)
  jx <- which(rownames(matrix)==g2)
  level <- levels(age)
  cor <- c()
  if (plot){
    par(mfrow=c(1, length(level)))
  }
  for (i in 1:length(level)){
    cor.i <- cor(as.numeric(matrix[ix, age==level[i]]), 
                 as.numeric(matrix[jx, age==level[i]]))
    cor <- c(cor, cor.i)
    if (plot){
      plot(matrix[ix, age==level[i]], matrix[jx, age==level[i]])
    }
  }
  return(cor)
}

# compute the whole correlation matrix
cor.matrix <- function(matrix, age){
  ngenes <- dim(matrix)[1]
  level <- levels(age)
  cormatrix <- array(0, c(ngenes, ngenes, length(level)))
  for (i in 1:ngenes-1){
    for (j in i+1:ngenes){
      for (k in 1:length(level)){
        cormatrix[i,j,k] <- cor(matrix[i, age==level[k]], matrix[j, age==level[k]])
      }
    }
    if (i%%1000==0){
      print(sprintf("%.1f genes have been processed (total %.1f)", i, ngenes))
    }
  }
  return(cormatrix)
}


# entropy
entropy.plugin <- function(freqs, unit=c("log", "log2", "log10"))
{
  unit = match.arg(unit)
  
  freqs = freqs/sum(freqs) # just to make sure ...
  
  H = -sum( ifelse(freqs > 0, freqs*log(freqs), 0) )
  
  if (unit == "log2")  H = H/log(2)  # change from log to log2 scale
  if (unit == "log10") H = H/log(10) # change from log to log10 scale
  
  return(H)
}

# bin number
DoaneBinNumber <- function(data){
  n <- length(data)
  k <- moments::kurtosis(data)
  nbin <- 1+log(n)+log(1+k*(n/6)**0.5)
  return(round(nbin))
}

RiceBinNumber <- function(data){
  return(round(2 * length(data)**(1/3))) 
}

entropy.cell.bin <- function(counts){

  nbin <- RiceBinNumber(counts[,1])
  
  pb = txtProgressBar(1, ncol(counts), style = 3)
  scEntropy <- c()
  for (i in 1:ncol(counts)){
    setTxtProgressBar(pb, i)
    cut <- cut(counts[,i], breaks=nbin)
    freq <- data.frame(table(cut)/length(cut))
    scEntropy <- c(scEntropy, entropy.plugin(freq$Freq))
  }
  close(pb)
  
  return(scEntropy)
}

calcEntropy <- function(counts, threshold = 0.05){
  keep.genes = Matrix::rowSums(counts>0) >= ncol(counts) * 0.05
  return (entropy.cell.bin(counts[keep.genes,]))
}


# Weighted Geomean 
geomean <- function(distributions, weights=NULL){
  if (is.null(weights)){
    weights = rep(1/dim(distributions)[2], dim(distributions)[2])
  }
  return(exp(log(distributions) %*% weights))
}

# Calculate the Barycenter of n distributions by Bregman's method
barycenter.bregman <- function(
  distributions, cost, reg, 
  weights=NULL, numItermax=1000, stopThr=1e-4, verbose=F, log=F){
  # distributions: d*n, n distributions of size d
  # cost: d*d, cost matrix
  # reg: regularization term > 0
  # weight: length = n
  # numItermax: max number of iterations
  # stopThr : Stop threshol on error (>0)
  # verbose : bool. Print information along iterations
  # log : record log if True
  # Returns
  if (is.null(weights)){
    weights = rep(1/dim(distributions)[2], dim(distributions)[2])
  }
  stopifnot(length(weights)==dim(distributions)[2])
  
  if (log){
    l =  list(err=c(),niter=c())
  }
  
  K = exp(-cost/reg)
  
  cpt = 0
  err = 1
  UKv = K %*% (distributions / rowSums(K))
  u = as.vector(geomean(UKv)) / UKv
  
  while (err > stopThr && cpt < numItermax){
    cpt = cpt + 1
    UKv = u * (K %*% (distributions / (K %*% u)))
    u = u * as.vector(geomean(UKv, weights = weights)) / UKv

    if (cpt%%10==1){
      err = sum(genefilter::rowSds(UKv), na.rm = T)
      if (log){
        l$err = append(l$err, err)
      }
      if (verbose){
        if (cpt%%200==0){
          cmd = sprintf("Number of iterations: %d   Error: %.2e", cpt, err)[1]
          print(eval(cmd))
        }
      }
    }
  }
  
  if (log){
    l$niter = cpt
    return(list(barycenter=geomean(UKv, weights), log=l))
  }
  return(geomean(UKv, weights))
}

findInterpolation <- function(
  feature, group, num_bin = 1000, reg=1e-3, replace_na_by = 1e-5, 
  weights=NULL, numItermax=1000, stopThr=1e-4, verbose=F, log=F){
  level = levels(group) # group must be ordered
  num_group = length(level)
  num_bin = num_bin
  
  xrange = seq(min(feature), max(feature), length.out = num_bin)
  cost = as.matrix(dist(xrange))**2
  
  results = matrix(0, nrow = num_bin, ncol = num_group*2 - 1)
  
  for (i in 1:(num_group - 1)){
    d1 = density(as.numeric( feature[group == level[i]] ))
    d2 = density(as.numeric( feature[group == level[i+1]] ))
    d.i1 = approx(d1$x, d1$y, xout=xrange)[2][[1]]
    d.i2 = approx(d2$x, d2$y, xout=xrange)[2][[1]]
    d.i1[is.na(d.i1)] = replace_na_by
    d.i2[is.na(d.i2)] = replace_na_by
    
    # Find the barycenter
    d = cbind(d.i1, d.i2)
    bc = barycenter.bregman(
      distributions = d, cost, reg, weights = weights, numItermax = numItermax,
      stopThr = stopThr, verbose = verbose, log = log)
    #bc = geomean(d)
    
    results[, 2*i - 1] = d.i1
    results[, 2*i] = bc
    results[, 2*i + 1] = d.i2
  }
  
  return(results)
}
