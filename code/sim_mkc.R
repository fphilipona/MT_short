source("code/utils.R")

library(MASS)
library(gamlss)

load("data/mkc_estimation.rdata")
cop0 = x[["cop0"]]
cop1 = x[["cop1"]]
cop2 = x[["cop2"]]
coef.mat = x[["coef.mat"]]
spot.coef = x[["spot.coef"]]
genes = rownames(x[["coef.mat"]])
coords = x[["coords"]]
regions = x[["regions"]]
n.cells = x[["n.cells"]]
gbar.mat = x[["gbar.mat"]]
decon.res = x[["decon"]]
h.mat = helper.mat1(regions)
spots = rownames(decon.res)

seed = 27
set.seed(seed)

# simulate from exactly same distribution as estimated.
# For each region first draw from a multivariate normal with mean all zeroes and covariance matrix
# the estimated copula correlation. Then, by applying standard normal cdfs, get uniformly distributed
# marginal distributions on [0,1]. Then, apply the estimated inverse cdfs to those uniforms to get the
# final gene expressions. (Details in actual master thesis.)

# simulate region 0
gaussian.draw.r0 = mvrnorm(length(which(regions==0)), mu=rep(0,dim(cop0)[1]), Sigma=cop0)
gaussian.draw.r0 = t(gaussian.draw.r0)
colnames(gaussian.draw.r0)=names(regions[which(regions==0)])
rownames(gaussian.draw.r0)=genes
# gaussian draw to PIT (probability integral transform. this gives marginal uniforms on [0,1])
pit.draw.r0 = matrix(nrow=nrow(gaussian.draw.r0), ncol=ncol(gaussian.draw.r0), 
                     dimnames=list(genes, names(regions[which(regions==0)])))
for (i in 1:nrow(pit.draw.r0)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  pit.draw.r0[i,] = sapply(1:ncol(gaussian.draw.r0),
                           function(i) pnorm(gaussian.draw.r0[gene,i]))
}
# PIT to simulated expressions (apply inverse estimated cdfs qNBI)
expr.mat.r0 = matrix(nrow=length(genes), ncol=length(which(regions==0)))
rownames(expr.mat.r0) = genes
colnames(expr.mat.r0) = names(regions[which(regions==0)])
for (i in 1:length(genes)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  disp.rate1 = get.rates.regions(gene, coef.mat, decon.res, gbar.mat, spot.coef, n.cells, h.mat)
  disp.rate1 = disp.rate1[which(regions==0),]
  expr.mat.r0[i,] = sapply(1:ncol(pit.draw.r0),
                           function(i) qNBI(pit.draw.r0[gene,i], mu=disp.rate1[i,2]+0.000001, sigma=disp.rate1[i,1]+0.000001))
  
}

# simulate region 1
gaussian.draw.r1 = mvrnorm(length(which(regions==1)), mu=rep(0,dim(cop1)[1]), Sigma=cop1)
gaussian.draw.r1 = t(gaussian.draw.r1)
colnames(gaussian.draw.r1)=names(regions[which(regions==1)])
rownames(gaussian.draw.r1)=genes
# gaussian draw to PIT
pit.draw.r1 = matrix(nrow=nrow(gaussian.draw.r1), ncol=ncol(gaussian.draw.r1), 
                     dimnames=list(genes, names(regions[which(regions==1)])))
for (i in 1:nrow(pit.draw.r1)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  pit.draw.r1[i,] = sapply(1:ncol(gaussian.draw.r1),
                           function(i) pnorm(gaussian.draw.r1[gene,i]))
}
# PIT to simulated expressions
expr.mat.r1 = matrix(nrow=length(genes), ncol=length(which(regions==1)))
rownames(expr.mat.r1) = genes
colnames(expr.mat.r1) = names(regions[which(regions==1)])
for (i in 1:length(genes)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  disp.rate1 = get.rates.regions(gene, coef.mat, decon.res, gbar.mat, spot.coef, n.cells, h.mat)
  disp.rate1 = disp.rate1[which(regions==1),]
  expr.mat.r1[i,] = sapply(1:ncol(pit.draw.r1),
                           function(i) qNBI(pit.draw.r1[gene,i], mu=disp.rate1[i,2]+0.000001, sigma=disp.rate1[i,1]+0.000001))
  
}

# simulate region 2
gaussian.draw.r2 = mvrnorm(length(which(regions==2)), mu=rep(0,dim(cop2)[1]), Sigma=cop2)
gaussian.draw.r2 = t(gaussian.draw.r2)
colnames(gaussian.draw.r2)=names(regions[which(regions==2)])
rownames(gaussian.draw.r2)=genes
# gaussian draw to PIT
pit.draw.r2 = matrix(nrow=nrow(gaussian.draw.r2), ncol=ncol(gaussian.draw.r2), 
                     dimnames=list(genes, names(regions[which(regions==2)])))
for (i in 1:nrow(pit.draw.r2)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  pit.draw.r2[i,] = sapply(1:ncol(gaussian.draw.r2),
                           function(i) pnorm(gaussian.draw.r2[gene,i]))
}
# PIT to simulated expressions
expr.mat.r2 = matrix(nrow=length(genes), ncol=length(which(regions==2)))
rownames(expr.mat.r2) = genes
colnames(expr.mat.r2) = names(regions[which(regions==2)])
for (i in 1:length(genes)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  disp.rate1 = get.rates.regions(gene, coef.mat, decon.res, gbar.mat, spot.coef, n.cells, h.mat)
  disp.rate1 = disp.rate1[which(regions==2),]
  expr.mat.r2[i,] = sapply(1:ncol(pit.draw.r2),
                           function(i) qNBI(pit.draw.r2[gene,i], mu=disp.rate1[i,2]+0.000001, sigma=disp.rate1[i,1]+0.000001))
  
}

expr.mat = cbind(expr.mat.r0, expr.mat.r1, expr.mat.r2)
expr.mat = expr.mat[,spots]

output = list("seed"=seed, "expr.mat"=expr.mat, "coords" = coords)

saveRDS(output, "data/mkc_sim_gt0.rds")



