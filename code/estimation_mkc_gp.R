source("code/utils.R")

library(randomForest)
library(entropy)
library(gamlss)
library(mgcv)
library(gamlss.add)
library(Rfast)

# read input
load("data/mkc_prep_est/mkc_prep.rdata")
#input = readRDS("data/mkc_prep.rds")
decon.res = x[["decon"]]
n.cells = x[["n.cells"]]
counts.spe = x[["counts.spe"]]
spe = x[["spe"]]
coords = x[["coords"]]
gbar.mat = x[["gbar.mat"]]
genes = rownames(counts.spe)[which(apply(counts.spe, 1, sum)>9)] # 27 dropped
regions = as.factor(x[["regions"]])

seed = 125
set.seed(seed)

# spot coef
entropies = apply(decon.res, 1, entropy)

spot.coef.iter = function(spots, data, ntree){
  rf1 = randomForest(libsizes~., data=data[!rownames(data) %in% spots,], ntree=ntree)
  return(data[spots,]$libsizes/predict(rf1, newdata=data[spots,]))
}
get.spot.coef = function(celltypes.spots, n.cells, counts.spe, ntree=500){
  n.spots = dim(counts.spe)[2]
  spots = sample(colnames(counts.spe))
  libsizes = apply(counts.spe, 2, sum)
  n.cells.types = sweep(celltypes.spots, 1, n.cells, FUN="*")
  df1 = data.frame(n.cells.types)
  df1$libsizes = libsizes; df1$entropies = entropies
  n.iter = n.spots%/%30
  spot.coef.vector = vector()
  for (i in 1:n.iter){
    print(i)
    spot.coef.vector = c(spot.coef.vector, spot.coef.iter(spots[(1+(i-1)*30):(i*30)], df1, ntree))
  }
  if (n.iter*30 < n.spots){
    spot.coef.vector = c(spot.coef.vector, spot.coef.iter(spots[(1+n.iter*30):n.spots], df1, ntree))
  }
  return(spot.coef.vector)
}
spot.coef = get.spot.coef(celltypes.spots = decon.res, n.cells = n.cells, counts.spe = counts.spe, ntree=200)

###################################### rest of marginal parameters #################################################
coef.mat = matrix(nrow=length(genes), ncol=3+ncol(counts.spe))
for (i in 1:length(genes)){
  gene = genes[i]
  if (i%%10==0){
    print(Sys.time())
    print(i)
  }
  df1 = makedf_gp(gene=gene, coords, celltypes.spots=decon.res, gbar.mat=gbar.mat,
                  counts.spe=counts.spe, n.cells=n.cells, spots.coef=spot.coef)
  reg1 = gamlss(formula=expr~ga(~s(coord1, coord2, bs="gp", k=20))+offset(log(off+0.01)), 
                sigma.formula = ~log(off+0.01), family=NBI(mu.link=log, sigma.link=log), data=df1, trace=FALSE)
  coef.mat[i,] = c(reg1$sigma.coefficients, coef(reg1)[1], reg1$mu.s)
}
rownames(coef.mat) = genes

######################################## copula estimation ####################################################
get.rates.gp = function(gene, coef.mat, celltypes.spots, gbar.mat, spots.coef, n.cells){
  offset = n.cells*celltypes.spots%*%gbar.mat[,gene]*spots.coef
  sigma = coef.mat[gene, 1] + coef.mat[gene, 2]*log(offset+0.01)
  intercept.and.spatial.effect = (offset+0.01)*exp(coef.mat[gene, 4:dim(coef.mat)[2]]+coef.mat[gene, 3])
  return(cbind(exp(sigma), intercept.and.spatial.effect))
}

# distribution transform
dist.transform = function(counts.spe, disp.rate, gene, i, u){
  ui = u[i]
  expr = counts.spe[gene,i]
  mu1 = disp.rate[i,2]+0.0001
  sigma1 = disp.rate[i,1]+0.0001
  lower = ifelse(expr==0, 0, ui*pNBI(expr-1, mu1, sigma1))
  return(lower+(1-ui)*pNBI(expr, mu1, sigma1))
}

# pit mat
pit.mat.disttrans = matrix(nrow=length(genes), ncol=ncol(counts.spe))
rownames(pit.mat.disttrans) = genes
colnames(pit.mat.disttrans) = colnames(counts.spe)
for (i in 1:length(genes)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  disp.rate = get.rates.gp(gene, coef.mat, decon.res, gbar.mat, spot.coef, n.cells)
  u = runif(counts.spe)
  pit.mat.disttrans[i,] = sapply(1:ncol(counts.spe),
                                 function(i) dist.transform(counts.spe, disp.rate, gene, i, u))
}
pit.mat.disttrans[which(pit.mat.disttrans==1)] = 0.999999999
gaussian.mat = matrix(qnorm(c(pit.mat.disttrans)), nrow=nrow(pit.mat.disttrans), ncol=ncol(pit.mat.disttrans))
cor.mat = cora(t(pit.mat.disttrans))

output = list("seed"=seed, "spot.coef"=spot.coef, "coef.mat"=coef.mat, "pit.mat.disttrans"=pit.mat.disttrans,
              "decon" = decon.res, "gbar.mat" = gbar.mat, "n.cells" = n.cells, "coords" = coords, "cop" = cor.mat)

saveRDS(output, "data/mkc_estimation_gp.rds")









