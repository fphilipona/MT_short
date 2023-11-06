source("code/utils.R")

library(randomForest)
library(entropy)
library(gamlss)
library(mgcv)


seed = 26
set.seed(seed)


# read input
x = readRDS("data/mkc_prep.rds")
#input = readRDS("data/mkc_prep.rds")
decon.res = x[["decon"]]
n.cells = x[["n.cells"]]
counts.spe = x[["counts.spe"]]
spe = x[["spe"]]
coords = x[["coords"]]
gbar.mat = x[["gbar.mat"]]
genes = rownames(counts.spe)[which(apply(counts.spe, 1, sum)>9)] # 27 dropped
regions = as.factor(x[["regions"]])

############################################## spot coefficients #############################################
# out of sample random forest predictions (based on number of cells, cell types and cell type vector entropy)
# for the librarysize (sum of all gene expressions) per spot are compared to the true librarysizes.
# the ratio of the two is the spot detection efficiency (for details see the actual master thesis).

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
# the gene detection efficiency, regional effects on the level of the mean and the dispersion parameters can now
# (that everything else has been estimated beforehand) be estimated with gamlss.

genes = genes[-c(1490,2640)] # 2640 does not work
coef.mat = matrix(nrow=length(genes),ncol=5) 
rownames(coef.mat) = genes

t0 = Sys.time()
for (i in 1:length(genes)){ # if error message slightly modify offset +0.01 term (maybe +0.001)
  # for this seed this happens at 2696. I ran the rest (2696:length(genes)) with +0.01 in both 
  # the formula for the mean and for sigma
  gene = genes[i] 
  print(i)
  if (i%%10==0){
    print(Sys.time()-t0)
  }
  df1 = makedf_regions(gene=gene, regions, celltypes.spots=decon.res, gbar.mat=gbar.mat,
                       counts.spe=counts.spe, n.cells=n.cells, spots.coef=spot.coef)
  reg1 = gamlss(formula=expr~offset(log(off+0.01))+region, 
                sigma.formula = ~log(off+0.01), family="NBI", data=df1, trace=FALSE)
  coef.mat[i,] = c(reg1$sigma.coefficients, coef(reg1))
}
# with seed 26, 15 times no convergence warning
colnames(coef.mat) = c("sigma.int", "sigma.off", "gde", "region2", "region3")

coef.mat = extend.coef.mat(coef.mat)

######################################## copula estimation ####################################################
# the copulas can now be estimated. first a distribution transform has to applied to the discrete marginal 
# distributions (which have now been estimated). Then probability integral transforms can be applied to get
# uniform marginal distributions. Taking the correlation of these new distributions gets us the correlation
# matrix defining the gaussian copulas. Correlations are taken for each region separately yielding one copula
# per region. The copulas give an estimated dependency structure per spot. (Details: see actual master thesis.)

h.mat = helper.mat1(regions)
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
  disp.rate = get.rates.regions(gene, coef.mat, decon.res, gbar.mat, spot.coef, n.cells, h.mat)
  u = runif(counts.spe)
  pit.mat.disttrans[i,] = sapply(1:ncol(counts.spe),
                                 function(i) dist.transform(counts.spe, disp.rate, gene, i, u))
}
pit.mat.disttrans[which(pit.mat.disttrans==1)] = 0.999999999
# copula region 0
pit.mat.r0 = pit.mat.disttrans[,which(regions==0)]
gaussian.mat.r0 = matrix(qnorm(c(pit.mat.r0)), nrow=nrow(pit.mat.r0), ncol=ncol(pit.mat.r0))
cor.mat.r0 = cor(t(gaussian.mat.r0))
# copula region 1
pit.mat.r1 = pit.mat.disttrans[,which(regions==1)]
gaussian.mat.r1 = matrix(qnorm(c(pit.mat.r1)), nrow=nrow(pit.mat.r1), ncol=ncol(pit.mat.r1))
cor.mat.r1 = cor(t(gaussian.mat.r1))
# copula region 2
pit.mat.r2 = pit.mat.disttrans[,which(regions==2)]
gaussian.mat.r2 = matrix(qnorm(c(pit.mat.r2)), nrow=nrow(pit.mat.r2), ncol=ncol(pit.mat.r2))
cor.mat.r2 = cor(t(gaussian.mat.r2))

output = list("seed"=seed, "spot.coef"=spot.coef, "coef.mat"=coef.mat, "pit.mat.disttrans"=pit.mat.disttrans,
              "cop0"=cor.mat.r0, "cop1"=cor.mat.r1, "cop2"=cor.mat.r2, "decon" = decon.res, "gbar.mat" = gbar.mat,
              "n.cells" = n.cells, "regions" = regions, "coords" = coords)

library(rlist)
list.save(output, "data/mkc_estimation.rdata")