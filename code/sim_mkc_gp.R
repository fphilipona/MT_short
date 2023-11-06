source("code/utils.R")

library(MASS)
library(gamlss)

input = readRDS("data/mkc_estimation_gp.rds")
cop = input[["cop"]]
coef.mat = input[["coef.mat"]]
spot.coef = input[["spot.coef"]]
genes = rownames(input[["coef.mat"]])
coords = input[["coords"]]
n.cells = input[["n.cells"]]
gbar.mat = input[["gbar.mat"]]
decon.res = input[["decon"]]
spots = rownames(decon.res)

seed = 27
set.seed(seed)

# simulate from exactly same distribution as estimated
gaussian.draw = mvrnorm(length(spots), mu=rep(0,dim(cop)[1]), Sigma=cop)
gaussian.draw = t(gaussian.draw)
colnames(gaussian.draw) = spots
rownames(gaussian.draw) = genes
# gaussian draw to PIT
pit.draw = matrix(nrow=nrow(gaussian.draw), ncol=ncol(gaussian.draw), 
                     dimnames=list(genes, spots))
for (i in 1:nrow(pit.draw)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  pit.draw[i,] = sapply(1:ncol(gaussian.draw),
                           function(i) pnorm(gaussian.draw[gene,i]))
}
# PIT to simulated expressions
expr.mat = matrix(nrow=length(genes), ncol=length(spots))
rownames(expr.mat) = genes
colnames(expr.mat) = spots
for (i in 1:length(genes)){
  if ((i%%10)==0){
    print(i)
  }
  gene = genes[i]
  disp.rate = get.rates.gp(gene, coef.mat, decon.res, gbar.mat, spot.coef, n.cells)
  expr.mat[i,] = sapply(1:ncol(pit.draw),
                           function(i) qNBI(pit.draw[gene,i], mu=disp.rate[i,2]+0.000001, sigma=disp.rate[i,1]+0.000001))
  
}

output = list("seed"=seed, "expr.mat"=expr.mat, "coords" = coords)
saveRDS(output, "data/mkc_sim_gp_gt0.rds")

# check estimation
# read input
load("data/mkc_prep_est/mkc_prep.rdata")
#input = readRDS("data/mkc_prep.rds")
decon.res = x[["decon"]]
n.cells = x[["n.cells"]]
counts.spe = x[["counts.spe"]]
spe = x[["spe"]]
coords = x[["coords"]]
gbar.mat = x[["gbar.mat"]]

var.sim.genes = apply(expr.mat, 1, var)
var.true.genes = apply(counts.spe[rownames(expr.mat),], 1, var)
boxplot(sqrt(var.sim.genes)[-which.max(var.sim.genes)], sqrt(var.true.genes)[-which.max(var.sim.genes)], main="standard deviations", names=c("sim", "true"))
plot(sqrt(var.sim.genes)[-which.max(var.sim.genes)], 
     sqrt(var.true.genes)[-which.max(var.sim.genes)],
     main="gene sds",
     ylab="true sds", xlab="sim sds")

which.max(var.sim.genes)
# Lyz1 has very high standard deviation in simulated data
compare.cell.plot(counts.spe["Lyz1",], expr.mat["Lyz1",], spe, titles=c("real expr", "sim expr"), main="Lyz1")
boxplot(counts.spe["Lyz1",], expr.mat["Lyz1",])
max(counts.spe["Lyz1",]); max(expr.mat["Lyz1",])
sp.effect.lyz1 = coef.mat["Lyz1", 4:ncol(coef.mat)]
cell.plot(exp(sp.effect.lyz1), spe)
offset.lyz1 = n.cells*decon.res%*%gbar.mat[,"Lyz1"]*spot.coef
sd.fit.lyz1 = coef.mat["Lyz1", 1] + coef.mat["Lyz1", 2]*log(offset.lyz1+0.01)
cell.plot(exp(sd.fit.lyz1), spe)

mean.sim.genes = apply(expr.mat, 1, mean)
mean.true.genes = apply(counts.spe[rownames(expr.mat),], 1, mean)
boxplot(mean.sim.genes, mean.true.genes, main="means", names=c("sim", "true"))
which.max(mean.sim.genes) # Lyz1
plot(sqrt(mean.sim.genes)[-which.max(mean.sim.genes)], 
     sqrt(mean.true.genes)[-which.max(mean.true.genes)],
     main="gene means",
     ylab="true means", xlab="sim means")

library(ape)
library(DR.SC)
library(Seurat)
library('stringr')

colnames(counts.spe) = sub("\\.", "-", colnames(counts.spe))
colnames(coords) = c("col", "row")
s.objc = CreateSeuratObject(counts = counts.spe, meta.data=data.frame(round(coords/100)))
adj = getAdj(s.objc, platform="Visium")
adj = as.matrix(adj)

#dists = as.matrix(dist(coords))
#dists.inv = 1/dists
#diag(dists.inv) = 0
I.vec.sim = vector()
skipped = vector()
for (i in 1:length(genes)){
  gene = genes[i]
  if(sum(expr.mat[gene,])==0){
    skipped = c(skipped, gene)
    next
  }
  if ((i %% 10) == 0){
    print(i)
  }
  I.vec.sim = c(I.vec.sim, Moran.I(expr.mat[gene,], adj)$observed)
}
names(I.vec.sim) = setdiff(genes, skipped)
I.vec.true = vector()
for (i in 1:length(genes)){
  gene = genes[i]
  if (gene %in% skipped){
    next
  }
  if ((i %% 10) == 0){
    print(i)
  }
  I.vec.true = c(I.vec.true, Moran.I(counts.spe[gene,], adj)$observed)
}
plot(I.vec.sim, I.vec.true, xlim=c(0,1), ylim=c(0,1))
abline(a=0, b=1)








