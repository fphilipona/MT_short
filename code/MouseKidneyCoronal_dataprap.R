source("code/utils.R")

library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(TENxVisiumData)
library(TabulaMurisSenisData)
library(mgcv)
library(ggpubr)
library(spatstat.geom)

# 1 sampling step; set seed
seed = 25
set.seed(seed)

###############################################################################################################
# FROM VIGNETTE https://bioconductor.org/packages/devel/bioc/vignettes/SPOTlight/inst/doc/SPOTlight_kidney.html
###############################################################################################################

# data
spe <- MouseKidneyCoronal()  #spatial experiment data
# Use symbols instead of Ensembl IDs as feature names
rownames(spe) <- rowData(spe)$symbol
sce <- TabulaMurisSenisDroplet(tissues = "Kidney")$Kidney # matching single scell experiment data
# Keep cells from 18m mice
sce <- sce[, sce$age == "18m"]
# Keep cells with clear cell type annotations
sce <- sce[, !sce$free_annotation %in% c("nan", "CD45")]

# first: get marker genes for each cell type
# library size normalization
sce <- logNormCounts(sce)

# get highly variable genes
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n = 3000)

# marker genes
colLabels(sce) <- colData(sce)$free_annotation
# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))
# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

# only keep relevant genes
mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

# downsample the cells
# split cell indices by identity
idx <- split(seq(ncol(sce)), sce$free_annotation)
# downsample to at most 50 per identity & subset
n_cells <- 50
cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]

res <- SPOTlight(
  x = sce,
  y = spe,
  groups = as.character(sce$free_annotation),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

###############################################################################################################
# END OF VIGNETTE
###############################################################################################################

# counts sce and counts spe
hvg = hvg[hvg %in% rownames(counts(spe))]
counts.sce = counts(sce)[hvg,]
counts.spe = counts(spe)[hvg,]
counts.sce = matrix(counts.sce, nrow = nrow(counts.sce), ncol = ncol(counts.sce), dimnames=list(rownames(counts.sce),
                                                                                                colnames(counts.sce)))
counts.spe = matrix(counts.spe, nrow = nrow(counts.spe), ncol = ncol(counts.spe), dimnames=list(rownames(counts.spe),
                                                                                               colnames(counts.spe)))

# Coordinates
coords = spatialCoords(spe)

# Gbar mat
idx.cidentity <- split(seq(ncol(sce)), sce$free_annotation) # 1, 19 are not relevant
gbar = lapply(idx.cidentity[-c(1,19)], function(x){ # 1, 19 are not relevant
  apply(counts.sce[,x], 1, mean)
}) # quite a few entries are zero -> later on log will have to be taken
gbar.mat = matrix(nrow = length(gbar), ncol=length(gbar[[1]]))
for (i in 1:length(gbar)){
  for (j in 1:length(gbar[[i]])){
    gbar.mat[i,j] = gbar[[i]][j]
  }
}
rownames(gbar.mat) = names(gbar)
colnames(gbar.mat) = names(gbar[[1]])

# number of cells using krigin
#### number of cells 
spe$x <- coords[, 1]
spe$y <- coords[, 2]
# 
libsize.spot = apply(counts.spe, 2, sum)
obs.libsize = data.frame("libsize" = libsize.spot, "coord1" = coords[, 1], "coord2" = coords[, 2])
reg.libsize = gam(libsize ~ s(coord1, coord2, bs="gp", k=100), family="nb", data=obs.libsize)
pred = reg.libsize$fitted.values
# divide smoothed libsize (means) by cell type proportions times mean contributions
theoretical.means = res$mat %*% gbar.mat # mean, if only cell types at spots would count
theoretical.libsizes = apply(theoretical.means, 1, sum)
hist(pred/theoretical.libsizes) # ~ numbers of cells not accounting for measurement error
n.cells.withnoise = pred/theoretical.libsizes

true.mean.ncells = 10
factor = 10/mean(n.cells.withnoise)
n.cells.notrounded = n.cells.withnoise*factor
n.cells = round(n.cells.notrounded)

n.cells[which(n.cells==0)]=1

# the two outliers
outlier.names = names(c(which(coords[,2]<2500),which(coords[,2]>8500)))
counts.spe = counts.spe[,!(colnames(counts.spe) %in% outlier.names)]
coords = coords[!(rownames(coords)%in%outlier.names),]
n.cells = n.cells[!(names(n.cells)%in%outlier.names)]
res.mat = res$mat[!(rownames(res$mat)%in%outlier.names),]
newspe = SpatialExperiment(assay = list(counts=counts.spe[,!(colnames(counts.spe) %in% outlier.names)]),
                           spatialCoords = coords[!(rownames(coords)%in%outlier.names),],
                           colData = colData(spe)[!(rownames(colData(spe))%in%outlier.names),])
newspe$x <- spatialCoords(newspe)[,1]
newspe$y <- spatialCoords(newspe)[,2]

# regions
# choose some genes that show "nice/desired" regional differences in their expressions
# Mettl7a2, Acadm, Arg2, Aqp1, 4930461G14Rik
regions = ifelse(counts.spe["Aqp1",]>60&counts.spe["Acadm",]>40&counts.spe["Mettl7a2",]>15, 1, 0)
names(regions) = colnames(counts.spe)

nn5.mat = make_nnmat(coords, k=5) # indices of 5 nearest neighbors for each spot
abc = sum_regions_neighbors(nn5.mat, regions) # sum up regions of neighbors
regions2 = regions
regions2[which(regions2!=1&abc>2)]=1

regions3 = regions2
regions3[which(regions2!=1&coords[,1]>5100&coords[,1]<8000&coords[,2]>4750&coords[,2]<6500)]=2
regions3[which(regions3==2&coords[,1]>5100&coords[,1]<5300&coords[,2]>4750&coords[,2]<5100)]=0
regions3[which(regions3==0&coords[,1]>6400&coords[,1]<7650&coords[,2]>4600&coords[,2]<4900)]=2
regions3[which(regions3==0&coords[,1]>8000&coords[,1]<8200&coords[,2]>5200&coords[,2]<5900)]=2
regions3[which(regions3==2&coords[,1]>7900&coords[,1]<8100&coords[,2]>5900&coords[,2]<6200)]=0
regions3[which(regions3==2&coords[,1]>7900&coords[,1]<8100&coords[,2]>5900&coords[,2]<6400)]=0
regions3[which(regions3==0&coords[,1]>6000&coords[,1]<7400&coords[,2]>6300&coords[,2]<6600)]=2
cell.plot(regions3,
          newspe)


# create output
output= list("seed"=seed, "counts.spe"=counts.spe, "counts.sce"=counts.sce,
             "coords"=coords, "decon"=res.mat, "gbar.mat"=gbar.mat, "n.cells"=n.cells,
             "spe"=newspe, "regions"=regions3)

saveRDS(output, "data/mkc_prep.rds")









