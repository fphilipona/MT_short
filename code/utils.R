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
library(patchwork)


# function to plot stati for the spatial experiment
cell.plot = function(stati, spe, main=" ", limits_scale = 0){
  e1 = stati
  if (sum(limits_scale)==0){
    limits_scale = c(min(stati), max(stati))
    ggcells(spe, aes(x, y, color = e1)) +
      geom_point() +
      scale_color_viridis_c(limits = limits_scale) +
      coord_fixed() +
      theme_bw() +
      labs(title=main)+
      theme(legend.title=element_blank())+
      theme(plot.title = element_text(size=10))
  }else{
    ggcells(spe, aes(x, y, color = e1)) +
      geom_point() +
      scale_color_viridis_c(limits = limits_scale) +
      coord_fixed() +
      theme_bw() +
      labs(title=main)+
      theme(legend.title=element_blank())+
      theme(plot.title = element_text(size=10))
  }
}

compare.cell.plot = function(expr, f.mean, spe1, titles=c("true expr","fitted mean"), main=" "){
  min.range = min(c(expr, f.mean))
  max.range = max(c(expr, f.mean))
  e1 = expr
  true = ggcells(spe1, aes(x, y, color = e1)) +
    geom_point() +
    scale_color_viridis_c(limits = range(min.range, max.range)) +
    coord_fixed() +
    theme_bw() +
    labs(title=paste(titles[1]))+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(size=10))
  e2 = f.mean
  predicted = ggcells(spe1, aes(x, y, color = e2)) +
    geom_point() +
    scale_color_viridis_c(limits = range(min.range, max.range)) +
    coord_fixed() +
    theme_bw() +
    labs(title=paste(titles[2]))+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(size=10))
  figure <- (true + predicted) + plot_annotation(title=main)&theme(plot.title = element_text(hjust = 0.5))
  figure
}

####################### for data prep ##################################
# get indices of the 5 nearest neighbors for each coordinate pair
make_nnmat = function(coords, k){
  for (i in 1:k){
    if (i == 1){
      nn.mat = nnwhich(coords)
    }
    else{
      nn.mat = cbind(nn.mat, nnwhich(coords, k=i))
    }
  }
  rownames(nn.mat) = rownames(coords)
  return(nn.mat)
}


sum_regions_neighbors = function(nn.mat, regions.vec){
  apply(nn.mat, 1, function(x) sum(regions.vec[x]))
}
#########################################################


######################################## estimation ###################################################

# function to make the dataframe used in the estimation with gamlss
makedf_regions = function(gene, regions, celltypes.spots, gbar.mat, n.cells, counts.spe, spots.coef){
  data.frame("expr"=counts.spe[gene,], "region" = as.factor(regions), 
             "off"= n.cells*celltypes.spots%*%gbar.mat[,gene]*spots.coef)
}

makedf_gp = function(gene, coords, celltypes.spots, gbar.mat, n.cells, counts.spe, spots.coef){
  data.frame("expr"=counts.spe[gene,], "coord1" = coords[,1], "coord2" = coords[,2], 
             "off"= n.cells*celltypes.spots%*%gbar.mat[,gene]*spots.coef)
}







######################################## check_estimation ######################################## 
helper.mat = function(regions){
  for(i in 1:length(regions)){
    if(regions[i]==0){
      row.iter = c(1,0,0)
    }
    if(regions[i]==1){
      row.iter = c(1,1,0)
    }
    if(regions[i]==2){
      row.iter = c(1,0,1)
    }
    if (i == 1){
      help.mat = row.iter
    }
    else{
      help.mat = rbind(help.mat, row.iter)
    }
  }
  return(help.mat)
}

helper.mat2 = function(regions){
  for(i in 1:length(regions)){
    if(regions[i]==1){
      row.iter = c(1,0,0,0,0,0,0)
    }
    if(regions[i]==2){
      row.iter = c(0,1,0,0,0,0,0)
    }
    if(regions[i]==3){
      row.iter = c(0,0,1,0,0,0,0)
    }
    if(regions[i]==4){
      row.iter = c(0,0,0,1,0,0,0)
    }
    if(regions[i]==5){
      row.iter = c(0,0,0,0,1,0,0)
    }
    if(regions[i]==6){
      row.iter = c(0,0,0,0,0,1,0)
    }
    if(regions[i]==7){
      row.iter = c(0,0,0,0,0,0,1)
    }
    if (i == 1){
      help.mat = row.iter
    }
    else{
      help.mat = rbind(help.mat, row.iter)
    }
  }
  return(cbind(rep(1, length(regions)), help.mat))
}

helper.mat1 = function(regions){
  for(i in 1:length(regions)){
    if(regions[i]==0){
      row.iter = c(1,0,0)
    }
    if(regions[i]==1){
      row.iter = c(0,1,0)
    }
    if(regions[i]==2){
      row.iter = c(0,0,1)
    }
    if (i == 1){
      help.mat = row.iter
    }
    else{
      help.mat = rbind(help.mat, row.iter)
    }
  }
  return(cbind(rep(1, length(regions)), help.mat))
}

extend.coef.mat = function(coef.mat){
  # function to get the regional effects to a mean of 0 (their exp = 1) and therefore a "true" gde
  genes = rownames(coef.mat)
  ncols = ncol(coef.mat)
  coef.mat.keep = coef.mat[,c(1,2)]
  regional.effects = cbind(rep(0, nrow(coef.mat)), coef.mat[,4:ncols])
  mean.regional.effects = apply(regional.effects, 1, function(x) mean(exp(x)))
  new.gde = log(exp(coef.mat[,3])*mean.regional.effects)
  new.regional.effects = sweep(regional.effects, 1, STATS=log(mean.regional.effects), FUN="-")
  new.coef.mat = cbind(coef.mat.keep, new.gde, new.regional.effects)
  rownames(new.coef.mat) = genes
  colnames(new.coef.mat) = c(colnames(coef.mat)[1:3], "region1", colnames(coef.mat)[4:ncols])
  return(new.coef.mat)
}



# estimated rates and dispersions for a gene from the coefficient matrix
get.rates.regions = function(gene, coef.mat, celltypes.spots, gbar.mat, spots.coef, n.cells, help.mat,
                             subset.spots){
  ncols = ncol(coef.mat)
  offset = n.cells*celltypes.spots%*%gbar.mat[,gene]*spots.coef
  sigma = coef.mat[gene, 1] + coef.mat[gene, 2]*log(offset+0.01)
  intercept.and.spatial.effect = help.mat%*%coef.mat[gene, 3:ncols]
  return(cbind(exp(sigma), (offset+0.01)*exp(intercept.and.spatial.effect)))
}

get.rates.gp = function(gene, coef.mat, celltypes.spots, gbar.mat, spots.coef, n.cells){
  offset = n.cells*celltypes.spots%*%gbar.mat[,gene]*spots.coef
  sigma = coef.mat[gene, 1] + coef.mat[gene, 2]*log(offset+0.01)
  intercept.and.spatial.effect = (offset+0.01)*exp(coef.mat[gene, 4:dim(coef.mat)[2]]+coef.mat[gene, 3])
  return(cbind(exp(sigma), intercept.and.spatial.effect))
}
################################################################################# 





