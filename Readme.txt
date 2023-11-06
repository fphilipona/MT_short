This R project contains excerpts of code used in my master thesis. In the thesis, I proposed a new way to simulate spatial transcriptomic data.
Briefly, I proposed a parametric model to simulate from (for details see the master thesis in document MT_fphilipona.pdf) which is a combination
of marginal gene expression distributions combined with a gaussian copula to joint distributions per spot. The model first has to be estimated,
before simulation can be done. The code contained in this project estimates the model based on one real world dataset, checks the estimated parameters,
simulates from the model and checks the simulated data.


Files in code folder:

1) MouseKidneyCoronal_dataprep.R
input: mouse kidney coronal dataset through download (R library: TENxVisiumData) and a single cell experiment through download (R library: TabulaMurisSenisData)
Some datapreparation is done on the mouse kidney coronal dataset, a SPOTLIGHT deconvolution is applied to get cell types per spots, the number of cells 
per spot is estimated, a mean contribution matrix for each cell type to the expression of each gene is constructed and the spatial dataset is divided into
3 regions.

2a) estimation_mkc.R
input: estimated parameters from MouseKidneyCoronal_dataprep.R
Spot detection efficiencies are estimated using random forests, regional effects and gene detection efficiencies are estimated via gamlss. A Gaussian copula
is estimated.

2b) estimation_mkc_gp.R
input: estimated parameters from MouseKidneyCoronal_dataprep.R
Alternative way to estimate the model. Here, there are no regional effects, instead a gaussian process smoother is used to estimate the geographic component of the model.

3a) sim_mkc.R
input: all estimated parameters (from estimation_mkc.R)
Data are simulated for the different regions separately by first drawing from the Gaussian copula and then transforming the result with the marginal gene expression
distributions. 

3b) sim_mkc_gp.R
input: all estimated parameters (from estimation_mkc_gp.R)
Simulate data from the alternative model using a gaussian process smoother for the geographic component instead of regional effects.

4) utils
contains utility functions


Data folder:
Contains datafiles created with some of the R scripts above. mkc_estimation.rdata not included because too large


Markdown files folder:

1) evaluation_mkc.pdf
file showing some evaluation statistics of the estimated model and the simulated data for the mouse kidney coronal data set. created with evaluation_mkc.rnw

2) MT_FPhilipona
the entire master thesis.

