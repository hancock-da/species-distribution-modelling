library(dplyr)
library(CoordinateCleaner)
library(dismo)
library(maptools)
library(raster)
library(biomod2)
library(remotes)
library(ENMTools)


gbif_download = readr::read_tsv("datasets/bank_vole_occurence.csv")

### LOADING AND CLEANING DATA ###
# Filter coordinates
gbif_download <- gbif_download %>%
  setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
  filter(occurrencestatus  == "PRESENT") %>% # only keep presence data
  filter(!is.na(decimallongitude)) %>% # remove null coordinate occurrences, if any
  filter(!is.na(decimallatitude)) %>% 
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")) %>%  
  filter(!establishmentmeans %in% c("MANAGED", "INTRODUCED", "INVASIVE", "NATURALISED")) %>%
  filter(year >= 1900) %>% # remove old records
  filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% # remove values with high uncertainty, missing values kept
  filter(coordinateuncertaintyinmeters < 10000 | is.na(coordinateuncertaintyinmeters)) %>% # remove values with high uncertainty, missing values kept
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% # known uncertainty values
  filter(!decimallatitude == 0 | !decimallongitude == 0) %>% # remove points along prime meridian or equator
  cc_cen(buffer = 2000) %>% # remove country centroids within 2km 
  cc_cap(buffer = 2000) %>% # remove capitals centroids within 2km
  cc_inst(buffer = 2000) %>% # remove zoo and herbaria within 2km 
  cc_sea() %>% # remove from ocean, as not marine species
  distinct(decimallongitude,decimallatitude,specieskey,datasetkey, .keep_all = TRUE) %>% # remove duplicates
  rename(lon = decimallongitude, lat = decimallatitude) %>%
  glimpse() # look at results of pipeline

# quick visual check of coordinates
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, xlim=c(-20,40), ylim=c(10,70), axes=TRUE, col="light yellow")
# restore the box around the map
box()
# add the points
points(gbif_download$lon, gbif_download$lat, col='orange', pch=20, cex=0.75)

### SAMPLING BIAS ###
# subsample occurrence records necessary as sampling highly biased to Western Europe 
acg <- gbif_download
coordinates(acg) <- ~lon+lat
crs(acg) <- crs(wrld_simpl)
# create a raster layer with the extent of acg
r <- raster(acg)
# set the resolution of the cells to (e.g.) 1 degree
res(r) <- 1
#expand the extent of the rasterlayer a little bit
r <- extend(r, extent(r)+10)
#sample for training set
acsel <- gridSample(acg, r, n=2, chess='white')
# use chessboard sampling to sample the same number of points for test set
actest <- gridSample(acg, r, n=2, chess='black')
# to illustrate the method and show the result
p <- rasterToPolygons(r)
plot(p, border='gray')
# plot all original points
points(acg)
# selected training points in red and test in blue
points(acsel, cex=1, col='red',pch='x')
points(actest, cex=1, col='blue',pch='x')

### ENVIRONMENTAL DATA ###
path <- "datasets/bioclim"
files <- list.files(path, pattern='tif$',full.names=TRUE)

# create a rasterStack of predictor variables and crop to same extent as species occurences
predictors <-  stack(files)
#cropped_predictors <- crop(predictors[[1]], extent(r))


# function to crop layers of a stack to a specified extent and keep the output as a RasterStack
# (as opposed to a RasterBrick if using crop() function) for use in biomod2
cropRasterStack <- function(stack, e) {
    cropped_predictors <- stack()
    for (i in 1:length(names(stack))) {
      cropped_layer <- crop(stack[[i]], extent(e))
      cropped_predictors <- addLayer(cropped_predictors, cropped_layer)
    }
    return(cropped_predictors)
}

cropped_predictors <- cropRasterStack(predictors, r)
plot(cropped_predictors)

# drop explanatory variables to remove collinearity and prevent overfitting
cor_mat <- raster.cor.matrix(cropped_predictors, method = 'pearson')
cor_mat > 0.7

cor_plot <- raster.cor.plot(cropped_predictors, method='pearson')
cor_plot

# keep most biologically easy to explain layers
test_predictors <- dropLayer(predictors, c('wc2.1_2.5m_bio_10','wc2.1_2.5m_bio_11','wc2.1_2.5m_bio_3','wc2.1_2.5m_bio_5','wc2.1_2.5m_bio_9', 'wc2.1_2.5m_bio_6','wc2.1_2.5m_bio_4',
                                           'wc2.1_2.5m_bio_13','wc2.1_2.5m_bio_14','wc2.1_2.5m_bio_16','wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_19'))
test_predictors_cropped <- cropRasterStack(test_predictors, r)
cor_mat <- raster.cor.matrix(test_predictors_cropped, method = 'pearson')
cor_mat > 0.7 # no more correlated explanatory variables


### MODELLING ###
# produce presence vector for training and test sets
myResp <- rep(1, nrow(acsel))
myEval <- rep(1, nrow(actest))
SpecName <- 'BankVole'


# for pseudo absences, using surface range envelop (sre) - pseudo absence points
# are selected in conditions that differ from a defined proportion (PA.sre.quant)
# of the presence data. Forces PAs to be selected outside of the broadly defined 
# environmental conditions of the species


# first use BIOMOD to generate pseudoabsences for training data
myBiomodPA <- BIOMOD_FormatingData(resp.var = myResp,
                                   expl.var = test_predictors_cropped,
                                   resp.xy = acsel,
                                   PA.nb.rep = 1,
                                   PA.nb.absences = length(myResp),
                                   PA.strategy = 'sre',
                                   PA.sre.quant = 0.1,
                                   resp.name = SpecName)

# grab pseudoabsence data from formatted dataset
train_sp <- SpatialPointsDataFrame(coords = myBiomodPA@coord, 
                                 data = data.frame(species = myBiomodPA@data.species), 
                                 proj = CRS(proj4string(test_predictors_cropped)))
# replace NAs with 0 for pseudobsences
train_sp$species[is.na(train_sp$species)] <- 0


# now use BIOMOD to generate pseudoabsences for test data
myBiomodData <- BIOMOD_FormatingData(resp.var = myEval,
                                     expl.var = test_predictors_cropped,
                                     resp.xy = actest,
                                     resp.name = SpecName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = nrow(actest),
                                     PA.strategy = 'sre',
                                     PA.sre.quant = 0.1)
# grab pseudoabsence data from formatted dataset
test_sp <- SpatialPointsDataFrame(coords = myBiomodPA@coord, 
                                   data = data.frame(species = myBiomodPA@data.species), 
                                   proj = CRS(proj4string(test_predictors_cropped)))
# replace NAs with 0 for pseudobsences
test_sp$species[is.na(test_sp$species)] <- 0


# check training datapoints
plot(test_predictors_cropped[[1]])
points(train_sp[train_sp$species == 1,], pch=19)
points(train_sp[train_sp$species == 0,], pch=24, col='red')

# check test datapoints
plot(test_predictors_cropped[[1]])
points(test_sp[test_sp$species == 1,], pch=19)
points(test_sp[test_sp$species == 0,], pch=24, col='red')

# input training and test datasets to BIOMOD formatted data
myBiomodData <- BIOMOD_FormatingData(resp.var = train_resp,
                                     expl.var = test_predictors_cropped,
                                     resp.xy = train_coords,
                                     resp.name = SpecName,
                                     eval.resp.var = test_resp,
                                     eval.expl.var = test_predictors_cropped,
                                     eval.resp.xy = test_coords)


# defining models options using default options
myBiomodOption <- BIOMOD_ModelingOptions(MAXENT = list(path_to_maxent.jar = "maxent"))

myBiomodModel0ut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM', 'RF', 'FDA', 'ANN', 'MARS', 'MAXENT.Phillips'),
  models.options = myBiomodOption,
  NbRunEval = 3,
  DataSplit = 80, # split 80% of data for training in cross-validation
  Prevalence=0.5,
  VarImport = 3, 
  models.eval.meth = c('ROC','TSS','ACCURACY'),
  SaveObj = TRUE,
  recal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(SpecName,"FirstModeling",sep="")
)

# get all models evaluation
myBiomodModelEval <- get_evaluations(myBiomodModel0ut)
myBiomodModelEval

# graph evalation scores by model
models_scores_graph(
  myBiomodModel0ut,
  by = "models",
  metrics=c("ROC","TSS"),
  xlim=c(0,1),
  ylim=c(0,1)
)

myBiomodModelEval["TSS", "Testing.data",,,]
myBiomodModelEval["TSS", "Evaluating.data",,,]
myBiomodModelEval["ROC", "Testing.data",,,]
myBiomodModelEval["ROC", "Evaluating.data",,,]

# retrieve variable importance
MyBiomodModelVarImp <- get_variables_importance(myBiomodModel0ut)
MyBiomodModelVarImp # bio_7 has high importance in best performing models

# get model names
get_built_models(myBiomodModel0ut)

# Ensemble Modeling, choose models with highest TSS
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModel0ut,
  chosen.models = 'all',
  em.by='PA_dataset+repet',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.8),
  prob.mean = T,
  prob.cv = F,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = F,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

myBiomodEMEval <- get_evaluations(myBiomodEM)
myBiomodEMEval

## Ensemble mode built with weighted means performs best.
# specifically random forests worked very well on test data 
# and the three RF runs formed the basis of the ensemble.
# RF on its own seems to outperform by ROC score.

# Projections
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModel0ut,
  new.env = test_predictors_cropped,
  proj.name = 'current',
  selected.models = grep('_RF', get_built_models(
    myBiomodModel0ut), value=TRUE),
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# get projected map
myCurrentProj <- get_predictions(myBiomodProj)

# ensemble forecasting
myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProj,
  selected.models = grep('_EMwmeanByTSS',get_built_models(
    myBiomodEM), value=TRUE)
)
# plot the 3 weighted mean ensemble models
plot(myBiomodEF) 
