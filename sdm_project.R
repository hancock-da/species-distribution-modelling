library(dplyr)
library(CoordinateCleaner)
library(dismo)
library(maptools)
library(raster)
library(biomod2)
library(remotes)
library(ENMTools)

setwd("./species-distribution-modelling")

SPECIES <- 'bank_vole'

gbif_download = readr::read_tsv("bank_vole_occurence.csv")

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
xlims <- c(min(gbif_download$lon)-10, max(gbif_download$lon)+10)
ylims <- c(min(gbif_download$lat)-10, max(gbif_download$lat)+10)
plot(wrld_simpl, xlim=xlims, ylim=ylims, axes=TRUE, col="light yellow")
# restore the box around the map
box()
# add the points
points(gbif_download$lon, gbif_download$lat, col='orange', pch=20, cex=0.75)

### SAMPLING BIAS ###
nrow(gbif_download)
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
# selected points in red
points(acsel, cex=1, col='red',pch='x')
points(actest, cex=1, col='blue',pch='x')


### ENVIRONMENTAL DATA ###
bioclim_path <- "map_files/bioclim"
bioclim_files <- list.files(bioclim_path, pattern='tif$',full.names=TRUE)
# create a rasterStack of predictor variables
bioclim_predictors <-  stack(bioclim_files)

# function to project and aggregate rasters to same as bioclim rasters
resizeRasterStack <- function(source_stack, target_stack, output_path, agg_factor=4) {
  print(names(source_stack))
  # if raster stacks share the same CRS, only need to aggregate.
  if (identical(crs(source_stack), crs(target_stack))) {
    resolution <- res(target_stack)
    fact <- resolution/res(source_stack)
    agg_raster <- aggregate(source_stack, fact, fun='mean')
  } else {
    gc()
    agg_raster <- aggregate(source_stack, fact=agg_factor) # this line is to make the projection manageable with small amounts of RAM but may need changing
    agg_raster <- projectRaster(agg_raster, res=res(target_stack[[1]]), crs=crs(target_stack[[1]])) # project to match resolution and CRS.
  }
  for (r in 1:nlayers(agg_raster)) {
    writeRaster(agg_raster[[r]], paste(output_path,names(source_stack)[[r]],'.tif',sep=""))
  }
  return(NULL)
}

# Process Human Footprint Maps to same CRS and resolution
# Uncomment if processing again.
# raw_hf_path <- "map_files/human_footprint/raw"
# raw_hf_files <- list.files(raw_hf_path, pattern='tif$', full.names=TRUE)
# raw_hf_predictors <- stack(raw_hf_files)
# 
# resizeRasterStack(raw_hf_predictors, bioclim_predictors, "map_files/human_footprint/aggregated")

# Read in aggregated human footprint maps
hf_path <- "map_files/human_footprint/aggregated"
hf_files <- list.files(hf_path, pattern='tif$',full.names=TRUE)
hf_predictors <- stack(hf_files)

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

# set extent to slightly larger than study area for cropping rasters
e <- extent(xlims, ylims)

# crop bioclim and human footprint maps to same extent
cropped_predictors <- cropRasterStack(bioclim_predictors, e)
cropped_hf <- resample(hf_predictors, cropped_predictors)

# stack all predictor rasters together
cropped_predictors <- addLayer(cropped_predictors, cropped_hf)

# drop explanatory variables to remove collinearity and prevent overfitting
sampled_preds <- extract(cropped_predictors, coordinates(acg), 'bilinear')
cor_mat <- raster.cor.matrix(cropped_predictors, method = 'pearson')
cor_mat > 0.7
cor_mat < -0.7

# cor_plot <- raster.cor.plot(cropped_predictors, method='pearson')
# cor_plot

# keep most biologically easy to explain layers. Also dropping elev which the other bioclim maps are interpolated from.
test_predictors_cropped <- dropLayer(cropped_predictors, c('wc2.1_2.5m_bio_10','wc2.1_2.5m_bio_11','wc2.1_2.5m_bio_2','wc2.1_2.5m_bio_3','wc2.1_2.5m_bio_5',
                                                           'wc2.1_2.5m_bio_6','wc2.1_2.5m_bio_7','wc2.1_2.5m_bio_9',
                                                           'wc2.1_2.5m_bio_13','wc2.1_2.5m_bio_14','wc2.1_2.5m_bio_16',
                                                           'wc2.1_2.5m_bio_17','wc2.1_2.5m_bio_18','wc2.1_2.5m_bio_19',
                                                           'Lights2009','wc2.1_2.5m_elev'))

# drop further features that add little importance to best performing models (GBM and RF) to minimise overfitting
test_predictors_cropped <- dropLayer(test_predictors_cropped, c('Navwater2009','Pasture2009','croplands2005','wc2.1_2.5m_bio_8'))

# check for no more correlated explanatory variables
cor_mat <- raster.cor.matrix(test_predictors_cropped, method = 'pearson')
cor_mat > 0.7 
cor_mat < -0.7


### MODELLING ###
# produce presence vector for training and test sets
myResp <- rep(1, nrow(acsel))
myEval <- rep(1, nrow(actest))

# for pseudo absences, using surface range envelop (sre) - pseudo absence points
# are selected in conditions that differ from a defined proportion (PA.sre.quant)
# of the presence data. Forces PAs to be selected outside of the broadly defined 
# environmental conditions of the species


# first use BIOMOD to generate pseudoabsences for training data
myBiomodPA_train <- BIOMOD_FormatingData(resp.var = myResp,
                                         expl.var = test_predictors_cropped,
                                         resp.xy = acsel,
                                         PA.nb.rep = 1,
                                         PA.nb.absences = length(myResp),
                                         PA.strategy = 'sre',
                                         PA.sre.quant = 0.1,
                                         resp.name = SPECIES)

# grab pseudoabsence data from formatted dataset
train_sp <- SpatialPointsDataFrame(coords = myBiomodPA_train@coord, 
                                   data = data.frame(species = myBiomodPA_train@data.species), 
                                   proj = CRS(proj4string(test_predictors_cropped)))
# replace NAs with 0 for pseudobsences
train_sp$species[is.na(train_sp$species)] <- 0


# now use BIOMOD to generate pseudoabsences for test data
myBiomodPA_test <- BIOMOD_FormatingData(resp.var = myEval,
                                        expl.var = test_predictors_cropped,
                                        resp.xy = actest,
                                        resp.name = SPECIES,
                                        PA.nb.rep = 1,
                                        PA.nb.absences = nrow(actest),
                                        PA.strategy = 'sre',
                                        PA.sre.quant = 0.1)
# grab pseudoabsence data from formatted dataset
test_sp <- SpatialPointsDataFrame(coords = myBiomodPA_test@coord, 
                                  data = data.frame(species = myBiomodPA_test@data.species), 
                                  proj = CRS(proj4string(test_predictors_cropped)))
# replace NAs with 0 for pseudobsences
test_sp$species[is.na(test_sp$species)] <- 0

train_sp
test_sp
# check training datapoints
plot(test_predictors_cropped[[1]])
points(train_sp[train_sp$species == 1,], pch=19)
points(train_sp[train_sp$species == 0,], pch=24, col='red')

# check test datapoints
plot(test_predictors_cropped[[1]])
points(test_sp[test_sp$species == 1,], pch=19)
points(test_sp[test_sp$species == 0,], pch=24, col='red')

# input training and test datasets to BIOMOD formatted data
myBiomodData <- BIOMOD_FormatingData(resp.var = train_sp$species,
                                     expl.var = test_predictors_cropped,
                                     resp.xy = train_sp@coords,
                                     resp.name = SPECIES,
                                     eval.resp.var = test_sp$species,
                                     eval.expl.var = test_predictors_cropped,
                                     eval.resp.xy = test_sp@coords)


# defining models options using default options
myBiomodOption <- BIOMOD_ModelingOptions(MAXENT = list(path_to_maxent.jar = "maxent"))

# build and train the individual models
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM', 'GAM', 'CTA', 'SRE', 'RF', 'FDA', 'ANN', 'MARS', 'MAXENT.Phillips'),
  models.options = myBiomodOption,
  NbRunEval = 3,
  DataSplit = 80, # split 80% of data for training in cross-validation
  Prevalence=0.5,
  VarImport = 3, 
  models.eval.meth = c('ROC','TSS','ACCURACY'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(SPECIES,"FirstModeling",sep="")
)

myBiomodModelOut
# get all models evaluation
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myBiomodModelEval

# graph evaluation scores by model
models_scores_graph(
  myBiomodModelOut,
  by = "models",
  metrics=c("ROC","ACCURACY"),
  xlim=c(0,1),
  ylim=c(0,1)
)

myBiomodModelEval["TSS", "Testing.data",,,]
myBiomodModelEval["TSS", "Evaluating.data",,,]
rowMeans(myBiomodModelEval["TSS", "Evaluating.data",,,])
myBiomodModelEval["ROC", "Testing.data",,,]
myBiomodModelEval["ROC", "Evaluating.data",,,]
rowMeans(myBiomodModelEval["ROC", "Evaluating.data",,,])
myBiomodModelEval["ACCURACY", "Testing.data",,,]
myBiomodModelEval["ACCURACY", "Evaluating.data",,,]


# retrieve variable importance
MyBiomodModelVarImp <- get_variables_importance(myBiomodModelOut)
MyBiomodModelVarImp
rowMeans(MyBiomodModelVarImp)


# get model names
get_built_models(myBiomodModelOut)
plot(myBiomodelOut)

# Ensemble Modeling, choose best 5 models with highest ROC
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.845),
  prob.mean = T,
  prob.cv = F,
  prob.ci = F,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = F,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

myBiomodEMEval <- get_evaluations(myBiomodEM, as.data.frame=TRUE)
myBiomodEMEval
myBiomodEMEval[myBiomodEMEval$Eval.metric=="ROC", c("Model.name","Testing.data", "Evaluating.data")]


# Projections
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = test_predictors_cropped,
  proj.name = 'current',
  selected.models = c('bank.vole_AllData_RUN1_GBM','bank.vole_AllData_RUN1_RF','bank.vole_AllData_RUN2_RF','bank.vole_AllData_RUN3_GBM',
                      'bank.vole_AllData_RUN3_RF'),
  #binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# get projected map
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj
# get projection of single best individual model 
bestIndiv <- myCurrentProj$bank.vole_AllData_RUN1_RF
plot(bestIndiv)

# ensemble forecasting
myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProj,
  selected.models = grep('_EMmedianByROC',get_built_models(
    myBiomodEM), value=TRUE)
)
# plot ensembles
plot(myBiomodEF)

# conversion to geotiff for circuitscape
mygrd <- raster("bank.vole/proj_current/proj_current_bank.vole_ensemble.grd")
mygrd
# replace nas (the sea) and 0s with the minimum value to allow some small chance of dispersal across oceans and seas.
mygrd[is.na(mygrd)] <- min(mygrd[mygrd>0])
mygrd[mygrd==0] <- min(mygrd[mygrd>0])

plot(mygrd)
mygrd
# write ensemble raster as ascii for use in circuitscape
myasc <- writeRaster(mygrd, paste("circuitscape/",SPECIES,"_ensemble.asc", sep=""),format="ascii",overwrite=TRUE)
