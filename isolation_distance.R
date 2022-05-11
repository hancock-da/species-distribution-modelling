library(pegas)
library(tabulizer)
library(dplyr)
library(plyr)
library(ade4)
library(raster)
library(stringr)
library(hierfstat)
library(vegan)
library(ecodist)

setwd("/species-distribution-modelling")

## Extracting locality and sample data from pdf files
# this function runs tabulizer on a given pdf file to extract tables between page numbers given as arguments.
# sometimes extra empty rows are given at the top of each table so there is an option to skip these.
extract_pdf_tables <- function(file, start_page, end_page, rows_to_skip) {
  final_df = tab = ldply(extract_tables(file,pages=start_page,output="data.frame"),data.frame)[FALSE,]
  for (page in start_page:end_page) {
    tab = extract_tables(file, pages=page, output="data.frame")
    df = ldply(tab,data.frame)
    df = df[-1:-rows_to_skip,]
    final_df <- rbind(final_df, df)
  }
  return(final_df)
}

# extract location data from pdf table
bv_extract <- extract_pdf_tables("supplementary_table.pdf",start_page=2,end_page=27,
                                 rows_to_skip = 4)
bv_df <- bv_extract %>% dplyr::select(X, X.1, X.3, X.4)
bv_df$X <- paste('CG_', bv_df$X, sep='')
names(bv_df) <- c('id','pop','lat','long')


# read in vcf file and preprocess
bv_vcf <- read.vcf("SNP.vcf")
bv_vcf$id <- rownames(bv_vcf)

# join snps with location data to match coordinates with individuals
bv_final <- merge(bv_df, bv_vcf)
bv_loci <- within(bv_final, rm('long','lat','id'))

loci_df <- cbind(bv_loci$pop, bv_loci[,2:ncol(bv_loci)] %>%
                   mutate_all(funs(str_replace_all(.,"C","1"))) %>%
                   mutate_all(funs(str_replace_all(.,"G","2"))) %>%
                   mutate_all(funs(str_replace_all(.,"A","3"))) %>%
                   mutate_all(funs(str_replace_all(.,"T","4"))) %>%
                   mutate_all(funs(str_replace_all(.,"/",""))) %>%
                   mutate_all(function(x) as.numeric(as.character(x))))

# extract just coordinates
coords <- unique(bv_final[,c('pop','long','lat')])
coords <- coords[!(coords$pop=='Wicken' & coords$lat == 52.28),] # remove duplicate coordinate for Wicken pop
coords <- coords[order(coords$pop),]
coords <- coords %>%
  mutate(pop = factor(pop)) %>%
  mutate_all(function(x) as.numeric(x))
coords
coords <- coords[,c('long','lat')]

## Calculate Genetic Distance ##
genet_dist <- genet.dist(loci_df,diploid=TRUE,method="Fst")
genet_dist
length(genet_dist)

## Calculate Geographic Distance ##
geog_dist <- as.dist(geod(coords))
geog_dist
length(geog_dist)


## Calculate environmental distance ##
path <- "map_files/bioclim"
bioclim_files <- list.files(path, pattern='tif$',full.names=TRUE)
# create a rasterStack of predictor variables and crop to same extent as species occurences
predictors <-  stack(bioclim_files)
# extract environmental values at coordinates
extraction <- extract(predictors, coords, 'bilinear')

# PCA for feature extraction to determine environmental distance
pca_mod <- rda(extraction, scale=TRUE)
summary(pca_mod)
biplot(pca_mod)
pca_scores <- scores(pca_mod,choices=c(1:4))$sites[,c(1:4)]
env_dist <- dist(pca_scores)
env_dist
length(env_dist)

## Calculate resistance distance ##
bvres <- read.table("circuitscape/circuitscape_resistances.out", 
                    header=TRUE, row.names=1)
max(bvres)
bvres
res_dist <- as.dist(bvres)
length(res_dist)
res_dist

## Multiple Regresion on Distance Matrices
mrm <- MRM(genet_dist ~ geog_dist + env_dist + res_dist, nperm=999, mrank=TRUE) # using spearmans rank as non-linear relationships between variables
mrm

plot(genet_dist, geog_dist)
plot(genet_dist, env_dist)
plot(genet_dist, res_dist)
mod <- lm(res_dist ~ genet_dist)
mod
abline(mod)