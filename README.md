# species-distribution-modelling
*This script is an excerpt from my Master's research project aimed at determining the drivers behind genetic differentiation in mammals and includes data for just a single species - the Bank Vole.*

For this project I built species distribution models which predict the possible distribution of a species by modelling the presence/absence of occurrence based on bioclimatic spatial data (average rainfall/temperature etc). To achieve this a number of different individual algorithms were employed including Random Forests, Artificial Neural Networks, Generalized Linear Models, Boosted Regression Trees and MAXENT using the biomod2 package in R.

I then evaluated these models against a separate test dataset using ROC and TSS scores and chose the best performing models to build an ensemble.

The performance of ensembles is often assumed to be better than individual models and comparisons are rarely made or reported (Hao et al. 2019). I therefore compared performances on the separate test dataset finding that in the case of the Bank Vole, the ensemble performed similarly to the best individual Random Forests model.

> Biomclimatic data from: *https://www.worldclim.org/data/bioclim.html.*

> Species occurrence records from: *GBIF.org (10 December 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.u9qgut*

> Hao, T., Elith, J., Guillera-Arroita, G. & Lahoz-Monfort, J.J. (2019) ‘A review of evidence about use and performance of species distribution modelling ensembles like BIOMOD’, Diversity and Distributions, 25(5), pp. 839–852.
