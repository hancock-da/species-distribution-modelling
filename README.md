# species-distribution-modelling
*This script is an excerpt from my Master's research project aimed at determining the drivers behind genetic differentiation in mammals and includes data for just a single species - the Bank Vole.*

For this project I built species distribution models which predict the possible distribution of a species by modelling the presence/absence of occurrence based on bioclimatic spatial data (average rainfall/temperature etc) and urbanisation data. To achieve this a number of different individual algorithms were employed including Random Forests, Artificial Neural Networks, Generalized Linear Models, Boosted Regression Trees and MAXENT using the biomod2 package in R.

I then evaluated these models against a separate test dataset using Receiver Operating Characteristic (ROC) scores and chose the best performing algorithm on average to build an ensemble of three separate runs with differing internal cross-validation.

The performance of ensembles is often assumed to be better than individual models and comparisons are rarely made or reported (Hao et al. 2019). I therefore compared performances in predicting occurrences in the test dataset finding that in the case of the Bank Vole, the ensemble performed worse than individual Boosted Regression Tree models. The ensembles will further be evaluated against individual models in their ability to best explain the genetic structure of the species in future work.

> Biomclimatic data from: *https://www.worldclim.org/data/bioclim.html.*

> Urbanisation data from *https://land.copernicus.eu/global/*

> Species occurrence records from: *GBIF.org (10 December 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.u9qgut*

> Hao, T., Elith, J., Guillera-Arroita, G. & Lahoz-Monfort, J.J. (2019) ‘A review of evidence about use and performance of species distribution modelling ensembles like BIOMOD’, Diversity and Distributions, 25(5), pp. 839–852.

# License

This project is licensed under the MIT license - see LICENSE.txt file for more details.
