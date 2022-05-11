# species-distribution-modelling
*This script is an excerpt from my Master's research project aimed at determining the drivers behind genetic differentiation in mammals and includes data for just a single species - the Bank Vole.*

For this project I built species distribution models which predict the possible distribution of a species by modelling the presence/absence of occurrence based on bioclimatic spatial data (average rainfall/temperature etc) and urbanisation data. To achieve this a number of different individual algorithms were tested including Random Forests, Artificial Neural Networks, Generalized Linear Models, Boosted Regression Trees and MAXENT using the biomod2 package in R.

I then evaluated these models against a separate test dataset using Receiver Operating Characteristic (ROC) scores and chose the best performing algorithms to build an ensemble with differing internal cross-validation.

The performance of ensembles is often assumed to be better than individual models and comparisons are rarely made or reported (Hao et al. 2019). I therefore compared performances in predicting occurrences in the test dataset finding that in the case of the Bank Vole, the ensemble performed worse than individual Boosted Regression Tree models and Random Forest models. The ensembles will further be evaluated against individual models in their ability to best explain the genetic structure of the species in future work.

> Bank Vole genetic data from Marková, S., Horníková, M., Lanier, H.C., Henttonen, H., Searle, J.B., Weider, L.J. & Kotlík, P. (2020) ‘High genomic diversity in the bank vole at the northern apex of a range expansion: The role of multiple colonizations and end-glacial refugia’, Molecular Ecology, 29(9), pp. 1730–1744.

> Biomclimatic data from: *https://www.worldclim.org/data/bioclim.html.*

> Human impact data from Venter, O., Sanderson, E.W., Magrach, A., Allan, J.R., Beher, J., Jones, K.R., Possingham, H.P., Laurance, W.F., Wood, P., Fekete, B.M., Levy, M.A. & Watson, J.E.M. (2016) ‘Global terrestrial Human Footprint maps for 1993 and 2009’, Scientific Data, 3(1), p. 160067.

> Species occurrence records from: *GBIF.org (10 December 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.u9qgut*

> Hao, T., Elith, J., Guillera-Arroita, G. & Lahoz-Monfort, J.J. (2019) ‘A review of evidence about use and performance of species distribution modelling ensembles like BIOMOD’, Diversity and Distributions, 25(5), pp. 839–852.

# Using Circuitscape

To run Circuitscape, the latest version of Julia should be installed. At the Julia prompt run the following commands to install Circuitscape:

    using Pkg
    Pkg.add("Circuitscape")

To run the job and reproduce the resistance maps, navigate to the circuitscape folder and run the preprepared INI file.

    cd("/species-distribution-modelling/circuitscape")

    using Circuitscape
    compute("circuitscape.ini")

# License

This project is licensed under the MIT license - see LICENSE.txt file for more details.
