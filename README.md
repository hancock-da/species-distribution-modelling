# species-distribution-modelling
*This is a small part of a wider research project investigating the drivers behind genetic differentiation in mammals and includes data for just a single species - the Bank Vole.*

This species distribution model predicts the possible distribution of bank voles in Europe by predicting the presence/absence of occurrence based on coordinates of known sightings with predictor variables including bioclimatic features (average rainfall/temperature etc) and human impact spatial data (population density, roads etc). To achieve this a number of different individual algorithms were tested including Random Forests, Artificial Neural Networks, Generalized Linear Models, Boosted Regression Trees and MAXENT using the biomod2 package in R (Thuiler et al. 2009). The best performing models according to the Receiver Operating Characteristic (ROC) score were merged to build an ensemble for presence prediction.

The resulting species distribution model is then fed into Circuitscape (Anantharaman et al., 2019), an open-source Julia program that uses electric circuit theory to model connectivity in heterogeneous landscapes. The output of pairwise comparisons between populations is a 'resistance distance' or the strength of resistance to dispersal between populations. 

The resistance distances are then evaluated in isolation_distance.R in their ability to explain genetic distance (Fst) between populations, also taking into account geographic distances between coordinates and environmental distance, calculated by performing principle component analysis on extracted values of the 19 bioclimatic variables at population coordinates. The strength of effect and significance of each distance metrics ability to explain genetic differentiation is assessed by multiple regression on distance matrices (Lichstein 2007).

# Using Circuitscape to recreate resistance maps

To run Circuitscape, the latest version of Julia should be installed. At the Julia prompt run the following commands to install Circuitscape:

    using Pkg
    Pkg.add("Circuitscape")

To run the job and reproduce the resistance maps, navigate to the circuitscape folder and run the preprepared INI file.

    cd("/species-distribution-modelling/circuitscape")

    using Circuitscape
    compute("circuitscape.ini")

# References used above

> Bank Vole genetic data from Marková, S., Horníková, M., Lanier, H.C., Henttonen, H., Searle, J.B., Weider, L.J. & Kotlík, P. (2020) ‘High genomic diversity in the bank vole at the northern apex of a range expansion: The role of multiple colonizations and end-glacial refugia’, Molecular Ecology, 29(9), pp. 1730–1744.

> Biomclimatic data from: *https://www.worldclim.org/data/bioclim.html.*

> Human impact data from Venter, O., Sanderson, E.W., Magrach, A., Allan, J.R., Beher, J., Jones, K.R., Possingham, H.P., Laurance, W.F., Wood, P., Fekete, B.M., Levy, M.A. & Watson, J.E.M. (2016) ‘Global terrestrial Human Footprint maps for 1993 and 2009’, Scientific Data, 3(1), p. 160067.

> Species occurrence records from: *GBIF.org (10 December 2021) GBIF Occurrence Download https://doi.org/10.15468/dl.u9qgut*

> Anantharaman, R., Hall, K., Shah, V. & Edelman, A. (2019) ‘Circuitscape in Julia: High Performance Connectivity Modelling to Support Conservation Decisions’, arXiv:1906.03542 [q-bio]

> Lichstein, J.W. (2007) ‘Multiple regression on distance matrices: a multivariate spatial analysis tool’, Plant Ecology, 188(2), pp. 117–131.

> Thuiller, W., Lafourcade, B., Engler, R. & Araújo, M.B. (2009) ‘BIOMOD – a platform for ensemble forecasting of species distributions’, Ecography, 32(3), pp. 369–373.



# License

This project is licensed under the MIT license - see LICENSE.txt file for more details.
