This folder contains all the files used in the [Applied Bayesian School 2025](https://abs25.imati.cnr.it/) on "_Spatio-temporal methods in environmental epidemiology with R_" by [Alexandra M. Schmidt (McGill University)](https://alex-schmidt.research.mcgill.ca/) and Carlo Zaccardi (Università degli Studi 'Gabriele d'Annunzio'​ di Chieti-Pescara).

## IMPORTANT NOTE
Please install the following software on your PC in advance to start your lessons smoothly:
- R (>= 4.0)
- **Windows users only**: [RTools](https://cran.r-project.org/bin/windows/Rtools/). Select the default options everywhere. To check the installation, the following commands must return `TRUE`:
```r
install.packages("devtools")
devtools::find_rtools()
```
- INLA: follow the instructions for your platform [here](https://www.r-inla.org/download-install)
- Stan: follow the instructions for your platform [here](https://mc-stan.org/install/). We will use two interfaces:
  - RStan: choose the CRAN Installer and follow the instructions.
  - CmdStanR: choose GitHub (Source) as the Installer and follow the instructions.

Please install also the following R packages (their installation may require some time, depending on your system):
```r
install.packages("tidyverse")
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("RcppEigen")
install.packages("sf")
install.packages("terra")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
devtools::install_github("ropensci/rnaturalearthhires")
install.packages("nimble")
install.packages("spdep")
install.packages("loo")
```

