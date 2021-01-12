# Watershed-Area-Model
## Code to run watershed-area model to estiamte LRPs for wCVI Chinook

-----  

**Authors: Carrie Holt, Kendra Holt, and Brooke Davis, adapted from Parken et al. (2006)**  
**Date: 2021-01-12 (ONGOING)**

-----



### Summary
This repository contains files to run the watershed-area model to estimate benchmarks and LRPs for West Coast Vancouver Island Chinook salmon. The principal file to run the model, `runWCVILRPs.R`, first runs the integrated watershed-area model within the function `runIWAM()` using a synpotic data set of 25 Chinook stocks across BC, Alaska, Washington, and Oregon. The data and model are adapted from Parken et al. (2006) and Liermann et al. (2010). The `runIWAM()` function also estiamtes SMSY and SREP for WCVI indicator stocks, inlets, and CUs based on the watershed-area regression. Then, the function `get.LRP()` is run to estimate Sgen and adjusted SMSY values for indicators stocks, inlets, and CUs. Finally,`get.LRP()` estiamtes a stock-manatement unit level LRP based on a logistic regression between the aggregate abundances and proportion of CUs above the red zone (i.e., CUs with all component inlets > inlet-level Sgen). Example plot outputs both including and excluding enhanced stocks are shown in `RmdReports/WCVI_LRPs.Rmd`.