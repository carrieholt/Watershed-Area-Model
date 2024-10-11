# Watershed-Area-Model
Code to run the integrated accessible watershed-area model for WCVI Chinook

Primary contact: Carrie Holt, carrie.holt@dfo-mpo.gc.ca
Date created: 2021-01-12, updates ongoing


### Summary
This repository contains files to run the accessible watershed-area model to estimate benchmarks and logistic regression based reference points for West Coast Vancouver Island  (WCVI) Chinook salmon. The accessible watershed-area model is adapted from Parken et al. (2006) and Liermann et al. (2011) and used to derive benchmarks and reference points for Holt, K. et al. (2023) and Brown et al. (in revision). Citations are provided below. Benchmarks are provided at the population (or stock), Conservation Unit (CU) and inlet scale. For WCVI Chinook inlets are nested within CUs. Logistic-regression reference points are provided at the Stock Management Unit (SMU) scale, which represent various probabilities of all component inlets or CUs being above their lower benchmark (see Holt, K. et al. (2023) for more details).

Benchmarks for all escapement indicators except those associated with major hatcheries (AllExMH) are reported in Table 4.1 of Brown et al. (in revision). These escapement indicators are: Artlish, Bedwell/Ursus, Burman, Cayeghle, Gold, Kaouk, Leiner, Marble, Megin, Moyeha, Nahmint, San Juan, Sarita, Tahsis, Tahsish, Tranquil and Zeballos.

Benchmarks for all extensive indicators (ExtInd) are reported Appendix B of Brown et al. (in revision). There are 55 extensive indicator systems.

### Contents of repository
Code and associated files are organized into the following sub-folders:  

1. **DataIn** Folder of input data required for analyses
2. **DataOut** Folder of output results generated from analyses
3. **R** Folder of R code used in analyses
4. **RmdReports** Folder of Rmd reports used for preliminary analyses for Holt, K. et al. (2023)
5. **stan_Files** Folder of stan files for estimating the integrated watershed-area model (Not currently used)
6. **TMB_Files** Folder of TMB files used to estimate the integrated watershed-area model
 
### Steps for running analyses:
**Step 1)** Run the integrated watershed-area model using the runIWAM() function. In this function, parameters of the watershed-area-model are estimated from a synoptic data set of spawner-recruitment time-series and watershed areas, and SMSY and SREP benchmarks for out-of-sample  stocks (default being WCVI Chinook) are predicted with prediction intervals. This step can be run on different sets of WCVI Chinook indicators, such as escapement indicators excluding major hatchery facilities (AllExMH) or all extensive indicators (ExtInd).

*File*: 'R/IWAM.R' Arguments for the runIWAM() function are described at the top of IWAM.R

*Input*: 

'DataIn/SRinputfile.csv' Synoptic data set of spawner-recruitment time-series used to estimate accessible watershed-area model. These data were used in Parken et al (2006) and Liermann et al (2011) and are available upon request.

'DataIn/Surv.csv' file of marine survival covariates for Cowichan and Harrison Chinook, not currently used

'DataIn/WatershedArea.csv' Accessible watershed areas for synoptic data sets, provided in Parken et al. (2006)

'DataIn/WCVIStocks.csv' Names and accessible watershed areas for escapement indicators for WCVI Chinook

'DataIn/WCVIStocks_ExtInd.csv' Names and accessible watershed areas for extensive indicators for WCVI Chinook


*Output*: 

'DataOut/SR_std.png' Plot of Ricker spawner-recruitment fits for synoptic data sets (if plot==TRUE)

'DataOut/SRL_std.png' Plot of linearized Ricker spawner-recruitment fits for synoptic data sets (if plot==TRUE)

'DataOut/SRResid_std.png' Plot of residuals of Ricker spawner-recruitment fits for synoptic data sets (if plot==TRUE)

'DataOut/ACF_std.png' Plot of autocorrelation in residuals of Ricker spawner-recruitment fits for synoptic data sets (if plot==TRUE)

'DataOut/WAregSMSY_std_wBC.png' Plot of linear fit of watershed area vs SMSY for synoptic data set (if plot==TRUE)

'DataOut/WAregSREP_std_wBC.png' Plot of linear fit of watershed area vs SREP for synoptic data set (if plot==TRUE)

'DataOut/WCVI_SMSY_AllExMH.csv' File of benchmarks SMSY and SREP for escapement indicators except those associated with major hatchery facilities, inlets and CUs, including 95% prediction intervals (lower limit, LL and upper limit, UL).

OR

'DataOut/WCVI_SMSY_ExtInd.csv' File of benchmarks SMSY and SREP for escapement indicators except those associated with major hatchery facilities, including 95% prediction intervals (lower limit, LL and upper limit, UL).

**Step 2)** Bootstrap benchmark estimates to derive uncertainty intervals for SMSY, SREP and SGEN, and estimate logistic regression reference point (optionally) using the function Get.LRP.bs(). This function has two sections. In the first section, for each stock (or indicator), benchmarks are calculated from stock-specific SREP values and inferences about productivity (log alpha). Stock-specific bootstrapped values of SREP are drawn from a distribution defined by the prediction intervals in step 1 for that stock, and from a distribution of productivities defined by a life-cycle model, a run reconstruction, or as inferred from the synoptic data sets used in step 1 (as in Parken et al. 2006).  In the second section, logistic-regression based reference points are estimated, as described in Holt, K. et al. (2023) and in the RMD file, 'RmdReports/WCVI_LRPs_Bernoulli.Rmd'. The second section is optional.

*File*: 'R/WCVILRPs_boostrapped.R' The function Get.LRP.bs() gets a single draw of the stock-specific benchmarks from distributions of SREP and productivity, which is iterated over a large number of bootstraps (8000, the number required to stabilize benchmark values to two significant digits). The distribution of bootstrapped benchmarks, Sgen, SMSY, and SREP, are then summarized by the median, and 2.5th and 97.5th percentiles.

*Input*: 

'DataOut/WCVI_SMSY_ExtInd.csv' File of predicted SMSY and SREP values for the accessible watershed area model

*Output*: 

'DataOut/wcviCK-BootstrappedRPs_ExtInd.csv' File of SGEN, SMSY and SREP values with uncertainty intervals from bootstrapped draws of SREP and productivity.

**Step 3)** Plot distribution of productivity for WCVI Chinook under three different assumptions: (1) from a life-cycle model for WCVI Chinook salmon (Luedke pers. comm), (2) from a run-reconstruction and CU-specific ricker stock-recruitment models (see repository: pacific-salmon-assess/SalmonLRP_wCVI_CK, file: runWCVIChinook_projLRP.r), and (3) as inferred from accessible watershed-area model, as in Parken et al. (2006).

*File*: 'R/Plot_lnalphas.R'

*Inputs*:

'DataIn/riclogArr.csv' estimates of productivity from run reconstruction

'DataOut/WCVI_SMSY_ExtInd.csv' estimates of SMSY And SREP from accessible watershed-area-model used to estimate inferred productivity from that model

*Outputs*:

'DataOut/alpha_plots.png' Distrubion of productivities, as in Fig. 4.1  in Brown et al. (in revision)

**Step 4)** Bar plot of benchmarks (SGEN, SMSY, SREP) WCVI Chinook escapement indicators under assumptions about productivity: (1) from life-cycle model for WCVI Chinoook salmon (Luedke pers. comm), from run-reconstruction and Ricker stock-recruitment modelling (see repository: pacific-salmon-assess/SalmonLRP_WCVI_CK, file: runWCVIChinook_projLRP.r, and (3) as inferred from accessible watershed-area model, as in Parken et al. (2006)

*File*: 'R/Plot_lnalphas.R'

*Inputs*:

'DataOut/wcviCK-BootstrappedRPs_ExtInd_LifeCycleModel.csv' Bootstrapped benchmarks from life-cycle model for all extensive indicators (see Step 2 above)

'DataOut/wcviCK-BootstrappedRPs_ExtInd_RunReconstruction.csv' Bootstrapped benchmarks from run reconstruction for all extensive indicators (see Step 2 above)

'DataOut/wcviCK-BootstrappedRPs_ExtInd_Parken.csv' Bootstrapped benchmarks from accessible watershed-area-model (as in Parken et al. (2006)) for all extensive indicators (see Step 2 above)

*Outputs*:

'DataOut/SGEN_plot_prod_assumption.png' Bar plots of SGEN estimates under three productivity assumptions, as in Fig. 4.2 in Brown et al. (in revision)

'DataOut/SMSY_plot_prod_assumption.png' Bar plots of SMSY estimates under three productivity assumptions, as in Appendix B in Brown et al. (in revision)

'DataOut/SREP_plot_prod_assumption.png' Bar plots of SREP estimates under three productivity assumptions, as in Appendix B in Brown et al. (in revision)

### Citations
Brown, N. et al. (in revision). West Coast of Vancouver Island Natural-Origin Chinook Salmon (Oncorhynchus tshawytscha) Stock Assessment. CSAS Working Paper20xx/nnn.

Holt, K.R., Holt, C.A., Warkentin, L., Wor, C., Davis, B., Arbeider, M., Bokvist, J., Crowley, S., Grant, S., Luedke, W., McHugh, D., Picco, C., and Van  Will, P. 2023. Case Study Applications of LRP Estimation Methods to Pacific  Salmon Stock Management Units. DFO Can. Sci. Advis. Sec. Res. Doc. 2023/010.  iv+129p.

Liermann, M.C., Sharma, R., and Parken, C.K. 2010. Using accessible watershed size to predict management parameters for Chinook salmon, Oncorhynchus tshawytscha, populations with little or no spawner-recruit data: A Bayesian hierarchical modelling approach. Fisheries Management and Ecology 17(1): 40–51.

Parken, C.K., McNicol, R.E., and Irvine, J.R. 2006. Habitat-based methods to estimate escapement goals for data limited Chinook salmon stocks in British Columbia, 2004. DFO Can. Sci. Advis. Res. Doc. 2006/83. vii + 67 p.
