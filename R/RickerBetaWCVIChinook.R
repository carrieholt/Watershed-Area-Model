# Calculation of CU/Inlet specific Ricker beta values for WCKI Chinook
# 6 May 2021, updated 21 Aug 2024
# See also: LRPs/CaseStudies/WCVIChinook/RunReconstructionTMBoutputs_BetaCalc.xlsx

# Consider moving this to SalmonLRP_WCVI_CK repository and then sourcing riclogArr.csv and riclogArr_nBC.csv files for productivities and downloading SREP files, as in line 114 of plot_lnalphas


# CU-specific alpha values derived from Diana Dobson in a run reconstruction
# WCVI_term_model_revisions_updated.xlsx
# Provided in 2019

remove.EnhStocks <- TRUE

cus <- c("WCVI_Nootka_Kyuquot","WCVI_North", "WCVI_South")
# lnalpha_cu <- c( 1.576, 1.433, 0.815)
lnalpha_cu <- c(1.578412, 1.529675, 1.138784)

# inlets <- c("Kyuquot", "Clayoquot", "Quatsino", "Barkley", "Nootka/Esperanza")
# lnalpha_inlet <- c(1.576, 0.815, 1.433, 0.815, 1.576)
inlets <- c("Kyuquot", "Clayoquot", "Quatsino", "Barkley", "Nootka/Esperanza", 
            "San Juan")
lnalpha_inlet_wbc <- c(1.578412, 1.138784, 1.529675, 1.138784, 1.578412, 
                       1.138784)
lnalpha_inlet_nbc <- c(1.343629, 0.817727, 1.297002, 0.817727, 1.343629, 
                       0.817727)
# See DataIn/riclogArr.csv and DataIn/riclogArr_nBC.csv in SalmonLRP_WCVI_CK 
# repository, created in R file, runWCVIChinook_projLRP.R, section 13 with bias 
# correction on or off, respectively.

dum <- data.frame(inlets=inlets, lnalpha_inlet_wbc=lnalpha_inlet_wbc, 
                  lnalpha_inlet_nbc=lnalpha_inlet_nbc)


# Read in SREP values from watershed-area model
# if(remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_noEnh.csv")
wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_AllExMH.csv")

dum2 <- wcviRPs_long %>%  filter (Param == "SREP") %>% 
  filter(Stock %in% inlets) %>% rename(inlets=Stock, SREP=Estimate) %>% select(-c(LL, UL, Param, X))
# Remove double instance of San Juan, as it's a indicator stock and an inlet
dum2 <- dum2[-1,]
dum2$SREP <- signif(dum2$SREP,2)

out <- dum2 %>% left_join(dum, by="inlets")
out <- out %>% mutate(beta=lnalpha_inlet_nbc/SREP)
out


# OLD- do not use
# These beta values are used in CUpars for WCVI Chinook

# SREP           inlets lnalpha_inlet         beta
#   649          Barkley         0.815 0.0012557781
#  7024        Clayoquot         0.815 0.0001160308
#  4802          Kyuquot         1.576 0.0003281966
#  1166 Nootka/Esperanza         1.576 0.0013516295
#  3155         Quatsino         1.433 0.0004541997