# Calculation of CU/Inlet specific Ricker beta values for WCKI Chinook
# 6 May 2021


# CU-specific alpha values derived from Diana Dobson in a run reconstruction
# WCVI_term_model_revisions_updated.xlsx
# Provided in 2019

remove.EnhStocks <- TRUE

cus <- c("WCVI_Nootka_Kyuquot","WCVI_North", "WCVI_South")
lnalpha_cu <- c( 1.576, 1.433, 0.815)

inlets <- c("Kyuquot", "Clayoquot", "Quatsino", "Barkley", "Nootka/Esperanza")
lnalpha_inlet <- c(1.576, 0.815, 1.433, 0.815, 1.576)
dum <- data.frame(inlets=inlets, lnalpha_inlet=lnalpha_inlet)


# Read in SREP values from watershed-area model
if(remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_noEnh.csv")

dum2 <- wcviRPs_long %>%  filter (Param == "SREP") %>% filter(Stock %in% inlets) %>% rename(inlets=Stock, SREP=Estimate) %>% select(-c(LL, UL, Param, X))

out <- dum2 %>% left_join(dum, by="inlets")
out <- out %>% mutate(beta=lnalpha_inlet/SREP)
out

# These beta values are used in CUpars for WCVI Chinook

# SREP           inlets lnalpha_inlet         beta
#   649          Barkley         0.815 0.0012557781
#  7024        Clayoquot         0.815 0.0001160308
#  4802          Kyuquot         1.576 0.0003281966
#  1166 Nootka/Esperanza         1.576 0.0013516295
#  3155         Quatsino         1.433 0.0004541997