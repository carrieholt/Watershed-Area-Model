#-------------------------------------------------------------------------------
# Code to plot various assumptions about productivity for WCVI Chinook and its 
# impacts on benchmark estimates for WCVI Chinook salmon, derived from the 
# accessible watershed-area model, adapted from Parken et al. (2006).
# Plots required for Brown et al. (in revision)
# 
# Brown, N. et al. (in revision). West Coast of Vancouver Island Natural-Origin 
# Chinook Salmon (Oncorhynchus tshawytscha) Stock Assessment. CSAS Working 
# Paper20xx/nnn.
# Parken, C.K., McNicol, R.E., and Irvine, J.R. 2006. Habitat-based methods to 
# estimate escapement goals for data limited Chinook salmon stocks in British 
# Columbia, 2004. DFO Can. Sci. Advis. Res. Doc. 2006/83. vii + 67 p.
#---------------------------------------------------------

#-------------------------------------------------------------------------------
# Libraries and source code
#-------------------------------------------------------------------------------
require(tidyverse)
require(ggplot2)
source(here::here('R/helperFunctions.R'))

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
# These functions are used to estimate inferred log(a) when SMSY and SREP are 
# estimated from the accessible watershed-area model (Parken et al. 2006)

# Function to minimize to estimate loga, where beta=logA/SREP
calc_loga <- function(loga, SMSY, SREP){
  abs( (1-(loga/SREP) * SMSY) * exp(loga - (loga/SREP) * SMSY) - 1)
  # Scheurell (2016) forms- not currently working
  # 1 - log( gsl::lambert_Wm1(1 - SMSY * loga / SREP)) - loga
  # (1 - gsl::lambert_W0 (exp( 1 - loga))) * (SREP / SMSY) - loga 
  
}

# Estimation of loga from SMSY and SREP given minimization function above
est_loga <- function(SMSY, SREP, shortloga=FALSE){
  
  loga <- nlminb(start = (0.5 - SMSY/SREP) / 0.07, 
                 objective = calc_loga, 
                 SMSY= SMSY, 
                 SREP=SREP)$par
  if(shortloga) loga <- (0.5 - SMSY/SREP) / 0.07
  beta <- loga/SREP
  return( list( loga = loga , beta = beta, SMSY = SMSY, SREP = SREP) )
}

# lnalpha_Parken_example <- nlminb(start = (0.5 - exampleSMSY/exampleSREP) /0.07, 
#                                  objective = calc_loga,  SMSY = exampleSMSY, 
#                                  SREP = exampleSREP)$par



#-------------------------------------------------------------------------------
# Code to identify standard deviation in loga that is aligned with uncertainty
# intervals (0-2, mean =1) derived from life-cycle model (Luedke pers. comm.)
#-------------------------------------------------------------------------------

test <- seq(0,4, len=40)
lna_sigma <- 0.51
plot(x=test, y=dnorm(test, 1,lna_sigma), type="l", xlab="LogA",
ylab="Probability Density", ylim=c(0,1))
# With this lna_sigma (0.51), 95% of probabilty density is within bounds mean
# +/- 1.0 (assuming range 0-2, mean=1). 0.51*1.96 = 1.0
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Code to get estimates of inferred productivity from accessible watershed-area
# model, i.e., using assumptions about productivity as in Parken et al. (2006)
#-------------------------------------------------------------------------------

# Indicators used in Holt, K. et al. (2023)
# if(remove.EnhStocks) wcviRPs_long<-read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv")
# if(!remove.EnhStocks) wcviRPs_long<-read.csv("DataOut/WCVI_SMSY_wEnh_wBC.csv")

# Indicators used in Brown et al. (in revision)- bar plot, Fig. 4.2, all 
# escapement indicators including enhanced stocks
wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_ExtInd.csv")
stocks <- unique(wcviRPs_long$Stock)

# Calculate inferred log(a) for each population in the the synoptic data set 
# (from their derived SMSY and SREP)
# First get SMSY and SREP values: 
wcvi_SMSY <- wcviRPs_long %>% filter(Param == "SMSY") %>% select(Estimate) 
wcvi_SREP <- wcviRPs_long %>% filter(Param == "SREP") %>% select(Estimate) 

## See that the ratio is ~ constant at about 1/3 SMSY:SREP
ratio <- wcvi_SMSY$Estimate/wcvi_SREP$Estimate

# This function returns inferred alpha values from the watershed-area, as used 
# for Parken et al. (2006). However, the underling watershed-area model is iwam, 
# which differs from Parken et al. (2006). The inferred productivity is ~ 
# constant across stocks at 2.2
lnalpha_Parken <- purrr::map2_dfr (wcvi_SMSY, wcvi_SREP, shortloga=FALSE, 
                                   est_loga)
hist(lnalpha_Parken$loga)
abline(v=median(lnalpha_Parken$loga), lty="dashed")


#-------------------------------------------------------------------------------
# Code to plot various assumptions about productivity: (1) from life-cycle model
# for WCVI Chinoook salmon (Luedke pers. comm), from run-reconstruction and
# ricker stock-recruitment modelling (see repository:
# pacific-salmon-assess/SalmonLRP_wCVI_CK, file: runWCVIChinook_projLRP.r, and
# (3) as inferred from accessible watershed-area model, as in Parken et al.
# (2006)
#-------------------------------------------------------------------------------

# Get MLE estimates of productivity from run reconstruction, assuming simple 
# Ricker model. See Pacific-salmon-assess/SalmonLRP_WCVI_CK repository, file
# runWCVIChinook_ProjLRPs.r, Section 13

download.file('https://raw.githubusercontent.com/Pacific-salmon-assess/SalmonLRP_WCVI_CK/main/WCVIChinookStudy/DataIn/riclogArr.csv', './DataIn/riclogArr.csv',  mode="wb")

# Read in productivity (logA) values
lnalpha_cu <- read.csv("DataIn/riclogArr.csv")
# order: WCVI South, WCVI N&K, WCVI North

# Plot alpha distributions
set.seed(123)
x <- c(rnorm(1000, 1, lna_sigma), 
       rnorm(1000, lnalpha_cu$logRicA[1], lna_sigma), 
       rnorm(1000, lnalpha_cu$logRicA[2], lna_sigma), 
       rnorm(1000, lnalpha_cu$logRicA[3], lna_sigma) )
group <- c(rep("Life-Cycle-Model", 1000), 
           rep("Run-Reconstruction", 1000*3))
cu <- c(rep("Life-Cycle-Model", 1000),
        rep("Run Reconstruction-South", 1000),
        rep("Run Reconstruction-Nootka Kyuquot", 1000),
        rep("Run Reconstruction-North", 1000))

df <- data.frame(Value = x, Productivity = group, cu=cu)
# df <- df %>% filter(Productivity != "Life-History-Model")
alpha_plot <- ggplot(df) + 
  aes(x = Value, fill = Productivity, alpha = cu, linetype = Productivity) + 
  geom_density(bw=0.2) + 
  scale_alpha_manual(values = c("Life-Cycle-Model" = 0.5,
                                "Run Reconstruction-South" = 0.2,
                                "Run Reconstruction-Nootka Kyuquot" = 0.2,
                                "Run Reconstruction-North" = 0.2), 
                     guide = "none") + 
  scale_linetype_manual(values = c("Life-Cycle-Model" = "solid", 
                                   "Run-Reconstruction" = "dashed")) +
  labs(x = "ln(alpha)", y = "Probability density") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = alpha('black', 0.1))) + 
  theme(legend.title = element_blank(), legend.position = "bottom") + 
  geom_vline (xintercept = median(lnalpha_Parken$loga), linetype = "dotted") +
  xlim(c(0,3))
ggsave(plot=alpha_plot,filename=here::here("DataOut/alpha_plots.png"))





# -----------------------------------------------------------------------------
# Plot Sgen, SMSY and SREP values under 3 assumptions about productivity for 
# WCVI Chinook escapement indicators
# Input: *.csv files generated in 'WCVILRPs_bootstrap.R', using one of the three 
# assumptions about productivity in the function, Get.LRP.bs()
# -----------------------------------------------------------------------------


# LC_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRPs.csv"))
# RR_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRPs_ProdRR.csv"))
# Pa_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRP-ParkenProd.csv"))
LC_sgens <- read.csv(here::here("DataOut",
                                "wcviCK-BootstrappedRPs_ExtInd_LifeCycleModel.csv"))
RR_sgens <- read.csv(here::here("DataOut",
                                "wcviCK-BootstrappedRPs_ExtInd_RunReconstruction.csv"))
Pa_sgens <- read.csv(here::here("DataOut",
                                "wcviCK-BootstrappedRPs_ExtInd_Parken.csv"))

core_ind <- read.csv(here::here("DataIn", "WCVIStocks_ExtInd.csv"))

LC_sgens <- LC_sgens %>% mutate(Prod="LifeCycle")
RR_sgens <- RR_sgens %>% mutate(Prod="RunReconstruction")
Pa_sgens <- Pa_sgens %>% mutate(Prod="Parken")

df <- rbind(LC_sgens, RR_sgens, Pa_sgens)
df <- left_join(df, core_ind, by="Stock")
df <- df %>% filter(CoreInd==1)

benchmark <- "SGEN"#"SMSY"#"SREP"#

df <- df %>% filter(RP==benchmark)
order <- c("LifeCycle", "RunReconstruction", "Parken")
df <-   df %>% mutate(Prod=factor(Prod, levels=order))
upper_y_lim <- max(df$upr)
Benchmark_plot_prod_assumption <-ggplot(df, 
                                   aes(x=Stock, y=Value,#ymin=lwr, ymax=upr, 
                                       fill=as.factor(Prod))) + 
       # aes(x=reorder(Prod <- factor(Prod, 
       #                                  levels=names("LifeHistory", "RunReconstruction", "Parken"))))) +
  geom_bar(aes( y=Value), stat="identity",
           position=position_dodge(), alpha=0.7) + 
  geom_errorbar(aes(x=Stock, y=Value,ymin=lwr, ymax=upr),
                width=0.1,
                alpha = 0.9,
                size=0.8,
                position=position_dodge(width=0.9), colour="grey") +
  coord_cartesian(ylim=c(0, upper_y_lim)) +
  theme(legend.position = "bottom", legend.text=element_text(size=8)) + 
  theme(strip.background = element_rect(fill = alpha('black', 0.1))) + 
  ylab(benchmark) +
  xlab(element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  guides(fill=guide_legend(title="Assumption about Productivity:")) +
  theme(legend.title=element_text(size=8))
  
  
ggsave(plot=Benchmark_plot_prod_assumption, 
       filename=here::here("DataOut", 
                           paste(benchmark,"_plot_prod_assumption.png", 
                                 sep="")))



#-------------------------------------------------------------------------------
# Test plots and code- not needed
#-------------------------------------------------------------------------------

lnalpha_cu <- read.csv("DataIn/CUPars_wBC.csv") %>%
  select(alpha,stkName) %>% rename(inlets=stkName, lnalpha=alpha)
lines(x=test, y=dnorm(test, lnalpha_cu$lnalpha[1], lna_sigma), col="grey")
lines(x=test, y=dnorm(test, lnalpha_cu$lnalpha[2], lna_sigma), col="grey")
lines(x=test, y=dnorm(test, lnalpha_cu$lnalpha[3], lna_sigma), col="grey")


if (remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv")
if (!remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_wEnh_wBC.csv")

stocks <- unique(wcviRPs_long$Stock)

# Pull an example SMSY and SREP to calculate the inferred productivity
# lnalpha_Parken
exampleSMSY <- wcviRPs_long %>% 
  filter(Stock=="Artlish"&Param=="SMSY") %>% 
  pull(Estimate)
exampleSMSY_LL <- wcviRPs_long %>% 
  filter(Stock=="Artlish"&Param=="SMSY") %>% 
  pull(LL)
exampleSMSY_UL <- wcviRPs_long %>% 
  filter(Stock=="Artlish"&Param=="SMSY") %>% 
  pull(UL)

exampleSREP <- wcviRPs_long %>% 
  filter(Stock=="Artlish"&Param=="SREP") %>% 
  pull(Estimate)
exampleSREP_LL <- wcviRPs_long %>% 
  filter(Stock=="Artlish"&Param=="SREP") %>% 
  pull(LL)
exampleSREP_UL <- wcviRPs_long %>% 
  filter(Stock=="Artlish"&Param=="SREP") %>% 
  pull(UL)

sd_logSMSY <- (log(exampleSMSY_UL)-log(exampleSMSY)) / 1.96
sd_logSREP <- (log(exampleSREP_UL)-log(exampleSREP)) / 1.96
SMSY_mc <- data.frame(SMSY = rlnorm(1000,log(exampleSMSY), sd_logSMSY) )
SREP_mc <- data.frame(SREP = rlnorm(1000,log(exampleSREP), sd_logSREP) )



# NEXT: Run this function over all SMSY and SREPs to get a distribution of ln_alpha, and Sgen (and LL and UL of Sgen)


# est_loga(SMSY = exampleSMSY, SREP = exampleSREP)

# Calculate inferred log(a) for each population in the the synoptic data set 
# (from their derived SMSY and SREP)
wcvi_SMSY <- wcviRPs_long %>% filter(Param == "SMSY") %>% select(Estimate) 
wcvi_SREP <- wcviRPs_long %>% filter(Param == "SREP") %>% select(Estimate) 
# See that the ratio is ~ constant at about 1/3 SMSY:SREP
ratio <- wcvi_SMSY$Estimate/wcvi_SREP$Estimate

# This function returns inferred alpha values from the watershed-area, as used for Parken. However, the underling watershed-area model is iwam, which differs from Parken
# The inferred productivities are ~ constant at 2.2
lnalpha_Parken <- purrr::map2_dfr (wcvi_SMSY, wcvi_SREP, shortloga=FALSE, est_loga)
hist(lnalpha_Parken$loga)
abline(v=median(lnalpha_Parken$loga), lty="dashed")


# These are extracted here to demonstrate the impact of using the
# Calculate distribution of inferred log(a) for one population in the the 
# synoptic data set (from the distributions of their derived SMSYs and SREPs)

## Cannot run this function on all randomly drawn SMSY and SREP, as the relationship between SMSY and SREP determine alpha
## lnalpha_Parken <- purrr::map2_dfr (SMSY_mc, SREP_mc, est_loga)


#-------------------------------------------------------------------------------
# Code to estimate Sgen benchmarks for WCVI Chinook escapement indicators  
# based on assumption about productivity from accessible watershed-area model
# as in Parken et al. (2006)
# NOT NEEDED, as boostrapped values already provided in WCVILRPs_bootstrap.R
#-------------------------------------------------------------------------------

# Use median lnalpha_Parken and draw random draws of SREP to estimate Sgen
# This give meanSgens, but is recalculated below for all stocks
# wcvi_SGEN <- purrr::map2_dfr (exp(median(lnalpha_Parken$loga)), 
#                               wcvi_SREP$Estimate, 
#                               Sgen.fn2)

# Use stock-specific lnalpha_Parken and draw random draws of SREP to estimate
# stock-specific Sgen values
SGEN_out <- data.frame(Value=double(), lwr=double(), upr=double(), 
                       Stock=character(), RP=character())

for(i in 1:length(stocks)){
  meanSREP <- wcviRPs_long %>% 
    filter(Stock==stocks[i]&Param=="SREP") %>% 
    pull(Estimate)
  LLSREP <- wcviRPs_long %>% 
    filter(Stock==stocks[i]&Param=="SREP") %>% 
    pull(LL)
  ULSREP<- wcviRPs_long %>% 
    filter(Stock==stocks[i]&Param=="SREP") %>% 
    pull(UL)
  
  sd_logSREP <- (log(ULSREP)-log(meanSREP)) / 1.96
  SREP_mc <- data.frame(SREP = rlnorm(1000,log(meanSREP) - 0.5*sd_logSREP^2, 
                                      sd_logSREP) )
  # Includes back-transformation correction here, but differences with and
  # without are relatively minor:
  # SREP_mc <- data.frame(SREP = rlnorm(1000,log(meanSREP), sd_logSREP) )
  
  SGEN_mc <- purrr::map2_dfr (exp(median(lnalpha_Parken$loga)), 
                              SREP_mc$SREP, 
                              Sgen.fn2)
  LLSGEN <- quantile(SGEN_mc$SGEN, p=c(0.025))
  ULSGEN <- quantile(SGEN_mc$SGEN, p=c(0.975))
  medSGEN <- quantile(SGEN_mc$SGEN, p=c(0.5))
  
  meanSGEN <- Sgen.fn2(exp(median(lnalpha_Parken$loga)),  meanSREP )$SGEN
  # LLSGEN
  # ULSGEN
  # medSGEN
  
  dum <- data.frame(Value = round(meanSGEN,0), 
                    lwr = round(LLSGEN,0), 
                    upr = round(ULSGEN,0), 
                    Stock = stocks[i], 
                    RP = "SGEN")
  SGEN_out <- rbind(SGEN_out, dum)
  SGEN_out <- SGEN_out %>% tibble::remove_rownames() 
  # round 
}

wcviCK_BootstrappedRPs_ParkenProd <- wcviRPs_long %>% 
  transmute(Value=Estimate, lwr=LL, upr=LL, Stock=Stock, RP=Param)
wcviCK_BootstrappedRPs_ParkenProd <- rbind(wcviCK_BootstrappedRPs_ParkenProd, 
                                           SGEN_out)
write.csv(wcviCK_BootstrappedRPs_ParkenProd, 
          here::here("DataOut", "wcviCK-BootstrappedRPs_Parken.csv") )


# Check that backcalculated SMSY is close 
# examplebeta <- lnalpha_Parken_example/exampleSREP
# backcalculatedSMSY <- 
# (1 - gsl::lambert_W0( exp( 1 - lnalpha_Parken_example))) / examplebeta

