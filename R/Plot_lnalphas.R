#-------------------------------------------------------------------------
# Code to plot productivity options
# Steps:
#---------------------------------------------------------

#-------------------------------------------------------------------------------
# Libaries and source code
#-------------------------------------------------------------------------------
require(tidyverse)
require(ggplot2)
source(here::here('R/helperFunctions.R'))
#-------------------------------------------------------------------------------
# Inputs
#-------------------------------------------------------------------------------
remove.EnhStocks <- TRUE

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

test <- seq(0,4, len=40)
lna_sigma <- 0.51
plot(x=test, y=dnorm(test, 1,lna_sigma), type="l", xlab="LogA",
ylab="Probability Density", ylim=c(0,1))
# With this lna_sigma, 95% of probabilty density is within bounds mean
# +/- 1.0 (assuming range 0-2, mean=1). 0.51*1.96 = 1.0

lnalpha_cu <- read.csv("DataIn/CUPars_wBC.csv") %>%
  select(alpha,stkName) %>% rename(inlets=stkName, lnalpha=alpha)
lines(x=test, y=dnorm(test, lnalpha_cu$lnalpha[1], lna_sigma), col="grey")
lines(x=test, y=dnorm(test, lnalpha_cu$lnalpha[2], lna_sigma), col="grey")
lines(x=test, y=dnorm(test, lnalpha_cu$lnalpha[3], lna_sigma), col="grey")
# 



if (remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv")
if (!remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_wEnh_wBC.csv")

stocks <- unique(wcviRPs_long$Stock)

# Pull an example SMSY and SREP to calculate the inferred productivity
# lnalpha_Parkin
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


# hist(out)
# abline(v=exampleSMSY)
# abline(v=exampleSMSY_LL, lty="dotted")
# abline(v=exampleSMSY_UL, lty="dotted")

# Function to estimate log(alpha) from SMSY and SREP, using Scheurell (2016)
# form equation for SMSY and SREP=ln(alpha)/beta

# NEXT: Run this function over all SMSY and SREPs to get a distribution of ln_alpha, and Sgen (and LL and UL of Sgen)

# Functions to estimate inferred log(a) from SMSY and SREP from IWAM or Parken et al
calc_loga <- function(loga, SMSY, SREP){
  # 1 - log( gsl::lambert_Wm1(1 - SMSY * loga / SREP)) - loga
  # (1 - gsl::lambert_W0 (exp( 1 - loga))) * (SREP / SMSY) - loga 
  abs( (1-(loga/SREP) * SMSY) * exp(loga - (loga/SREP) * SMSY) - 1)
  
}
# lnalpha_Parkin_example <- nlminb(start = (0.5 - exampleSMSY/exampleSREP) / 0.07, objective = calc_loga,  SMSY = exampleSMSY, SREP = exampleSREP)$par

est_loga <- function(SMSY, SREP, shortloga=FALSE){
  
  loga <- nlminb(start = (0.5 - SMSY/SREP) / 0.07, 
                 objective = calc_loga, 
                 SMSY= SMSY, 
                 SREP=SREP)$par
  if(shortloga) loga <- (0.5 - SMSY/SREP) / 0.07
  beta <- loga/SREP
  return( list( loga = loga , beta = beta, SMSY = SMSY, SREP = SREP) )
}

# est_loga(SMSY = exampleSMSY, SREP = exampleSREP)

# Calculate inferred log(a) for each population in the the synoptic data set 
# (from their derived SMSY and SREP)
wcvi_SMSY <- wcviRPs_long %>% filter(Param == "SMSY") %>% select(Estimate) 
wcvi_SREP <- wcviRPs_long %>% filter(Param == "SREP") %>% select(Estimate) 
# See that the ratio is ~ constant at about 1/3 SMSY:SREP
ratio <- wcvi_SMSY$Estimate/wcvi_SREP$Estimate

# This function returns inferred alpha values from the watershed-area, as used for Parken. However, the underling watershed-area model is iwam, which differs from Parken
# The inferred productivities are ~ constant at 2.2
lnalpha_Parkin <- purrr::map2_dfr (wcvi_SMSY, wcvi_SREP, shortloga=FALSE, est_loga)
hist(lnalpha_Parkin$loga)
abline(v=median(lnalpha_Parkin$loga), lty="dashed")


# These are extracted here to demonstrate the impact of using the
# Calculate distribution of inferred log(a) for one population in the the 
# synoptic data set (from the distributions of their derived SMSYs and SREPs)

## Cannot run this function on all randomly drawn SMSY and SREP, as the relationship between SMSY and SREP determine alpha
## lnalpha_Parkin <- purrr::map2_dfr (SMSY_mc, SREP_mc, est_loga)

# Instead assume median lnalpha_Parkin and random draws of SREP to estimate Sgen
# This give meanSgens, but is recalculated below.
wcvi_SGEN <- purrr::map2_dfr (exp(median(lnalpha_Parkin$loga)), 
                              wcvi_SREP$Estimate, 
                              Sgen.fn2)

SGEN_out <- data.frame(Value=double(), lwr=double(), upr=double(), Stock=character(), RP=character())


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
  SREP_mc <- data.frame(SREP = rlnorm(1000,log(meanSREP), sd_logSREP) )
  # Could add back-transformation correction here: actual differences are relatively minor
  # SREP_mc <- data.frame(SREP = rlnorm(1000,log(meanSREP)-sd_logSREP^2/2, sd_logSREP) )
  
  SGEN_mc <- purrr::map2_dfr (exp(median(lnalpha_Parkin$loga)), 
                              SREP_mc$SREP, 
                              Sgen.fn2)
  LLSGEN <- quantile(SGEN_mc$SGEN, p=c(0.025))
  ULSGEN <- quantile(SGEN_mc$SGEN, p=c(0.975))
  medSGEN <- quantile(SGEN_mc$SGEN, p=c(0.5))
  
  meanSGEN <- Sgen.fn2(exp(median(lnalpha_Parkin$loga)),  meanSREP )$SGEN
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

wcviCK_BootstrappedRPs_ParkenProd <- wcviRPs_long %>% transmute(Value=Estimate, lwr=LL, upr=LL, Stock=Stock, RP=Param)
wcviCK_BootstrappedRPs_ParkenProd <- rbind(wcviCK_BootstrappedRPs_ParkenProd, SGEN_out)
write.csv(wcviCK_BootstrappedRPs_ParkenProd, here::here("DataOut", "wcviCK-BootstrappedRP-ParkenProd.csv") )


# # Check that backcalculated SMSY is close 
# examplebeta <- lnalpha_Parkin_example/exampleSREP
# backcalculatedSMSY <- (1 - gsl::lambert_W0( exp( 1 - lnalpha_Parkin_example))) / examplebeta


# Plot alpha distributions
set.seed(123)
x <- c(rnorm(1000, 1, lna_sigma), 
       rnorm(1000, lnalpha_cu$lnalpha[1], lna_sigma), 
       rnorm(1000, lnalpha_cu$lnalpha[2], lna_sigma), 
       rnorm(1000, lnalpha_cu$lnalpha[3], lna_sigma) )
group <- c(rep("Life-History-Model", 1000), 
           rep("Run-Reconstruction", 1000*3))
cu <- c(rep("Life-History-Model", 1000),
        rep("Run Reconstruction-South", 1000),
        rep("Run Reconstruction-Nootka Kyuquot", 1000),
        rep("Run Reconstruction-North", 1000))

df <- data.frame(Value = x, Productivity = group, cu=cu)
# df <- df %>% filter(Productivity != "Life-History-Model")
alpha_plot <- ggplot(df) + 
  aes(x = Value, fill = Productivity, alpha = cu, linetype = Productivity) + 
  geom_density(bw=0.2) + 
  scale_alpha_manual(values = c("Life-History-Model" = 0.5,
                                "Run Reconstruction-South" = 0.2,
                                "Run Reconstruction-Nootka Kyuquot" = 0.2,
                                "Run Reconstruction-North" = 0.2), 
                     guide = "none") + 
  scale_linetype_manual(values = c("Life-History-Model" = "solid", 
                                   "Run-Reconstruction" = "dashed")) +
  labs(x = "ln(alpha)", y = "Density") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = alpha('black', 0.1))) + 
  theme(legend.title = element_blank(), legend.position = "bottom") + 
  geom_vline (xintercept = median(lnalpha_Parkin$loga), linetype = "dotted") +
  xlim(c(0,3))
ggsave(plot=alpha_plot,filename=here::here("DataOut/alpha_plots_RR.png"))

# plot Sgen values
# 

# LH_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRPs.csv"))
# RR_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRPs_ProdRR.csv"))
# Pa_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRP-ParkenProd.csv"))
LH_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRPs_ExtInd.csv"))
RR_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRPs_ExtInd_RunReconstruction.csv"))
Pa_sgens <- read.csv(here::here("DataOut", "wcviCK-BootstrappedRPs_ExtInd_Parken.csv"))

core_ind <- read.csv(here::here("DataIn", "WCVIStocks_ExtInd.csv"))

LH_sgens <- LH_sgens %>% mutate(Prod="LifeHistory")
RR_sgens <- RR_sgens %>% mutate(Prod="RunReconstruction")
Pa_sgens <- Pa_sgens %>% mutate(Prod="Parken")

df <- rbind(LH_sgens, RR_sgens, Pa_sgens)
df <- left_join(df, core_ind, by="Stock")
df <- df %>% filter(CoreInd==1)

benchmark <- "SGEN"#"SMSY"#"SREP"

df <- df %>% filter(RP==benchmark)
order <- c("LifeHistory", "RunReconstruction", "Parken")
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
       filename=here::here("DataOut", paste(benchmark,"_plot_prod_assumption.png", sep="")))




