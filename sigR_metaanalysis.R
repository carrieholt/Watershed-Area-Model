#Code to estimate sigR from PSE data

#---------------------------------------------------------
# Libaries

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)
library(viridis)
library(hrbrthemes)


SRDatwNA <- read.csv("DataIn/PSE_ChinookData.csv")
#need to sort out missing years, and mayabe add year column to csv