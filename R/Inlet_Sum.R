
# Create Inlet_Sums.csv

remove.EnhStocks <- FALSE
CoreInd <- TRUE#FALSE#TRUE
AllExMH <- FALSE


#----------------------------------------------------------------------------
# Sum escapements across indicators within inlets
#----------------------------------------------------------------------------

WCVIEsc <- data.frame(read.csv("DataIn/WCVIEsc.csv", row.names="Yr")) %>% 
  dplyr::select (-"Little.Zeballos")

# Take "." out of name as in escapement data
WCVIEsc_names <- sapply(colnames(WCVIEsc), 
                        function(x) (gsub(".", " ", x, fixed=TRUE) ) )
WCVIEsc_names <- sapply(WCVIEsc_names, function(x) 
  (gsub("Bedwell Ursus", "Bedwell/Ursus", x, fixed=TRUE) ) )
WCVIEsc_names <- sapply(WCVIEsc_names, function(x) 
  (gsub("Nootka Esperanza", "Nootka/Esperanza", x, fixed=TRUE) ) )
colnames(WCVIEsc) <- WCVIEsc_names 

EnhStocks <- data.frame(read.csv("DataIn/WCVIstocks.csv")) %>% filter (Enh==1) %>%
  pull(Stock)
EnhStocks <- as.character(EnhStocks)

#EnhStocks <- c("Burman",  "Conuma", "Leiner", "Nitinat", "Sarita",  
#               "Somass",  "Zeballos", "San Juan", "Tranquil")
# Artlish removed from Enhanced stocks 23 Dec. 2020
# Tranquil added 18 Jan 2021


if (remove.EnhStocks) {WCVIEsc <- WCVIEsc %>% dplyr::select(-EnhStocks) }



Years <- rownames(WCVIEsc)

# Get stock information for WCVI Chinook & Remove Cypre as it's not an 
# indicator stocks
WCVIStocks <- read.csv("DataIn/WCVIStocks.csv") %>% 
  filter (Stock != "Cypre")
if(!CoreInd & !AllExMH) {
  if (remove.EnhStocks) WCVIStocks <- WCVIStocks %>% 
      filter(Stock %not in% EnhStocks) 
}
if(CoreInd){
  WCVIStocks <- WCVIStocks %>% 
    filter(CoreInd == 1)
}
if(AllExMH){
  WCVIStocks <- WCVIStocks %>% 
    filter(highEnh == 0)
}
Inlet_Names <- unique(WCVIStocks$Inlet)
Inlet_Sum <- matrix(NA, nrow=length(Years), ncol=length(Inlet_Names))
colnames(Inlet_Sum) <- Inlet_Names
# CU_Names <- unique(WCVIStocks$CU)


# Sum escapements across stocks within inlets
for (i in 1:length(Inlet_Names)) {
  # For each inlet, which are the component indicator stocks
  Ins <- WCVIStocks %>% filter(Inlet==Inlet_Names[i]) %>% pull(Stock)
  WCVIEsc_Inlets <- matrix(NA, nrow= length(Years), ncol= length(Ins))
  
  #  Make a matrix of escapements of component indicators
  for (j in 1:length(Ins)){
    WCVIEsc_Inlets[,j] <- WCVIEsc %>% 
      dplyr::select(as.character(Ins[j])) %>% pull()
    
  }
  
  # Sum the escapement of component indicators, setting sum=NA for years 
  # where there are any NAs
  Inlet_Sum[,i] <- apply(WCVIEsc_Inlets, 1, sum, na.rm=F)
}

# Create csv files of inlet-level escapements for projection-based LRP
Inlet_Sum_csv <- as.data.frame(Inlet_Sum) %>% 
  add_column(Years=as.numeric(Years))
Inlet_Sum_csv <- Inlet_Sum_csv %>%  pivot_longer(cols=Inlet_Names) %>% 
  rename(BroodYear=Years, Inlet=name, Escapement=value)

n_inlets <- length(WCVIStocks %>% select(CU, Inlet) %>% distinct() %>% pull(Inlet))

WCVIStocks_inlets <- WCVIStocks %>% select(CU, Inlet) %>% distinct() %>% 
  add_column(Inlet_ID=1:n_inlets)


Inlet_Sum_csv <- Inlet_Sum_csv %>% left_join(WCVIStocks_inlets, by="Inlet")

Inlet_Sum_csv <-  Inlet_Sum_csv %>% rename(Inlet_Name=Inlet, CU_Name=CU, 
                                           Spawners=Escapement) %>% 
  add_column(Recruits = NA)

# Inlet_Sum_csv <- Inlet_Sum_csv %>% add_column(row.names(WCVIEsc))

if(remove.EnhStocks & !CoreInd & !AllExMH) write.csv(Inlet_Sum_csv, file="DataOut/Inlet_Sum.csv", row.names=F)
if(!remove.EnhStocks & !CoreInd & !AllExMH) write.csv(Inlet_Sum_csv, file="DataOut/Inlet_Sum_wEnh.csv", row.names=F)
# Copy this file over to SalmonLRP_RetroEval/WCVICaseStudy/DataIn/


if(CoreInd) write.csv(Inlet_Sum_csv, "DataOut/Inlet_Sum_CoreInd.csv") 
if(AllExMH) write.csv(Inlet_Sum_csv, "DataOut/Inlet_Sum_AllExMH.csv") 


