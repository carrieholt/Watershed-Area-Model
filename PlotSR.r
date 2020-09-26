# Plot SR models: curves and linearized mddel

# Functions
# Functions
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  t.col
  #invisible(t.col)
}

# Plot SR curves

PlotSRCurve <- function(SRDat, All_Est, SMSY_std, stksNum_ar, stksNum_surv, stks_surv, r2, removeSkagit, mod) {
  Stks <- unique(SRDat$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1)
  
  for (i in Stks){
    names <- All_Est %>% select ("Name", "Stocknumber") %>% distinct()
    name <- All_Est %>% filter (Stocknumber==i) %>% select ("Name") %>% distinct()
    
    R <- SRDat %>% filter (Stocknumber==i) %>% select(Rec) 
    S <- SRDat %>% filter (Stocknumber==i) %>% select(Sp) 
    # what is the scale of Ricker b estimate?
    Sc <- SRDat %>% filter (Stocknumber==i) %>% select(Scale) %>% distinct() %>% as.numeric()
    if(name$Name != "Skagit") plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(R$Rec) ) )
    if(name$Name == "Skagit") plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(0,max(R$Rec) ) )
    
    a <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
      summarise(A=exp(Estimate)) %>% as.numeric()
    # Divide b by scale
    b <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logB") %>% 
      summarise(B=exp(Estimate)/Sc) %>% as.numeric()
    
    
    if(mod!="IWAM_FixedSep_RicStd" & mod!="Liermann" & mod!="Liermann_SepRicA" & mod!="Liermann_HalfNormRicVar"& mod!="Liermann_HalfCauchyRicVar" & mod!="Liermann_HalfNormRicVar_NormDeltaSig" & mod!="Liermann_HalfNormRicVar_FixedDelta"){
      if (i %in% stksNum_surv) {  #stocknumber 0 and either 22 or 23, depending on if Skagit is removed
        surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name==name$Name)
        #if(i==0) surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name==name$Name) 
        #if(i==22|i==23)surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name=="Cowichan") 
        #if(i==22|i==23) { surv.dat <- surv.dat %>% filter(Yr >= 1985 & Yr !=1986 & Yr != 1987) }
        if(name$Name == "Cowichan") { surv.dat <- surv.dat %>% filter(Yr >= 1985 & Yr !=1986 & Yr != 1987) }
        mean.log.surv <- surv.dat %>% summarize(mean = mean(log(Surv)))
        gamma <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="gamma") %>% select(Estimate)
        a <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
          summarise(A=exp(Estimate + gamma$Estimate*mean.log.surv$mean)) %>% as.numeric()
      }
      
    }
    #Parken values for skagit
    skagit_alpha <- 7.74
    skagit_beta <- 0.0000657
    RR_skagit <- NA
    SS <- RR<- NA
    #RR_std <- NA
    
    if(removeSkagit==FALSE){
      for (j in 1:100){
        if (i!=22) SS[j] <- j*(max(S$Sp)/100)
        if (i==22) SS[j] <- j*(max(S$Sp*3)/100)
        RR[j] <- a * SS[j] * exp(-b * SS[j])
        if(i==22) {RR_skagit[j] <- skagit_alpha * SS[j] * exp(-skagit_beta * SS[j])}
        #if (i %in% stks_ar) {RR_std[j] <- A_std$A[which(A_std$Stocknumber==i)] * SS[j] *  exp(-B_std$B[which(B_std$Stocknumber==i)] * SS[j])}
      }
    } 
    if(removeSkagit==TRUE){
      for (j in 1:100){
        SS[j] <- j*(max(S$Sp)/100)
        RR[j] <- a * SS[j] * exp(-b * SS[j])
      }
    } 
    
    if (i %not in% c(stksNum_ar, stksNum_surv)) col.use <- "black"
    if (i %in% stksNum_ar) col.use <- "red"
    if (i %in% stksNum_surv) col.use <- "blue"
    if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") col.use <- "black"
    lines(x=SS, y=RR, col=col.use) 
    
    #For Skagit, add Parken et al. 2006 model curve
    if(removeSkagit==FALSE){if(i==22) lines(x=SS, y=RR_skagit, lty="dashed")}
    
    mtext(name$Name, side=3, cex=0.8)
    
    # Plot SMSYs (black for std, red for AR(1), and dashed for Parken et al. 2006)
    smsy <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY = Estimate * Sc) %>% as.numeric()
    smsy_ul <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY_ul = Estimate * Sc + 1.96 * Std..Error * Sc ) %>% as.numeric()
    smsy_ll <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY_ul = Estimate * Sc - 1.96 * Std..Error * Sc ) %>% as.numeric()
    
    
    abline(v = smsy, col=col.use)

    if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="Ricker_AllMod"){
      if (i %in% stksNum_ar) polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(1,0,0, alpha=0.1), border=NA ) 
      if (i %in% stksNum_surv) polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(0,0,1, alpha=0.1), border=NA ) 
      if (i %not in% c(stksNum_ar, stksNum_surv))  polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
    }
    
    if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta")  polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
    #else polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(0,0,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
    
    SMSY_std <- SMSY_std %>% right_join(names) %>% filter(Name==name$Name)#filter(Stocknumber != 22)
    
    if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="Ricker_AllMod"){
      if(i %in% stksNum_ar) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*Scale.stock[i+1] , col="black")
      if(i %in% stksNum_surv) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*Scale.stock[i+1] , col="black")
    }
    
    ParkenSMSY <- read.csv("DataIn/ParkenSMSY.csv")
    #if (removeSkagit==TRUE) ParkenSMSY <- ParkenSMSY %>% filter(Name != "Skagit")
    ParkenSMSY <- ParkenSMSY %>% filter(Name==as.character(name$Name)) %>% select (SMSY) %>% as.numeric()
    abline(v=ParkenSMSY, lty="dashed")
    if(is.data.frame(r2)==TRUE) {
      lab <-  r2 %>% filter(Stocknumber==i) %>% select(r2) %>% as.numeric() %>% round(2)
      legend("topright", legend = "", title= paste0("r2=",lab), bty="n")
    }
  }
  
}


# Plot SR linearized model

PlotSRLinear <- function(SRDat, All_Est, SMSY_std, stksNum_ar, stksNum_surv, r2, removeSkagit) {
  Stks <- unique(SRDat$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)
  
  for (i in Stks){
    R <- SRDat %>% filter (Stocknumber==i) %>% select(Rec) 
    S <- SRDat %>% filter (Stocknumber==i) %>% select(Sp) 
    LogRS <- log(R$Rec/S$Sp)
    
    # what is the scale of Ricker b estimate?
    Sc <- SRDat %>% filter (Stocknumber==i) %>% select(Scale) %>% distinct() %>% as.numeric()
    plot(x=S$Sp, y=LogRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(LogRS) ) )
    #if(i !=22 & i!=0) plot(x=S$Sp, y=LogRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(LogRS) ) )
    #if(i ==22) plot(x=S$Sp, y=LogRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(0,max(LogRS) ) )
    #if(i ==0) plot(x=S$Sp, y=LogRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(min(LogRS),max(LogRS) ) )
    
    name <- All_Est %>% filter (Stocknumber==i) %>% select ("Name") %>% distinct()
    mtext(name$Name, side=3)
    
     
    LogA <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
      summarise(A=Estimate) %>% as.numeric()
    # Divide b by scale
    B <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logB") %>% 
      summarise(B=exp(Estimate)/Sc) %>% as.numeric()

    if(mod!="IWAM_FixedSep_RicStd" & mod!="Liermann" & mod!="Liermann_SepRicA" & mod!="Liermann_HalfNormRicVar" & mod!="Liermann_HalfCauchyRicVar" & mod!="Liermann_HalfNormRicVar_NormDeltaSig" & mod!="Liermann_HalfNormRicVar_FixedDelta"){
      if (i %in% stksNum_surv) {  #stocknumber 0 and either 22 or 23, depending on if Skagit is removed
        
        if(i==0) surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name=="Harrison") 
        if(i==22|i==23)surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name=="Cowichan") 
        if(i==22|i==23) { surv.dat <- surv.dat %>% filter(Yr >= 1985 & Yr !=1986 & Yr != 1987) }
        mean.log.surv <- surv.dat %>% summarize(mean = mean(log(Surv)))
        gamma <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="gamma") %>% select(Estimate)
        LogA <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
          summarise(LogA.adj=(Estimate + gamma$Estimate*mean.log.surv$mean)) %>% as.numeric()
      }
    }
    
    if (i %in% stksNum_ar) col.use <- "red"
    if (i %in% stksNum_surv) col.use <- "blue"
    if (i %not in% c(stksNum_ar, stksNum_surv)) col.use <- "black"
    if (mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") col.use <- "black"
    
    abline(a=LogA, b=-B, col=col.use)
    
    # For Skagit, add the best fit based on Parken et al. 2006 Ricker estimates
    skagit_alpha <- 7.74
    skagit_beta <- 0.0000657
    if(removeSkagit==FALSE) {if (i==22)     abline(a=log(skagit_alpha), b=-skagit_beta, col="black", lty="dashed")}

    lab <-  r2 %>% filter(Stocknumber==i) %>% select(r2) %>% as.numeric() %>% round(2)
    legend("topright", legend = "", title= paste0("r2=",lab), bty="n")

  }
  
  
}

#-------------------------------------------------------------------------------------------------
#Plot standardized residuals

PlotStdResid <- function(SRes){
  
  Stks <- unique(SRes$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)
  
  for (i in Stks){
    SR <- SRes %>% filter (Stocknumber==i) %>% select(Res) 
    plot(x=1:length(SR$Res), y=SR$Res, xlab="", ylab="", pch=20, ylim=c(min(-4, min(SR$Res)),max(4, max(SR$Res)) ) ) 
    abline(h=c(3,-3), col="red")
    abline(h=0, col="black")
    name <- SRes %>% filter (Stocknumber==i) %>% select ("Name") %>% distinct()
    mtext(name$Name, side=3)
  }
}

Plotacf <- function(Preds){
  Stks <- unique(Preds$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(2, 2, 1.5, 1) + 0.1)
  for (i in Stks){
    Res <- Preds %>% filter (Stocknumber==i) %>% select(Res)
    acf(Res$Res, plot=T)
    acf1 <- acf(Res$Res, plot=F)$acf[2]
    len <- length(Res$Res)
    acf1.ci <- qnorm((1 + 0.95)/2)/sqrt(len)
    if (abs(acf1)>acf1.ci) col="red" else col="black"
    name <- Preds %>% filter (Stocknumber==i) %>% select ("Name") %>% distinct()
    mtext(name$Name, side=3, col=col)
  }
  
}
# What are the upper and lower bounds of KSR, Stikine and Cowichan?

# KSR.SMSY <- All_Est %>% filter(Name=="KSR") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
# KSR.SMSY.ul <- KSR.SMSY + 1.96 * (All_Est %>% filter(Name=="KSR") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# KSR.SMSY.ul <- KSR.SMSY.ul * SRDat %>% filter (Name=="KSR") %>% select(Scale) %>% distinct() %>% as.numeric()
# 
# Stikine.SMSY <- All_Est %>% filter(Name=="Stikine") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
# Stikine.SMSY.ul <- Stikine.SMSY + 1.96 * (All_Est %>% filter(Name=="Stikine") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# Stikine.SMSY.ul <- Stikine.SMSY.ul * SRDat %>% filter (Name=="Stikine") %>% select(Scale) %>% distinct() %>% as.numeric()
# 
# Cow.SMSY <- All_Est %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
# Cow.SMSY.ul <- Cow.SMSY + 1.96 * (All_Est %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# Cow.SMSY.ul <- Cow.SMSY.ul * SRDat %>% filter (Name=="Cowichan") %>% select(Scale) %>% distinct() %>% as.numeric()
# Cow.SMSY.ll <- Cow.SMSY - 1.96 * (All_Est %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# Cow.SMSY.ll <- Cow.SMSY.ll * SRDat %>% filter (Name=="Cowichan") %>% select(Scale) %>% distinct() %>% as.numeric()


#------------------------------------------------------------------
# Plot  WA regression

plotWAregressionSMSY <- function (All_Est, All_Deltas, SRDat, Stream, WA,  PredlnSMSY=NA, PredlnWA, title1, mod) {

  SMSY <- All_Est %>% filter(Param=="SMSY") %>% mutate(ModelOrder=0:(length(unique(All_Est$Stocknumber))-1))
  # what is scale of SMSY?
  Sc <- SRDat %>% select(Stocknumber, Scale) %>% distinct()
  SMSY <- SMSY %>% left_join(Sc, by="Stocknumber") %>% mutate(rawSMSY=Estimate*Scale)
  #if (mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_Constm"|mod=="IWAM_FixedSep_Constyi") SMSY <- SMSY %>% left_join(Stream, by=c("Stocknumber","ModelOrder"))
  #if (mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar") SMSY <- SMSY %>% left_join(Stream, by=c("Stocknumber"))
  lnSMSY <- log(SMSY$rawSMSY)
  lnWA <- log(WA$WA)
  
  par(cex=1.5)
  col.use <- NA
  for(i in 1:length(SMSY$lh)) {if (SMSY$lh[i]==0) col.use[i] <- "forestgreen" else col.use[i] <- "dodgerblue3"}
  plot(y=lnSMSY, x=lnWA, pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SMSY)")
  points(y=lnSMSY, x=lnWA, pch=20, col=col.use, cex=1.5)
  #points(y=lnSMSY[18:length(SMSY$lh)], x=lnWA[18:length(SMSY$lh)], pch=3, col=col.use[18:length(SMSY$lh)], cex=1.5)
  logD1 <- All_Deltas %>% filter(Param=="logDelta1") %>% select(Estimate) %>% pull()
  logD2 <- All_Deltas %>% filter(Param=="logDelta2") %>% select(Estimate) %>% pull()
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="IWAM_FixedSep_Constm"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") {
    logD1o <- All_Deltas %>% filter(Param=="logDelta1ocean") %>% select(Estimate) %>% pull() + logD1}
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="IWAM_FixedSep_Constyi"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") {
    D2o <- exp(All_Deltas %>% filter(Param=="logDelta2ocean") %>% select(Estimate) %>% pull() ) + exp(logD2)
    if(nrow(All_Deltas %>% filter(Param=="Delta2ocean"))>=1)  D2o <- (All_Deltas %>% filter(Param=="Delta2ocean") %>% select(Estimate) %>% pull() ) + exp(logD2)
  }
  #if (mod=="IWAM_FixedCombined") D2 <- All_Deltas %>% filter(Param=="Delta2_bounded") %>% select(Estimate) %>% pull()
  if(length(logD1)==1&length(logD2)==2){
    abline(a=logD1[1], b=exp(logD2[1]), col="forestgreen", lwd=2)
    abline(a=logD1[1], b=exp(logD2[2]), col="dodgerblue3", lwd=2)
  }
  if(length(logD1)==2&length(logD2)==2){
    abline(a=logD1[1], b=exp(logD2[1]), col="forestgreen", lwd=2)
    abline(a=logD1[2], b=exp(logD2[2]), col="dodgerblue3", lwd=2)
  }
  if(mod=="IWAM_FixedCombined") abline(a=logD1, b=exp(logD2), col="maroon", lwd=2)#Actually pulls Delta2_bounded, so no need to log
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") {
    abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
    abline(a=logD1o, b=D2o, col="dodgerblue3", lwd=2)
  }
  if(mod=="IWAM_FixedSep_Constm") {
    abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
    abline(a=logD1o, b=exp(logD2), col="dodgerblue3", lwd=2)
  }
  if(mod=="IWAM_FixedSep_Constyi") {
    abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
    abline(a=logD1, b=D2o, col="dodgerblue3", lwd=2)
  }
  
  
  if(exists("PredlnSMSY")){
    PredlnSMSY <- PredlnSMSY %>% mutate (up = Estimate + 1.96 * Std..Error, lo=Estimate - 1.96*Std..Error) 
    #up_S <- PredlnSMSY %>% filter(Param== "PredlnSMSY_S") %>% select(up) %>% pull()
    #lo_S <- PredlnSMSY %>% filter(Param== "PredlnSMSY_S") %>% select(lo) %>% pull()
    #up_O <- PredlnSMSY %>% filter(Param== "PredlnSMSY_O") %>% select(up) %>% pull()
    #lo_O <- PredlnSMSY %>% filter(Param== "PredlnSMSY_O") %>% select(lo) %>% pull()
    up_S <- PredlnSMSY %>% filter(Param== "PredlnSMSYs_CI") %>% select(up) %>% pull()
    lo_S <- PredlnSMSY %>% filter(Param== "PredlnSMSYs_CI") %>% select(lo) %>% pull()
    up_O <- PredlnSMSY %>% filter(Param== "PredlnSMSYo_CI") %>% select(up) %>% pull()
    lo_O <- PredlnSMSY %>% filter(Param== "PredlnSMSYo_CI") %>% select(lo) %>% pull()
    up <- PredlnSMSY %>% filter(Param== "PredlnSMSY_CI") %>% select(up) %>% pull()
    lo <- PredlnSMSY %>% filter(Param== "PredlnSMSY_CI") %>% select(lo) %>% pull()
    if(is.na(up_S[1])==FALSE) polygon(x=c(PredlnWA, rev(PredlnWA)), y=c(up_S, rev(lo_S)), col=rgb(0,0.4,0, alpha=0.2), border=NA)
    if(is.na(up_O[1])==FALSE) polygon(x=c(PredlnWA, rev(PredlnWA)), y=c(up_O, rev(lo_O)), col=rgb(0,0.2,0.4, alpha=0.2), border=NA)
    if(is.na(up[1])==FALSE) polygon(x=c(PredlnWA, rev(PredlnWA)), y=c(up, rev(lo)), col=rgb(0.6,0.2,0.4, alpha=0.2), border=NA)
  }

  if(length(logD1)==2&length(logD2)==2){
    text(x=6, y=10.5,labels= paste0("Ocean-type\nlog(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("Stream-type\nlog(Delta1)=",round(logD1[3],2), ", \nDelta2=", round(exp(logD2[2]),2)), col="forestgreen", cex=0.8)
  }
  if(length(logD1)==1&length(logD2)==2){
    text(x=6, y=10.5,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[2]),2)), col="forestgreen", cex=0.8)
  }
  if(mod=="IWAM_FixedCombined"){
    text(x=6, y=10.5,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="maroon", cex=0.8)
  
  }
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta"){
    text(x=9, y=7,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Delta1)=",round(logD1o[1],2), ", \nDelta2=", round(D2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  if(mod=="IWAM_FixedSep_Constm"){
    text(x=9, y=7,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Delta1)=",round(logD1o[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="dodgerblue3", cex=0.8)
  }
  if(mod=="IWAM_FixedSep_Constyi"){
    text(x=9, y=7,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(D2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  
  
  title(title1, cex.main=0.9)
  
}


plotWAregressionSREP <- function (All_Est, All_Deltas, SRDat, Stream, WA,  PredlnSREP=NA, PredlnWA, title1, mod) {
  
  SREP <- All_Est %>% filter(Param=="SREP") %>% mutate(ModelOrder=0:(length(unique(All_Est$Stocknumber))-1))
  # what is scale of SREP?
  Sc <- SRDat %>% select(Stocknumber, Scale) %>% distinct()
  SREP <- SREP %>% left_join(Sc) %>% mutate(rawSREP=Estimate*Scale)
  #if (mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_Constm"|mod=="IWAM_FixedSep_Constyi") SREP <- SREP %>% left_join(Stream, by=c("Stocknumber","ModelOrder"))
  #if (mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar") SREP <- SREP %>% left_join(Stream, by=c("Stocknumber"))
  lnSREP <- log(SREP$rawSREP)
  lnWA <- log(WA$WA)
  
  par(cex=1.5)
  col.use <- NA
  for(i in 1:length(SREP$lh)) {if (SREP$lh[i]==0) col.use[i] <- "forestgreen" else col.use[i] <- "dodgerblue3"}
  plot(y=lnSREP, x=lnWA, pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SREP)")
  points(y=lnSREP, x=lnWA, pch=20, col=col.use, cex=1.5)
  #points(y=lnSREP[18:length(SREP$lh)], x=lnWA[18:length(SREP$lh)], pch=3, col=col.use[18:length(SREP$lh)], cex=1.5)
  logN1 <- All_Deltas %>% filter(Param=="logNu1") %>% select(Estimate) %>% pull()
  logN2 <- All_Deltas %>% filter(Param=="logNu2") %>% select(Estimate) %>% pull()
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="IWAM_FixedSep_Constm"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") logN1o <- All_Deltas %>% filter(Param=="logNu1ocean") %>% select(Estimate) %>% pull() + logN1
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="IWAM_FixedSep_Constyi"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") {
    N2o <- exp(All_Deltas %>% filter(Param=="logNu2ocean") %>% select(Estimate) %>% pull() ) + exp(logN2)
    if(nrow(All_Deltas %>% filter(Param=="Nu2ocean"))>=1)  N2o <- (All_Deltas %>% filter(Param=="Nu2ocean") %>% select(Estimate) %>% pull() ) + exp(logN2)
  }
  #if (mod=="IWAM_FixedCombined") D2 <- All_Deltas %>% filter(Param=="Delta2_bounded") %>% select(Estimate) %>% pull()
  if(length(logN1)==1&length(logN2)==2){
    abline(a=logN1[1], b=exp(logN2[1]), col="forestgreen", lwd=2)
    abline(a=logN1[1], b=exp(logN2[2]), col="dodgerblue3", lwd=2)
  }
  if(length(logN1)==2&length(logN2)==2){
    abline(a=logN1[1], b=exp(logN2[1]), col="forestgreen", lwd=2)
    abline(a=logN1[2], b=exp(logN2[2]), col="dodgerblue3", lwd=2)
  }
  if(mod=="IWAM_FixedCombined") abline(a=logN1, b=exp(logN2), col="maroon", lwd=2)#Actually pulls Delta2_bounded, so no need to log
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") {
    abline(a=logN1, b=exp(logN2), col="forestgreen", lwd=2)
    abline(a=logN1o, b=N2o, col="dodgerblue3", lwd=2)
  }
  if(mod=="IWAM_FixedSep_Constm") {
    abline(a=logN1, b=exp(logN2), col="forestgreen", lwd=2)
    abline(a=logN1o, b=exp(logN2), col="dodgerblue3", lwd=2)
  }
  if(mod=="IWAM_FixedSep_Constyi") {
    abline(a=logN1, b=exp(logN2), col="forestgreen", lwd=2)
    abline(a=logN1, b=N2o, col="dodgerblue3", lwd=2)
  }
  
  
  if(exists("PredlnSREP")){
    PredlnSREP <- PredlnSREP %>% mutate (up = Estimate + 1.96 * Std..Error, lo=Estimate - 1.96*Std..Error) 
    #up_S <- PredlnSREP %>% filter(Param== "PredlnSREP_S") %>% select(up) %>% pull()
    #lo_S <- PredlnSREP %>% filter(Param== "PredlnSREP_S") %>% select(lo) %>% pull()
    #up_O <- PredlnSREP %>% filter(Param== "PredlnSREP_O") %>% select(up) %>% pull()
    #lo_O <- PredlnSREP %>% filter(Param== "PredlnSREP_O") %>% select(lo) %>% pull()
    up_S <- PredlnSREP %>% filter(Param== "PredlnSREPs_CI") %>% select(up) %>% pull()
    lo_S <- PredlnSREP %>% filter(Param== "PredlnSREPs_CI") %>% select(lo) %>% pull()
    up_O <- PredlnSREP %>% filter(Param== "PredlnSREPo_CI") %>% select(up) %>% pull()
    lo_O <- PredlnSREP %>% filter(Param== "PredlnSREPo_CI") %>% select(lo) %>% pull()
    up <- PredlnSREP %>% filter(Param== "PredlnSREP_CI") %>% select(up) %>% pull()
    lo <- PredlnSREP %>% filter(Param== "PredlnSREP_CI") %>% select(lo) %>% pull()
    if(is.na(up_S[1])==FALSE) polygon(x=c(PredlnWA, rev(PredlnWA)), y=c(up_S, rev(lo_S)), col=rgb(0,0.4,0, alpha=0.2), border=NA)
    if(is.na(up_O[1])==FALSE) polygon(x=c(PredlnWA, rev(PredlnWA)), y=c(up_O, rev(lo_O)), col=rgb(0,0.2,0.4, alpha=0.2), border=NA)
    if(is.na(up[1])==FALSE) polygon(x=c(PredlnWA, rev(PredlnWA)), y=c(up, rev(lo)), col=rgb(0.6,0.2,0.4, alpha=0.2), border=NA)
  }
  
  if(length(logN1)==2&length(logN2)==2){
    text(x=6, y=10.5,labels= paste0("Ocean-type\nlog(Nu1)=",round(logN1[1],2), ", \nNu2=", round(exp(logN2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("Stream-type\nlog(Nu1)=",round(logN1[3],2), ", \nNu2=", round(exp(logN2[2]),2)), col="forestgreen", cex=0.8)
  }
  if(length(logN1)==1&length(logN2)==2){
    text(x=6, y=10.5,labels= paste0("log(Nu1)=",round(logN1[1],2), ", \nNu2=", round(exp(logN2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("log(Nu1)=",round(logN1[1],2), ", \nNu2=", round(exp(logN2[2]),2)), col="forestgreen", cex=0.8)
  }
  if(mod=="IWAM_FixedCombined"){
    text(x=6, y=10.5,labels= paste0("log(Nu1)=",round(logN1[1],2), ", \nNu2=", round(exp(logN2[1]),2)), col="maroon", cex=0.8)
    
  }
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_SepRicA"|mod=="Liermann_HalfNormRicVar"|mod=="Liermann_HalfCauchyRicVar"|mod=="Liermann_HalfNormRicVar_NormDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta"){
    text(x=9, y=7,labels= paste0("log(Nu1)=",round(logN1[1],2), ", \nNu2=", round(exp(logN2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=10.5,labels= paste0("log(Nu1)=",round(logN1o[1],2), ", \nNu2=", round(N2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  if(mod=="IWAM_FixedSep_Constm"){
    text(x=9, y=7,labels= paste0("log(Nu1)=",round(logN1[1],2), ", \nNu2=", round(exp(logN2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Nu1)=",round(logN1o[1],2), ", \nNu2=", round(exp(logN2[1]),2)), col="dodgerblue3", cex=0.8)
  }
  if(mod=="IWAM_FixedSep_Constyi"){
    text(x=9, y=7,labels= paste0("log(Nu1)=",round(logN1[1],2), ", \nNu2=", round(exp(logN2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Nu1)=",round(logN1[1],2), ", \nNu2=", round(N2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  
  
  title(title1, cex.main=0.9)
  
}


plotWAregression_Parken <- function(data, All_Deltas){
  par(cex=1.5)
  col.use <- NA
  Parken_data <- as.data.frame(read.csv("DataIn/WA_Parken.csv"))
  for(i in 1:length(Parken_data$lh)) {if (Parken_data$lh[i]==0) col.use[i] <- "forestgreen" else col.use[i] <- "dodgerblue3"}
  plot(x=log(data$WA), y=log(data$SMSY*data$Scale),pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SMSY)")
  points(x=log(data$WA), y=log(data$SMSY*data$Scale), pch=20, col=col.use, cex=1.5)
  logD1 <- All_Deltas %>% filter(Param=="logDelta1") %>% select(Estimate) %>% pull()
  D2 <- All_Deltas %>% filter(Param=="Delta2_bounded") %>% select(Estimate) %>% pull()
  abline(a=logD1, b=D2, col="maroon", lwd=2)
  
  #abline(lm(log(data$SMSY*data$Scale) ~ log(data$WA)))
  #text(x=9, y=8,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
  text(x=6, y=10.5,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(D2[1],2)), col="maroon", cex=0.8)
  
  title("Parken Data: Combined")
  
}

plotWAregression_ParkenSep <- function(data, All_Deltas){
  par(cex=1.5)
  col.use <- NA
  Parken_data <- as.data.frame(read.csv("DataIn/WA_Parken.csv"))
  for(i in 1:length(Parken_data$lh)) {if (Parken_data$lh[i]==0) col.use[i] <- "forestgreen" else col.use[i] <- "dodgerblue3"}
  plot(x=log(data$WA), y=log(data$SMSY*data$Scale),pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SMSY)")
  points(x=log(data$WA), y=log(data$SMSY*data$Scale), pch=20, col=col.use, cex=1.5)
  logD1 <- All_Deltas %>% filter(Param=="logDelta1") %>% select(Estimate) %>% pull()
  logD1o <- All_Deltas %>% filter(Param=="logDelta1ocean") %>% select(Estimate) %>% pull() + logD1
  logD2 <- All_Deltas %>% filter(Param=="logDelta2") %>% select(Estimate) %>% pull()
  D2o <- exp(All_Deltas %>% filter(Param=="logDelta2ocean") %>% select(Estimate) %>% pull() ) + exp(logD2)
  abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
  abline(a=logD1o, b=D2o, col="dodgerblue3", lwd=2)
  
  #abline(lm(log(data$SMSY*data$Scale) ~ log(data$WA)))
  text(x=9, y=7,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
  text(x=6, y=10.5,labels= paste0("log(Delta1)=",round(logD1o[1],2), ", \nDelta2=", round(D2o[1],2)), col="dodgerblue3", cex=0.8)
  
  title("Parken Data: Separated")
  
}



plotRicA <- function (){#All_Est_Liermann, All_Est_Ricker_AllMod, All_Est_Liermann_SepRicA=NA){
  #par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
  All_Est_Liermann <- readRDS("DataOut/All_Est_Liermann.RDS")
  All_Est_Ricker_AllMod <- readRDS("DataOut/All_Est_Ricker_AllMod.RDS")
  All_Est_Ricker_std <- readRDS("DataOut/All_Est_Ricker_std.RDS")
  All_Est_Liermann_SepRicA <- readRDS("DataOut/All_Est_Liermann_SepRicA.RDS")
  All_Est_Liermann_HalfNormRicVar <- readRDS("DataOut/All_Est_Liermann_HalfNormRicVar.RDS")
  All_Est_Liermann_HalfCauchyRicVar <- readRDS("DataOut/All_Est_Liermann_HalfCauchyRicVar.RDS")
  All_Est_Liermann_SepRicA_uniformSigmaPrior <- readRDS("DataOut/All_Est_Liermann_SepRicA_uniformSigmaPrior.RDS")
  All_Est_Liermann_SepRicA_noSigmaPrior <- readRDS("DataOut/All_Est_Liermann_SepRicA_noSigmaPrior.RDS")
  All_Est_Liermann_SepRicA_invGamma0.001 <- readRDS("DataOut/All_Est_Liermann_SepRicA_invGamma0.001.RDS")
  All_Est_Liermann_SepRicA_invGamma0.1 <- readRDS("DataOut/All_Est_Liermann_SepRicA_invGamma0.1.RDS")
  All_Est_Liermann_SepRicA_invGamma0.01_invGammaA0.001 <- readRDS("DataOut/All_Est_Liermann_SepRicA_invGamma0.01_invGammaA0.001.RDS")
  All_Est_Liermann_SepRicA_invGamma0.001_invGammaA0.01 <- readRDS("DataOut/All_Est_Liermann_SepRicA_invGamma0.001_invGammaA0.01.RDS")
  
  RicAL <- All_Est_Liermann %>% filter (Param=="logA")
  RicALsep <- All_Est_Liermann_SepRicA %>% filter (Param=="logA")
  RicALhNRV <- All_Est_Liermann_HalfNormRicVar %>% filter (Param=="logA")
  RicALhCRV <- All_Est_Liermann_HalfCauchyRicVar %>% filter (Param=="logA")
  RicALsep_uni <- All_Est_Liermann_SepRicA_uniformSigmaPrior %>% filter (Param=="logA")
  RicALsep_0.1 <- All_Est_Liermann_SepRicA_invGamma0.1 %>% filter (Param=="logA")
  RicALsep_0.001 <- All_Est_Liermann_SepRicA_invGamma0.001 %>% filter (Param=="logA")
  RicALsep_0.001_0.01 <- All_Est_Liermann_SepRicA_invGamma0.001_invGammaA0.01 %>% filter (Param=="logA")
  RicALsep_0.01_0.001 <- All_Est_Liermann_SepRicA_invGamma0.01_invGammaA0.001 %>% filter (Param=="logA")
  RicALsep_none <- All_Est_Liermann_SepRicA_noSigmaPrior %>% filter (Param=="logA")
  RicAR <- All_Est_Ricker_AllMod %>% filter (Param=="logA")
  RicARstd <- All_Est_Ricker_std %>% filter (Param=="logA")#(Param=="logA_std")
  nRicAL <- length(RicALsep$Estimate)
  nRicAstd <- length(RicARstd$Estimate)

  
  #box.data <- data.frame(RicLogA=c(RicAL$Estimate,RicALsep$Estimate, RicAR$Estimate), Model=c(rep("Hierarchical_logA_combined",nRicA),rep("Hierarchical_logA_sep",nRicA), rep("Fixed_logA",nRicA)))
  box.data <- data.frame(RicLogA = c(RicARstd$Estimate, RicALsep_0.1$Estimate, RicALsep$Estimate, RicALsep_0.001$Estimate, 
                                   RicALhNRV$Estimate, RicALhCRV$Estimate, RicALsep_uni$Estimate, RicALsep_none$Estimate), 
                         Model = c(rep("Fixed\neffects",nRicAstd), rep("InvGamma\n(0.1, 0.1)", nRicAL), 
                                 rep("InvGamma\n(0.01,0.01)",nRicAL), rep("Inv\nGamma\n(0.001, 0.001)", nRicAL),
                                 rep("Half\nNormal\n(0,1)",nRicAL), 
                                 rep("Half\nCauchy\n(0,1)",nRicAL),  
                                 rep("Uniform\n(0,2)",nRicAL), rep("No priors\nSigmas",nRicAL)), 
                         lh = c(RicARstd$lh, RicALsep_0.1$lh, RicALsep$lh, RicALsep_0.001$lh, 
                                      RicALhNRV$lh, RicALhCRV$lh, RicALsep_uni$lh, RicALsep_none$lh) )
  
  order <- data.frame(Model=unique(box.data$Model), order=1:8)
  box.data <- box.data %>% left_join(order, by="Model")
  #Add reorder to as.factor(Model), I think, but need to specify order as above 1:8. Add column to dataframe to do this, of length 200.
  # see my old dplyr code on adding new columns with string of numbers aligned with a factor (IWAM.r)
  cols<-viridis(4, alpha=0.9)
  cols.order <- c(grey(0.5), cols[3], cols[2], t_col(color=cols[1], percent=70), t_col(cols[1], percent=50), t_col(cols[1], percent=30), "white", cols[4])
  #print cols.order and then cut and past into line: scale_fill_manual below. These are the colours I used in the plot of priors. can add a legend to quickly see how ggplot reorders this befor plotting
  
  ggplot(box.data, aes(x=reorder(as.factor(Model),order), y=RicLogA, fill=as.factor(Model))) +  
    geom_boxplot() + 
    scale_fill_manual(values=c("#808080", "#35B779E6", "#31688EE6", "#4401544C", "#4401547F", "#440154B2", "white", "#FDE725E6" )) +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) + 
    geom_jitter(color="black", shape=as.factor(box.data$lh), size=2, alpha=0.9) + #,aes(shape=lh)?. I like shape 16,17
    scale_shape_manual(values=c(16,17)) + #Does not work 
    ggtitle("Distribution of Ricker LogA") +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=20,face="bold"),
          plot.title = element_text(size = 20)) + 
    xlab("Model")
  #boxplot(list(Hierarchical=RicAL$Estimate, Fixed=RicAR$Estimate), ylab="Ricker log(a)", outline=FALSE, col="forestgreen")
  #points(x=jitter(rep(1,nRicA)), y=RicAL$Estimate, pch=20)
  #points(x=jitter(rep(2,nRicA)), y=RicAR$Estimate, pch=20)
}

#png(paste("DataOut/RicADist_ComparePriors.png", sep=""), width=7, height=7, units="in", res=500)
#plotRicA()
#dev.off()
