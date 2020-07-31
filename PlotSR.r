# Plot SR models: curves and linearized mddel

# Plot SR curves

PlotSRCurve <- function(SRDat, All_Est, SMSY_std, stksNum_ar, stksNum_surv, stks_surv, r2, removeSkagit) {
  Stks <- unique(SRDat$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)
  
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
    lines(x=SS, y=RR, col=col.use) 
    
    #For Skagit, add Parken et al. 2006 model curve
    if(removeSkagit==FALSE){if(i==22) lines(x=SS, y=RR_skagit, lty="dashed")}
    
    mtext(name$Name, side=3)
    
    # Plot SMSYs (black for std, red for AR(1), and dashed for Parken et al. 2006)
    smsy <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY = Estimate * Sc) %>% as.numeric()
    smsy_ul <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY_ul = Estimate * Sc + 1.96 * Std..Error * Sc ) %>% as.numeric()
    smsy_ll <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY_ul = Estimate * Sc - 1.96 * Std..Error * Sc ) %>% as.numeric()
    
    
    abline(v = smsy, col=col.use)

    if (i %in% stksNum_ar) polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(1,0,0, alpha=0.1), border=NA ) 
    if (i %in% stksNum_surv) polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(0,0,1, alpha=0.1), border=NA ) 
    if (i %not in% c(stksNum_ar, stksNum_surv))  polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
    #else polygon(x=c(smsy_ul, smsy_ll, smsy_ll, smsy_ul), y=c(0,0,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
    
    SMSY_std <- SMSY_std %>% right_join(names) %>% filter(Name==name$Name)#filter(Stocknumber != 22)
    
    if(i %in% stksNum_ar) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*Scale.stock[i+1] , col="black")
    if(i %in% stksNum_surv) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*Scale.stock[i+1] , col="black")
    
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
    if(i !=22 & i!=0) plot(x=S$Sp, y=LogRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(LogRS) ) )
    if(i ==22) plot(x=S$Sp, y=LogRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(0,max(LogRS) ) )
    if(i ==0) plot(x=S$Sp, y=LogRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(min(LogRS),max(LogRS) ) )
    
    name <- All_Est %>% filter (Stocknumber==i) %>% select ("Name") %>% distinct()
    mtext(name$Name, side=3)
    
     
    LogA <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
      summarise(A=Estimate) %>% as.numeric()
    # Divide b by scale
    B <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logB") %>% 
      summarise(B=exp(Estimate)/Sc) %>% as.numeric()

    if (i %in% stksNum_surv) {  #stocknumber 0 and either 22 or 23, depending on if Skagit is removed
      
      if(i==0) surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name=="Harrison") 
      if(i==22|i==23)surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name=="Cowichan") 
      if(i==22|i==23) { surv.dat <- surv.dat %>% filter(Yr >= 1985 & Yr !=1986 & Yr != 1987) }
      mean.log.surv <- surv.dat %>% summarize(mean = mean(log(Surv)))
      gamma <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="gamma") %>% select(Estimate)
      LogA <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
        summarise(LogA.adj=(Estimate + gamma$Estimate*mean.log.surv$mean)) %>% as.numeric()
    }
    
    if (i %in% stksNum_ar) col.use <- "red"
    if (i %in% stksNum_surv) col.use <- "blue"
    if (i %not in% c(stksNum_ar, stksNum_surv)) col.use <- "black"
    
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
  par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)
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
