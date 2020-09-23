initial<-list()
initial$mean.a<-1
initial$var.a<-1
initial$b <- 2
initial$ln.sig.e<-1
initial$ln.sig.w<-1
initial$Ts<-0
initial$EstB<-TRUE

data<-read.csv("data/Stellako.csv")
x <- data$ETS
y <- log(data$Rec/data$ETS)


sdKFa <- sd(kf.rw(initial=initial, x=data$ETS, y=log(data$Rec/data$ETS))$smoothe.mean.a, na.rm=T)
kfAICc <- kf.rw(initial=initial, x=data$ETS, y=log(data$Rec/data$ETS))$AICc # uses concentrated nLL
kfAICc


linearnLL <- -logLik( lm ( log (data$Rec/data$ETS) ~ data$ETS ))[1] #full nLL
rss <- sum(lm ( log (data$Rec/data$ETS) ~ data$ETS )$resid^2)
sigEst <- sd(lm ( log (data$Rec/data$ETS) ~ data$ETS )$resid)
linearnLLconc <- - ( -(N/2)*log(sigEst^2) - (1/(2*sigEst^2)) * rss ) #Concentrated nLL
lin_param <- 3
N <- length(data$ETS) #assuming complete data
linAICc <- 2 * linearnLLconc + 2 * lin_param * (N/(N - lin_param - 1))
linAICc

if( linAICc > (kfAIC-2) ) use <- 1

