#library(tidyverse)
#library(lubridate)
library(ggplot2)
library(mgcv)
library(readr)
library(data.table)
library(gridExtra)
setwd("D:/Dissertation")
plotDir <- file.path(getwd(),"Plots")

load("dengue2013April15thEpiWeek.RData")
load("dengue2013April15thEpiWeekFull.RData")

dengue.data <- data.table( dengue.data )
dengue.data.full <- data.table( dengue.data.full )

# rename columns to include "delay"
dengue.data.full[,Obs := NULL]
index <- -c(1:5)
names(dengue.data.full)[index] <- paste0("delay",names(dengue.data.full)[index])
names(dengue.data.full)[names(dengue.data.full)=="Total"] <- "y"


# only consider up to N delays
Ndelay <- 10
delay_index <- grep("delay",names(dengue.data.full))
col_index <- c(1:5,delay_index[1:(Ndelay+1)])
dat <- dengue.data.full[,..col_index]
delay_indexHere <- delay_index[1:(Ndelay+1)]
dat[,x := rowSums(.SD), .SDcols = delay_indexHere]


## Loop over nowcast periods starting at the last outbreak
plotList <- list()
jj <- 1
for(k in 159:170){

  # introduce NAs for nowcasting and forecasting
  datS <- copy(dat)
  for(i in 1:Ndelay){
    ind <- delay_indexHere[i] + 1
    datS[(k-i+1):k,ind] <- NA
  }
  # Now set everyhing to NA during forecasting
  for(i in (k+1):(nrow(datS))){
    datS[i,delay_indexHere] <- NA
  }
  # this will set x to NA at all nowcasting
  datS[,x := rowSums(.SD), .SDcols = delay_indexHere]


  # Model for X
  #Xmodel <- gam(x ~ s(Time, k=15, bs="tp") + s(EpiWeek, k=10, bs="cc"), data=datS, 
  #            family=nb(link="log"),method="REML",knots=list(EpiWeek=c(0,52)))
  #k.check(Xmodel)
  # Deterministic preds of totals X
  #Xpred <- predict(Xmodel,type="response",newdata = datS)
  # check how well it predicts Y
  #plot(y~Time,data=datS,type="l",lwd=3)
  #lines(datS$Time,Xpred,col="blue")

  # Joint model for Z,X. First need to put in long format
  Xind <- which(names(datS)=="x")
  longDat <- melt(datS,id.vars=c("DATE","EpiWeek","Time","x"),measure.vars = c(grep("delay", colnames(datS) ),Xind)  )

  ## now fit gfam GAM
  # matrix response
  longDat[,Nvariable := as.numeric(variable)]
  longDat[,omega := factor(Nvariable)]
  longDat[,y := value]
  longDat <- as.data.frame(longDat)
  y <- longDat$y
  fin <- longDat$Nvariable
  longDat$y <- cbind(y,fin)


  XZmodel <- gam(y~omega+s(Time,k=20,bs="tp",m=1) +
              s(Time,omega,bs="sz",xt=list(bs="tp"),k=5,m=1) +
              s(EpiWeek, k=15, bs="cc"), 
              #s(EpiWeek,omega,bs="sz",xt=list(bs="cc")),
              family=gfam(list(nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb)),method="REML",
              data=longDat,knots=list(EpiWeek=c(0,52))) 
  #k.check(XZmodel)

  # predict:
  predDat <- data.table(longDat)
  predDat[,y.y := NULL]
  predDat[,y.fin := NULL]
  predDat[,y := Nvariable]
  #
  predDat$preds <- predict(XZmodel,newdata=predDat,type="response")
  modelMatrix <- predict(XZmodel,newdata=predDat,type="lpmatrix")
  n.sims <- 1000
  b <- rmvn(n.sims,coef(XZmodel),XZmodel$Vp)
  LP <- tcrossprod(modelMatrix,b)
  Mean <- exp(LP)
  myThetas <- exp( XZmodel$family$getTheta() )
  Y <- apply(Mean,2,function(x){rnbinom(length(x),mu=x,size=myThetas[12])})

  predDat$preds <- apply(Y,1,median)
  predDat$lwr <- apply(Y,1,quantile,probs=0.025)
  predDat$upr <- apply(Y,1,quantile,probs=0.975)

#plot(y~Time,data=datS,type="l",lwd=3)
#lines(datS$Time,predDat[Nvariable==12,preds],col="blue",lwd=2)
#lines(datS$Time,predDat[Nvariable==12,lwr],col="blue",lwd=2,lty=2)
#lines(datS$Time,predDat[Nvariable==12,upr],col="blue",lwd=2,lty=2)
#abline(v=k)
#abline(v=k-Ndelay+1)
## ggplot
  plotHere <- ggplot(data=dat,aes(x=Time)) + theme_bw() + geom_line(aes(y=y)) + 
    geom_ribbon(data=predDat[Nvariable==12,],aes(x=Time,y=preds,ymin=lwr,ymax=upr),alpha=0.2) + 
    geom_line(data=predDat[Nvariable==12,],aes(x=Time,y=preds),col="blue") + 
    ylim(0,12000) + geom_vline(xintercept = c(k+1,k-Ndelay+1)) + 
    labs(title=paste0("Nowcast weeks: ",k-Ndelay+1,"-",k,", Forecast weeks: ",k+1,"-",nrow(datS)))
  plotList[[jj]] <- plotHere
  jj <- jj+1
}

test <- do.call(grid.arrange,plotList)
test
ggsave("joint_model.pdf",test,width=15,height=10, path=plotDir)



## Longer loop
plotList <- list()
jj <- 1
for(k in 125:170){
  
  # introduce NAs for nowcasting and forecasting
  datS <- copy(dat)
  for(i in 1:Ndelay){
    ind <- delay_indexHere[i] + 1
    datS[(k-i+1):k,ind] <- NA
  }
  # Now set everyhing to NA during forecasting
  for(i in (k+1):(nrow(datS))){
    datS[i,delay_indexHere] <- NA
  }
  # this will set x to NA at all nowcasting
  datS[,x := rowSums(.SD), .SDcols = delay_indexHere]
  
  
  # Model for X
  #Xmodel <- gam(x ~ s(Time, k=15, bs="tp") + s(EpiWeek, k=10, bs="cc"), data=datS, 
  #            family=nb(link="log"),method="REML",knots=list(EpiWeek=c(0,52)))
  #k.check(Xmodel)
  # Deterministic preds of totals X
  #Xpred <- predict(Xmodel,type="response",newdata = datS)
  # check how well it predicts Y
  #plot(y~Time,data=datS,type="l",lwd=3)
  #lines(datS$Time,Xpred,col="blue")
  
  # Joint model for Z,X. First need to put in long format
  Xind <- which(names(datS)=="x")
  longDat <- melt(datS,id.vars=c("DATE","EpiWeek","Time","x"),measure.vars = c(grep("delay", colnames(datS) ),Xind)  )
  
  ## now fit gfam GAM
  # matrix response
  longDat[,Nvariable := as.numeric(variable)]
  longDat[,omega := factor(Nvariable)]
  longDat[,y := value]
  longDat <- as.data.frame(longDat)
  y <- longDat$y
  fin <- longDat$Nvariable
  longDat$y <- cbind(y,fin)
  
  
  XZmodel <- gam(y~omega+s(Time,k=20,bs="tp",m=1) +
                   s(Time,omega,bs="sz",xt=list(bs="tp"),k=5,m=1) +
                   s(EpiWeek, k=15, bs="cc"), 
                 #s(EpiWeek,omega,bs="sz",xt=list(bs="cc")),
                 family=gfam(list(nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb)),method="REML",
                 data=longDat,knots=list(EpiWeek=c(0,52))) 
  #k.check(XZmodel)
  
  # predict:
  predDat <- data.table(longDat)
  predDat[,y.y := NULL]
  predDat[,y.fin := NULL]
  predDat[,y := Nvariable]
  #
  predDat$preds <- predict(XZmodel,newdata=predDat,type="response")
  modelMatrix <- predict(XZmodel,newdata=predDat,type="lpmatrix")
  n.sims <- 1000
  b <- rmvn(n.sims,coef(XZmodel),XZmodel$Vp)
  LP <- tcrossprod(modelMatrix,b)
  Mean <- exp(LP)
  myThetas <- exp( XZmodel$family$getTheta() )
  Y <- apply(Mean,2,function(x){rnbinom(length(x),mu=x,size=myThetas[12])})
  
  predDat$preds <- apply(Y,1,median)
  predDat$lwr <- apply(Y,1,quantile,probs=0.025)
  predDat$upr <- apply(Y,1,quantile,probs=0.975)
  
  #plot(y~Time,data=datS,type="l",lwd=3)
  #lines(datS$Time,predDat[Nvariable==12,preds],col="blue",lwd=2)
  #lines(datS$Time,predDat[Nvariable==12,lwr],col="blue",lwd=2,lty=2)
  #lines(datS$Time,predDat[Nvariable==12,upr],col="blue",lwd=2,lty=2)
  #abline(v=k)
  #abline(v=k-Ndelay+1)
  ## ggplot
  plotHere <- ggplot(data=dat,aes(x=Time)) + theme_bw() + geom_line(aes(y=y)) + 
    geom_ribbon(data=predDat[Nvariable==12,],aes(x=Time,y=preds,ymin=lwr,ymax=upr),alpha=0.2) + 
    geom_line(data=predDat[Nvariable==12,],aes(x=Time,y=preds),col="blue") + 
    ylim(0,12000) + geom_vline(xintercept = c(k+1,k-Ndelay+1)) + 
    labs(title=paste0("Nowcast weeks: ",k-Ndelay+1,"-",k,", Forecast weeks: ",k+1,"-",nrow(datS)))
  plotList[[jj]] <- plotHere
  jj <- jj+1
}

test <- do.call(grid.arrange,plotList)
test
ggsave("joint_model.pdf",test,width=15,height=10, path=plotDir)









### For talk
k <- 160
datS <- copy(dat)
for(i in 1:Ndelay){
  ind <- delay_indexHere[i] + 1
  datS[(k-i+1):k,ind] <- NA
}
# Now set everyhing to NA during forecasting
for(i in (k+1):(nrow(datS))){
  datS[i,delay_indexHere] <- NA
}
# this will set x to NA at all nowcasting
datS[,x := rowSums(.SD), .SDcols = delay_indexHere]

# Joint model for Z,X. First need to put in long format
Xind <- which(names(datS)=="x")
longDat <- melt(datS,id.vars=c("DATE","EpiWeek","Time","x"),measure.vars = c(grep("delay", colnames(datS) ),Xind)  )
longDat[,Nvariable := as.numeric(variable)]
longDat[,omega := factor(Nvariable)]
longDat[,y := value]
longDat <- as.data.frame(longDat)
y <- longDat$y
fin <- longDat$Nvariable
longDat$y <- cbind(y,fin)
XZmodel <- gam(y~omega+s(Time,k=20,bs="tp",m=1) +
                 s(Time,omega,bs="sz",xt=list(bs="tp"),k=5,m=1) +
                 s(EpiWeek, k=15, bs="cc"), 
               #s(EpiWeek,omega,bs="sz",xt=list(bs="cc")),
               family=gfam(list(nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb)),method="REML",
               data=longDat,knots=list(EpiWeek=c(0,52))) 

# predict:
predDat <- data.table(longDat)
predDat[,y.y := NULL]
predDat[,y.fin := NULL]
predDat[,y := Nvariable]
#
predDat$preds <- predict(XZmodel,newdata=predDat,type="response")
modelMatrix <- predict(XZmodel,newdata=predDat,type="lpmatrix")
n.sims <- 1000
b <- rmvn(n.sims,coef(XZmodel),XZmodel$Vp)
LP <- tcrossprod(modelMatrix,b)
Mean <- exp(LP)
myThetas <- exp( XZmodel$family$getTheta() )
Y <- apply(Mean,2,function(x){rnbinom(length(x),mu=x,size=myThetas[12])})
predDat$preds <- apply(Y,1,median)
predDat$lwr <- apply(Y,1,quantile,probs=0.025)
predDat$upr <- apply(Y,1,quantile,probs=0.975)
## ggplot
datS[,x_miss := rowSums(.SD,na.rm=T), .SDcols = delay_indexHere]
plotHere <- ggplot(data=datS[-c(1:35)],aes(x=DATE)) + theme_bw() + geom_line(aes(y=x_miss),col="red")+
  ylim(0,10000) + geom_vline(xintercept = c(predDat$DATE[k+1],predDat$DATE[k-Ndelay+1]))+
  annotate("text", x = predDat$DATE[k-5], y = 9500, label = "Nowcast",size = 3) + 
  annotate("text", x = predDat$DATE[k+10], y = 9500, label = "Forecast",size = 3) +
  labs(title=paste0("Nowcast weeks: ",k-Ndelay+1,"-",k,", Forecast weeks: ",k+1,"-",nrow(datS)))+
  ylab("Number of Dengue cases")
plotHere
ggsave("plot1.pdf",plotHere,width=10,height=5, path=plotDir)
plotHere <- ggplot(data=datS[-c(1:35)],aes(x=DATE)) + theme_bw() + geom_line(aes(y=x_miss),col="red")+
  geom_line(data=dat[-c(1:35)],aes(x=DATE,y=y)) +
  ylim(0,10000) + geom_vline(xintercept = c(predDat$DATE[k+1],predDat$DATE[k-Ndelay+1]))+
  annotate("text", x = predDat$DATE[k-5], y = 9500, label = "Nowcast",size = 3) + 
  annotate("text", x = predDat$DATE[k+10], y = 9500, label = "Forecast",size = 3) +
  labs(title=paste0("Nowcast weeks: ",k-Ndelay+1,"-",k,", Forecast weeks: ",k+1,"-",nrow(datS)))+
  ylab("Number of Dengue cases")
plotHere
ggsave("plot2.pdf",plotHere,width=10,height=5, path=plotDir)
plotHere <- ggplot(data=datS[-c(1:35)],aes(x=DATE)) + theme_bw() + 
  geom_ribbon(data=predDat[Nvariable==12,][-c(1:35)],aes(x=DATE,y=preds,ymin=lwr,ymax=upr),alpha=0.2) + 
  geom_line(data=predDat[Nvariable==12,][-c(1:35)],aes(x=DATE,y=preds),col="blue") +
  geom_line(aes(y=x_miss),col="red")+
  geom_line(data=dat[-c(1:35)],aes(x=DATE,y=y)) +
  ylim(0,10000) + geom_vline(xintercept = c(predDat$DATE[k+1],predDat$DATE[k-Ndelay+1]))+
  annotate("text", x = predDat$DATE[k-5], y = 9500, label = "Nowcast",size = 3) + 
  annotate("text", x = predDat$DATE[k+10], y = 9500, label = "Forecast",size = 3) +
  labs(title=paste0("Nowcast weeks: ",k-Ndelay+1,"-",k,", Forecast weeks: ",k+1,"-",nrow(datS)))+
  ylab("Number of Dengue cases")
plotHere
ggsave("plot3.pdf",plotHere,width=10,height=5, path=plotDir)


















gercutoff <- 160
dummy <- dengue.data[1:cutoff,]
dengue.data.full <- dengue.data.full[1:cutoff,]
full.length <- nrow(dengue.data)
dummy.length <- nrow(dummy)
for(i in 0:86){
  ind <- which(is.na(dengue.data[full.length-i,]))
  dummy[dummy.length-i,ind] <- NA
}
dengue.data <- dummy
rm(dummy)

# Same as delay.data.obs.train but with all the dates
delay.data.obs.trian2 <- data.frame(dengue.data[,6:16])

# adding NA values
make.df.trian <- function(M){
  Time <- nrow(M)
  Delay <- ncol(M)
  aux.df <- data.frame(Y = as.vector(as.matrix(M)), 
                       Time = rep(x = 1:Time, times = Delay),
                       
                       Delay = rep(x = 0:(Delay-1), each=Time)
  )
  aux.df
}

# creating my data
num.weeks <- nrow(delay.data.obs.trian2)
weeks <- rep( 1:52, 5 )
time_df2 <- 1:num.weeks
delay_df2 <- 0:(ncol(delay.data.obs.trian2)-1)
n_df2 <- as.vector(as.matrix(delay.data.obs.trian2))
eachtime2 <- length(n_df2)/length(time_df2)
eachdelay2 <- length(n_df2)/length(delay_df2)
CC <- rep(weeks[1:num.weeks], times=11)
mydata_Z_miss <- data.frame(Cases=n_df2, Time=rep(time_df2), Delay=rep(delay_df2,each=eachdelay2), Week=CC) # df3

# total number of cases
X <- 0
for (i in 1:num.weeks) {
  T2 <- filter(mydata_Z_miss, Time == i)
  X[i] <- sum(T2$Cases)
}
weeks <- weeks[1:num.weeks]
mydata_X_miss <- data.frame(Cases=X, Time=1:num.weeks, Week=weeks) # df4

# repeat with not missing data
delay.data.obs.trian.full <- data.frame(dengue.data.full[,6:16])
names(delay.data.obs.trian.full) <- names(delay.data.obs.trian2)
rownames(delay.data.obs.trian.full) <- dengue.data.full$DATE

# creating my data
fn_df2 <- as.vector(as.matrix(delay.data.obs.trian.full))
mydata_Z_full <- mydata_Z_miss 
mydata_Z_full$Cases <- fn_df2 #  df5

# total number of cases
X <- 0
for (i in 1:num.weeks) {
  fT2 <- filter(mydata_Z_full, Time == i)
  X[i] <- sum(fT2$Cases)
}
mydata_X_full <- data.frame(Cases=X, Time=1:num.weeks, Week=weeks) # df6



# Model for X
Xmodel <- gam(Cases ~ s(Time, k=7, bs="tp") + s(Week, k=5, bs="cc"), data=mydata_X_miss, family=nb(link="log"))
# Deterministic preds of totals X
Xpred <- predict(Xmodel,type="response",newdata = mydata_X_miss)

mydata_X_miss$predictedtotals <- mydata_X_miss$Cases
mydata_X_miss$predictedtotals[is.na(mydata_X_miss$Cases)] <- Xpred[is.na(mydata_X_miss$Cases)]
mydata_Z_miss$predictedLtotals <- log(rep(mydata_X_miss$predictedtotals,11))
mydata_Z_miss$Totals <- rep(mydata_X_miss$Cases, 11)
mydata_Z_miss$Ltotals <- log(mydata_Z_miss$Totals)
mydata_Z_full$Totals <- rep(mydata_X_full$Cases, 11)
mydata_Z_full$Ltotals <- log(mydata_Z_full$Totals)

Zmodel <- gam(Cases ~ s(Week, bs="cc", k=15) + ti(Delay, k=10, bs="tp") 
                +ti(Time, k=15, bs="tp") + ti(Time,Delay, k=5, bs="tp"), data=mydata_Z_miss, 
                family=nb(link="log"), method="REML",offset=mydata_Z_miss$Ltotals)

Zmodel <- gam(Cases ~ s(Week, bs="cc", k=15) + te(Time,Delay, k=10, bs="tp"), data=mydata_Z_miss, 
              family=nb(link="log"), method="REML",offset=mydata_Z_miss$Ltotals)

n <- 1000
betas <- rmvn(n,coef(Zmodel),Zmodel$Vp) ## 1000 values from the posterior of the coefficients
ModMatrix <- predict(Zmodel,newdata=mydata_Z_full,type="lpmatrix") # model matrix / predict out of sample?
logMUpreds <- ModMatrix %*% t(betas) ## 1000 values from the posterior of the linear predictor
MUpreds <- exp(logMUpreds + mydata_Z_miss$predictedLtotals) # add the predicted offset
phi <- Zmodel$family$getTheta(TRUE)
Zpreds <- apply(MUpreds,1,function(x){rnbinom(length(x),mu=x,size=phi)})
ZpredsArray <- array(dim=c(num.weeks,11,1000))
for (i in 1:n){
  ZpredsArray[,,i] <- matrix(Zpreds[i,],ncol = 11)
}
TotalsPreds <- apply(ZpredsArray,c(1,3),sum)

# plot fitted
Xobs <- dengue.data.full$Obs #mydata_X_full$Cases
plot(Xobs,type="l", xlab="Time (Weeks)", ylab="Cases")
lines(apply(TotalsPreds,1,mean),col="blue")
lines(apply(TotalsPreds,1,quantile,probs=0.025),col="blue",lty=2)
lines(apply(TotalsPreds,1,quantile,probs=0.975),col="blue",lty=2)
# nowcast
z_vector <- as.vector(as.matrix(delay.data.obs.trian2))
na.index <- which(!is.na(z_vector))
Znow <- Zpreds
for(i in na.index){ Znow[,i] <- z_vector[i] }
ZnowArray <- array(dim=c(num.weeks,11,1000))
for (i in 1:n){   ZnowArray[,,i] <- matrix(Znow[i,],ncol = 11) }
TotalsNow <- apply(ZnowArray,c(1,3),sum)

# plot
plot((num.weeks-10):(num.weeks),Xobs[(num.weeks-10):(num.weeks)],type="l", xlab="Time (Weeks)", ylab="Cases",xlim=c(num.weeks-10,num.weeks),ylim=c(0,max(apply(TotalsNow,1,quantile,probs=0.975)[(num.weeks-9):(num.weeks)])))
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,mean)[(num.weeks-9):(num.weeks)],col="blue")
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,quantile,probs=0.975)[(num.weeks-9):(num.weeks)],col="blue",lty=2)
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,quantile,probs=0.025)[(num.weeks-9):(num.weeks)],col="blue",lty=2)




# Bayesian + monte carlo stuff
# marginal p(Z) <- stays the same... outside for loop
n <- 1000
nY <-20000
phiX <- Xmodel$family$getTheta(TRUE)

Y_Z_nowcast <- matrix(nrow=10,ncol=4)
delays <- 10:1
betas <- rmvn(n,coef(Xmodel),Xmodel$Vp)

for(k in 1:10){
  print(k)
  pY_Z <- rep(0,nY)
  Time <- num.weeks-delays[k]+1
  Week <- weeks[num.weeks-delays[k]+1]
  ModMatrix <- predict(Xmodel,newdata=data.frame(Time=Time, Week=Week),type="lpmatrix") 
  logMUpreds <- ModMatrix %*% t(betas)
  MUpreds <- exp(logMUpreds)
  Xpreds <- rnbinom(1000, mu=MUpreds, size=phiX) ### 1000 values from post. pred. of corresponding X
  muZ <- predict(Zmodel,newdata=data.frame(Time=Time, Week=Week, Delay=0:10),type="response")[1:delays[k]]#*Xpreds
  Zs <- mydata_Z_miss$Cases[mydata_Z_miss$Delay%in%c(0:(delays[k]-1)) & mydata_Z_miss$Time==Time]
  pZ <- Zs
  for(kk in 1:length(Zs)){
    pZ[kk] <- log( mean( dnbinom(Zs[kk],mu=muZ[kk]*Xpreds,size=phi) ) )
  }
  pZ <- sum(pZ)
  #pZ <- log( mean( dnbinom(mydata_Z_miss$Cases[mydata_Z_miss$Delay==0 & mydata_Z_miss$Time==num.weeks],mu=muZ,size=phi) ) )

  # p(X|Z)
  for (Yt in 1:nY){
    # p(y) first need 1000 values of y_t
    #mu1 <- predict(Xmodel,newdata=data.frame(Time=171, Week=15),type="response")
    pY <- dnbinom(Yt, mu=predict(Xmodel,newdata=data.frame(Time=Time, Week=Week),type="response"), size=phiX,log=T) #mean( dnbinom(Yt, mu=MUpreds, size=phiX,log=T) )
    # p(Zt | Yt)
    muZ <-predict(Zmodel,newdata=data.frame(Time=num.weeks, Week=weeks[num.weeks], Delay=0:10),type="response")[1:delays[k]] * Yt
    pZ_Y <- sum( dnbinom(Zs,mu=muZ,size=phi,log=T) )
    pY_Z[Yt] <- exp(pY + pZ_Y - pZ)
  }
  YgivenZ <- sample(1:nY,100000,replace=T,prob=pY_Z)
  Y_Z_nowcast[k,1] <- mean(YgivenZ)
  Y_Z_nowcast[k,2] <- median(YgivenZ)
  Y_Z_nowcast[k,3:4] <- quantile(YgivenZ,probs=c(0.025,0.975))
}

# plot
x11()
par(mfrow=c(2,1))
plot(Xobs,type="l", xlab="Time (Weeks)", ylab="Cases")#,xlim=c(150,num.weeks))
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,median)[(num.weeks-9):(num.weeks)],col="blue")
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,quantile,probs=0.975)[(num.weeks-9):(num.weeks)],col="blue",lty=2)
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,quantile,probs=0.025)[(num.weeks-9):(num.weeks)],col="blue",lty=2)
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,2],col="red")
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,3],col="red",lty=2)
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,4],col="red",lty=2)
legend("topleft",c("GAM-Bayes","GAM"),lty=c(1,1),col=c("red","blue"))
plot((num.weeks-10):(num.weeks),Xobs[(num.weeks-10):(num.weeks)],type="l", xlab="Time (Weeks)", ylab="Cases",xlim=c(num.weeks-10,num.weeks),ylim=c(0,max(Y_Z_nowcast[,4])))
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,median)[(num.weeks-9):(num.weeks)],col="blue")
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,quantile,probs=0.975)[(num.weeks-9):(num.weeks)],col="blue",lty=2)
lines((num.weeks-9):(num.weeks),apply(TotalsNow,1,quantile,probs=0.025)[(num.weeks-9):(num.weeks)],col="blue",lty=2)
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,2],col="red")
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,3],col="red",lty=2)
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,4],col="red",lty=2)
legend("topleft",c("GAM-Bayes","GAM"),lty=c(1,1),col=c("red","blue"))

# now inla model as per paper
half_normal_sd <- function(sigma){
  return(
    paste("expression:
              sigma = ",sigma,";
              precision = exp(log_precision);
              logdens = -0.5*log(2*pi*sigma^2)-1.5*log_precision-1/(2*precision*sigma^2);
              log_jacobian = log_precision;
              return(logdens+log_jacobian);",sep='')
  )
} 
model <- Cases ~ 1 +
  f(Delay2, model = "rw1", constr=T, hyper = list("prec" = list(prior = half_normal_sd(1)))) +
  f(Week, model = "rw2", constr=T, hyper = list("prec" = list(prior = half_normal_sd(1))),cyclic=TRUE) +
  f(Time, model = "rw1", constr= T, hyper = list("prec" = list(prior = half_normal_sd(0.1)))) +
  f(Time2, model = "rw1", constr= T, hyper = list("prec" = list(prior = half_normal_sd(0.1))),replicate=Delay2) 
mydata_Z_miss$Time2 <- mydata_Z_miss$Time
mydata_Z_miss$Delay2 <- mydata_Z_miss$Delay+1
output <- inla(model, family = "nbinomial", data = mydata_Z_miss,
               control.predictor = list(link = 1, compute = T),
               control.compute = list( config = T, waic=TRUE, dic=TRUE),
               control.family = list( 
                 hyper = list("theta" = list(prior = "loggamma", param = c(1, 0.1))))
)
# Sampling from the posterior distribution using inla.posterior.sample 
n.sim <- 1000 # Number of posterior samples
delay.samples.list <- inla.posterior.sample(n = n.sim, output) # Posterior samples list
# Sampling z_t
vector.samples <- sapply(X = delay.samples.list, 
                         FUN = function(x, idx = 1:nrow(mydata_Z_miss)) 
                           rnbinom(n = idx, 
                                   mu = exp(x$latent[idx]), 
                                   size = x$hyperpar[1])) 
# Estimating Y_t as sum of z_t
Yt.fitted <- as_tibble(vector.samples) %>% 
  bind_cols(Time = mydata_Z_miss$Time) %>% 
  group_by(Time) %>% 
  summarise_all(sum) %>% gather(key = V, value = Ynew, -Time) %>% 
  group_by(Time) %>% summarise(
    Mean = mean(Ynew),
    Median = quantile(Ynew, probs = 0.5),
    Lower = quantile(Ynew, probs = 0.025),
    Upper = quantile(Ynew, probs = 0.975)
  )
dengue.fit <- select(dengue.data, EpiYear, EpiWeek, DATE, Time, Obs, Total) %>% left_join(Yt.fitted, by = "Time") 
# nowcast
index.missing <- which(is.na(mydata_Z_miss$Cases))
Yt.nowcast <- as_tibble(vector.samples[index.missing,]) %>% 
  bind_cols(Time = mydata_Z_miss$Time[index.missing]) %>% 
  group_by(Time) %>% 
  summarise_all(sum) %>% gather(key = V, value = Ynew, -Time) %>% 
  group_by(Time) %>% summarise(
    Mean = mean(Ynew),
    Median = quantile(Ynew, probs = 0.5),
    Lower = quantile(Ynew, probs = 0.025),
    Upper = quantile(Ynew, probs = 0.975)
  )
dengue.now <- select(dengue.data, EpiYear, EpiWeek, DATE, `10`, Time, Obs, Total) %>% 
  left_join(Yt.nowcast, by = "Time") %>%
  mutate(
    Mean = Obs + Mean,
    Median = Median + Obs,
    Lower = Lower + Obs,
    Upper = Upper + Obs
  )
x11()
par(mfrow=c(2,1))
plot(Xobs,type="l", xlab="Time (Weeks)", ylab="Cases")
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,2],col="red")
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,3],col="red",lty=2)
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,4],col="red",lty=2)
lines((num.weeks-9):(num.weeks),unlist(dengue.now[(num.weeks-9):(num.weeks),'Median']),col="blue")
lines((num.weeks-9):(num.weeks),unlist(dengue.now[(num.weeks-9):(num.weeks),'Lower']),col="blue",lty=2)
lines((num.weeks-9):(num.weeks),unlist(dengue.now[(num.weeks-9):(num.weeks),'Upper']),col="blue",lty=2)
legend("topleft",c("GAM-Bayes","INLA"),lty=c(1,1),col=c("red","blue"))
x11(width=10,height=6)
par(mfrow=c(1,1),mar = c(4, 4, 1, 3),cex=1.3,lwd=1)
plot((num.weeks-10):(num.weeks),Xobs[(num.weeks-10):(num.weeks)],type="l", xlab="Time (Weeks)", ylab="Cases",xlim=c(num.weeks-10,num.weeks),ylim=c(0,max(Y_Z_nowcast[,4])),lwd=2)
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,2],col="red")
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,3],col="red",lty=2)
lines((num.weeks-9):(num.weeks),Y_Z_nowcast[,4],col="red",lty=2)
lines((num.weeks-9):(num.weeks),unlist(dengue.now[(num.weeks-9):(num.weeks),'Median']),col="blue")
lines((num.weeks-9):(num.weeks),unlist(dengue.now[(num.weeks-9):(num.weeks),'Lower']),col="blue",lty=2)
lines((num.weeks-9):(num.weeks),unlist(dengue.now[(num.weeks-9):(num.weeks),'Upper']),col="blue",lty=2)
legend("topleft",c("Obs","GAM-Bayes","INLA"),lty=c(2,1,1),col=c("black","red","blue"))
dev.print(pdf,"GAMvsINLA_151.pdf")






