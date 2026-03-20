library(ggplot2)
library(mgcv)
library(data.table)
library(gridExtra)
library(readr)
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
stopLabel <- 163
maxT      <- max(dat$Time)

plotList <- list()
jj       <- 1
weeks    <- 159:170

# Ensure this exists before the loop
diag_table <- data.table()

for(k in weeks){

  #k <- 160
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


  # # Model for X
  # Xmodel <- gam(x ~ s(Time, k=15, bs="tp") + s(EpiWeek, k=10, bs="cc"), data=datS, 
  #             family=nb(link="log"),method="REML",knots=list(EpiWeek=c(0,52)))
  # k.check(Xmodel)
  # # Deterministic preds of totals X
  # Xpred <- predict(Xmodel,type="response",newdata = datS)
  # check how well it predicts Y
  #plot(y~Time,data=datS,type="l",lwd=3)
  #lines(datS$Time,Xpred,col="blue")

  # Model for Z. First need to put in long format
  # Xind <- which(names(datS)=="x")
  longDat <- melt(datS,id.vars=c("DATE","EpiWeek","Time"),measure.vars = grep("delay", colnames(datS) )  )
  longDat
  # delay covariate
  longDat[,delay := parse_number(as.character( variable ))]
  

  ## now fit GAM
  Zmodel <- gam(value ~ s(Time,k=20,bs="tp") +
              s(delay,k=8,bs="tp") + 
              ti(Time,delay,bs="tp",k=5) +
              s(EpiWeek, k=15, bs="cc"),
              family=nb,method="REML",
            data=longDat,knots=list(EpiWeek=c(0,52))) 
 k.check(Zmodel)
  
 gam.check(Zmodel)
 
 par(mfrow = c(2,2))
 plot(Zmodel, pages = 1, residuals = TRUE, shade = TRUE, seWithMean = TRUE)
 par(mfrow = c(1,1))
 
 vis.gam(Zmodel, view = c("Time","delay"),
         plot.type = "contour", too.far = 0.1)

 
 vis.gam(Zmodel, view = c("Time","delay"),
         plot.type = "persp", theta = 40, phi = 25, ticktype = "detailed")
 

 library(data.table)
 
 # ---- Predict from Zmodel ----
 pred_Z <- predict(Zmodel, type = "link", se.fit = TRUE)
 fit_Z  <- exp(pred_Z$fit)
 lwr_Z  <- exp(pred_Z$fit - 1.96 * pred_Z$se.fit)
 upr_Z  <- exp(pred_Z$fit + 1.96 * pred_Z$se.fit)
 
 truth_Z <- Zmodel$model$value  # actual counts used in fit
 
 # ---- Predict from XZmodel ----
 # pick the "total" margin
 mf_XZ <- model.frame(XZmodel)
 pred_XZ <- predict(XZmodel, type = "link", se.fit = TRUE)
 fit_XZ  <- exp(pred_XZ$fit[mf_XZ$omega == "total"])  # adjust label
 se_XZ   <- pred_XZ$se.fit[mf_XZ$omega == "total"]
 lwr_XZ  <- exp(pred_XZ$fit[mf_XZ$omega == "total"] - 1.96 * se_XZ)
 upr_XZ  <- exp(pred_XZ$fit[mf_XZ$omega == "total"] + 1.96 * se_XZ)
 
 truth_XZ <- mf_XZ$y[mf_XZ$omega == "total"]
 
 
 library(gratia)
 draw(Zmodel, select = "s(Time)")
 draw(Zmodel, select = "s(delay)")
 draw(Zmodel, select = "s(EpiWeek)")
 draw(Zmodel, select = "ti(Time,delay)")  
 
 graphics.off()
 dev.new()  # open a fresh plotting window
 
 
  # ---- Try to extract diagnostics ----
  mf <- model.frame(Zmodel)
  res_all <- residuals(Zmodel, type = "deviance")
  
  # figure out which level is "x" total
  x_idx <- unique(longDat$Nvariable[longDat$variable == "x"])
  x_level <- factor(x_idx, levels = levels(mf$omega))
  keep_x  <- (mf$omega == x_level) & is.finite(res_all)
  res_x   <- res_all[keep_x][order(mf$Time[keep_x])]
  res_x   <- res_x[is.finite(res_x)]
  
  # GOF ratio
  gof_ratio <- summary(Zmodel)$deviance / Zmodel$df.residual
  
  # Collect into a row
  met_row <- data.table(
    k_window_start = k - Ndelay + 1,
    k_window_end   = k,
    aic            = AIC(Zmodel),
    dev_expl       = summary(Zmodel)$dev.expl,
    gof_ratio      = gof_ratio
  )
  
  # Append
  diag_table <- rbind(diag_table, met_row, fill = TRUE)
  
  # Optional: print as you go
  print(met_row)
  
  # ---- Try to make plots but don’t stop if they fail ----
  try({
    png(sprintf("diag_gamcheck_k%03d.png", k), width = 1200, height = 900, res = 140)
    par(mfrow = c(2,2))
    gam.check(Zmodel, rep = 300, k.sample = 5000)
    dev.off()
  }, silent = TRUE)
  
  try({
    if (length(res_x) > 2) {
      png(sprintf("diag_acf_x_k%03d.png", k), width = 900, height = 700, res = 140)
      acf(res_x, main = sprintf("Residual ACF (x margin) weeks %d–%d", k - Ndelay + 1, k))
      dev.off()
    }
  }, silent = TRUE)
}

print(diag_table) 

  # predict:
  modelMatrix <- predict(Zmodel,newdata=longDat,type="lpmatrix")
  n.sims <- 1000
  b <- rmvn(n.sims,coef(Zmodel),Zmodel$Vp)
  LP <- tcrossprod(modelMatrix,b)
  Mean <- exp(LP)
  myThetas <- exp( Zmodel$family$getTheta() )
  Y <- apply(Mean,2,function(x){rnbinom(length(x),mu=x,size=myThetas)})
  # merge with longdat
  longDatPred <- cbind(longDat,Y)
  shortDatPred <- longDatPred[,lapply(.SD,sum),by=DATE,.SDcols=V1:V1000]
  shortDatPred[,preds := apply(shortDatPred[,V1:V1000],1,median)]
  shortDatPred[,upr := apply(shortDatPred[,V1:V1000],1,quantile,probs=0.025)]
  shortDatPred[,lwr := apply(shortDatPred[,V1:V1000],1,quantile,probs=0.975)]
  
  # shortDatPred[,preds := do.call(mean,.SD),.SDcols=V1:V1000,by=seq_len(nrow(shortDatPred))]
  # shortDatPred[,lwr := do.call(quantile,.SD,probs=0.025),.SDcols=V1:V1000,by=seq_len(nrow(shortDatPred))]
  # shortDatPred[,upr := do.call(quantile,.SD,probs=0.975),.SDcols=V1:V1000,by=seq_len(nrow(shortDatPred))]
  
  predDat <- merge(dat,shortDatPred[,c("DATE","preds","upr","lwr")],by="DATE")
  
# ## ggplot
# plotHere <- ggplot(data=dat,aes(x=Time)) + theme_bw() + geom_line(aes(y=y)) + 
#   geom_ribbon(data=predDat,aes(x=Time,y=preds,ymin=lwr,ymax=upr),alpha=0.2) + 
#   geom_line(data=predDat,aes(x=Time,y=preds),col="blue") + 
#   ylim(0,12000) + geom_vline(xintercept = c(k,k-Ndelay+1)) + 
#   labs(title=paste0("Nowcast weeks: ",k-Ndelay+1,"-",k,", Forecast weeks: ",k+1,"-",nrow(datS)))
# plotList[[jj]] <- plotHere
# jj <- jj+1


# Base ggplot with the yellow Nowcasting band
p <- ggplot(dat, aes(x = Time)) +
  theme_bw() +
  annotate("rect",
           xmin  = k - Ndelay + 1,
           xmax  = k + 1,
           ymin  = -Inf,
           ymax  = Inf,
           fill  = "yellow", alpha = 0.2
  ) +
  annotate("text",
           x     = (k - Ndelay + 1 + k + 1) / 2,
           y     = Inf,
           label = "NC",
           vjust = 2, fontface = "bold"
  ) +
  
  # ALWAYS draw the orange Forecasting band
  annotate("rect",
           xmin  = k + 1,
           xmax  = maxT,
           ymin  = -Inf,
           ymax  = Inf,
           fill  = "orange", alpha = 0.2
  )

# ONLY draw the Forecasting label for k ≤ stopLabel
if (k <= stopLabel) {
  p <- p +
    annotate("text",
             x     = (k + 1 + maxT) / 2,
             y     = Inf,
             label = "FC",
             vjust = 2, fontface = "bold"
    )
}

# then add your data layers, legend, and title
p <- p +
  geom_line(aes(y = y, colour = "Observed"), size = 1) +
  geom_ribbon(
    data    = predDat[Nvariable == 12, ],
    aes(x = Time, ymin = lwr, ymax = upr, fill = "95% PI"),
    alpha   = 0.2
  ) +
  geom_line(
    data    = predDat[Nvariable == 12, ],
    aes(x = Time, y = preds, colour = "Predicted"),
    size    = 1
  ) +
  scale_colour_manual(
    name   = NULL,
    values = c("Observed" = "black", "Predicted" = "blue")
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("95% PI" = "blue")
  ) +
  ylim(0, 12000) +
  labs(
    title = paste0(
      "Nowcast weeks: ",  k - Ndelay + 1, "–", k,
      ", Forecast weeks: ", k + 1, "–", maxT
    ),
    x = "Time",
    y = "Total cases"
  )

plotList[[jj]] <- p
jj <- jj + 1
}

# save each plot in plotList as a numbered PNG
for (i in seq_along(plotList)) {
  ggsave(
    filename = sprintf("outputs_og/plot_%02d.png", i),
    plot     = plotList[[i]],
    width    = 6,
    height   = 4,
    units    = "in",
    dpi      = 300
  )
}

#Make an animated plot
install.packages("magick")
library(magick)
library(dplyr)

# read in and sort the files
imgs <- list.files("outputs_og", pattern = "plot_\\d+\\.png$", full.names = TRUE) %>%
  sort() %>%
  image_read()

# animate at, say, 1 frame per second (fps)
anim <- image_animate(imgs, fps = 4)

# write to disk
image_write(anim, "nowcast_animation_og.gif")


library(ggplot2)

ggplot(dat, aes(x = DATE)) +
  # shaded area between current and eventual cases
  geom_ribbon(
    aes(ymin = delay0, ymax = x),
    fill  = "yellow",
    alpha = 0.4
  ) +
  # solid black for total cases
  geom_line(
    aes(y = x,
        colour   = "Total reported cases (account for delays)",
        linetype = "Total reported cases (account for delays)"),
    size = 1
  ) +
  # dashed red for delay-0 cases
  geom_line(
    aes(y = delay0,
        colour   = "Reported cases at Week 0",
        linetype = "Reported cases at Week 0"),
    size = 1
  ) +
  # manual legend entries
  scale_colour_manual(
    name   = NULL,
    values = c(
      "Total reported cases (account for delays)"   = "black",
      "Reported cases at Week 0"    = "red"
    )
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c(
      "Total reported cases (account for delays)"   = "solid",
      "Reported cases at Week 0"    = "11"
    )
  ) +
  # start x-axis around July 2010, keep up to the latest date
  scale_x_date(
    limits     = c(as.Date("2010-07-01"), max(dat$DATE, na.rm = TRUE)),
    date_breaks = "12 months",
    date_labels = "%Y"
  ) +
  labs(
    x = "Time (in epidemic weeks)",
    y = "Number of dengue cases"
  ) +
  theme_classic() +
  theme(
    legend.position     = c(0.25, 0.8),
    legend.background   = element_rect(fill = "white", colour = "black"),
    legend.key.size     = unit(0.4, "cm"),          # smaller key boxes
    legend.text         = element_text(size = 9)    # smaller text
  )




test <- do.call(grid.arrange,plotList)
test
ggsave("orig_model.pdf",test,width=15,height=10, path=plotDir)


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






