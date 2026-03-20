#library(tidyverse)
#library(lubridate)
install.packages("gratia")
library(ggplot2)
library(mgcv)
library(readr)
library(data.table)
library(gridExtra)
library(gratia)   # for appraise()
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

# Ensure this exists before the loop
diag_table <- data.table()

# Loop
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
  
  #Joint model 
  XZmodel <- gam(y~omega+s(Time,k=20,bs="tp",m=1) +
                   s(Time,omega,bs="sz",xt=list(bs="tp"),k=5,m=1) +
                   s(EpiWeek, k=15, bs="cc"), 
                 #s(EpiWeek,omega,bs="sz",xt=list(bs="cc")),
                 family=gfam(list(nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb)),method="REML", #12 NBs
                 data=longDat,knots=list(EpiWeek=c(0,52))) 
  k.check(XZmodel)
  
  gam.check(XZmodel)
  
  par(mfrow = c(2,2))
  plot(XZmodel, pages = 1, residuals = TRUE, shade = TRUE, seWithMean = TRUE)
  par(mfrow = c(1,1))
  
  # ---- Try to extract diagnostics ----
  mf <- model.frame(XZmodel)
  res_all <- residuals(XZmodel, type = "deviance")
  
  # figure out which level is "x" total
  x_idx <- unique(longDat$Nvariable[longDat$variable == "x"])
  x_level <- factor(x_idx, levels = levels(mf$omega))
  keep_x  <- (mf$omega == x_level) & is.finite(res_all)
  res_x   <- res_all[keep_x][order(mf$Time[keep_x])]
  res_x   <- res_x[is.finite(res_x)]
  
  # GOF ratio
  gof_ratio <- summary(XZmodel)$deviance / XZmodel$df.residual
  
  # Collect into a row
  met_row <- data.table(
    k_window_start = k - Ndelay + 1,
    k_window_end   = k,
    aic            = AIC(XZmodel),
    dev_expl       = summary(XZmodel)$dev.expl,
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
    gam.check(XZmodel, rep = 300, k.sample = 5000)
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
# }



# fwrite(diag_table, "diagnostics_summary.csv")


## amended code with legend
stopLabel <- 163
maxT      <- max(dat$Time)

plotList <- list()
jj       <- 1
k <- max(159:170)   # or simply k <- 170

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
  
  # predict:
  predDatXZ <- data.table(longDat)
  predDatXZ[,y.y := NULL]
  predDatXZ[,y.fin := NULL]
  predDatXZ[,y := Nvariable]
  #
  predDatXZ$preds <- predict(XZmodel,newdata=predDatXZ,type="response")
  modelMatrix <- predict(XZmodel,newdata=predDatXZ,type="lpmatrix")
  n.sims <- 1000
  b <- rmvn(n.sims,coef(XZmodel),XZmodel$Vp)
  LP <- tcrossprod(modelMatrix,b)
  Mean <- exp(LP)
  myThetas <- exp( XZmodel$family$getTheta() )
  Y <- apply(Mean,2,function(x){rnbinom(length(x),mu=x,size=myThetas[12])})
  
  predDatXZ$preds <- apply(Y,1,median)
  predDatXZ$lwr <- apply(Y,1,quantile,probs=0.025)
  predDatXZ$upr <- apply(Y,1,quantile,probs=0.975)
  
  colnames(predDatXZ)
  
  # Define metrics
  rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
  cov95 <- function(obs, lo, hi) mean(obs >= lo & obs <= hi, na.rm = TRUE)
  pi_width <- function(lo, hi) mean(hi - lo, na.rm = TRUE)
  
  
  # For Zmodel
  metrics_Z <- data.frame(
    Model = "Model Z",
    AIC = AIC(Zmodel),
    "Deviance explained" = summary(Zmodel)$dev.expl,
    RMSE = rmse(predDatZ$y, predDatZ$preds),
    MAE  = mae (predDatZ$y, predDatZ$preds),
    "95% IC" = cov95(predDatZ$y, predDatZ$lwr, predDatZ$upr),
    "PI Width" = pi_width(predDatZ$lwr, predDatZ$upr)
  )
  
  # For XZmodel
  metrics_XZ <- data.frame(
    Model = "Model XZ",
    AIC = AIC(XZmodel),
    "Deviance explained" = summary(XZmodel)$dev.expl,
    RMSE = rmse(predDatXZ$y, predDatXZ$preds),
    MAE  = mae (predDatXZ$y, predDatXZ$preds),
    "95% IC" = cov95(predDatXZ$y, predDatXZ$lwr, predDatXZ$upr),
    "PI Width" = pi_width(predDatXZ$lwr, predDatXZ$upr)
  )
  
  # Combine
 comp_table <- rbind(metrics_Z, metrics_XZ)
  # print(comp_table)
  
# Export to upload on Overleaf 
  # library(knitr)
  # library(kableExtra)
  
  comp_table %>%
    kable(format = "latex", booktabs = TRUE, digits = 3) %>%
    kable_styling(latex_options = c("hold_position")) %>%
    cat(file = "performance_comparison_table.tex")
  
  
  
  library(dplyr)
  
# metrics by delayed week  
  perf_XZ_by_margin <- predDatXZ %>%
    group_by(Nvariable) %>%
    summarise(
      rmse = rmse(y, preds),
      mae  = mae(y, preds),
      cov95 = cov95(y, lwr, upr),
      pi_width = pi_width(lwr, upr)
    )
  
summary(XZmodel)

plot(XZmodel)
  
  #plot(y~Time,data=datS,type="l",lwd=3)
  #lines(datS$Time,predDat[Nvariable==12,preds],col="blue",lwd=2)
  #lines(datS$Time,predDat[Nvariable==12,lwr],col="blue",lwd=2,lty=2)
  #lines(datS$Time,predDat[Nvariable==12,upr],col="blue",lwd=2,lty=2)
  #abline(v=k)
  #abline(v=k-Ndelay+1)
  ## ggplot
  # plotHere <- ggplot(data=dat,aes(x=Time)) + theme_bw() + geom_line(aes(y=y)) + 
  #   geom_ribbon(data=predDat[Nvariable==12,],
  #               aes(x=Time,y=preds,ymin=lwr,ymax=upr),alpha=0.2) + 
  #   geom_line(data=predDat[Nvariable==12,],
  #             aes(x=Time,y=preds),col="blue") + 
  #   ylim(0,12000) 
  # + geom_vline(xintercept = c(k+1,k-Ndelay+1)) + 
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

p

k.check(XZmodel)
gam.check(XZmodel)


# create an “outputs” folder
# if (!dir.exists("outputs")) dir.create("outputs")

# save each plot in plotList as a numbered PNG
for (i in seq_along(plotList)) {
  ggsave(
    filename = sprintf("outputs_joint/plot_%02d.png", i),
    plot     = plotList[[i]],
    width    = 6,
    height   = 4,
    units    = "in",
    dpi      = 300
  )
}

#Make an animated plot
# install.packages("magick")
library(magick)
library(dplyr)

# read in and sort the files
imgs <- list.files("outputs_joint", pattern = "plot_\\d+\\.png$", full.names = TRUE) %>%
  sort() %>%
  image_read()

# animate at, say, 1 frame per second (fps)
anim <- image_animate(imgs, fps = 4)

# write to disk
image_write(anim, "nowcast_animation_joint.gif")




test <- do.call(grid.arrange,plotList)
test
ggsave("joint_model.pdf",test,width=15,height=10, path=plotDir)


# 
# ## Longer loop
# plotList <- list()
# jj <- 1
# for(k in 125:170){
#   
#   # introduce NAs for nowcasting and forecasting
#   datS <- copy(dat)
#   for(i in 1:Ndelay){
#     ind <- delay_indexHere[i] + 1
#     datS[(k-i+1):k,ind] <- NA
#   }
#   # Now set everyhing to NA during forecasting
#   for(i in (k+1):(nrow(datS))){
#     datS[i,delay_indexHere] <- NA
#   }
#   # this will set x to NA at all nowcasting
#   datS[,x := rowSums(.SD), .SDcols = delay_indexHere]
#   
#   
#   # Model for X
#   #Xmodel <- gam(x ~ s(Time, k=15, bs="tp") + s(EpiWeek, k=10, bs="cc"), data=datS, 
#   #            family=nb(link="log"),method="REML",knots=list(EpiWeek=c(0,52)))
#   #k.check(Xmodel)
#   # Deterministic preds of totals X
#   #Xpred <- predict(Xmodel,type="response",newdata = datS)
#   # check how well it predicts Y
#   #plot(y~Time,data=datS,type="l",lwd=3)
#   #lines(datS$Time,Xpred,col="blue")
#   
#   # Joint model for Z,X. First need to put in long format
#   Xind <- which(names(datS)=="x")
#   longDat <- melt(datS,id.vars=c("DATE","EpiWeek","Time","x"),measure.vars = c(grep("delay", colnames(datS) ),Xind)  )
#   
#   ## now fit gfam GAM
#   # matrix response
#   longDat[,Nvariable := as.numeric(variable)]
#   longDat[,omega := factor(Nvariable)]
#   longDat[,y := value]
#   longDat <- as.data.frame(longDat)
#   y <- longDat$y
#   fin <- longDat$Nvariable
#   longDat$y <- cbind(y,fin)
#   
#   
#   XZmodel <- gam(y~omega+s(Time,k=20,bs="tp",m=1) +
#                    s(Time,omega,bs="sz",xt=list(bs="tp"),k=5,m=1) +
#                    s(EpiWeek, k=15, bs="cc"), 
#                  #s(EpiWeek,omega,bs="sz",xt=list(bs="cc")),
#                  family=gfam(list(nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb,nb)),method="REML",
#                  data=longDat,knots=list(EpiWeek=c(0,52))) 
#   #k.check(XZmodel)
#   
#   # predict:
#   predDat <- data.table(longDat)
#   predDat[,y.y := NULL]
#   predDat[,y.fin := NULL]
#   predDat[,y := Nvariable]
#   #
#   predDat$preds <- predict(XZmodel,newdata=predDat,type="response")
#   modelMatrix <- predict(XZmodel,newdata=predDat,type="lpmatrix")
#   n.sims <- 1000
#   b <- rmvn(n.sims,coef(XZmodel),XZmodel$Vp)
#   LP <- tcrossprod(modelMatrix,b)
#   Mean <- exp(LP)
#   myThetas <- exp( XZmodel$family$getTheta() )
#   Y <- apply(Mean,2,function(x){rnbinom(length(x),mu=x,size=myThetas[12])})
#   
#   predDat$preds <- apply(Y,1,median)
#   predDat$lwr <- apply(Y,1,quantile,probs=0.025)
#   predDat$upr <- apply(Y,1,quantile,probs=0.975)
#   
#   #plot(y~Time,data=datS,type="l",lwd=3)
#   #lines(datS$Time,predDat[Nvariable==12,preds],col="blue",lwd=2)
#   #lines(datS$Time,predDat[Nvariable==12,lwr],col="blue",lwd=2,lty=2)
#   #lines(datS$Time,predDat[Nvariable==12,upr],col="blue",lwd=2,lty=2)
#   #abline(v=k)
#   #abline(v=k-Ndelay+1)
#   ## ggplot
#   plotHere <- ggplot(data=dat,aes(x=Time)) + theme_bw() + geom_line(aes(y=y)) + 
#     geom_ribbon(data=predDat[Nvariable==12,],aes(x=Time,y=preds,ymin=lwr,ymax=upr),alpha=0.2) + 
#     geom_line(data=predDat[Nvariable==12,],aes(x=Time,y=preds),col="blue") + 
#     ylim(0,12000) + geom_vline(xintercept = c(k+1,k-Ndelay+1)) + 
#     labs(title=paste0("Nowcast weeks: ",k-Ndelay+1,"-",k,", Forecast weeks: ",k+1,"-",nrow(datS)))
#   plotList[[jj]] <- plotHere
#   jj <- jj+1
# }
# 
# test <- do.call(grid.arrange,plotList)
# ggsave("joint_model.pdf",test,width=15,height=10, path=plotDir)





