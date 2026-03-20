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
k        <- 170

# Ensure this exists before the loop
diag_table <- data.table()

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
  
  # 
  # library(gratia)
  # draw(Zmodel, select = "s(Time)")
  # draw(Zmodel, select = "s(delay)")
  # draw(Zmodel, select = "s(EpiWeek)")
  # draw(Zmodel, select = "ti(Time,delay)")  
  # 
  # graphics.off()
  # dev.new()  # open a fresh plotting window
  # 
  # 
  # # ---- Try to extract diagnostics ----
  # mf <- model.frame(Zmodel)
  # res_all <- residuals(Zmodel, type = "deviance")
  # 
  # # figure out which level is "x" total
  # x_idx <- unique(longDat$Nvariable[longDat$variable == "x"])
  # x_level <- factor(x_idx, levels = levels(mf$omega))
  # keep_x  <- (mf$omega == x_level) & is.finite(res_all)
  # res_x   <- res_all[keep_x][order(mf$Time[keep_x])]
  # res_x   <- res_x[is.finite(res_x)]
  # 
  # # GOF ratio
  # gof_ratio <- summary(Zmodel)$deviance / Zmodel$df.residual
  # 
  # # Collect into a row
  # met_row <- data.table(
  #   k_window_start = k - Ndelay + 1,
  #   k_window_end   = k,
  #   aic            = AIC(Zmodel),
  #   dev_expl       = summary(Zmodel)$dev.expl,
  #   gof_ratio      = gof_ratio
  # )
  # 
  # # Append
  # diag_table <- rbind(diag_table, met_row, fill = TRUE)
  # 
  # # Optional: print as you go
  # print(met_row)

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
shortDatPred[,upr := apply(shortDatPred[,V1:V1000],1,quantile,probs=0.975)]
shortDatPred[,lwr := apply(shortDatPred[,V1:V1000],1,quantile,probs=0.025)]

# shortDatPred[,preds := do.call(mean,.SD),.SDcols=V1:V1000,by=seq_len(nrow(shortDatPred))]
# shortDatPred[,lwr := do.call(quantile,.SD,probs=0.025),.SDcols=V1:V1000,by=seq_len(nrow(shortDatPred))]
# shortDatPred[,upr := do.call(quantile,.SD,probs=0.975),.SDcols=V1:V1000,by=seq_len(nrow(shortDatPred))]

predDatZ <- merge(dat,shortDatPred[,c("DATE","preds","upr","lwr")],by="DATE")
colnames(predDatZ)

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

p

#Model checking
# Draw just the QQ plot for the model
mgcv::qq.gam(Zmodel, rep = 100, level = 0.95)  # envelope
title("Original model")         # add the header manually

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




