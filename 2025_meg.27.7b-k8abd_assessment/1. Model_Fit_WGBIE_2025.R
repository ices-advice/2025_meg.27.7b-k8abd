## Run analysis, write model results

library(devtools)
library(FLCore)
library(icesTAF)
library(FLa4a)
library(FLEDA)
library(ggplotFL)
library(gridExtra)
library(icesAdvice)
library(FLFishery)
library(FLasher)
library(msy)
library(a4adiags)
library(devtools)
library(ggplot2)
library(data.table)
library(readr)
library(tidyverse)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/6.4.Assessment_mgw78_cat1_FINAL")

# utilities
SavePlot <- function(plotname, width = 10, height = 7) {
  file <- paste0(getwd(), "/report/", "meg78_data_", plotname, ".png")
  dev.print(png, file, width = width, height = height, units = "in", res = 300)
}

mkdir("model")
set.seed(1234)

# load index, stock, the fit and the submodels
# load("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/6.4.meg78_review_FINAL/model/MegFit_FINAL.Rdata")
# load("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/6.4.meg78_review_FINAL/model/SubModels.Rdata")

# load index and stock
load("bootstrap/data/stock/meg78_stock.RData")
load("bootstrap/data/indices/meg78_indices.RData")

## ---------------------------------------------------------------------------------------------------------------------------------
# Fit the assessment model
## ---------------------------------------------------------------------------------------------------------------------------------

## Surveys

index<- tun.sel[c("SP_PORC","CPUE.IRLFRsurvey")]

plot(index[[1]]@index[1,])+geom_point()+ggtitle("index age 1")+ylab("Ind./30min haul")
plot(index[[2]]@index[1,])+geom_point()+ggtitle("index age 1")+ylab("numbers per 10 hours fished")

index[[2]]@index[ac(1:3),ac(2015:2021)] <- NA # IRLFR
index[[2]]@index['1','2011'] <- NA # strange outlier, set to NA in IRLFRsurvey 2011

index[[1]]@index[ac(1:3),ac(2015:2021)] <- NA # Porcupine

## Catch

plot(stock@catch.n['1',])

stock@catch.n['1',as.character(1984:2000)] <- NA # We do not really believe that the increase in 1-year-olds in the catch is real so we shouldn't formulate a model that treats it as real.

## Define the submodels

qmod <- list(~ I(1 / (1 + exp(-age))),~ I(1 / (1 + exp(-age))))

srmod <- ~ factor(year) 

fmod <- ~s(ifelse(year <= 2013 & age > 7, 7, 
                  ifelse(year >= 2014 & age > 6, 6, age)), k = 3, by = breakpts(year, 2013)) + factor(year)

n1mod <- ~s(age, k = 3)

vmod <- list(~s(age, k = 3), ~1, ~1)

## Fit
fit1 <- sca(stock, index, fmodel = fmod, qmodel = qmod, srmodel = srmod, n1model = n1mod, vmodel = vmod)

stk1 <- stock + fit1

plot(stk1)

submodels(fit1)

AIC(fit1)
BIC(fit1)

save(stk1, stock, index, fit1, file = "model/MegFit_FINAL.Rdata")
stk_est <- stk1
save(stk_est, file = "model/Meg78_WGMIX_stock_estimated.Rdata")
stk_obs <- stock
save(stk_obs, file = "model/Meg78_WGMIX_stock_input.Rdata")

## ---------------------------------------------------------------------------------------------------------------------------------
# Observed and predicted catches
## ---------------------------------------------------------------------------------------------------------------------------------                      

df <- tibble(
  Year = as.numeric(dimnames(catch(stock))$year),  
  Observed = as.numeric(catch(stock)),
  Predicted = as.numeric(catch(stk1))
) %>%
  pivot_longer(cols = c(Observed, Predicted), names_to = "Type", values_to = "Catch")

ggplot(df, aes(x = Year, y = Catch, color = Type)) +
  geom_line(size = 1.2) +
  labs(x = "Year",
       y = "Catch (tonnes)",
       color = "Catch type") +
  theme_bw()

Obs_pred_catch <- df %>%
  pivot_wider(names_from = Type, values_from = Catch) %>%
  mutate(
    Observed = round(Observed),
    Predicted = round(Predicted)
  )

write.table(Obs_pred_catch,'report/Obs_pred_catch.csv', row.names = FALSE, sep = ";")
## ---------------------------------------------------------------------------------------------------------------------------------
# Retrospective analysis
## ---------------------------------------------------------------------------------------------------------------------------------                      

## Run the retro without smoothers
back <- 5
retro <- split(1:back, 1:back)
retro <- lapply(retro, function(x){
  yr <- range(stock)["maxyear"] - x
  stk <- window(stock, end=yr)
  idx <- window(index, end=yr)
  stk + sca(stk, idx,fmodel=fmod, qmodel=qmod, srmodel=srmod)
})

retro$"0" <- stock+fit1


new_names <- c("Rec" = "Recruitment", "SB" = "SSB", "C" = "Catch", "F" = "F")
plot(FLStocks(retro),col=1,lwd=1)+ facet_wrap(~qname,scales='free_y', labeller = labeller(qname = new_names)) + theme_bw() +
  labs(color = "N years removed")

save(retro, file = "model/MegRetroFit_FINAL.Rdata")

SavePlot('Retro_Fit_FINAL')

## ERROR in Retro value calculation

error <- function(years,x,se){
  #  polygon(c(years,rev(years)),c(x-se,rev(c(x+se))),col='#0000FF20',border=NA)
  polygon(c(years,rev(years)),c(x-1.96*se,rev(c(x+1.96*se))),col='#0000FF20',border=NA)
}

pal <- c("#7570B3","#1B9E77","#E6AB02","#D95F02","#E7298A")
years <- stk1@range[4]:stk1@range[5]
nyears <-length(years)
fits <- simulate(fit1, 1000)
stks <- stock + fits

# F
Retro_F <- data.frame(Y0=c(fbar(retro$`0`))
                      ,Y1=c(fbar(retro$`1`),NA)
                      ,Y2=c(fbar(retro$`2`),NA,NA)
                      ,Y3=c(fbar(retro$`3`),NA,NA,NA)
                      ,Y4=c(fbar(retro$`4`),NA,NA,NA,NA)                
                      ,Y5=c(fbar(retro$`5`),NA,NA,NA,NA,NA)
)

mohn(Retro_F,details=T) 
mohn(Retro_F,plot=T) 

fbar <- fbar(stk1)
fbarse <- apply(fbar(stks),2,sd) # AI: considering 1sd
fbar2sd <- 2*apply(fbar(stks),2,sd) # AI: considering 2sd
ylim <- c(0,max(fbar+1.96*fbarse))
ylim <- c(0,max(fbar+1.96*fbar2sd))
xlim <-range(years)
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='F',main=paste0('Fbar ',stk1@range[6],'-',stk1@range[7]))
error(years,fbar,fbarse)
error(years,fbar,fbar2sd)
for(i in 1:5) lines(years[1:(nyears-i)],fbar(retro[[i]]),lwd=1,col=pal[i])
lines(years,fbar,lwd=2,col=4)
#abline(h=0.28,lty=2)
#abline(h=0.181,lty=3)
#abline(h=0.39,lty=3)
#lines(1986:2005,c(0.3486,0.3176,0.3423,0.3854,0.3835,0.3433,0.271,0.2041,0.2105,0.2664,0.3087,0.3065,0.2512,0.155,0.1261,0.182,0.1981,0.2328,0.2376,0.1622))
#legend('topleft',c('Fbar','Fmsy','Frange'),lty=c(1,2,3),col=c(4,1,1),lwd=c(2,1,1),bty='n',y.intersp=0.8)
SavePlot('Retro_Fit_Fbar')

# SSB
Retro_SSB <- data.frame(Y0=c(ssb(retro$`0`))
                        ,Y1=c(ssb(retro$`1`),NA)
                        ,Y2=c(ssb(retro$`2`),NA,NA)
                        ,Y3=c(ssb(retro$`3`),NA,NA,NA)
                        ,Y4=c(ssb(retro$`4`),NA,NA,NA,NA)                
                        ,Y5=c(ssb(retro$`5`),NA,NA,NA,NA,NA)
)

mohn(Retro_SSB) 
mohn(Retro_SSB,details=T) 
mohn(Retro_SSB,plot=T) 

ssb <- ssb(stk1)/1000
ssbse <- apply(ssb(stks),2,sd)/1000 # AI: considering 1sd
ssb2sd <- 2*apply(ssb(stks),2,sd)/1000 # AI: considering 2sd
ylim <- c(0,max(ssb+1.96*ssbse))# AI: considering 1sd
ylim <- c(0,max(ssb+ssb2sd))# AI: considering 2sd
xlim <-range(years)
plot(NA,xlim=xlim,ylim=ylim,xlab='Year',ylab='SSB',main=paste0('SSB '))
error(years,ssb,ssbse) #1sd
error(years,ssb,ssb2sd) #2sd
for(i in 1:5) lines(years[1:(nyears-i)],ssb(retro[[i]])/1000,lwd=1,col=pal[i])
lines(years,ssb,lwd=2,col=4)
SavePlot('Retro_Fit_SSB')

# Rec
recr <- function(x) x@stock.n[1,]
Retro_R <- data.frame(Y0=c(recr(retro$`0`))
                      ,Y1=c(recr(retro$`1`),NA)
                      ,Y2=c(recr(retro$`2`),NA,NA)
                      ,Y3=c(recr(retro$`3`),NA,NA,NA)
                      ,Y4=c(recr(retro$`4`),NA,NA,NA,NA)                
                      ,Y5=c(recr(retro$`5`),NA,NA,NA,NA,NA)
)
mohn(Retro_R) 
mohn(Retro_R,details=T) 
mohn(Retro_R,plot=T) 

## --------------------------------------------------------------------------------------------------------------------------------- 
# Plot residuals
## --------------------------------------------------------------------------------------------------------------------------------- 

mkdir("report")

# load fit, index and stock
load("model/MegFit_FINAL.Rdata")

# load retro run
load("model/MegRetroFit_FINAL.Rdata")

# hacked function from a4aFitresiduals-class
plotr <- function(x, y=missing, auxline="smooth",...){
  args <- list()
  args$data <- as.data.frame(x)
  args$x <- as.formula("data~year|factor(age)*qname")
  args$type=c("p", auxline)
  args$groups <- quote(qname)
  args$cex=0.3
  args$lwd=0
  args$ylab="standardized residuals"
  args$xlab=""
  args$panel=function(x,y,...){
    panel.abline(h=0, col.line="gray80")
    panel.xyplot(x,y,...)
  }
  args$par.settings=list(
    superpose.symbol=list(col=1, pch=19, cex=0.2),
    superpose.line=list(col="gray75", lty=1, lwd=0, col=NA),
    strip.background=list(col="gray90"),
    strip.border=list(col="black"),
    box.rectangle=list(col="gray90"))
  args$main="log residuals of catch and abundance indices by age"
  if(is(latticeExtra::useOuterStrips, "function")) latticeExtra::useOuterStrips(do.call("xyplot", args)) else do.call("xyplot", args)
}

# Res1
res <- residuals(fit1, stock, index)
print(plotr(res, auxline = "none"))
SavePlot('Residuals1',10,6)

# Res2
res$catch.n[is.na(res$catch.n)] <- 0 #hack
#res$FR_EVHOE[is.na(res$FR_EVHOE)] <- 0 #hack
res$SP_PORC[is.na(res$SP_PORC)] <- 0 #hack
#res$CPUE.Vigo84[is.na(res$CPUE.Vigo84)] <- 0 #hack
#res$CPUE.Vigo99[is.na(res$CPUE.Vigo99)] <- 0 #hack
#res$LPUE.ITBB[is.na(res$LPUE.ITBB)] <- 0 #hack
#res$CPUE.IAM[is.na(res$CPUE.IAM)] <- 0 #hack
res$CPUE.IRLFRsurvey[is.na(res$CPUE.IRLFRsurvey)] <- 0 #hack
bubbles(res)
SavePlot('Residuals2',10,6)

# Res3
res <- residuals(fit1, stock, index)
qqmath(res)
SavePlot('Residuals3')

wireframe(fit1@harvest)
#a <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,0.9),ylab='F',main='Fishing mortality')
a <- xyplot(data~age,groups=year,stk1@harvest,type='b',ylim=c(0,1),ylab='F',main='Fishing mortality')
a
SavePlot('F')


## --------------------------------------------------------------------------------------------------------------------------------- 
# Summary plots
## --------------------------------------------------------------------------------------------------------------------------------- 

fitted <- predict(pars(fit1))

#a1 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchabilty',main='Q FR_EVHOE',ylim=c(0,4e-5))
#a1
a2 <- xyplot(data~age,groups=year,data=fitted$qmodel[1],type='b',ylab='Catchabilty',main=names(index)[1])
a2
#a3 <- xyplot(data~age,groups=year,data=fitted$qmodel[3],type='b',ylab='Catchabilty',main='Q CPUE.Vigo84',ylim=c(0,0.005))
#a3
#a4 <- xyplot(data~age,groups=year,data=fitted$qmodel[4],type='b',ylab='Catchabilty',main='Q CPUE.Vigo99',ylim=c(0,0.0015))
#a4
#a5 <- xyplot(data~age,groups=year,data=fitted$qmodel[5],type='b',ylab='Catchabilty',main='Q LPUE.ITBB',ylim=c(0,0.02))
#a5
#a6 <- xyplot(data~age,groups=year,data=fitted$qmodel[6],type='b',ylab='Catchabilty',main='CPUE.IAM',ylim=c(0,0.03))
#a6
a7 <- xyplot(data~age,groups=year,data=fitted$qmodel[2],type='b',ylab='Catchabilty',main=names(index)[2])
a7

grid.arrange(a,a2,a7,ncol=2)
SavePlot('F and Q',8,6)

# simulate

fits <- simulate(fit1, 1000)
flqs <- FLQuants(sim=iterMedians(stock.n(fits)), det=stock.n(fit1))
keylst <- list(points=FALSE, lines=TRUE, space="right")
xyplot(data~year|factor(age), groups=qname, data=flqs, type="l", main="Median simulations VS fit", scales=list(y=list(relation="free")), auto.key=keylst)
#fits@stock.n[,1:2] <- fit1@stock.n[,1:2]
#fits@catch.n[,1:2] <- fit1@catch.n[,1:2]
stks <- stock + fits
plot(stks) + facet_wrap(~qname,scales='free_y', labeller = labeller(qname = new_names)) +theme_bw()
SavePlot('summary0')

## --------------------------------------------------------------------------------------------------------------------------------- 
# Summary table for the report
## --------------------------------------------------------------------------------------------------------------------------------- 

# Recruitment intermediate year and advice year:
years <- stk1@range[4]:stk1@range[5]
nyears <- length(years)

Rassumption <- "lastwoyr"
if(Rassumption == "gmean"){
  GM <- round(exp(mean(log(c(stk1@stock.n[1,1:(nyears-2)])))),0) # thousands. GM time series (without the 	last 2 years)
}
if(Rassumption == "lastwoyr"){
  GM <- round(exp(mean(log(c(stk1@stock.n[1,(nyears-1):nyears])))),0) # thousands. GM last two years
}

GM

### F status quo (Fsq)

# Extract fishing mortality (F) and plot it
f_values <- harvest(stk1)
f_mean_per_year <- apply(f_values[3:6,], 2, mean)
plot(f_mean_per_year)

# F setting
Fassumption <- "Funscaled" # use the definitive assumption not both!! "Fscaled" or "Funscaled"

if (Fassumption == "Fscaled"){
  fsq <- fbar(stk1)[,nyears] # F SCALED to Fbar of final assessment year (if there is a decreasing trend of F in the resutls of the assessment time-series)
}
if(Fassumption == "Funscaled"){
  fsq <-mean(fbar(stk1)[,nyears-2:0]) # F UNSCALED (Average F last 3 years) 
}

# AI: 08/06/2022 in the ADG the 2sd is requested and it is included in all xxse variables for example: tsbse <- 2*apply(tsb(stks),2,sd).

# # Summary table using coeficient interval of 1sd
# tsb <- tsb(stk1)
# tsbse <- apply(tsb(stks),2,sd)
# ssb <- ssb(stk1)
# ssbse <- apply(ssb(stks),2,sd)
# ssbInt <- NA #get this from stf
# catchobs <- stock@catch
# discardsobs <- stock@discards
# catch <- stk1@catch
# recr <- stk1@stock.n[1,]
# recrse <- apply(stks@stock.n[1,],2,sd)
# fbar <- fbar(stk1)
# fbarse <- apply(fbar(stks),2,sd)

# # Summary table using coeficient interval of 2sd
# 
# tsb <- tsb(stk1)
# tsb2sd <- 2*apply(tsb(stks),2,sd)
# ssb <- ssb(stk1)
# ssb2sd <- 2*apply(ssb(stks),2,sd)
# ssbInt <- NA #get this from stf
# catchobs <- stock@catch
# discardsobs <- stock@discards
# catch <- stk1@catch
# recr <- stk1@stock.n[1,]
# recr2sd <- 2*apply(stks@stock.n[1,],2,sd)
# fbar <- fbar(stk1)
# fbar2sd <- 2*apply(fbar(stks),2,sd)
# 
# 
# 
# sum1 <-
#   data.frame(
#     Year = c(years, paste0(max(years) + 1, "*")),
#     Lan = c(catchobs - discardsobs, NA),
#     Dis = c(discardsobs, NA),
#     Cat = c(catchobs, NA),
#     CatEst = c(catch, NA),
#     Tsb = c(tsb, NA),
#     Ssb = c(ssb, ssbInt),
#     SsbCv = c(ssbse / ssb, NA),
#     Recr = c(recr, GM),
#     RecrCv = c(recrse / recr, NA),
#     Fbar = c(fbar, fsq),
#     FbarCv = c(fbarse / fbar, NA)
#   )
# 
# knitr::kable(sum1,row.names=F,digits=c(rep(3,12)))
# write.csv(sum1,'report/Summary_FINAL.csv',row.names=F)
# 
# plot(years,tsb,type='l',ylim=c(0,120))
# lines(years,ssb,type='l',col=2)
# plot(years,ssb/tsb,type='l',ylim=c(0,1))


# # Summary for standard graphs template with coeficient interval of 1sd
# sum2 <- data.frame(year=c(years,max(years+1))
#                    ,recrlo=c(recr-recrse,NA)
#                    ,recr=c(recr,GM)
#                    ,recrhi=c(recr+recrse,NA)
#                    ,tsblo=c(tsb-tsbse,NA)
#                    ,tsb=c(tsb,NA)
#                    ,tsbhi=c(tsb+tsbse,NA)
#                    ,ssblo=c(ssb-ssbse,NA)
#                    ,ssb=c(ssb,ssbInt)
#                    ,ssbhi=c(ssb+ssbse,NA)
#                    ,catch=c(catchobs,NA)
#                    ,lan=c(catchobs-discardsobs,NA)
#                    ,dis=c(discardsobs,NA)
#                    ,ibc=NA
#                    ,ur=NA
#                    ,yssb=NA
#                    ,flo=c(fbar-fbarse,NA)
#                    ,f=c(fbar,fsq)
#                    ,fhi=c(fbar+fbarse,NA)
# )
# knitr::kable(sum2,row.names=F,digits=c(rep(0,16),rep(3,3)))
# write.csv(sum2,'report/Summary_2_1sd.csv', row.names=F)
# 
# 
# 
# # Summary for standard graphs template with coeficient interval of 2sd
# sum3 <- data.frame(year=c(years,max(years+1))
#                    ,recrlo=c(recr-recr2sd,NA)
#                    ,recr=c(recr,GM)
#                    ,recrhi=c(recr+recr2sd,NA)
#                    ,tsblo=c(tsb-tsb2sd,NA)
#                    ,tsb=c(tsb,NA)
#                    ,tsbhi=c(tsb+tsb2sd,NA)
#                    ,ssblo=c(ssb-ssb2sd,NA)
#                    ,ssb=c(ssb,ssbInt)
#                    ,ssbhi=c(ssb+ssb2sd,NA)
#                    ,catch=c(catchobs,NA)
#                    ,lan=c(catchobs-discardsobs,NA)
#                    ,dis=c(discardsobs,NA)
#                    ,ibc=NA
#                    ,ur=NA
#                    ,yssb=NA
#                    ,flo=c(fbar-fbar2sd,NA)
#                    ,f=c(fbar,fsq)
#                    ,fhi=c(fbar+fbar2sd,NA)
# )
# knitr::kable(sum3,row.names=F,digits=c(rep(0,16),rep(3,3)))
# sum3 <- as.data.table(sum3)
# 
# sum3[, c(2:13) := lapply(.SD, round, digits=0), .SDcols = 2:13]
# sum3[, c(17:19) := lapply(.SD, icesRound), .SDcols = 17:19]
# write.csv(sum3,'report/Summary_3_2sd.csv', row.names=F)

# Summary table using %90 confidence interval
# Lognormal-based confidence intervals (90%) for key variables

tsb <- tsb(stk1)
ssb <- ssb(stk1)
catchobs <- stock@catch
discardsobs <- stock@discards
catch <- stk1@catch
recr <- stk1@stock.n[1,]
fbar <- fbar(stk1)

# Parámetro z para 90% CI
z <- 1.6449

# Valores medios
tsb <- tsb(stk1)
ssb <- ssb(stk1)
recr <- stk1@stock.n[1,]
fbar <- fbar(stk1)

# Desviaciones estándar (estimadas a partir de las iteraciones en stks)
tsbse <- apply(tsb(stks), 2, sd)
ssbse <- apply(ssb(stks), 2, sd)
recrse <- apply(stks@stock.n[1,], 2, sd)
fbarse <- apply(fbar(stks), 2, sd)

# Convertir a vectores
tsb_mean <- as.numeric(tsb)
ssb_mean <- as.numeric(ssb)
recr_mean <- as.numeric(recr)
fbar_mean <- as.numeric(fbar)

# Catch, landings y discards
catch <- as.numeric(stk1@catch)
discards <- as.numeric(stock@discards)
landings <- catch - discards

# Función para calcular IC log-normal
calc_lognorm_CI <- function(mean, se, z = 1.6449) {
  log_sd <- sqrt(log(1 + (z * se / mean)^2))
  lower <- exp(log(mean) - log_sd)
  upper <- exp(log(mean) + log_sd)
  list(lower = lower, upper = upper)
}

# Aplicar a los indicadores
tsb_CI <- calc_lognorm_CI(tsb_mean, tsbse)
ssb_CI <- calc_lognorm_CI(ssb_mean, ssbse)
recr_CI <- calc_lognorm_CI(recr_mean, recrse)
fbar_CI <- calc_lognorm_CI(fbar_mean, fbarse)

# Añadir valores para 2025
tsb_mean <- c(tsb_mean, NA)
tsb_CI$lower <- c(tsb_CI$lower, NA)
tsb_CI$upper <- c(tsb_CI$upper, NA)

ssb_mean <- c(ssb_mean, NA)
ssb_CI$lower <- c(ssb_CI$lower, NA)
ssb_CI$upper <- c(ssb_CI$upper, NA)

recr_mean <- c(recr_mean, GM)
recr_CI$lower <- c(recr_CI$lower, NA)
recr_CI$upper <- c(recr_CI$upper, NA)

fbar_mean <- c(fbar_mean, fsq)
fbar_CI$lower <- c(fbar_CI$lower, NA)
fbar_CI$upper <- c(fbar_CI$upper, NA)

catch <- c(catch, NA)
landings <- c(landings, NA)
discards <- c(discards, NA)

# Crear tabla resumen final
summary_df <- data.frame(
  year = c(years, max(years) + 1),
  recr_lo = recr_CI$lower,
  recr = recr_mean,
  recr_hi = recr_CI$upper,
  tsb_lo = tsb_CI$lower,
  tsb = tsb_mean,
  tsb_hi = tsb_CI$upper,
  ssb_lo = ssb_CI$lower,
  ssb = ssb_mean,
  ssb_hi = ssb_CI$upper,
  fbar_lo = fbar_CI$lower,
  fbar = fbar_mean,
  fbar_hi = fbar_CI$upper,
  catch=c(catchobs,NA),
  lan=c(catchobs-discardsobs,NA),
  dis=c(discardsobs,NA)
)
summary_df <- as.data.table(summary_df)
summary_df[, c(2:10,14:16) := lapply(.SD, round, digits=0), .SDcols = c(2:10,14:16)]
summary_df[, c(11:13) := lapply(.SD, icesRound), .SDcols = c(11:13)]
# Exportar CSV
write.csv(summary_df, "report/Summary_SAG_lognorm.csv", row.names = FALSE)

# Sensitivity for template - not needed anymore?
ages <- stk1@range[1]:stk1@range[2]
stk3y <- window(stk1, start = max(years) - 2)
p <- apply(landings.n(stk3y) / catch.n(stk3y), 1, mean, na.rm = T)
p <- c(ifelse(is.na(p), 0, p))
sen <- data.frame(
  Age = ages,
  M = c(apply(m(stk3y), 1, mean)),
  Mat = c(apply(mat(stk3y), 1, mean)),
  PF = c(apply(harvest.spwn(stk3y), 1, mean)),
  PM = c(apply(m.spwn(stk3y), 1, mean)),
  Sel = c(apply(harvest(stk3y), 1, mean)) * p,
  WeCa = c(apply(landings.wt(stk3y), 1, mean)),
  Fd = c(apply(harvest(stk3y), 1, mean)) * (1 - p),
  WeCad = c(apply(discards.wt(stk3y), 1, mean)),
  Fi = 0,
  WeCai = 0
)
knitr::kable(sen, row.names = F, digits = c(0, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3))
write.csv(sen, "report/Sen.csv", row.names = F)



plot(fit1, stock)
SavePlot("Fit1")
plot(fit1, index[1])
SavePlot("Fit2")
plot(fit1, index[2])
SavePlot("Fit3")

fitSumm(fit1)

# iters
# 1
# nopar       9.800000e+01
# nlogl       5.262983e+02
# maxgrad     1.220970e-05
# nobs        8.030000e+02
# gcv         1.999128e-01
# convergence 0.000000e+00
# accrate               NA
# nlogl_comp1 8.359910e+01
# nlogl_comp2 1.944460e+02
# nlogl_comp3 2.482530e+02

#### KOBE PLOT (MAE: updated RP for megrim 78 ICES. 2024 SA.)

dev.off()

kobe <- data.frame(Btrig=43947,Fmsy=0.27,B=rev(c(ssb(stk1))),F=rev(c(fbar(stk1))))
with(kobe[1,],plot(B/Btrig,F/Fmsy,xlim=c(0,3.5),ylim=c(0,3),main='meg78'))
rect(-1,-1,5,5,col='yellow',border=NA)
rect(-1,1,1,5,col='red',border=NA)
rect(1,-1,5,1,col='green',border=NA)
with(kobe[1,],text(B/Btrig,F/Fmsy,max(years),pos=1))
with(kobe[1,],points(B/Btrig,F/Fmsy,pch=16))
with(kobe,lines(B/Btrig,F/Fmsy))
#plot(kobe)
SavePlot('Kobe',6,6)

dev.off()
