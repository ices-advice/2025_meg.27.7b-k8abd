# model.R - Run forecast
# 2020_sol.27.4_forecast/model.R

# Copyright Iago MOSQUEIRA (WMR), 2020
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# install.packages("devtools")
# install.packages("FLCore", repo = "http://flr-project.org/R")
# install_github("ices-tools-prod/msy")
# install.packages(c("FLa4a", "FLasher", "ggplotFL"), repos="https://flr-project.org/R")
# install.packages(c("ggplot2", "snpar", "foreach", "data.table"))
# devtools::install_github("flr/a4adiags")

library(FLCore)
library(icesTAF)
library(ggplot2)
library(ggplotFL)
library(FLFishery)
library(FLasher)
library(FLa4a)
library(a4adiags)
library(icesAdvice)

setwd("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/6.4.Assessment_mgw78_cat1_FINAL")
mkdir("model")
set.seed(1234)

# setting
Fassumption <- "Funscaled" #"Fscaled" or "Funscaled"
Rassumption <- "lastwoyr" # "gmean" or "meanblowmean" or "lastwoyr"

# utilities
SavePlot <- function(plotname, width = 10, height = 7) {
  file <- paste0(getwd(), "/model/", "meg78_fwd_", plotname, "_", Fassumption, "_", Rassumption, ".png")
  dev.print(png, file, width = width, height = height, units = "in", res = 300)
}

# Load assessment

# load("C:/Users/maltuna/OneDrive - AZTI/Oilarra/WGBIE/WGBIE_2025/6.5.Reference points/2024_meg.27.7b-k8abd_assessment_dgarcia/RefpPoints_Sept25/Megrim_2024_EqSim_Workspace.Rdata")
# load('Refpoints_WKMEGRIM_megrim78_run3/Megrim_2020_EqSim_Workspace_CPUEIRLFR10_corrected.RData')
#AI: importante cargar primero el archivo de Refpoints y luego el Fit, porque en este ultimo tenemos todos los años del assessment.

refPts <- data.frame(
  MSYBtrigger = 43946.96,
  Bpa = 43946.96,
  Blim = 35795.40,
  Fpa = 0.4894419,
  Fp05 = 0.4894419,
  Fmsy_unconstr = 0.2692693,
  Fmsy = 0.2692693,
  Fmsyupper_unconstr = 0.4564565,
  Fmsyupper = 0.4564565,
  Fmsylower_unconstr = 0.1631632,
  Fmsylower = 0.1631632
)

load('model/MegFit_FINAL.RData')

run <- stk1 
years <- run@range[4]:run@range[5]
nyears <-length(years)
ages <- run@range[1]:run@range[2]
nages <- length(ages)

# --- SETUP

# refpts: Bpa, Fpa, Blim, BmsyTrig, Fmsy, FmsyLow, FmsyUpp

# MODEL, ADVICE and FORECAST year

dy <- dims(run)$maxyear #last year of assessment input data
ay <- dy + 1
fy <- ay + 1

# TAC & advice current year

advice <- FLQuant(21144, dimnames=list(age='all', year=ay), units="tonnes") #en WGBIE 2024, advice para el año 2025 son 21144 tonas.
tac <- advice


### RECRUITMENT
if(Rassumption == "meanblowmean"){
# Mean of the R below historical mean (historical mean using all years minus last 2)
  gmean <- exp(mean(log(window(stock.n(run)["1",], end=-2))))
  recsq <- mean(window(stock.n(run)["1",], end=-2)[window(stock.n(run)["1",], end=-2) < gmean])
}
if(Rassumption == "gmean"){
# GEOMEAN but last year
#recsq <- exp(mean(log(window(stock.n(run)["1",], end=-1))))

# GEOMEAN using all years minus last 2 (distintas formas de obtener "recsq")
# specity the number of years removing from the end of the time-series
recsq <- exp(mean(log(window(stock.n(run)["1",], end=-2))))
## specity the start and final years
# recsq <- exp(mean(log(window(stock.n(run)["1",], start=1984, end=2022))))

# GEOMEAN: para cambiar el reclutamiento del ultimo año del assessment por geomean si consideramos que tenemos mucha incertidumbre.
#plot(run@stock.n['1',])
#run@stock.n['1','2022'] <- recsq #AI: Para la evalucion con datos hasta 2022.
}
if(Rassumption == "lastwoyr"){
  # GEOMEAN but lasttwo year
  recsq <- exp(mean(log(window(stock.n(run)["1",], start=2023, end=2024))))
}

# --- SETUP future

# 3 years, 5 years wts/selex, 3 years discards
fut <- stf(run, nyears=3, wts.nyears = 5, fbar.nyears=5, disc.nyears=3)

### F status quo (Fsq)

if (Fassumption == "Fscaled"){
  Fsq <- expand(fbar(fut)[, ac(dy)], year=ay) # F SCALED to Fbar of final assessment year
}
if(Fassumption == "Funscaled"){
  Fsq <-expand(yearMeans(fbar(fut)[, ac(seq(dy - 2, dy))]), year= ay) # F UNSCALED (Average F last 3 years)
  # Fsq <-FLCore::expand(yearMeans(fbar(fut)[, ac(seq(dy - 2, dy))]), year= ay) # Be careful not to have package (tidyr)instaled. In case of error use: FLCore::expand
  
}

########################################################################################

# SET geomean SRR
gmsrr <- predictModel(model=rec~a, params=FLPar(c(recsq), units="thousands",
  dimnames=list(params="a", year=seq(ay, length=3), iter=1)))

# > gmsrr   # AI: reemplaza el reclutamiento con recsq del año intermedio (2025) y los dos siguientes.
# An object of class "FLQuants": EMPTY
# model:  
#   rec ~ a
# 
# params:  
#   An object of class "FLPar"
# year
# params    2025      2026      2027 
# a        111160      111160    111160
# units:  thousands 

# GENERATE targets from refpts
refpts <- FLPar(refPts[1,])
targets <- expand(as(refpts, 'FLQuant'), year=fy)

# ##############################################################################
# #AI: added to include a new variable in the output
# # GENERATE targets from refpts
# # need to run the forecast once to see what the SSB at the start of the advice year is: as.numeric(ssb(runs[[1]])[,as.character(fy)])
# 
# as.numeric(ssb(runs[[1]])[,as.character(fy)])
# # 55163.19 unscaled
# 
# 
# # Then enter it manually
# refPts <- cbind(refPts,  55163.19)
# dimnames(refPts)[[2]][ncol(refPts)] <- "constSSB"
# refpts <- FLPar(refPts[1,])
# targets <- expand(as(refpts, 'FLQuant'), year=fy)
# 
# #AI: to complete in the advice the assumptions for interim year and forecast.
# 
# icesRound(fbar(runs$Fmsy)[, '2025']) # 0.24
# icesRound(fbar(runs$Fmsy)[, '2026']) # 0.27
# ssb(runs$Fmsy)[, '2026'] # 55163
# rec(runs$Fmsy)[, '2026'] # 111160
# catch(runs$Fmsy)[, '2025'] # 13684
# landings(runs$Fmsy)[, '2025'] # 11426
# discards(runs$Fmsy)[, '2025'] # 2258
# 
# ssb(runs$Fmsy)[, '2025'] # La SSB del año intermedio (2025) hay que incluirla en la tabla del advice donde se da la summary table pero este ultimo a?o no hay valor (NA).
# # 57378
# catch(runs$Fmsy)[, '2024'] # fitted catch for the last assessment year to put in the report Tab 5.1.1.2
# # 9691
##############################################################################


# --- PROJECT catch options

# Targets
C0 <- FLQuant(0, dimnames=list(age='all', year=fy))
Fiy <- FLQuants(fbar=append(Fsq, C0))
Ciy <- FLQuants(catch=append(tac, C0))

# RUN for Fsq
Fsqrun <- fwd(fut, sr=gmsrr, control=as(Fiy, "fwdControl"))

# TEST if Cay < TACay
if(catch(Fsqrun)[, ac(ay)] <= tac) {
  itarget <- Fiy
} else {
  itarget <- Ciy
  Fsqrun <- fwd(fut, sr=gmsrr, control=as(Ciy, "fwdControl"))
}

# DEFINE catch options
catch_options <- list(

  # FMSY
  Fmsy=FLQuants(fbar=targets["Fmsy",]),

  # lowFMSY
  lFmsy=FLQuants(fbar=targets["Fmsylower_unconstr",]),

  # uppFMSY
  uFmsy=FLQuants(fbar=targets["Fmsyupper",]),

  # F0
  F0=FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=fy))),

  # Fpa
  Fpa=FLQuants(fbar=targets["Fpa",]),

  # # Flim
  # Flim=FLQuants(fbar=targets["Flim",]),

  # Bpa
  Bpa=FLQuants(ssb_flash=targets["Bpa",]),

  # Blim
  Blim=FLQuants(ssb_flash=targets["Blim",]),

  # MSYBtrigger
  MSYBtrigger=FLQuants(ssb_flash=targets["MSYBtrigger",]),

  # F sq
  Fsq=FLQuants(fbar=expand(fbar(Fsqrun)[, ac(ay)], year=fy)),
  
  # SB sq
  SBsq=FLQuants(ssb_flash=expand(ssb(Fsqrun)[, ac(ay+1)], year=fy)),
  
  # TAC
  rotac=FLQuants(catch=expand(tac, year=fy))#,
  # 
  # # 35% quantile rule
  # Quantile=FLQuants(catch=catchQR) # calculated in 7.SSB_vs_Blim.R
)

# C0
F0 <- FLQuants(fbar=FLQuant(0, dimnames=list(age='all', year=fy + 1)))

# CONVERT to fwdControl

fctls <- lapply(catch_options, function(x) {
  as(FLQuants(c(window(itarget, end=ay), x, F0)), "fwdControl")
}
)

# RUN!

runs <- FLStocks(lapply(fctls, function(x) fwd(fut, sr=gmsrr, control=x)))

# SSB




# COMPARE

Map(compare, runs, fctls)


# --- PROJECT F levels

flevels <- seq(0, 0.50, 0.01)

control <- as(as(c(lapply(window(itarget, end=ay), propagate, length(flevels)),
  FLQuants(fbar=FLQuant(flevels,  dimnames=list(year=fy, iter=seq(length(flevels))))),
  lapply(F0, propagate, length(flevels))), "FLQuants"), "fwdControl")

f_runs <- divide(fwd(fut, sr=gmsrr, control=control), names=flevels)

# Contribution of each age to the landing and SSB
par(mfrow=c(1,2),mar=c(5,8,4,1),cex=0.8)
yield <- c(runs$Fmsy@landings.n[,nyears+2])
prop <- paste0(round(100*yield/sum(yield)),'%')
labels <- paste0(max(years)-ages+2,rep(c(' (GM)',' (a4a)'),c(2,nages-2))) #if you dont replace last year recr only 2 ages GM
b <- barplot(yield,horiz=T,names=labels,las=1,xlab='Tonnes',main=paste('Landings yield',max(years)+2),xlim=c(0,max(yield)*1.25))
text(yield,b,prop,adj=-0.2)
mtext('Cohort',2,6)

ssb <- c((runs$Fmsy@stock.n*runs$Fmsy@stock.wt*runs$Fmsy@mat)[,nyears+3])
prop <- paste0(round(100*ssb/sum(ssb)),'%')
labels <- max(years)-ages+3
b <- barplot(ssb,horiz=T,names=labels,las=1,xlab='Tonnes',main=paste('SSB',max(years)+3),xlim=c(0,max(ssb)*1.25))
text(ssb,b,prop,adj=-0.2)
mtext('Cohort',2,4)

SavePlot('Contrib')

# SAVE
save(runs, f_runs, recsq, tac, advice, file=paste0("model/runs_", Fassumption, "_", Rassumption, ".RData"))
stk_est_stf <- runs[["Fmsy"]]
save(stk_est_stf, file = "model/Meg78_WGMIX_stock_estimated_stf.Rdata")
