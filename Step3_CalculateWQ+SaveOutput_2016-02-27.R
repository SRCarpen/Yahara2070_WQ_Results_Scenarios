# Calculate total loads and water quality for four lakes given direct drainage loads
# Include hypereutrophy estimate based on GLM
# Adapted for Factorial results
# SRC 2016-02-27

rm(list = ls())
graphics.off()

# Functions +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Estimate Chl from TP using the Filstrup model
# See Filstrup_Demo+MendotaData_2015-02-08.R for validation
# Paper is:
# Filstrup, C.T., T. Wagner, P.A. Soranno, E.H. Stanley,
# C.A. Stow, K.E. Webster and J.A. Downing. 2014.
# Regional variability among nonlinear chlorophyll-phosphorus
# relationships in lakes. Limnol. Oceanogr.59: 1691-1703.
filstrup = function(TP) {
  # Parameters from Table 3 of Filstrup et al.
  A = 2.14 # upper asymptote
  B = 0.49 # inflection point
  C = 3.21 # rate of increase
  D = 0.08 # lower asymptote
  TP.GM = 1.1287 # grand mean log10(TP) from Filstrup 2015-02-08
  lTP = log10(1000*TP) # log10 of TP in ug/L
  lChl = D + ( (A-D)/(1 + exp(-C*(lTP-TP.GM-B))) )  # Equation 1
  Chl = 10^lChl  # Chl in ug/L
  return(Chl)
}

# Forward transform
# Parameter must lie between 0 and Pmax; transform to unbounded parameter for fitting
FT = function(P,Pmax) {
  Ps = tan( (pi*P/Pmax) - (pi/2))
  return(Ps)
}

# Back transform
# Unbounded parameter to a value between 0 and Pmax
BT = function(Ps,Pmax) {
  P = (Pmax/pi)*( (pi/2) + atan(Ps))
}

# Function to Simulate P and export through one year
# This version does book-keeping on P0 each year, rather
# than starting over at a new value each year.
Ppred.v2 = function(par,P0,L,Pmax,dt,Nstep,NY) {
  S = BT(par[1],Pmax)
  H = BT(par[2],Pmax)
  W = BT(par[3],Pmax)
  P1hat = rep(0,NY)
  Ehat = rep(0,NY)
  # Compute estimated P1 and Export for each year
  for(iy in 1:NY) {
    Pstart = ifelse(iy==1,P0,P1hat[iy-1])
    Pt = rep(Pstart,Nstep)
    Et = rep(0,Nstep)
    # Integrate by Euler
    for(i in 2:Nstep){
      #Xt = H*(Pt[i-1] + W*L[iy])
      Xt = H*Pt[i-1] + W*L[iy]
      Et[i] = Et[i-1] + Xt*dt
      Pt[i] = Pt[i-1] + (L[iy]-S*Pt[i-1]- Xt)*dt
    } # end integration loop within one year
    P1hat[iy] = Pt[Nstep]
    Ehat[iy] = Et[Nstep]
  }  # end loop over years
  outlist = list(P1hat,Ehat)
  return(outlist)
}

# END OF FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Constants for all scenarios 

# Set Daphnia pulicaria to absent (0) or present (1) for Mendota and Monona
# Notes: (1) Presence/absence of Daphnia pulicaria will depend on time and WSC scenario
#        (2) Effect of Daphnia pulicaria has not been measured in Waubesa & Kegonsa
Dpul.Me = 1 # Daphnia pulicaria in Mendota
Dpul.Mo = 1 # Daphnia pulicaria in Monona

# Toggle stochastic within-year effects off (0) or on (1)
stoch = 0

# Constants for within-year simulation of P dynamics
Pmax = 1 # maximum value of balance model parameters, for the transform function
Nstep = 30 # number of time steps per year
dt = 1/Nstep  # time step for within-year dynamics


# <<<<<<<<<<<<<<<<
# Read Factorial data from Step 2
# AllVars is all data arranged by lake, scenario and year 
# The lake code is 1 = Mendota, 2 = Monona, 3 = Waubesa, 4 = Kegonsa
# The driver scenario code is 1 = AI, 2 = AR, 3 = CC, 4 = NW
# For definitions of the columns see the code for Mendota in Step 2
#
# Useful submatrices are as follows:
# LandUse is the land use variables only
# LandMgmt is the P management variables only (Manure + Fert are NOT areal)
# Snow is the snow cover variables only
# Runoff is the water runoff variables only
# SoilP is soil P content
# Pyield is the P runoff variables
# DDL is direct drainage load
# DDLareal is direct drainage load per ha
# LakeID is lake ID number (1 for Mendota through 4 for Kegonsa)
# Precip is the precip by scenario and year for the entire Yahara (not by lake)
# Precip4 is the precip matrix stacked 4X to account for all 4 lakes
# DDarea is the area in ha contributing to direct drainage
#save(AllVars,LandUse,LandMgmt,Snow,Runoff,SoilP,Pyield,DDL,DDLareal,LakeID,Precip,Precip4,DDarea,file='AllVarsByLake.Rdata')
load(file='AllVarsByLake.Rdata')

# Remove matrices that will be rebuilt after WQ calculation
rm(LandUse,LandMgmt,Snow,Runoff,SoilP,Pyield)

# Reorganize AllVars into 4 lake sub-matrices
DDL.Me = subset(DDL,subset=(LakeID==1))
DDL.Mo = subset(DDL,subset=(LakeID==2))
DDL.Wa = subset(DDL,subset=(LakeID==3))
DDL.Ke = subset(DDL,subset=(LakeID==4))

# >>>>>>>>>>>>>>>>

# Load the data files needed to complete the calculations for the four lakes ****************

# Load Water Quality Regressions
# Predictors for Log summer TP:
#   Mendota input and Dpul
#   Monona output and Dpul
#   Waubesa and Kegonsa: outputs
# Predictors for Log summer Secchi transparency:
#   Mendota input and Dpul
#   Monona output and Dpul
#   Waubesa and Kegonsa P1 (NOvember 1 P after the growing season)
#save(LTPreg.Me,LTPreg.Mo,LTPreg.Wa,LTPreg.Ke,
#     LSecreg.Me,LSecreg.Mo,LSecreg.Wa,LSecreg.Ke,
#     file = 'WaterQualityRegressions.Rdata')
load(file = 'WaterQualityRegressions.Rdata')

# Parameters from 'BalanceModelFit_Mendota+2sigmas.R'
# par.t.Me is the transformed parameters needed to simulate the model
# P0bar.Me is the mean P0
# Shat.Me is the sedimentation coefficient
# Hhat.Me is the hydrologic outflow coefficient for P0
# What.Me is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Me,P0bar.Me,Shat.Me,Hhat.Me,What.Me,sigmaP.Me,sigmaX.Me, file='Pbalance.Me.Rdata')
load(file='Pbalance.Me.Rdata')

# Parameters from 'BalanceModelFit_Monona+2sigmas.R'
# par.t.Mo is the transformed parameters needed to simulate the model
# P0bar.Mo is the mean P0
# Shat.Mo is the sedimentation coefficient
# Hhat.Mo is the hydrologic outflow coefficient for P0
# What.Mo is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Mo,P0bar.Mo,Shat.Mo,Hhat.Mo,What.Mo,sigmaP.Mo,sigmaX.Mo, file='Pbalance.Mo.Rdata')
load(file='Pbalance.Mo.Rdata')

# Parameters from 'BalanceModelFit_Waubesa+1s+2sigmas.R'
# par.t.Wa is the transformed parameters needed to simulate the model
# P0bar.Wa is the mean P0
# Shat.Wa is the sedimentation coefficient
# Hhat.Wa is the hydrologic outflow coefficient for P0
# What.Wa is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Wa,P0bar.Wa,Shat.Wa,Hhat.Wa,What.Wa,sigmaP.Wa,sigmaX.Wa, file='Pbalance.Wa.Rdata')
load(file='Pbalance.Wa.Rdata')

# Parameters from 'BalanceModelFit_Kegonsa+1s+2sigmas.R'
# par.t.Ke is the transformed parameters needed to simulate the model
# P0bar.Ke is the mean P0
# Shat.Ke is the sedimentation coefficient
# Hhat.Ke is the hydrologic outflow coefficient for P0
# What.Ke is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Ke,P0bar.Ke,Shat.Ke,Hhat.Ke,What.Ke,sigmaP.Ke,sigmaX.Ke, file='Pbalance.Ke.Rdata')
load(file='Pbalance.Ke.Rdata')

# Load GLM model for probability of hypereutrophy based on log total P
# See DRP_GLM-vs-TP_AllLakesTogether_2016-01-13.R
# save(GLM.DRP,file='GLM.TPtoDRP.Rdata')
load(file='GLM.TPtoDRP.Rdata')
HypEut.b0 = GLM.DRP$coefficients[1]
HypEut.b1 = GLM.DRP$coefficients[2]

# End loading of the additional data files ********************************************************

# Start loop over all 16 simulations  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# DDloads is NY*16 matrix of direct drainage loads
# Scenario sequence AR, AI, NW, CC; lake sequence Mo, Me, Wa, Ke
# Names in simulation code:
# AnnLoad.Me is load to Mendota (both total and direct drainage)
# DDload.Mo is direct drainage load to Monona
# DDload.Wa is direct drainage load to Waubesa
# DDload.Ke is direct drainage load to Kegonsa

# Matrices for output
# Columns for each lake will be TP, Chl, Secchi, Hyper
NY.per.cell = 57
NY.all = NY.per.cell*16
cnames=c('TotLoad','TP','Chl','Secchi','Hyper','Pexport')
WQ.Me = matrix(0,nr=NY.all,nc=6,dimnames=list(NULL,cnames))
WQ.Mo = matrix(0,nr=NY.all,nc=6,dimnames=list(NULL,cnames))
WQ.Wa = matrix(0,nr=NY.all,nc=6,dimnames=list(NULL,cnames))
WQ.Ke = matrix(0,nr=NY.all,nc=6,dimnames=list(NULL,cnames))

# Note NY in original WQ program is the same as NY.per.cell
NY = NY.per.cell

for(iscen in 1:16)  { # start loop over all 4 scenarios
  # lake indices for direct drainage loads
  i0 = NY.per.cell*(iscen-1)+1 # start
  i1 = i0+NY.per.cell-1
  # Direct drainage loads for scenario
  AnnLoad.Me = DDL.Me[i0:i1]
  DDload.Mo = DDL.Mo[i0:i1]
  DDload.Wa = DDL.Wa[i0:i1]
  DDload.Ke = DDL.Ke[i0:i1]

# Compute Mendota Budget and WQ *******************************************************************

# Simulate within-year dynamics to get Nov P and Export
Sim.Me = Ppred.v2(par.t.Me,P0bar.Me,AnnLoad.Me,Pmax,dt,Nstep,NY)
PNovsim.Me = Sim.Me[[1]]
XPsim.Me = Sim.Me[[2]]
# Save annual load
WQ.Me[i0:i1,1] = AnnLoad.Me
# Save Export
WQ.Me[i0:i1,6] = XPsim.Me

# Compute summer TP
LTPb0 = LTPreg.Me$coefficients[1] # intercept
LTPb1 = LTPreg.Me$coefficients[2] # effect of input
LTPb2 = LTPreg.Me$coefficients[3] # effect of D. pul.
LTPdf = LTPreg.Me$df.residual
LTPse = sd(LTPreg.Me$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
LTPdet.Me = LTPb0 + LTPb1*AnnLoad.Me + LTPb2*Dpul.Me
TPsummer.Me = exp( LTPb0 + LTPb1*AnnLoad.Me + LTPb2*Dpul.Me + noise )
WQ.Me[i0:i1,2] = TPsummer.Me

# Compute summer Chl
Chlsummer.Me = filstrup(TPsummer.Me)
WQ.Me[i0:i1,3] = Chlsummer.Me

# Compute summer Secchi transparency
LSecb0 = LSecreg.Me$coefficients[1] # intercept
LSecb1 = LSecreg.Me$coefficients[2] # effect of input
LSecb2 = LSecreg.Me$coefficients[3] # effect of D. pul.
LSecdf = LSecreg.Me$df.residual
LSecse = sd(LSecreg.Me$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Me = exp( LSecb0 + LSecb1*AnnLoad.Me + LSecb2*Dpul.Me + noise)
WQ.Me[i0:i1,4] = Secsummer.Me

# Compute probability of hypereutrophy
HypEut.Me = 1/(1 + exp(-(HypEut.b0 + HypEut.b1*LTPdet.Me)))
WQ.Me[i0:i1,5] = HypEut.Me

# End Mendota, Compute Monona Budget and WQ ***********************************

# Annual load is Direct Drainage plus river transfer
AnnLoad.Mo = DDload.Mo + XPsim.Me
WQ.Mo[i0:i1,1] = AnnLoad.Mo # save annual load

# Run within-year simulations to get P in Nov and Export
# Simulate within-year dynamics to get Nov P and Export
Sim.Mo = Ppred.v2(par.t.Mo,P0bar.Mo,AnnLoad.Mo,Pmax,dt,Nstep,NY)
PNovsim.Mo = Sim.Mo[[1]]
XPsim.Mo = Sim.Mo[[2]]

# Save Export
WQ.Mo[i0:i1,6] = XPsim.Mo

# Calculate summer TP using regression
LTPb0 = LTPreg.Mo$coefficients[1] # intercept
LTPb1 = LTPreg.Mo$coefficients[2] # effect of output
LTPb2 = LTPreg.Mo$coefficients[3] # effect of D. pul.
LTPdf = LTPreg.Mo$df.residual
LTPse = sd(LTPreg.Mo$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
LTPdet.Mo = LTPb0 + LTPb1*XPsim.Mo + LTPb2*Dpul.Mo
TPsummer.Mo = exp( LTPb0 + LTPb1*XPsim.Mo + LTPb2*Dpul.Mo + noise )
WQ.Mo[i0:i1,2] =TPsummer.Mo

# Compute summer Chl
Chlsummer.Mo = filstrup(TPsummer.Mo)
WQ.Mo[i0:i1,3] = Chlsummer.Mo

# Compute summer Secchi transparency
LSecb0 = LSecreg.Mo$coefficients[1] # intercept
LSecb1 = LSecreg.Mo$coefficients[2] # effect of input
LSecb2 = LSecreg.Mo$coefficients[3] # effect of D. pul.
LSecdf = LSecreg.Mo$df.residual
LSecse = sd(LSecreg.Mo$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Mo = exp( LSecb0 + LSecb1*XPsim.Mo + LSecb2*Dpul.Mo + noise)
WQ.Mo[i0:i1,4] = Secsummer.Mo

# Compute probability of hypereutrophy
HypEut.Mo = 1/(1 + exp(-(HypEut.b0 + HypEut.b1*LTPdet.Mo)))
WQ.Mo[i0:i1,5] = HypEut.Mo

# End Monona, Compute Waubesa Budget and WQ ***********************************

# Annual load is Direct Drainage plus river transfer
AnnLoad.Wa = DDload.Wa + XPsim.Mo
WQ.Wa[i0:i1,1] = AnnLoad.Wa  # save annual load

# Run within-year simulations to get P in Nov and Export
# Simulate within-year dynamics to get Nov P and Export
Sim.Wa = Ppred.v2(par.t.Wa,P0bar.Wa,AnnLoad.Wa,Pmax,dt,Nstep,NY)
PNovsim.Wa = Sim.Wa[[1]]
XPsim.Wa = Sim.Wa[[2]]

# Save Export
WQ.Wa[i0:i1,6] = XPsim.Wa

# Calculate summer TP using regression
LTPb0 = LTPreg.Wa$coefficients[1] # intercept
LTPb1 = LTPreg.Wa$coefficients[2] # effect of output
LTPdf = LTPreg.Wa$df.residual
LTPse = sd(LTPreg.Wa$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
LTPdet.Wa = LTPb0 + LTPb1*XPsim.Wa
TPsummer.Wa = exp( LTPb0 + LTPb1*XPsim.Wa + noise)
WQ.Wa[i0:i1,2] = TPsummer.Wa

# Compute summer Chl
Chlsummer.Wa = filstrup(TPsummer.Wa)
WQ.Wa[i0:i1,3] = Chlsummer.Wa

# Calculate summer Secchi transparency
LSecb0 = LSecreg.Wa$coefficients[1] # intercept
LSecb1 = LSecreg.Wa$coefficients[2] # effect of P1 (Nov. P mass)
LSecdf = LSecreg.Wa$df.residual
LSecse = sd(LSecreg.Wa$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Wa = exp( LSecb0 + LSecb1*PNovsim.Wa + noise)
WQ.Wa[i0:i1,4] = Secsummer.Wa

# Compute probability of hypereutrophy
HypEut.Wa = 1/(1 + exp(-(HypEut.b0 + HypEut.b1*LTPdet.Wa)))
WQ.Wa[i0:i1,5] = HypEut.Wa

# End Waubesa, Compute Kegonsa Budget and TP ***********************************

# Annual load is Direct Drainage plus river transfer
AnnLoad.Ke = DDload.Ke + XPsim.Wa
WQ.Ke[i0:i1,1] = AnnLoad.Ke  # save annual load

# Run within-year simulations to get P in Nov and Export
# Simulate within-year dynamics to get Nov P and Export
Sim.Ke = Ppred.v2(par.t.Ke,P0bar.Ke,AnnLoad.Ke,Pmax,dt,Nstep,NY)
PNovsim.Ke = Sim.Ke[[1]]
XPsim.Ke = Sim.Ke[[2]]

# Save Export
WQ.Ke[i0:i1,6] = XPsim.Ke

# Calculate summer TP using regression
LTPb0 = LTPreg.Ke$coefficients[1]
LTPb1 = LTPreg.Ke$coefficients[2]
LTPdf = LTPreg.Ke$df.residual
LTPse = sd(LTPreg.Ke$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
LTPdet.Ke = LTPb0 + LTPb1*XPsim.Ke
TPsummer.Ke = exp( LTPb0 + LTPb1*XPsim.Ke + noise)
WQ.Ke[i0:i1,2] = TPsummer.Ke

# Compute summer Chl
Chlsummer.Ke = filstrup(TPsummer.Ke)
WQ.Ke[i0:i1,3] = Chlsummer.Ke

# Calculate summer Secchi transparency
LSecb0 = LSecreg.Ke$coefficients[1] # intercept
LSecb1 = LSecreg.Ke$coefficients[2] # effect of P1 (Nov. P mass)
LSecdf = LSecreg.Ke$df.residual
LSecse = sd(LSecreg.Ke$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Ke = exp( LSecb0 + LSecb1*PNovsim.Ke + noise)
WQ.Ke[i0:i1,4] = Secsummer.Ke

# Compute probability of hypereutrophy
HypEut.Ke = 1/(1 + exp(-(HypEut.b0 + HypEut.b1*LTPdet.Ke)))
WQ.Ke[i0:i1,5] = HypEut.Ke

# END COMPUTATION OF BUDGETS AND WATER QUALITY FOR THE FOUR LAKES ************************
} # End loop over scenarios >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Merge WQ with lake-specific driver matrices ---------------------------------------------

Driver.Me = subset(AllVars,subset=(LakeID==1))
Driver.Mo = subset(AllVars,subset=(LakeID==2))
Driver.Wa = subset(AllVars,subset=(LakeID==3))
Driver.Ke = subset(AllVars,subset=(LakeID==4))

DriveWQ.Me = cbind(Driver.Me,WQ.Me)
DriveWQ.Mo = cbind(Driver.Mo,WQ.Mo)
DriveWQ.Wa = cbind(Driver.Wa,WQ.Wa)
DriveWQ.Ke = cbind(Driver.Ke,WQ.Ke)

AllDriverWQ = rbind(DriveWQ.Me,DriveWQ.Mo,DriveWQ.Wa,DriveWQ.Ke)

# Make some useful submatrices
Year = AllDriverWQ[,1]
LandUse = AllDriverWQ[,2:9]
LandMgmt = AllDriverWQ[,10:21]
Snow=AllDriverWQ[,22:24]
Runoff=AllDriverWQ[,25:28]
Pyield = AllDriverWQ[,30:32]

# WQ outputs
DDL = AllDriverWQ$DDL
DDLareal = AllDriverWQ$DDLareal
LakeID = AllDriverWQ$LakeID
LULM = AllDriverWQ$LULM
CLIM = AllDriverWQ$CLIM
TotLoad = AllDriverWQ$TotLoad
TP = AllDriverWQ$TP
Chl = AllDriverWQ$Chl
Secchi = AllDriverWQ$Secchi
Hyper = AllDriverWQ$Hyper
Pexport = AllDriverWQ$Pexport
WQoutput = cbind(LULM,CLIM,LakeID,Year,DDL,DDLareal,TotLoad,TP,Chl,Secchi,Hyper,Pexport)

# Make intuitive column names 

colnames(LandUse) = c('alfalfa','corn','soy','haypasture','grass','forest','wetland','devel')

colnames(LandMgmt) = c('Manure','Fert','MeanMan','MedMan','Man25','Man75','Man90',
                       'MeanFert','MedFert','Fert25','Fert75','Fert90')
  
colnames(Snow) = c('AnnSno','AvSnoMo','MaxSnoMo')

SoilP = AllDriverWQ$SoilP

# Save output created by Step 3
# AllDriverWQ is all data arranged by lake, scenario factorial cell, and year 
# The lake code is 1 = Mendota, 2 = Monona, 3 = Waubesa, 4 = Kegonsa
# The driver scenario code for LULM and CLIM is 1 = AI, 2 = AR, 3 = CC, 4 = NW
# For definitions of the columns see the code
#
# Useful submatrices for regressions:
# Precip is the precipitation variables only
# Precip4 is the precipitation variables stacked to conform to all 4 lakes
# LandUse is the land use variables only
# LandMgmt is the P management variables only (Manure + Fert are NOT areal)
# Snow is the snow cover variables only
# Runoff is the water runoff variables only
# SoilP is soil P content
# Pyield is the P runoff variables
# DDL is direct drainage load
# DDLareal is direct drainage load per ha
# WQ outputs are Total Load, TP mg/L, Chl ug/L, Secchi m, prob(measurable DRP), and P export
# Note Hyper is prob(DRP present) i.e. hypereutrophy and 1-Hyper is prob(DRP absent) i.e. eutrophy
# LakeID is lake ID number (1 for Mendota through 4 for Kegonsa)
# DDarea is the area in ha contributing to direct drainage
#
# Remember to rename data file as - or + Daphnia
save(AllDriverWQ,Precip,Precip4,LandUse,LandMgmt,Snow,Runoff,SoilP,Pyield,
     DDL,DDLareal,LakeID,LULM,CLIM,DDarea,TotLoad,TP,Chl,Secchi,Hyper,WQoutput,
     file='AllDrivers+WQ.Rdata')

# Save WQ in R plus Excel / csv formats
# WQoutput = cbind(LULM,CLIM,LakeID,Year,DDL,DDLareal,TotLoad,TP,Chl,Secchi,Hyper,Pexport)
# LULM and CLIM are scenarios for land use / management and climate drivers, respectively
# The driver scenario code for LULM and CLIM is 1 = AI, 2 = AR, 3 = CC, 4 = NW
# The LakeID code is 1 = Mendota, 2 = Monona, 3 = Waubesa, 4 = Kegonsa
# DDL is direct drainage load kg/y
# DDLareal is direct drainage load per unit area, kg/(ha*y)
# TotLoad is total (direct drainage + upstream lakes) load, kg/y
# TP is total P in summer mg/L
# Chl is chlorophyll a in summer ug/L
# Secchi is transparency in m
# Hyper is probability of hypereutrophy
# Pexport is outflow of P in kg/y
#
# Remember to rename data files as - or + Daphnia!
# R data file
save(WQoutput,file='FactorialWQ.Rdata') 
# CSV file
write.csv(WQoutput,file='FactorialWQ.csv')
