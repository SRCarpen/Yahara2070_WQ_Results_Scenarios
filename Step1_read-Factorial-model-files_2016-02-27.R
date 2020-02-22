# R code to read the water quality output files
# SRC 2016-02-27

rm(list = ls())
graphics.off()

# Function to make a matrix for one scenario combination
# Argument 'names' contains csv file name, land use scenario, climate scenario
# If you want a numeric matrix, rather than a data frame, save Xraw separately
MakeFrame = function(names) {
  cname=(1:171)
  Raw = read.csv(names[1],header=F,col.names=cname)
  Xraw = as.matrix(Raw)
  # Get size of input file
  dimRaw=dim(Raw)
  nr = dimRaw[1]
  nc = dimRaw[2]
  # Land use and climate columns
  Luse=rep(names[2],nr)
  Clim=rep(names[3],nr)
  X.all = data.frame(Xraw,Luse,Clim)
  return(X.all) # or return(Xraw) if you want only the numeric data
}


# <<<<<<<<<<<<<<<<
# Read the regressor for each scenario combination, from Missy Motew
# e.g. AINW is LU variables of AI combined with climate of NW
# see readme file for a list of variables

names=c('AIAI.csv','AI','AI')
X.AIAI = MakeFrame(names)

names=c('AIAR.csv','AI','AR')
X.AIAR = MakeFrame(names)

names=c('AICC.csv','AI','CC')
X.AICC = MakeFrame(names)

names=c('AINW.csv','AI','NW')
X.AINW = MakeFrame(names)

# End AI series 

names=c('ARAI.csv','AR','AI')
X.ARAI = MakeFrame(names)

names=c('ARAR.csv','AR','AR')
X.ARAR = MakeFrame(names)

names=c('ARCC.csv','AR','CC')
X.ARCC = MakeFrame(names)

names=c('ARNW.csv','AR','NW')
X.ARNW = MakeFrame(names)

# End AR series

names=c('CCAI.csv','CC','AI')
X.CCAI = MakeFrame(names)

names=c('CCAR.csv','CC','AR')
X.CCAR = MakeFrame(names)

names=c('CCCC.csv','CC','CC')
X.CCCC = MakeFrame(names)

names=c('CCNW.csv','CC','NW')
X.CCNW = MakeFrame(names)

# End CC series

names=c('NWAI.csv','NW','AI')
X.NWAI = MakeFrame(names)

names=c('NWAR.csv','NW','AR')
X.NWAR = MakeFrame(names)

names=c('NWCC.csv','NW','CC')
X.NWCC = MakeFrame(names)

names=c('NWNW.csv','NW','NW')
X.NWNW = MakeFrame(names)

# >>>>>>>>>>>>>>>>

# Concatenate all the data frames
XC = rbind(X.AIAI,X.AIAR,X.AICC,X.AINW,
           X.ARAI,X.ARAR,X.ARCC,X.ARNW,
           X.CCAI,X.CCAR,X.CCCC,X.CCNW,
           X.NWAI,X.NWAR,X.NWCC,X.NWNW)

# Save the composite dataframe for later analysis
save(XC,file='Factorial_All_V2.Rdata')

