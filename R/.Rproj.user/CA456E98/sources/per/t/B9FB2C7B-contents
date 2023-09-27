# setup script for the jittering project
# install required packages if necessary
if(FALSE) {
  install.packages("Matrix")
  install.packages("spam")
  install.packages("fields")
  install.packages("latex2exp")
  install.packages("xtable")
  install.packages("profvis")
  install.packages("viridis")
  install.packages("numDeriv")
  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  install.packages("data.table")
  install.packages("abind")
  install.packages("devtools")
  install.packages("ggplot2")
  install.packages("colorspace")
  install.packages("TMB")
  install.packages("plyr")
  install.packages("sf")
  install.packages("ncdf4")
  install.packages("pracma")
  install.packages("concaveman")
  install.packages("fmesher")
}

# load required packages and R scripts
library(Matrix)
library(spam)
library(fields)
library(latex2exp)
library(xtable)
library(profvis)
library(viridis)
library(numDeriv)
library(data.table)
library(abind)
library(devtools)
library(ggplot2)
library(colorspace)
library(TMB)
library(plyr)
#library(INLA)
library(sf)
library(ncdf4)
library(pracma)
library(concaveman)
library(devtools)
#library(fmesher)

#inf = sessionInfo()
#if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
#  setwd("~/git/jittering/")
#  options(error=recover)
#} else if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)" && inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
#  INLA:::inla.dynload.workaround()
  # avoid setting too many threads and thereby using too much memory
#  inla.setOption(num.threads=1)
#  options(error=traceback)
#  setwd("~/git/csz/")
#} else if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
#  setwd("~/git/csz/")
#  options(error=recover)
#} else {
#  setwd("~/git/csz/")
#  inla.setOption(num.threads=1) # consider raising
 # options(error=recover)
#}
#options(warn=1)

#setwd("~/git/csz")
#source("code/loadSubDat.R")
#source("code/slab.R")
#source("code/okada.R")
#source("code/plotter.R")
#source("code/genericSpatialPlottingFunctions.R")
#source("code/utilityFuns.R")
#source("code/test.R")

## load in global variables made from the following script: 
#if(FALSE) {
#  source('code/preprocess.R')
#}

# load("savedOutput/global/datMICS.RData")
#out = loadSubsidenceDat()
#allSubDat <- out$allSubDat
#majorSubDat <- out$majorSubDat

# determine version of PROJ for projections
#ver = terra::gdal(lib="proj")
#PROJ6 <- as.numeric(substr(ver, 1, 1)) >= 6



