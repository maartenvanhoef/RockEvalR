#Clear environment if needed
rm(list=ls())

#Install packages if needed
install.packages("tidyverse")
devtools::install_github("maartenvanhoef/RockEvalR")

#Load packages into library
library(tidyverse)
library(RockEvalR)

#Location of your RE data
dataloc<-"e.g. C:/.../BULK ROCK"

#Load all data files in folder
REdata<-RE_read(dataloc)

#Plot first data file
RE_plot(REdata,1)

## Choices of interpretation and adaptation ----

#1 Convert raw RE data to values with proper units
#  This step should always be done first to work with the proper units
REdata.adj<-RE_convert(REdata)

# Remove unconverted data (if no longer needed)
rm(REdata)

# After this everything could be done/ added to the data individually

#2 Calculate RE zones (S1, S2, ...)
REdata.adj<-RE_surfaces(REdata.adj)

#3 Compute the RE metrics (TOC, OI, ...)
REdata.adj<-RE_metrics(REdata.adj)

#4 Determine the combined C flux in each oven
REdata.adj<-RE_Ccombined(REdata.adj)

#5 Determine Sebag et al. (2016) parameters
REdata.adj<-RE_SebagIR(REdata.adj)

#6 Calculate T percentiles
REdata.adj<-RE_Tpercentiles(REdata.adj)

# Standard is T25, T50 and T75 but any custom vector can be used as:
REdata.adj<-RE_Tpercentiles(REdata.adj, c(10,30,50,70,90))

## Extracting data from the list ----

# All data can be extracted at once with
table<-RE_extracttable(REdata.adj)

# Specific data can be extracted with apply (or a for loop), e.g.:

#Extract zones
REzones<-as.data.frame(t(sapply(REdata.adj,
                                function(sample)sample[["Zones_C"]])))

#Extract data of single pyrogram
REpyro<-REdata.adj[[1]][["Pyrolysis"]]
par(mfrow = c(1,1))
plot(REpyro$t,REpyro$CO2,
     xlab="Time (s)", ylab="Pyro CO2 mg/g/s", main="example plot", type="l")


