# Clear environment if needed
rm(list=ls())

# Install packages if needed
install.packages("tidyverse")
install.packages("devtools")

library(devtools)
devtools::install_github("maartenvanhoef/RockEvalR")

# Load packages into library
library(tidyverse)
library(RockEvalR)

# Access help page
??RockEvalR

#1 Location of your RE data
dataloc<-"e.g. C:/.../BULK ROCK"

#2 Load all data files in folder
REdata<-RE_read(dataloc)

#3.1 Plot first data file
RE_plot(REdata,1)

#3.2 Plot specific thermogram
RE_plot(REdata,"name.B00","Pyrolysis","CO2")


## Choices of interpretation and adaptation ----

#4 Convert raw RE data to values with proper units
#  This step should always be done first to work with the proper units!
REdata.adj<-RE_convert(REdata)

#5 Adjust the pyrolysis CO2 cursor 2 to the valley between OC and IC peaks (only for SOIL and SOIL TS cycles)
#  By default this cursor is fixed at 400C, this is an optional adjustment to the minimum
REdata.adj<-RE_cursadjust(REdata, plot = F)

#6 Remove unconverted data (if no longer needed)
rm(REdata)

# After this everything could be done/ added to the data individually

#7 Calculate RE zones (S1, S2, ...)
REdata.adj<-RE_surfaces(REdata.adj)

#8 Compute the RE metrics (TOC, OI, ...)
REdata.adj<-RE_metrics(REdata.adj)

#9 Determine the combined C flux in each oven
REdata.adj<-RE_Ccombined(REdata.adj)

#10 Determine Sebag et al. (2016) parameters
REdata.adj<-RE_SebagIR(REdata.adj)

#11.1 Calculate T percentiles
REdata.adj<-RE_Tpercentiles(REdata.adj)

#11.2 Standard percentiles are T25, T50 and T75 but any custom vector can be used as:
REdata.adj<-RE_Tpercentiles(REdata.adj, c(10,30,50,70,90))

## Extracting data from the list ----

#12 All data can be extracted into a neat table at once with
table<-RE_extracttable(REdata.adj)



#13 Specific data can be extracted manually with apply (or a for loop), e.g.:

#Extract zones
REzones<-as.data.frame(t(sapply(REdata.adj,
                                function(sample)sample[["Zones_C"]])))

#Extract data of single pyrogram
REpyro<-REdata.adj[[1]][["Pyrolysis"]]
par(mfrow = c(1,1))
plot(REpyro$t,REpyro$CO2,
     xlab="Time (s)", ylab="Pyro CO2 mg/g/s", main="example plot", type="l")


# For additional information see help pages or Github
