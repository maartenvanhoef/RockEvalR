#' Trapezoidal integration
#'
#' This function does a simplified trapezoidal integration.
#'
#' @param x Independent variable for integral (e.g. time)
#' @param y Dependent variable of same length as x
#' @return Trapezoidal area under time series data
#' @export
RE_traparea<-function(x,y){
  # Area is calculated manually (alternatives did not perform better).

  #1.1 Determine time steps
  dx<-diff(x)

  #1.2 Determine mean values between each step
  my<-(head(y, -1) + tail(y, -1)) / 2

  #1.3 Integral is approximated by the sum of (y_n + y_n+1)/2 x dx
  area<-sum(my*dx)

area
}


#' Cumulative trapezoidal integration
#'
#' This function does a step wise trapezoidal integration.
#'
#' @param x Independent variable for integral (e.g. time)
#' @param y Dependent variable of same length as x
#' @return Trapezoidal area under time series data
#' @export
RE_cumtraparea<-function(x,y){
  # Area is calculated manually (alternatives did not perform better).

  #1.1 Determine time steps
  dx<-diff(x)

  #1.2 Determine mean values between each step
  my<-(head(y, -1) + tail(y, -1)) / 2

  #1.3 Integral is approximated by the sum of (y_n + y_n+1)/2 x dx
  area<-cumsum(my*dx)

  area
}


#' Convert Rock-Eval signal to right units
#'
#' This function derives the elemental amount using calibration, baseline and sample weight.
#'
#' @param list List with raw Rock-Eval data from RE_read
#' @return Input list with Rock-Eval thermograms in proper units (mg/g)
#' @export
RE_convert<-function(list){
  #1.1 Create new converted list as copy of input list
  list.converted<-list

  #2 Pyrolysis CH thermogram

  #2.1 Take baseline from measurement values using lapply
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["CH"]]<-sample[["Pyrolysis"]][["CH"]]-sample[["Cursors"]]["base1"]
    sample})

  #2.2 Apply the Geoworks calibration for the FID
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["CH"]]<-
      sample[["Pyrolysis"]][["CH"]]*sample[["Calibration"]]["cal_FID"]
    sample})

  #2.3 Sample weight normalisation for CH (mg/g)
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["CH"]]<-
      sample[["Pyrolysis"]][["CH"]]/as.numeric(sample[["Parameters"]]["Quant"])*1000
    sample})

  #3 Pyrolysis CO, CO2 and SO2 thermogram

  #3.1 Baseline correction for CO
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["CO"]]<-sample[["Pyrolysis"]][["CO"]]-sample[["Cursors"]]["base2"]
    sample})

  #3.2 Sample weight normalisation for CO (mg/g)
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["CO"]]<-
      sample[["Pyrolysis"]][["CO"]]/as.numeric(sample[["Parameters"]]["Quant"])/1000
    sample})

  #3.3 Baseline correction for CO2
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["CO2"]]<-sample[["Pyrolysis"]][["CO2"]]-sample[["Cursors"]]["base3"]
    sample})

  #3.4 Sample weight normalisation for CO2 (mg/g)
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["CO2"]]<-
      sample[["Pyrolysis"]][["CO2"]]/as.numeric(sample[["Parameters"]]["Quant"])/1000
    sample})

  #3.5 Baseline correction for SO2
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["SO2"]]<-sample[["Pyrolysis"]][["SO2"]]-sample[["Cursors"]]["base4"]
    sample})

  #3.6 Sample weight normalisation for SO2 (mg/g)
  list.converted<-lapply(list.converted,
    function(sample){
    sample[["Pyrolysis"]][["SO2"]]<-
      sample[["Pyrolysis"]][["SO2"]]/as.numeric(sample[["Parameters"]]["Quant"])/1000
    sample})

  #4 Oxidation CO, CO2 and SO2 thermogram

  #4.1 Baseline correction for CO
  list.converted<-lapply(list.converted,
                         function(sample){
                           sample[["Oxidation"]][["CO"]]<-sample[["Oxidation"]][["CO"]]-sample[["Cursors"]]["base5"]
                           sample})

  #4.2 Sample weight normalisation for CO (mg/g)
  list.converted<-lapply(list.converted,
                         function(sample){
                           sample[["Oxidation"]][["CO"]]<-
                             sample[["Oxidation"]][["CO"]]/as.numeric(sample[["Parameters"]]["Quant"])/1000
                           sample})

  #4.3 Baseline correction for CO2
  list.converted<-lapply(list.converted,
                         function(sample){
                           sample[["Oxidation"]][["CO2"]]<-sample[["Oxidation"]][["CO2"]]-sample[["Cursors"]]["base6"]
                           sample})

  #4.4 Sample weight normalisation for CO2 (mg/g)
  list.converted<-lapply(list.converted,
                         function(sample){
                           sample[["Oxidation"]][["CO2"]]<-
                             sample[["Oxidation"]][["CO2"]]/as.numeric(sample[["Parameters"]]["Quant"])/1000
                           sample})

  #4.5 Baseline correction for SO2
  list.converted<-lapply(list.converted,
                         function(sample){
                           sample[["Oxidation"]][["SO2"]]<-sample[["Oxidation"]][["SO2"]]-sample[["Cursors"]]["base7"]
                           sample})

  #4.6 Sample weight normalisation for SO2 (mg/g)
  list.converted<-lapply(list.converted,
                         function(sample){
                           sample[["Oxidation"]][["SO2"]]<-
                             sample[["Oxidation"]][["SO2"]]/as.numeric(sample[["Parameters"]]["Quant"])/1000
                           sample})

list.converted
}




#' Cursor adjustment for soil cycle
#'
#' This function adjusts the unaltered pyrolysis CO2 cursor which separates organic and inorganic C. By default this value is fixed at 400 C in Geoworks, however, for soils the valley between the two peaks may be more appropriate.
#' The function alter only the samples that are run with the SOIL or SOIL TS cycle and only samples that have their cursor not manually adjusted already.
#'
#' @param list List with raw Rock-Eval data from RE_read
#' @param plot Create plot of cursor determination TRUE or FALSE
#' @return List with adjusted pyrolysis CO2 cursor
#' @export
RE_cursadjust<-function(list, plot=FALSE){

  #1.1 Create new converted list as copy of input list
  list.converted<-list

  #1.2 define safe sequence in case of NAs
  s.seq<-function(x,y){
    if (all(is.na(x))|all(is.na(y))) {NA} else {
      seq(x,y,by=1)}
  }

  #2 Determine the adjsuted cursors (separated from actual changing for flexibility)

  #2.1 Apply function to converted list
  list.converted<-lapply(list.converted, function(sample){

    #2.1.1 Determine original cursor
    curs.Geo<-sample[["Cursors"]][["curs3.2"]] #cursor predetermined in Geoworks
    cursor.adj<-c(curs3.2=curs.Geo, certainty=9)

    #2.1.2 Check if right cycle for adjustment
    if (sample[["Parameters"]]["CyclN"]=="SOIL" | sample[["Parameters"]]["CyclN"]=="SOIL TS") {
      temp.curs<-ifelse(is.na(curs.Geo),NA,sample[["Pyrolysis"]][["T"]][curs.Geo])

    #2.1.3 Check if original cursor is not already altered
      if (temp.curs>=398 & temp.curs<=402 & !is.na(temp.curs)) {

    #2.2 Determine adjusted cursor

        #2.2.1 Determine pseudo-derivatives of the curve with varying lag (to exlcude random noise)
        deriv1 <-diff(sample[["Pyrolysis"]][["CO2"]], lag = 1, differences=1)
        deriv4 <-diff(sample[["Pyrolysis"]][["CO2"]], lag = 4, differences=1)
        deriv10<-diff(sample[["Pyrolysis"]][["CO2"]], lag = 10, differences=1)

        #2.2.2 Time steps where derivative is 0
        ts.deriv1 <-which(rowSums(embed(sign(deriv1) ,1)) == 0)
        ts.deriv4 <-which(rowSums(embed(sign(deriv4) ,1)) == 0)+2
        ts.deriv10<-which(rowSums(embed(sign(deriv10),1)) == 0)+5

        #2.2.3 Determine the three lowest values within the time range where valley would be
        tr.min<-s.seq(sample[["Cursors"]]["curs3.2"],sample[["Cursors"]]["curs3.3"])
        ts.min<-head(order(sample[["Pyrolysis"]][["CO2"]][tr.min]),3)+sample[["Cursors"]]["curs3.2"]

        #2.2.4 Optional plotting of cursor determination on curves
        if(plot == TRUE){
        plot(sample[["Pyrolysis"]][["CO2"]][300:1200],type="l",
             xlab="Time", ylab="RE Signal CO2")

        #points(ts.deriv1-300,  sample[["Pyrolysis"]][["CO2"]][ts.deriv1],  col = "grey")
        points(ts.deriv4-300,  sample[["Pyrolysis"]][["CO2"]][ts.deriv4],  col="grey")
        points(ts.deriv10-300, sample[["Pyrolysis"]][["CO2"]][ts.deriv10], col="red3")

        points(ts.min-300,sample[["Pyrolysis"]][["CO2"]][ts.min], col="blue3", pch=15)
        } # end of 3rd if

        #2.2.5 Check the certainty of the minimum to be a good fit
        val1.check<- ts.min[1] %in% ts.deriv1 + ts.min[1] %in% ts.deriv4 + ts.min[1] %in% ts.deriv10
        val2.check<- ts.min[2] %in% ts.deriv1 + ts.min[2] %in% ts.deriv4 + ts.min[2] %in% ts.deriv10
        val3.check<- ts.min[3] %in% ts.deriv1 + ts.min[3] %in% ts.deriv4 + ts.min[3] %in% ts.deriv10

        val.check<-c(val1.check,val2.check,val3.check)
        val.select<-which.max(val.check)

        ts.min<-ts.min[val.select]
        cr.min<-val.check[val.select]

        #2.2.6 Combine the adjsuted cursor and its certainty for output
        cursor.adj<-c(curs3.2=ts.min, certainty=cr.min)

            } # end of 2nd if
    } # end of if

    #2.2.7 Add updated cursor information to the list
    sample[["cursor.adj"]]<-cursor.adj
    sample

  }) # end of apply

  #3 Take the updated cursor info and give warning message if uncertain cursors are included.
  list.values<-as.data.frame(
    do.call(rbind, lapply(list.converted, function(sample)sample[["cursor.adj"]])))

  uncertain.val<-which(list.values["certainty"]<2)

  if(length(uncertain.val)>0){
    warning(sprintf("Was unable to automatically determine a valley for minimum value with high certainty in item %i: %s. \n",
                    uncertain.val,
                    row.names(list.values)[uncertain.val]))
  }

  #4 Overwrite the old cursors with the new and remove the unnecessary cursor information from the list.
  list.converted<-lapply(list.converted, function(sample){
  sample[["Cursors"]][["curs3.2"]]<-sample[["cursor.adj"]][["curs3.2"]]
  sample[["cursor.adj"]]<-NULL
  sample
  })

  list.converted
} # end of function


