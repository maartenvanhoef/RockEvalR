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



#' Calculate the Rock-Eval surfaces
#'
#' This function calculates the areas of the different Rock-Eval zones. Input should be converted already.
#'
#' @param list List with converted Rock-Eval data from RE_read
#' @return Input list with included areas of Rock-Eval zones
#' @export
RE_surfaces<-function(list){
  # Uses the weight converted values which is slightly different from Geoworks

  #1.1 take list to extend
  list.extended<-list

  #1.2 define safe sequence in case of NAs
  s.seq<-function(x,y){
  if (all(is.na(x))|all(is.na(y))) {NA} else {
    seq(x,y,by=1)}
  }

  #2 Determine zones for C
  list.extended<-lapply(list.extended, function(sample){

      #2.1.1 Determine time ranges between Rock-Eval cursors
      tr.s1<-s.seq(min(sample[["Pyrolysis"]]["t"]),sample[["Cursors"]]["curs1.1"])
      tr.s2<-s.seq(sample[["Cursors"]]["curs1.1"],length(sample[["Pyrolysis"]][["t"]]))

      tr.s3CO<-s.seq(min(sample[["Pyrolysis"]]["t"]),sample[["Cursors"]]["curs2.2"])
      tr.s3COi<-s.seq(sample[["Cursors"]]["curs2.2"],length(sample[["Pyrolysis"]][["t"]]))

      tr.s3CO2<-s.seq(min(sample[["Pyrolysis"]]["t"]),sample[["Cursors"]]["curs3.2"])
      tr.s3CO2i<-s.seq(sample[["Cursors"]]["curs3.2"],length(sample[["Pyrolysis"]][["t"]]))

      tr.s4CO<-s.seq(min(sample[["Oxidation"]]["t"]),sample[["Cursors"]]["curs5.2"])
      tr.s4COi<-s.seq(sample[["Cursors"]]["curs5.2"],length(sample[["Oxidation"]][["t"]]))

      tr.s4CO2<-s.seq(min(sample[["Oxidation"]]["t"]),sample[["Cursors"]]["curs6.2"])
      tr.s5<-s.seq(sample[["Cursors"]]["curs6.2"],length(sample[["Oxidation"]][["t"]]))

      #2.1.2 Compute the area between these cursors
      S1<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s1],sample[["Pyrolysis"]][["CH"]][tr.s1])
      S2<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s2],sample[["Pyrolysis"]][["CH"]][tr.s2])

      S3CO<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3CO],sample[["Pyrolysis"]][["CO"]][tr.s3CO])
      S3COi<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3COi],sample[["Pyrolysis"]][["CO"]][tr.s3COi])

      S3CO2<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3CO2],sample[["Pyrolysis"]][["CO2"]][tr.s3CO2])
      S3CO2i<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3CO2i],sample[["Pyrolysis"]][["CO2"]][tr.s3CO2i])

      S4CO<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s4CO],sample[["Oxidation"]][["CO"]][tr.s4CO])
      S4COi<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s4COi],sample[["Oxidation"]][["CO"]][tr.s4COi])

      S4CO2<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s4CO2],sample[["Oxidation"]][["CO2"]][tr.s4CO2])
      S5<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s5],sample[["Oxidation"]][["CO2"]][tr.s5])

      #2.1.3 Return the values as an addition to the original list
      zones<-c(S1=S1, S2=S2,
               S3CO=S3CO, S3COi=S3COi, S3CO2=S3CO2, S3CO2i=S3CO2i,
               S4CO=S4CO, S4COi=S4COi, S4CO2=S4CO2, S5=S5)
      sample[["Zones_C"]]<-zones
      sample
    })

list.extended
}


# to do
# -automated cursor adjust for pyr CO2
# -but only for SOIL cycle
