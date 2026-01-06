#' Calculate the Rock-Eval Sulphur surfaces
#'
#' This function calculates the areas of the different Rock-Eval zones for Sulphur. Input should be converted already. As an option the zones can also be included as a time series.
#'
#' @param list List with converted Rock-Eval data
#' @param time.include Include the curves as individual time series TRUE/FALSE
#' @return Input list with included areas of Rock-Eval zones
#'
#' @section How to use:
#' It is important that first the data have been converted to the right units (by calibration values) and that the baseline is set to 0.
#' Converting the raw measurements to proper units and baselined values can be done with the RE_convert() function.
#' This function integrates the different thermograms for S. Integration happens between the cursors that were determined automatically or manually adjusted.
#' Here the areas and nomenclature of Cohen-Sadon et al. (2022) are used instead of those in Geoworks. However, the cut-off temperatures are not fixed an based on the cursors and local minima.
#'
#' @section Parameter info:
#' \itemize{
#' \item S1S, Labile pyrolysable organic S (mg S/g)
#' \item S2S, Pyrolysable organic S (mg S/g)
#' \item S3S, Pyrolysable inorganic S (mg S/g)
#' \item S4S, Oxidised S (mg S/g)
#' \item S5S, Residual oxidised S (mg S/g)
#'}
#'
#' @section Cut-off info:
#' \itemize{
#' \item S1S, pyrolysed before start of the heating ramp
#' \item S2S, pyrolysed from start of heating until minimum between 475 and 550 C
#' \item S3S, pyrolysed after minimum between 475 and 550 C
#' \item S4S, oxidised below minimum between 550 and 650 C
#' \item S5S, oxidised above minimum between 550 and 650 C
#'}
#'
#' @references
#' Cohen-Sadon, et al. (2022). A new empirical approach for rapid quantification of organic and pyritic sulphur in sedimentary rocks using the Rock-Eval 7S. Organic Geochemistry, 166, 104350.
#'
#' Aboussou (2018). New Rock-Eval method for Pyritic and Organic Sulphur quantification: Application to study Organic matter preservation in Jurassic sediments. Earth Sciences, Sorbonne Université.
#' @export
RE_Ssurfaces<-function(list, time.include=FALSE){
  # Uses the weight converted values which is slightly different from Geoworks

  #1.1 take list to extend
  list.extended<-list

  #1.2 define safe sequence in case of NAs
  s.seq<-function(x,y){
    if (all(is.na(x))|all(is.na(y))) {NA} else {
      seq(x,y,by=1)}
  }

  #2 Determine zones for S
  list.extended<-lapply(list.extended, function(sample){

    #2.1.1 Determine the lowest value within the time range where minima need to be determined
    tr.min.s2S<-which(sample[["Pyrolysis"]]["T"] >= 400 &
                      sample[["Pyrolysis"]]["T"] <= 550 &
                      sample[["Pyrolysis"]]["t"] < sample[["Cursors"]]["curs4.3"])

    ts.min.s2S<-which.min(sample[["Pyrolysis"]][["SO2"]][tr.min.s2S])+which.max(sample[["Pyrolysis"]]["T"] >= 475)


    tr.min.s4S<-which(sample[["Oxidation"]]["T"] >= 550 &
                      sample[["Oxidation"]]["T"] <= 650)

    ts.min.s4S<-which.min(sample[["Oxidation"]][["SO2"]][tr.min.s4S])+which.max(sample[["Pyrolysis"]]["T"] >= 550)


    #2.1.2 Determine time ranges between Rock-Eval temperature ranges
    tr.s1S<-s.seq(1,sample[["Cursors"]]["curs4.1"])

    tr.s2S<-which(sample[["Pyrolysis"]]["t"] >= sample[["Cursors"]]["curs4.1"] &
                  sample[["Pyrolysis"]]["t"] <= ts.min.s2S)

    tr.s3S<-which(sample[["Pyrolysis"]]["t"] >= ts.min.s2S)


    tr.s4S<-which(sample[["Oxidation"]]["t"] <= ts.min.s4S)

    tr.s5S<-which(sample[["Oxidation"]]["t"] >= ts.min.s4S)


    #2.1.3 Compute the area between these cut-offs
    S1S<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s1S],sample[["Pyrolysis"]][["SO2"]][tr.s1S])*32.07/64.07
    S2S<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s2S],sample[["Pyrolysis"]][["SO2"]][tr.s2S])*32.07/64.07

    S3S<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3S],sample[["Pyrolysis"]][["SO2"]][tr.s3S])*32.07/64.07

    S4S<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s4S],sample[["Oxidation"]][["SO2"]][tr.s4S])*32.07/64.07
    S5S<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s5S],sample[["Oxidation"]][["SO2"]][tr.s5S])*32.07/64.07



    #2.1.3 Return the values as an addition to the original list
    zones<-c(S1S=S1S, S2S=S2S, S3S=S3S,
             S4S=S4S, S5S=S5S)
    sample[["Zones_S"]]<-zones



    #3 Include the surfaces as a time series in the data
    if (time.include==TRUE) {

      #3.1 Take values as time series instead of areas
      S1S.t<-sample[["Pyrolysis"]][["SO2"]][tr.s1S]
      S2S.t<-sample[["Pyrolysis"]][["SO2"]][tr.s2S]

      S3S.t<-sample[["Pyrolysis"]][["SO2"]][tr.s3S]

      S4S.t<-sample[["Oxidation"]][["SO2"]][tr.s4S]
      S5S.t<-sample[["Oxidation"]][["SO2"]][tr.s5S]

      #3.2 Make time series continuous over entire time range
      S1S<-c(S1S.t,S2S.t*0,S3S.t*0)
      S2S<-c(S1S.t*0,S2S.t,S3S.t*0)
      S3S<-c(S1S.t*0,S2S.t*0,S3S.t)

      S4S<-c(S4S.t,S5S.t*0)
      S5S<-c(S4S.t*0,S5S.t)

      #3.3 Combine and return
      values.P<-data.frame(S1S=S1S, S2S=S2S, S3S=S3S)
      values.O<-data.frame(S4S=S4S, S5S=S5S)

      sample[["Pyrolysis"]]<-cbind(sample[["Pyrolysis"]],values.P)
      sample[["Oxidation"]]<-cbind(sample[["Oxidation"]],values.O)

    }
    sample


  }) # end of lapply

  list.extended
}



#' Calculate Rock-Eval Sulphur metrics
#'
#' This function calculates main Rock-Eval Sulphur metrics from zone areas.
#'
#' @param list List with converted Rock-Eval data and zone areas
#' @param org.conversion.a Slope of pyrolysis to total organic S empirical conversion
#' @param org.conversion.b Intercept of pyrolysis to total organic S empirical conversion
#' @return Computed Rock-Eval metrics
#'
#' @section How to use:
#' It is necessary that the RE Sulphur surfaces have been integrated first.
#' This can be done with the RE_Ssurfaces() function.
#' Additionally, for the sulphur index (SI) the total C has to be determined with the RE_surfaces() function.
#' From these areas the basic Rock-Eval metrics are then computed.
#' The determination of most sulphur-species depends on empirical correction of matrix effects.
#' Based on type of sample (e.g., Kerogen type), and data availability different empirical adjustments may lead to better results.
#'
#' In this case total organic S is calculated from a relationship with the Tmax_CH during pyrolysis:
#'
#' (Pyrolysed org S) / (Total org S) = (org.conversion.a) x (Tmax_CH) + (org.conversion.b)
#'
#' (Total inorg S) = (Total S) - (Total org S)
#'
#' @section Parameter info:
#' \itemize{
#' \item Tpeak_S, RE7 temperature at organic SO2 peak during pyrolysis (C)
#' \item POS, pyrolysis organic S (%)
#' \item PIS, pyrolysis inorganic S (%)
#' \item RTS, total oxidised S (%)
#' \item TS, total measured S (%)
#' \item PyOS.emp, emperically determined ratio of pyrolysed organic S (%)
#' \item ROS.cor, empirically corrected residual organic S (%)
#' \item TOS.cor, empirically corrected total organic S (%)
#' \item RIS.cor, empirically corrected residual inorganic S (%)
#' \item TIS.cor, empirically corrected total inorganic S (%)
#' \item SI, sulphur index (mg total org. S/ g TOC)
#' \item SI.pyr, labile sulphur index (mg pyrolysed org. S/ g TOC)
#'}
#'
#' @seealso [RE_surfaces()]
#'
#' @references
#' Cohen-Sadon, et al. (2022). A new empirical approach for rapid quantification of organic and pyritic sulphur in sedimentary rocks using the Rock-Eval 7S. Organic Geochemistry, 166, 104350.
#'
#' Aboussou (2018). New Rock-Eval method for Pyritic and Organic Sulphur quantification: Application to study Organic matter preservation in Jurassic sediments. Earth Sciences, Sorbonne Université.
#'
#' Cohen-Sadon, et al. (2025). Tmax-S: A new proxy for the role of sulphur on sedimentary organic matter preservation and thermal maturation. Journal of Analytical and Applied Pyrolysis 190 - 107115.
#'
#' @export
RE_Smetrics<-function(list, org.conversion.a=-1.90, org.conversion.b=843){
  #1.1 take list to extend
  list.extended<-list

  #1.2 check if needed zones are present in list
  if(is.null(list.extended[[1]][["Zones_S"]])){
    stop("First the Rock-Eval zones (S1S, S2S, ...) need to be defined.")
  }

  #2 calculate the parameters
  list.extended<-lapply(list.extended, function(sample){

    #2.1 Temperature parameters
    tr.min.s2S<-which(sample[["Pyrolysis"]]["T"] >= 400 &
                      sample[["Pyrolysis"]]["T"] <= 550 &
                      sample[["Pyrolysis"]]["t"] < sample[["Cursors"]]["curs4.3"])

    ts.min.s2S<-which.min(sample[["Pyrolysis"]][["SO2"]][tr.min.s2S])+which.max(sample[["Pyrolysis"]]["T"] >= 475)

    tr.Tpeak_S<-which(sample[["Pyrolysis"]]["t"] <= ts.min.s2S)

    ts.Tpeak_S<-which.max(sample[["Pyrolysis"]][["SO2"]][tr.Tpeak_S])

    Tpeak_S<-sample[["Pyrolysis"]][["T"]][ts.Tpeak_S]
    Tpeak_S<-ifelse(is.null(Tpeak_S),NA,Tpeak_S)


    #2.2 Pyrolysable sulphur
    POS<-sample[["Zones_S"]][["S1S"]]/10+
         sample[["Zones_S"]][["S2S"]]/10
    PIS<-sample[["Zones_S"]][["S3S"]]/10


    #2.3 Residual (oxidised) total sulphur
    RTS<-sample[["Zones_S"]][["S4S"]]/10+
         sample[["Zones_S"]][["S5S"]]/10

    #2.4 Total sulphur (both ovens)
    TS<-POS+PIS+RTS

    #2.5 Empirically corrected organic sulphur species
    Tmax_CH<-sample[["Metrics_C"]][["Tpeak"]]-39 #using empirical conversion from Tpeak of S2 to Tmax

    PyOS.emp<-org.conversion.a*Tmax_CH + org.conversion.b
    PyOS.emp<-ifelse(PyOS.emp>100,NA,PyOS.emp)
    PyOS.emp<-ifelse(PyOS.emp<1,NA,PyOS.emp)

    TOS.cor<-POS*100/(PyOS.emp)
    ROS.cor<-TOS.cor-POS
    ROS.cor<-ifelse(ROS.cor<0,0,ROS.cor)

    #2.6 Emperically corrected inorganic sulphur species
    TIS.cor<-TS-TOS.cor
    RIS.cor<-TIS.cor-PIS


    #2.7 Sulphur index
    SI<-TOS.cor/sample[["Metrics_C"]][["TOC"]]*100
    SI.pyr<-POS/sample[["Metrics_C"]][["TOC"]]*100


    #3 Grouping all parameters and returning to the original list
    metrics<-c(Tpeak_S=Tpeak_S,
               POS=POS, PIS=PIS, RTS=RTS, TS=TS,
               PyOS.emp=PyOS.emp, ROS.cor=ROS.cor, TOS.cor=TOS.cor,
               RIS.cor=RIS.cor, TIS.cor=TIS.cor,
               SI=SI, SI.pyr=SI.pyr)
    sample[["Metrics_S"]]<-metrics
    sample
  })

  list.extended
}
