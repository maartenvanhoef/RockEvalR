#' Calculate the Rock-Eval surfaces
#'
#' This function calculates the areas of the different Rock-Eval zones. Input should be converted already. As an option the zones can also be included as a timeseries.
#'
#' @param list List with converted Rock-Eval data
#' @param time.include Include the curves as individual time series TRUE/FALSE
#' @return Input list with included areas of Rock-Eval zones
#' @export
RE_surfaces<-function(list, time.include=FALSE){
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
    tr.s1<-s.seq(1,sample[["Cursors"]]["curs1.1"])
    tr.s2<-s.seq(sample[["Cursors"]]["curs1.1"]+1,length(sample[["Pyrolysis"]][["t"]]))

    tr.s3CO<-s.seq(1,sample[["Cursors"]]["curs2.2"])
    tr.s3COi<-s.seq(sample[["Cursors"]]["curs2.2"]+1,length(sample[["Pyrolysis"]][["t"]]))

    tr.s3CO2<-s.seq(1,sample[["Cursors"]]["curs3.2"])
    tr.s3CO2i<-s.seq(sample[["Cursors"]]["curs3.2"]+1,length(sample[["Pyrolysis"]][["t"]]))

    tr.s4CO<-s.seq(1,sample[["Cursors"]]["curs5.2"])
    tr.s4COi<-s.seq(sample[["Cursors"]]["curs5.2"]+1,length(sample[["Oxidation"]][["t"]]))

    tr.s4CO2<-s.seq(1,sample[["Cursors"]]["curs6.2"])
    tr.s5<-s.seq(sample[["Cursors"]]["curs6.2"]+1,length(sample[["Oxidation"]][["t"]]))

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


    #3 Include the surfaces as a time series in the data
    if (time.include==TRUE) {

    #3.1 Take values as time series instead of areas
    S1.t<-sample[["Pyrolysis"]][["CH"]][tr.s1]
    S2.t<-sample[["Pyrolysis"]][["CH"]][tr.s2]

    S3CO.t<-sample[["Pyrolysis"]][["CO"]][tr.s3CO]
    S3COi.t<-sample[["Pyrolysis"]][["CO"]][tr.s3COi]

    S3CO2.t<-sample[["Pyrolysis"]][["CO2"]][tr.s3CO2]
    S3CO2i.t<-sample[["Pyrolysis"]][["CO2"]][tr.s3CO2i]

    S4CO.t<-sample[["Oxidation"]][["CO"]][tr.s4CO]
    S4COi.t<-sample[["Oxidation"]][["CO"]][tr.s4COi]

    S4CO2.t<-sample[["Oxidation"]][["CO2"]][tr.s4CO2]
    S5.t<-sample[["Oxidation"]][["CO2"]][tr.s5]

    #3.2 Make time series continuous over entire time range
    S1<-c(S1.t,S2.t*0)
    S2<-c(S1.t*0,S2.t)
    S3CO<-c(S3CO.t,S3COi.t/2)
    S3COi<-c(S3CO.t*0,S3COi.t/2)

    S3CO2<-c(S3CO2.t,S3CO2i.t*0)
    S3CO2i<-c(S3CO2.t*0,S3CO2i.t)

    S4CO<-c(S4CO.t,S4COi.t*0)
    S4COi<-c(S4CO.t*0,S4COi.t)

    S4CO2<-c(S4CO2.t,S5.t*0)
    S5<-c(S4CO2.t*0,S5.t)

    #3.3 Combine and return
    values.P<-data.frame(S1=S1, S2=S2,
                         S3CO=S3CO, S3COi=S3COi,
                         S3CO2=S3CO2, S3CO2i=S3CO2i)
    values.O<-data.frame(S4CO=S4CO, S4COi=S4COi,
                         S4CO2=S4CO2, S5=S5)

    sample[["Pyrolysis"]]<-cbind(sample[["Pyrolysis"]],values.P)
    sample[["Oxidation"]]<-cbind(sample[["Oxidation"]],values.O)

    }
    sample


  }) # end of lapply

  list.extended
}



#' Calculate Rock-Eval metrics
#'
#' This function calculates main Rock-Eval metrics from zone areas.
#'
#' @param list List with converted Rock-Eval data and zone areas
#' @return Computed Rock-Eval metrics
#'
#' @section Parameter info:
#' \itemize{
#' \item Tpeak, RE7 temperature at S2 peak (C)
#' \item POC, pyrolysis organic C (%)
#' \item ROC, residual organic C (%)
#' \item TOC, total organic C (%)
#' \item PIC, pyrolysis inorganic C (%)
#' \item RIC, residual inorganic C (%)
#' \item TIC, total inorganic C (%)
#' \item HI, hydrogen index (mg CH/ g TOC)
#' \item OI_CO, CO oxygen index (mg CO/ g TOC)
#' \item OI_CO2, main oxygen index (mg CO2/ g TOC)
#' \item OI_combi, combined oxygen index/ OIRE6 (mg O/ g TOC)
#'}
#'
#' @export
RE_metrics<-function(list){
  #1.1 take list to extend
  list.extended<-list

  #1.2 check if needed zones are present in list
  if(is.null(list.extended[[1]][["Zones_C"]])){
    stop("First the Rock-Eval zones (S1, S2, ...) need to be defined.")
  }

  #2 calculate the parameters
  list.extended<-lapply(list.extended, function(sample){

    #2.1 Temperature parameters
    Tpeak.t<-which(sample[["Pyrolysis"]][["t"]]==sample[["Cursors"]]["curs1.2"])
    Tpeak<-sample[["Pyrolysis"]][["T"]][Tpeak.t]
    Tpeak<-ifelse(is.null(Tpeak),NA,Tpeak)

    # Something about the Tmax seems to be different from Geoworks value...
    # Excluded in output for now
    Tmax<-Tpeak+(sample[["Calibration"]][["cal_Tmax"]]-sample[["Calibration"]][["cal_Tpeak"]])

    #2.2 Pyrolysable, residual and total organic C
    POC<-sample[["Zones_C"]][["S1"]]*0.83/10+
         sample[["Zones_C"]][["S2"]]*0.83/10+
         sample[["Zones_C"]][["S3CO"]]*12/(28*10)+
         sample[["Zones_C"]][["S3COi"]]*12/(2*28*10)+
         sample[["Zones_C"]][["S3CO2"]]*12/(44*10)
    ROC<-sample[["Zones_C"]][["S4CO"]]*12/(28*10)+
         sample[["Zones_C"]][["S4CO2"]]*12/(44*10)
    TOC<-POC+ROC

    #2.2 Pyrolysable, residual and total inorganic C
    PIC<-sample[["Zones_C"]][["S3COi"]]*12/(2*28*10)+
         sample[["Zones_C"]][["S3CO2i"]]*12/(44*10)
    RIC<-sample[["Zones_C"]][["S5"]]*12/(44*10)
    TIC<-PIC+RIC

    #Hydrogen and oxygen index (plus extra)
    HI<-sample[["Zones_C"]][["S2"]]/TOC*100

    OI_CO<-sample[["Zones_C"]][["S3CO"]]/TOC*100
    OI_CO2<-sample[["Zones_C"]][["S3CO2"]]/TOC*100
    OI_combi<-OI_CO*16/28+
              OI_CO2*32/44

  metrics<-c(Tpeak=Tpeak,
             POC=POC, ROC=ROC, TOC=TOC,
             PIC=PIC, RIC=RIC, TIC=TIC,
             HI=HI,
             OI_CO=OI_CO, OI_CO2=OI_CO2, OI_combi=OI_combi)
  sample[["Metrics_C"]]<-metrics
  sample
  })

list.extended
}

