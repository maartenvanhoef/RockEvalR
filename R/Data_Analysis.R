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

