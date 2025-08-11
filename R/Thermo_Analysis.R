#' Sebag (2016) areas and I/R indices
#'
#' Calculates the areas A1, A2, A3, A4 plus the I and R indices.
#'
#' @param list List with converted Rock-Eval data and zone areas
#' @return Computed indices from Sebag (2016)
#' @references Sebag, et al. (2016). Dynamics of soil organic matter based on new Rock-Eval indices. Geoderma, 284, 185-203.
#' @export
RE_SebagIR<-function(list){
  #1.1 take list to extend
  list.extended<-list

  #1.2 define safe sequence in case of NAs
  s.seq<-function(x,y){
    if (all(is.na(x))|all(is.na(y))) {NA} else {
      seq(x,y,by=1)}
  }

  #2 Determine areas for C
  list.extended<-lapply(list.extended, function(sample){

    #2.1.1 Determine time ranges for the areas S2, A1, A2, A3, A4
    tr.T1<-which(sample[["Pyrolysis"]]["T"] >= 200 &
                sample[["Pyrolysis"]]["T"] <= 340 &
                sample[["Pyrolysis"]]["t"] >= sample[["Cursors"]]["curs1.1"])
    tr.T2<-which(sample[["Pyrolysis"]]["T"] >= 340 &
                sample[["Pyrolysis"]]["T"] <= 400)
    tr.T3<-which(sample[["Pyrolysis"]]["T"] >= 400 &
                sample[["Pyrolysis"]]["T"] <= 460 &
                sample[["Pyrolysis"]]["t"] <= sample[["Cursors"]]["curs1.5"])
    tr.T4<-s.seq(which.max(sample[["Pyrolysis"]]["T"] >= 460),
                length(sample[["Pyrolysis"]][["T"]]))

    tr.R <-s.seq(which.max(sample[["Pyrolysis"]]["T"] >= 400),length(sample[["Pyrolysis"]][["t"]]))
    tr.s2<-s.seq(sample[["Cursors"]]["curs1.1"],length(sample[["Pyrolysis"]][["t"]]))

    #2.2.2 Compute the area for these ranges
    A1 <-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.T1],sample[["Pyrolysis"]][["CH"]][tr.T1])
    A2 <-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.T2],sample[["Pyrolysis"]][["CH"]][tr.T2])
    A3 <-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.T3],sample[["Pyrolysis"]][["CH"]][tr.T3])
    A4 <-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.T4],sample[["Pyrolysis"]][["CH"]][tr.T4])

    AR<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.R],sample[["Pyrolysis"]][["CH"]][tr.R])
    S2<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s2],sample[["Pyrolysis"]][["CH"]][tr.s2])

    #2.2.3 Calculate the indices with the ratios
    I<-log10((A1+A2)/A3)
    R<-AR/S2

    #2.2.4 Combine the areas into the list and return it
    Sebag_IR<-c(A1=A1, A2=A2, A3=A3, A4=A4,
                AR=AR, AS2=S2,
                I=I, R=R)
    sample[["Sebag_IR"]]<-Sebag_IR
    sample
  })

list.extended
}



#' Determine relative T percentiles
#'
#' Calculates the temperature percentiles of the thermograms
#'
#' @param list List with converted Rock-Eval data and zone areas
#' @param vector Vector with the percentiles to be calculated
#' @return Temperature percentiles
#' @export
RE_Tpercentiles<-function(list,vector=c(25,50,75)){

  #1.1 take list to extend
  list.extended<-list

  #1.2 name of zones and the T's based on the input vector
  zones<-c("S2","S3CO","S3CO2","S4CO","S4CO2","POC","ROC")
  names.T<-paste0(rep(zones,length(vector)),"_T",rep(vector, each=length(zones)))

  #1.3 define safe sequence in case of NAs
  s.seq<-function(x,y){
    if (all(is.na(x))|all(is.na(y))) {NA} else {
      seq(x,y,by=1)}
  }

  #2.1 determine the cumulative area at each time step for all zones
  list.extended<-lapply(list.extended, function(sample){

    #2.1.1 Determine time ranges between Rock-Eval cursors
    tr.s1<-s.seq(1,sample[["Cursors"]]["curs1.1"])
    tr.s2<-s.seq(sample[["Cursors"]]["curs1.1"],length(sample[["Pyrolysis"]][["t"]]))

    tr.s3CO<-s.seq(1,sample[["Cursors"]]["curs2.2"])
    tr.s3COi<-s.seq(sample[["Cursors"]]["curs2.2"],length(sample[["Pyrolysis"]][["t"]]))
    tr.s3CO2<-s.seq(1,sample[["Cursors"]]["curs3.2"])

    tr.s4CO<-s.seq(1,sample[["Cursors"]]["curs5.2"])
    tr.s4CO2<-s.seq(1,sample[["Cursors"]]["curs6.2"])

    #2.1.2 Compute the area between these cursors
    S1<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s1],sample[["Pyrolysis"]][["CH"]][tr.s1])
    S1.c<-RE_cumtraparea(sample[["Pyrolysis"]][["t"]][tr.s1],sample[["Pyrolysis"]][["CH"]][tr.s1])

    S2<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s2],sample[["Pyrolysis"]][["CH"]][tr.s2])
    S2.c<-RE_cumtraparea(sample[["Pyrolysis"]][["t"]][tr.s2],sample[["Pyrolysis"]][["CH"]][tr.s2])

    S3CO<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3CO],sample[["Pyrolysis"]][["CO"]][tr.s3CO])
    S3CO.c<-RE_cumtraparea(sample[["Pyrolysis"]][["t"]][tr.s3CO],sample[["Pyrolysis"]][["CO"]][tr.s3CO])

    S3COi<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3COi],sample[["Pyrolysis"]][["CO"]][tr.s3COi])
    S3COi.c<-RE_cumtraparea(sample[["Pyrolysis"]][["t"]][tr.s3COi],sample[["Pyrolysis"]][["CO"]][tr.s3COi])

    S3CO2<-RE_traparea(sample[["Pyrolysis"]][["t"]][tr.s3CO2],sample[["Pyrolysis"]][["CO2"]][tr.s3CO2])
    S3CO2.c<-RE_cumtraparea(sample[["Pyrolysis"]][["t"]][tr.s3CO2],sample[["Pyrolysis"]][["CO2"]][tr.s3CO2])

    S4CO<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s4CO],sample[["Oxidation"]][["CO"]][tr.s4CO])
    S4CO.c<-RE_cumtraparea(sample[["Oxidation"]][["t"]][tr.s4CO],sample[["Oxidation"]][["CO"]][tr.s4CO])

    S4CO2<-RE_traparea(sample[["Oxidation"]][["t"]][tr.s4CO2],sample[["Oxidation"]][["CO2"]][tr.s4CO2])
    S4CO2.c<-RE_cumtraparea(sample[["Oxidation"]][["t"]][tr.s4CO2],sample[["Oxidation"]][["CO2"]][tr.s4CO2])

    #2.1.3 Calculate the (cumulative) POC and ROC from these as well
    POC<-S1*0.83/10+S2*0.83/10+S3CO*12/(28*10)+S3COi*12/(2*28*10)+S3CO2*12/(44*10)
    ROC<-S4CO*12/(28*10)+S4CO2*12/(44*10)

    P.CH<-c(S1.c*0.83/10,
            tail(S1.c,1)*0.83/10+S2.c*0.83/10)
    P.CO<-c(S3CO.c*12/(28*10),
            tail(S3CO.c,1)*12/(28*10)+S3COi.c*12/(2*28*10))
    P.CO2<-c(S3CO2.c*12/(44*10),
             tail(S3CO2.c,1)*12/(44*10)+
             rep(0,length(sample[["Pyrolysis"]][["t"]])-1-length(S3CO2.c))) # keep same value untill end
    POC.c<-P.CH+P.CO+P.CO2

    R.CO<-c(S4CO.c*12/(28*10),
            tail(S4CO.c,1)*12/(28*10)+
            rep(0,length(sample[["Oxidation"]][["t"]])-1-length(S4CO.c)))
    R.CO2<-c(S4CO2.c*12/(44*10),
            tail(S4CO2.c,1)*12/(44*10)+
            rep(0,length(sample[["Oxidation"]][["t"]])-1-length(S4CO2.c)))
    ROC.c<-R.CO+R.CO2


  #2.2 for loop to apply variable number of T's to list
    value.v<-c(rep(NA,length(names.T)))

    for (i in 1:length(vector)) {
      # 2.2.1 determine time step at which each percentile is reached
      T.i<-vector[i]

      ts.S2.i<-which.max(S2.c>=S2*T.i/100)
      ts.S2.i<-ifelse(length(ts.S2.i)==0, NA, ts.S2.i)

      ts.S3CO.i<-which.max(S3CO.c>=S3CO*T.i/100)
      ts.S3CO.i<-ifelse(length(ts.S3CO.i)==0, NA, ts.S3CO.i)

      ts.S3CO2.i<-which.max(S3CO2.c>=S3CO2*T.i/100)
      ts.S3CO2.i<-ifelse(length(ts.S3CO2.i)==0, NA, ts.S3CO2.i)

      ts.S4CO.i<-which.max(S4CO.c>=S4CO*T.i/100)
      ts.S4CO.i<-ifelse(length(ts.S4CO.i)==0, NA, ts.S4CO.i)

      ts.S4CO2.i<-which.max(S4CO2.c>=S4CO2*T.i/100)
      ts.S4CO2.i<-ifelse(length(ts.S4CO2.i)==0, NA, ts.S4CO2.i)

      ts.POC.i<-which.max(POC.c>=POC*T.i/100)
      ts.POC.i<-ifelse(length(ts.POC.i)==0, NA, ts.POC.i)

      ts.ROC.i<-which.max(ROC.c>=ROC*T.i/100)
      ts.ROC.i<-ifelse(length(ts.ROC.i)==0, NA, ts.ROC.i)

      # 2.2.2 determine temperature at each percentile time step
      T.S2.i<-sample[["Pyrolysis"]][["T"]][min(tr.s2)+ts.S2.i]

      T.S3CO.i<-sample[["Pyrolysis"]][["T"]][min(tr.s3CO)+ts.S3CO.i]
      T.S3CO2.i<-sample[["Pyrolysis"]][["T"]][min(tr.s3CO2)+ts.S3CO2.i]

      T.S4CO.i<-sample[["Oxidation"]][["T"]][min(tr.s4CO)+ts.S4CO.i]
      T.S4CO2.i<-sample[["Oxidation"]][["T"]][min(tr.s4CO2)+ts.S4CO2.i]

      T.POC.i<-sample[["Pyrolysis"]][["T"]][ts.POC.i+0] # +0 so only one NA can be returned
      T.ROC.i<-sample[["Oxidation"]][["T"]][ts.ROC.i+0] # which is a lazy fix for now


      # 2.2.2 get the temperature value at this time step into vector
      value.v[1:length(zones)+length(zones)*(i-1)]<-
        c(T.S2.i,T.S3CO.i,T.S3CO2.i,T.S4CO.i,T.S4CO2.i,T.POC.i,T.ROC.i)
    }

      #2.2.3 name the values and return
      names(value.v)<-names.T
      value.v<-value.v[order(names(value.v))]
      sample[["Tpercentiles"]]<-value.v
      sample
  })

list.extended
}



#' Determine the total C released at each time step
#'
#' Combines the amount of C released of the different thermograms to calculate the flux of PC and RC over time.
#'
#' @param list List with converted Rock-Eval data
#' @return Combined C flux at each time step
#' @export
RE_Ccombined<-function(list){
  #1.1 make new list to be adapted from input
  list.adapted<-list

  #1.2 define safe sequence in case of NAs
  s.seq<-function(x,y){
    if (all(is.na(x))|all(is.na(y))) {NA} else {
      seq(x,y,by=1)}
  }

  #2 the total C in each oven is a simple addition
  list.adapted<-lapply(list.adapted, function(sample){

    PC.t<-
      sample[["Pyrolysis"]][["CH"]]*0.83/10+
      sample[["Pyrolysis"]][["CO"]]*12/(28*10)+
      sample[["Pyrolysis"]][["CO2"]]*12/(44*10)

   RC.t<-
      sample[["Oxidation"]][["CO"]]*12/(28*10)+
      sample[["Oxidation"]][["CO2"]]*12/(44*10)

  #3 the separation between organic and inorganic varies for the different curves (S1, S2, S3, ...)

    #3.1.1 Determine time ranges between Rock-Eval cursors
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

    #3.1.2 Compute the C flux between these cursors
    S1.t<-sample[["Pyrolysis"]][["CH"]][tr.s1]*0.83/10
    S2.t<-sample[["Pyrolysis"]][["CH"]][tr.s2]*0.83/10

    S3CO.t<-sample[["Pyrolysis"]][["CO"]][tr.s3CO]*12/(28*10)
    S3COi.t<-sample[["Pyrolysis"]][["CO"]][tr.s3COi]*12/(28*10)

    S3CO2.t<-sample[["Pyrolysis"]][["CO2"]][tr.s3CO2]*12/(44*10)
    S3CO2i.t<-sample[["Pyrolysis"]][["CO2"]][tr.s3CO2i]*12/(44*10)

    S4CO.t<-sample[["Oxidation"]][["CO"]][tr.s4CO]*12/(28*10)
    S4COi.t<-sample[["Oxidation"]][["CO"]][tr.s4COi]*12/(28*10)

    S4CO2.t<-sample[["Oxidation"]][["CO2"]][tr.s4CO2]*12/(44*10)
    S5.t<-sample[["Oxidation"]][["CO2"]][tr.s5]*12/(44*10)

    #3.1.3 Combine the separated curves over the whole time range
    POC.t<-c(S1.t,S2.t)+c(S3CO.t,S3COi.t/2)+c(S3CO2.t,S3CO2i.t*0)
    PIC.t<-c(S1.t*0,S2.t*0)+c(S3CO.t*0,S3COi.t/2)+c(S3CO2.t*0,S3CO2i.t)

    ROC.t<-c(S4CO.t,S4COi.t*0)+c(S4CO2.t,S5.t*0)
    RIC.t<-c(S4CO.t*0,S4COi.t)+c(S4CO2.t*0,S5.t)

    values.P<-data.frame(POC=POC.t, PIC=PIC.t, PC=PC.t)
    values.O<-data.frame(ROC=ROC.t, RIC=RIC.t, RC=RC.t)

    sample[["Pyrolysis"]]<-cbind(sample[["Pyrolysis"]],values.P)
    sample[["Oxidation"]]<-cbind(sample[["Oxidation"]],values.O)

    sample
  })

list.adapted
}



