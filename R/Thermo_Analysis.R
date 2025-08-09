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
                AR=AR, S2=S2,
                I=I, R=R)
    sample[["Sebag_IR"]]<-Sebag_IR
    sample
  })

list.extended
}



# to do
# -determination of T-zones
# -cursors and selective integration


