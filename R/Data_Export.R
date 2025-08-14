#' Extract all time independent Rock-Eval metrics
#'
#' Creates a table with all calculated metrics added to the data.
#'
#' @param list List with converted and extended Rock-Eval data
#' @return Simplest possible data type containing requested data
#' @export
RE_extracttable<-function(list){

  #1 Definition of data in the list that is to be exported.
  in.parameters<-c("Sample","Quant","CyclT")
  nin.lists<-c("Parameters","Pyrolysis","Oxidation","Cursors","Calibration","Coefficients")

  #2 Extraction of data to be exported
  names<-data.frame(Analysis=gsub(".B00","",names(list)))

  parameters.all<-as.data.frame(t(sapply(list,function(sample)sample[["Parameters"]])))
  parameters<-parameters.all[,in.parameters]

  metrics.all<- as.data.frame(t(
    sapply(list,function(sample){
      do.call(c,unname(sample[!(names(sample) %in% nin.lists)]))})
    ))

  # Combine all data into a single table and return.
  table<-cbind(names,parameters,metrics.all)
  rownames(table)<-NULL
  table
}



#' Do all non-custom calculations at once
#'
#' Combines all main functions in one go. Not recommended for first time use.
#'
#' @param list List with raw Rock-Eval data
#' @return Converted list with all basic metrics
#' @export
RE_calculateall<-function(list){

  list.converted<-RE_convert(list)
  list.converted<-RE_surfaces(list.converted)
  list.converted<-RE_metrics(list.converted)
  list.converted<-RE_SebagIR(list.converted)
  list.converted<-RE_Tpercentiles(list.converted)

  list.converted
}
