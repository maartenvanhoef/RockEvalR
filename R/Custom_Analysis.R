# to do
# -cursors and selective integration
# -automated cursor fix for soil program

#' Calculate custom Rock-Eval surfaces
#'
#' This function calculates the areas between a set of temperatures. Input should be converted already.
#' The areas between all temperatures is given plus the area between 0 and the first temperature and the area after the highest temperature is reached.
#'
#' @param list List with converted Rock-Eval data from RE_read
#' @param oven Pyrolysis or Oxidation oven
#' @param thermogram Must be one of the columns of the oven data, e.g., CO2
#' @param vector Vector containing temperature intervals (excluding 0 and max)
#' @return Input list with included areas of Rock-Eval zones
#' @export
RE_customarea<-function(list,oven,thermogram,temperatures){
  #1.1 take list to extend
  list.extended<-list

  #1.2 define safe sequence in case of NAs
  s.seq<-function(x,y){
    if (all(is.na(x))|all(is.na(y))) {NA} else {
      seq(x,y,by=1)}
  }

  #2.1 Determine zones for C
  table<-as.data.frame(t(sapply(list.extended, function(sample){

    #2.2 check if temperature vector is appropriate

    if(temperatures[length(temperatures)-1]>tail(sample[[oven]][["T"]],1)){
      warning("The highest given interval is above the final oven temperature (including cool down), which may give unwanted results. The input range should exclude your min and max temperature.")
    }

    if(temperatures[1]<=head(sample[[oven]][["T"]],1)){
      warning("The lowest given interval is below the starting oven temperature, which may give unwanted results. The input range should exclude your min and max temperature.")
    }

  invector<-c(0,sort(temperatures))
  areas<-c()

  #2.3 Loop over input vector to determine the zones
  for (i in 2:length(invector)) {
  tr.i<-which(sample[[oven]][["T"]]>=invector[i-1] &
              sample[[oven]][["T"]]<invector[i])

  a.tr<-sample[[oven]][["t"]][tr.i]
  a.val<-sample[[oven]][[thermogram]][tr.i]
  area.i<-RE_traparea(a.tr,a.val)

  names(area.i)<-paste0("A",invector[i-1],"-",invector[i])

  areas<-c(areas,area.i)
  }

  #2.4 Calulcate the part left after the highest temperature
  tr.end<-s.seq(which.max(sample[[oven]][["T"]]>=invector[length(invector)]),
                length(sample[[oven]][["T"]]))

  a.tr<-sample[[oven]][["t"]][tr.end]
  a.val<-sample[[oven]][[thermogram]][tr.end]
  area.end<-RE_traparea(a.tr,a.val)

  names(area.end)<-paste0("A",invector[length(invector)],"+")

  areas<-c(areas,area.end)


  }))) # end of sapply and data.frame
} # end of function


