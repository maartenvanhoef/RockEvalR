# Read all the compatible files in the given data location
# and convert these files to a list with organised data.

#' Load Rock-Eval data files
#'
#' This function loads all Rock-Eval raw data files into a list.
#'
#' @param dataloc Path to the input file
#' @return A list of the files in dataloc.
#' @export
read_RE <- function(dataloc){

  #1.1 list all files in data location
  files<-list.files(dataloc)

  #1.2 only take files with the right extension and exclude raw data back-up
  files<-files[stringr::str_detect(files,".B00")]
  files<-files[!stringr::str_detect(files,"~")]

  #2.1 create empty list for all data to be collected into
  file.list<-list()

  #2.2 for loop that reads the relevant data from the files into lists one by one
  for (i in 1:length(files)) {

    #2.2.1 read files
    file.i<-readLines(paste0(dataloc,"/",files[i]))

    #2.2.2 check if data is of the right type
    if (length(file.i)>1) {

      #2.2.3 description of the lines to cut
      file.i.rp1<-which(str_detect(file.i,"Curves Pyro"))+1 #row pyro start
      file.i.rp2<-which(str_detect(file.i,"Curves Oxi"))-1 #row pyro end
      file.i.ro1<-which(str_detect(file.i,"Curves Oxi"))+1 #row pyro start
      file.i.ro2<-length(file.i) #row pyro end


      #2.2.4 cut pyro part
      file.isp<-file.i[c(file.i.rp1:file.i.rp2)]
      file.isp<-file.isp[file.isp!=""]
      table.isp<-read.table(textConnection(file.isp))
      colnames(table.isp)<-c("t","T","CH","CO","CO2","SO2")

      #2.2.5 cut oxi part
      file.iso<-file.i[c(file.i.ro1:file.i.ro2)]
      file.iso<-file.iso[file.iso!=""]
      table.iso<-read.table(textConnection(file.iso))
      colnames(table.iso)<-c("t","T","CO","CO2","SO2")


    #2.2.6 check for Geoworks derived cursors and calibration info
    if (any(str_detect(file.i,"Geoworks"))){

      #2.2.7 description of the lines to cut
      file.i.rc1<-which(str_detect(file.i,"Curs manu_1")) #row Curs start
      file.i.rc2<-which(str_detect(file.i,"Calibration:Geoworks"))-2 #row Curs end

      file.i.curs<-file.i[file.i.rc1:file.i.rc2]
      rcv<-c(2:8,11:14,17:20,23:27,30:34,37:42,45:50) #rows with actual data
      file.i.curs<-str_extract(file.i.curs[rcv],"(?<==).*")
      names(file.i.curs)<-c("c1.1","c1.2","c1.3","c1.4","c1.5","c1.6","b1")

    }

    #2.2.6 escape when file has not been through Geoworks
    else{
      #2.2.8 combine pyro and oxi in list
      file.list[[i]]<-list(Pyro=table.isp,Oxi=table.iso)
      warning("No Geoworks cursors or calibration found.")
    }
    }

    #2.2.2 escape when data is incorrect
    else{
      file.list[[i]]<-"INCORRECT DATA ERROR"
      warning("Incorrect data found in one of the files!")
    }

  }
  names(file.list)<-files

  file.list
}





