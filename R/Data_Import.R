# Read all the compatible files in the given data location
# and convert these files to a list with organised data.

#' Load Rock-Eval data files
#'
#' This function loads all Rock-Eval raw data files into a list.
#'
#' @param dataloc Path to the input file
#' @return A list of the files in dataloc.
#' @export
RE_read <- function(dataloc){

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

      #2.2.7 description of the lines to cut for cursors
      file.i.rc1<-which(str_detect(file.i,"Curs manu_1")) #row Curs start
      file.i.rc2<-which(str_detect(file.i,"Calibration:Geoworks"))-2 #row Curs end

      #2.2.8 extract cursors
      file.i.curs<-file.i[file.i.rc1:file.i.rc2]
      rcv<-c(2:8,11:14,17:20,23:27,30:34,37:42,45:50) #rows with actual data
      file.i.curs<-as.numeric(str_extract(file.i.curs[rcv],"(?<==).*"))
      names(file.i.curs)<-c(
        paste0("curs1.",1:6),
        "base1",
        paste0("curs2.",1:3),
        "base2",
        paste0("curs3.",1:3),
        "base3",
        paste0("curs4.",1:4),
        "base4",
        paste0("curs5.",1:4),
        "base5",
        paste0("curs6.",1:5),
        "base6",
        paste0("curs7.",1:5),
        "base7")

      #2.2.9 description of the lines to cut for the calibration
      # currently only KID and Tpeak are not pre-calibrated
      file.i.rcal1<-which(str_detect(file.i,"P_HC_K1=")) #row calibration 1
      file.i.rcal2<-which(str_detect(file.i,"TempSTD=")) #row calibration 2
      file.i.rcal3<-which(str_detect(file.i,"TempFID=")) #row calibration 3

      #2.2.10 extract calibration
      file.i.cals<-file.i[c(file.i.rcal1,file.i.rcal2,file.i.rcal3)]
      file.i.cals<-as.numeric(str_extract(file.i.cals,"(?<==).*"))
      names(file.i.cals)<-c("cal_FID","cal_Tmax","cal_Tpeak")

    }    else{
      #2.2.6b escape when file has not been modified by Geoworks

      warning(paste0("No Geoworks cursors or calibration found for: ",
                     files[i]," setting values to NA."))

      #2.2.8b set cursors to NA
      file.i.curs<-rep(NA,37)
      names(file.i.curs)<-c(
        paste0("curs1.",1:6),
        "base1",
        paste0("curs2.",1:3),
        "base2",
        paste0("curs3.",1:3),
        "base3",
        paste0("curs4.",1:4),
        "base4",
        paste0("curs5.",1:4),
        "base5",
        paste0("curs6.",1:5),
        "base6",
        paste0("curs7.",1:5),
        "base7")

      #2.2.10b set calibration to NA
      file.i.cals<-rep(NA,3)
      names(file.i.cals)<-c("cal_FID","cal_Tmax","cal_Tpeak")
    }

      #2.2.11 read sample parameter info
      file.i.ri1<-which(str_detect(file.i,"Param"))+1 #row parameter info start
      file.i.ri2<-which(str_detect(file.i,"Param"))+12 #row parameter info end

      file.parameters<-str_extract(file.i[file.i.ri1:file.i.ri2],"(?<==).*")
      names(file.parameters)<-str_extract(file.i[file.i.ri1:file.i.ri2],".*(?==)")

      #2.2.12 read coefficient info
      file.i.rf1<-which(str_detect(file.i,"Coeff"))+1 #row coefficient info start
      file.i.rf2<-which(str_detect(file.i,"Coeff"))+33 #row coefficient info end

      file.coefficients<-as.numeric(
        str_extract(file.i[file.i.rf1:file.i.rf2],"(?<==).*"))
      names(file.coefficients)<-str_extract(file.i[file.i.rf1:file.i.rf2],".*(?==)")
      # remove 'empty gas' coefficients
      file.coefficients<-file.coefficients[!str_detect(names(file.coefficients),
                                                             "GAS")]

      #2.2.13 combine the thermograms, cursors and calibration into a list
      file.list[[i]]<-list(Parameters=file.parameters,
                           Pyrolysis=table.isp,Oxidation=table.iso,
                           Cursors=file.i.curs, Calibration=file.i.cals,
                           Coefficients=file.coefficients)

    }else{
      #2.2.2b escape when data is incorrect

      file.list[[i]]<-"INCORRECT DATA ERROR"
      warning(paste0("Incorrect data found in file: ",files[i],"."))
    }

  }
  names(file.list)<-files

  file.list
}





