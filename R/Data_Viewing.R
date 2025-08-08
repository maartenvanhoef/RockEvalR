#' Plot Rock-Eval thermograms
#'
#' This function previews Rock-Eval thermograms and cursors. Input contains a single Rock-Eval analysis.
#'
#' @param datalist List containing data of one or more Rock-Eval analyses
#' @param analysis Name of file containing data of one Rock-Eval analysis
#' @param oven Choice of pyrolysis or oxidation oven
#' @param thermogram Choice of selected thermogram
#' @return A plot of thermogram(s)
#' @export
RE_plot <- function(datalist, analysis, oven="both", thermogram="all"){

    #1.1 Check names

    # Error when impossible combination is made.
    if (oven=="Oxidation" & thermogram=="CH"){
      stop("A CH pyrogram doesn't exist in the oxidation oven.")
    }

    # Error when inappropriate names are given
    if (!(oven %in% c("both","Pyrolysis","Oxidation"))){
      stop("Only appropriate values for oven are both, Pyrolysis, and Oxidation.")
    }

    # Error when inappropriate names are given
    if (!(thermogram %in% c("all","CH","CO","CO2","SO2"))){
      stop("Only appropriate values for oven are all, CH, CO, CO2, and SO2.")
    }

  #1.2 Get selected data
  plotdata<-datalist[[analysis]]

  #1.3 determine dimensions
    # oven selection
    if (oven=="both"){
      oven.s<-c("Pyrolysis","Oxidation")
    } else{
      oven.s<-oven
    }

    # thermogram selection
    if (thermogram=="all"){
      thermo.s<-c("CH","CO","CO2","SO2")
      cross.s<-expand.grid(thermo.s,oven.s)

    } else{
      thermo.s<-thermogram
      cross.s<-expand.grid(thermo.s,oven.s)
    }

    # Par
    par(mfrow=c(length(oven.s),length(thermo.s)))

  #1.4 plot RE thermogram data

  for (i in 1:nrow(cross.s)) {
    #1.4.1 main plotting
    thermo.i<-as.character(cross.s[i,1])
    oven.i<-as.character(cross.s[i,2])

    plot(plotdata[[oven.i]][["t"]],plotdata[[oven.i]][[thermo.i]],
         type="l", lwd=2, col="gray11",
         xlab = "Time", ylab="RE Signal")

    #1.4.2 determine what cursors to include from LU table
    curs.s<-cursLU$curs[cursLU$oven %in% oven.i & cursLU$thermo %in% thermo.i & cursLU$no!=0]
    base.s<-cursLU$curs[cursLU$oven %in% oven.i & cursLU$thermo %in% thermo.i & cursLU$no==0]

    select.curs<-plotdata$Cursors[curs.s]
    select.base<-plotdata$Cursors[base.s]

    #1.4.3 plot cursors and baseline
    invisible(mapply(function(curs) abline(v=curs, col="red3"), curs=select.curs))
    abline(h=select.base, col="royalblue4")

    #1.4.4 add description text
    mtext(paste0(oven.i,"-",thermo.i),side=3)
  }

  }




