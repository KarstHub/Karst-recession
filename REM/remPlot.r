# Author: Tunde Olarinoye, University of Freiburg
# plot function for hydrograph, conduit and matrix recession curves
# after extracting the recession segment of conduit and matrix
# this function is used to plot the Q time series and recession segments

remPlot <- function(Qts, segList, plot.name, save.plot){
  
  if(!is.null(save.plot)){
    pdf(paste0("./Results_study2/Reces_segments/",save.plot, ".pdf"), height=3.0, width=4.5)
	par(mgp=c(1.0,0.1,0), cex.main=0.5,cex.lab=0.5, cex.axis=0.5, las=1, tck=0.05)
    xlim <- c(min(Qts$date), max(Qts$date))
    ylim <- c(min(na.omit(Qts$Qavg)), max(na.omit(Qts$Qavg)))
    
    plot(Qts$date, Qts$Qavg, xlim = xlim, ylim = ylim, main = save.plot, type = "l", lwd = 0.5, xlab = "Date", ylab = "Q (m3/day)")
    
    for(i in 1:length(segList)){
      seg <- segList[[i]] 
      
      if(nrow(seg) > 0){
        lines(as.Date(seg$Date), seg$Qavg, lwd=1.0, col="powderblue")
        segQc = seg[which(!is.na(seg$Qc)), ] 
        lines(as.Date(segQc$Date), segQc$Qc, lwd=1.0, col="purple")
        segQm = seg[which(!is.na(seg$Qm)), ] 
        lines(as.Date(segQm$Date), segQm$Qm, lwd=1.0, col="orange")
        legend("topright", legend=c("Measured spring discharge","Entire recession segment","Conduit recession","Matrix recession")
               ,lty=1, col=c("black","powderblue","purple","orange"), cex=0.5, bty="n")
      }
      else next
      
    }
	dev.off()
  }
  
  
  dev.new()
  # pdf(paste0("./Results_study2/Reces_segments/",plot.name, ".pdf"), height=4.0, width=6.5)
  xlim <- c(min(Qts$date), max(Qts$date))
  ylim <- c(min(na.omit(Qts$Qavg)), max(na.omit(Qts$Qavg)))
  
  plot(Qts$date, Qts$Qavg, xlim = xlim, ylim = ylim, main = plot.name, type = "l", lwd = 1.0, xlab = "Date", ylab = "Q (m3/day)")
  
  for(i in 1:length(segList)){
    seg <- segList[[i]] 
    
    if(nrow(seg) > 0){
      lines(as.Date(seg$Date), seg$Qavg, lwd=3, col="powderblue")
      segQc = seg[which(!is.na(seg$Qc)), ] 
      lines(as.Date(segQc$Date), segQc$Qc, lwd=1.5, col="purple")
      segQm = seg[which(!is.na(seg$Qm)), ] 
      lines(as.Date(segQm$Date), segQm$Qm, lwd=1.5, col="orange")
      legend("topright", legend=c("Measured spring discharge","Entire recession segment","Conduit recession","Matrix recession")
             ,lty=1, col=c("black","powderblue","purple","orange"), cex=0.5, bty="n")
    }
    else next
    
  }
  
}