# plot function for hydrograph, conduit and matrix recession curves
# after extracting the recession segment of conduit and matrix
# this function is used to plot the Q time series and recession segments

recPlot <- function(x, cumRP, plot.name){
      dev.new()
      xlim <- c(min(as.Date(x$Date,format = "%Y-%m-%d")), max(as.Date(x$Date,format = "%Y-%m-%d")))
      ylim <- c(min(na.omit(x$Q)), max(na.omit(x$Q)))
      plot(as.Date(x$Date,format = "%Y-%m-%d"), x$Q, xlim = xlim, ylim = ylim, main = plot.name, type = "l", lwd = 2, xlab = "Date", ylab = "Q")
      
      for(i in 1:length(cumRP)){
        R <- cumRP[[i]] 
        R$Date <- format(R$Date, format = "%Y-%m-%d")
        
        if(isTRUE(nrow(R) > 0)){
          par(new = T)
          #RC <- R[!is.na(R$Qc), ]
		  RC <- R[, c(1,3)]
          plot(as.Date(RC$Date, format = "%Y-%m-%d"), RC$Qc, xlim = xlim, ylim = ylim, type = "l", lwd = 2, col = "red", xlab = "",ylab = "", axes = F)
          par(new = T)
          #RM <- R[!is.na(R$Qm), ]
		  RM <- R[, c(1,4)]
          plot(as.Date(RM$Date, format = "%Y-%m-%d"), RM$Qm, xlim = xlim, ylim = ylim, type = "l", lwd = 2, col = "blue", xlab ="", ylab = "", axes = F)
        }
        else next
        
      }
    }