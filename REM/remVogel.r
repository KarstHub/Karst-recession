# Author: Tunde Olarinoye, University of Freiburg
# Modified Vogel recession extraction method for karst spring recession
# Reference: Vogel and Kroll(1992), Regional geoghydrologic- geomorphic relationships
# for the estimation of low-flow statistics. Water Res. Research.

# Methodology
# Decreasing spring discharge Q based on a 3 days (default) moving average
# exclusion: if Qi - Qi+1 / Qi+1 > 30% (default) - for low-flow selection
# filter: removes first 30% data points - for low flow selection

# Qts = Qts; MA = 5; dQ = 0.3; lambda = 0.3; len = 7; plot = T

# define function
  remVogel <- function(Qts, MA, dQ, lambda, len, plot, save.plot){
    
    # Input parameters
    # Qts = A cell array file in the format [date, Q, Qavg, dQdt] for each spring hydrograph; Date format = "yyyy-mm-dd"
    # MA = Number of moving average days for smoothing the hydrograph (default is 3 days)
    # dQ = Maximum percentage difference between consecutive discharge value, from 0 to 100
    # lambda = Fraction of recession segment influenced by spurious flow to be removed, 0.3 by default
    # len = Minimum length of recession event to retain for output
    # plot = logical input to plot extracted recession segments
    
    # Output dataframe
    # segList = List of dataframes of extracted recession segments ["Date","Q","Qc","Qm","ti","index","Qro"]
    #-------------

    # smooth discharge time series with moving day average
    colnames(Qts) <- c("date","discharge","Qavg","dQdt")
    Qts$date <- as.Date(Qts$date, "%Y-%m-%d")
    Qts$Qma <- filter(Qts[,3], rep(1/MA, MA), sides = 1, circular = F)
    Qts <- Qts[-seq(1,MA,1), ]
    
    # define variables
    seg <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(Qts))), c("Date","Q","Qavg","dQdt"))	#df for each extraxted recession segment
    segList <- list()
    
    a = 0
    
    # extract all recession segments of hydrograph
    # overview: select all recession points after smoothing with moving days
    reces <- subset(Qts, subset = diff(Qts$Qma) <= 0)
    reces[,1] <- as.Date(reces[,1], "%Y-%m-%d")
    colnames(reces) <- c("date","discharge","Qavg","subset_dQdt")
    nQts <- merge(Qts, reces, by="date", all.x = T)
    nQts <- nQts[, -c(5,6)]
	
	  N <- nrow(nQts)
    
    for(i in 1:(N-2)){
      # i=1
      if(!is.na(nQts[i,5]) || (is.na(nQts[i,5]) && !is.na(nQts[i+1,5]))
	      || (is.na(nQts[i,5]) && is.na(nQts[i+1,5]) && (!is.na(nQts[i+1,3]) && !is.na(nQts[i+2,3]))
	          && abs((nQts[i+1,3] - nQts[i+2,3])/nQts[i+1,3]) < 0.05)   # additonal condtion to extract longer recession
	  ){
        seg$Date[i] <- format(nQts[i, 1], "%Y-%m-%d")
        seg$Q[i] <- nQts[i, 2]
		seg$Qavg[i] <- nQts[i, 3]
        seg$dQdt[i] <- nQts[i, 4]
		
      }else{
        seg <- seg[seg$Q != 0, ]          # select only non-zero discharge value
        seg <- seg[!is.na(seg[,3]), ]     # remove NA values
        if(nrow(seg) > len){              # select only segments with specify number of dates
          a <- a+1
          segList[[a]] <- seg
          seg <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(Qts))), c("Date","Q","Qavg","dQdt"))
        }else{
          seg <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(Qts))), c("Date","Q","Qavg","dQdt"))
          next
        }
      }
    }
	
	#---------
    segList <- lapply(segList, function(x){
      for(i in 1:nrow(x)){
        if(x[i,3] < x[i+1, 3]){
          x <- x[-i, ]
        }else{
          break
        }
      }
      return(x)
    })
  
    #---------
   
    # separate conduit and matrix segment
    segList <- lapply(segList, function(x){
      
      # use restriction rules to filter conduit and matrix
      x$Qc <- NA; x$Qm <- NA
      N <- round(lambda * nrow(x))
      x$Qc[1:N] <- x[1:N, 3]
      
      for(i in (N+1):nrow(x)){
        if((x[i-1, 3] - x[i, 3])/x[i-1, 3] <= dQ && (x[i-1, 3] - x[i, 3])/x[i-1, 3] > -0.1){
          x$Qm[i] <- x[i, 3]
        }else{
          x$Qc[i] <- x[i, 3]
        }
      }
      
      # clean-up: define start of continuous conduit and matrix recession
      for(i in 1:nrow(x)){
        if(!is.na(x$Qm[i] && x$Qm[i+1] && x$Qm[i+2])){
          x$Qm[i:nrow(x)] <- x[i:nrow(x), 3]
          x$Qc[i:nrow(x)] <- NA
          break
        }else{
          x$Qm[i] <- NA
          x$Qc[i] <- x[i, 3]
        }
      }
      
      # get when conduit drainage stops (days)
      ti <- which(!is.na(x$Qm))[1]
      x$ti <- ti
      return(x)
    })
    
    #---------------
    # fit lm to matrix segment to get intercept, Qro - required for Mangin model
    segList <- lapply(segList, function(x) {
      # x <- segList[[2]]
      
      x$index <- seq_along(x$Date)-1
      ti <- x$ti[1]; N <- nrow(x)
      if(!is.na(ti)){
        mat.lm <- summary(lm(Qavg~index, data = x[ti:N, ]))
        x$Qro <- mat.lm$coefficients[1,1]
      }else{
        x$Qro <- NA
      }
      return(x)
    })

    # create plot
    if(plot == T){
	   plot.name <- "Recession segment extraction by Vogel"
	   remPlot(Qts=Qts, segList=segList, plot.name=plot.name, save.plot=save.plot)
	
    }
    
    return(segList)
  }
  
  