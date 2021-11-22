# Author: Tunde Olarinoye, University of Freiburg
# Modified Brutseart recession extraction method for karst spring recession
# Reference: W. Brutseart (2008). Long-term groundwater storage trends estimated 
# from streamflow records: Climatic perspective. Water Res. Research.

# Methodology
# Decreasing spring discharge Q values
# exclude 3 days after peak discharge and 4 days in case of major events
# exclude subsequent recession points with difference more than 30% (default)

# dQ = 0.3; len = 7; plot=T; save.plot=NULL
# Qf = quantile(Qts$discharge, probs = 0.7, na.rm = T)
# Qts <- hydrogr[[2]]

# define function
remBrut <- function(Qts, Qf, dQ, len, plot, save.plot){
  
  # Input parameters
  # Qts = A cell array file in the format [date, Q, Qavg, dQdt] for each spring hydrograph; Date format = "yyyy-mm-dd"
  # Qf = Streamflow percentile for major events, default value set to 30%
  # dQ = Maximum percentage difference (-dQ/dt) for unspurious baseflow condition, default to 30%
  # len = Minimum length of recession event to retain for output
  # plot = logical input to plot extracted recession segments
  
  # Output dataframe
  # segList = List of dataframes of extracted recession segments ["Date","Q",Qavg, dQdt,"Qc","Qm","ti","index","Qro"]
  #-------------
  
  # smooth discharge time series with moving day average
  colnames(Qts) <- c("date","discharge","Qavg","dQdt")
  Qts$date <- as.Date(Qts$date, "%Y-%m-%d")
  
  # define variables
  seg <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(Qts))), c("Date","Q","Qavg","dQdt"))	#df for each extraxted recession segment
  segList <- list()
  
  a = 0
  
  # extract all recession segments of hydrograph
  # overview: select all recession points with dQdt < 0
  reces <- subset(Qts, subset = dQdt <= 0)
  reces[,1] <- as.Date(reces[,1], "%Y-%m-%d")
  colnames(reces) <- c("date","discharge","Qavg","subset_dQdt")
  nQts <- merge(Qts, reces, by="date", all.x = T)
  nQts <- nQts[, -c(5,6)]
  
  N <- nrow(nQts)
  
  for(i in 2:(N-2)){
    # i=35
    if(!is.na(nQts[i,5]) || (is.na(nQts[i,5]) && !is.na(nQts[i+1,5]) && 
                             !is.na(nQts[i,3] && !is.na(nQts[i-1,3] &&
                             abs((nQts[i-1,3] - nQts[i,3])/nQts[i-1,3]) < 0.05)) || 
                             (is.na(nQts[i,5]) && is.na(nQts[i+1,5]) && 
                             (!is.na(nQts[i,3]) && !is.na(nQts[i-1,3])) &&
                             abs((nQts[i-1,3] - nQts[i,3])/nQts[i-1,3]) < 0.05))  # additonal condtion to extract longer recession
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
    
    # x <- segList[[130]]
    # use restriction rules to filter conduit and matrix
    x$Qc <- NA; x$Qm <- NA
    N <- nrow(x)
    
    # detect major events
    flag <- which(x[,3] < Qf)[1]
    if(is.na(flag)){
      k = N
    }else{
      if(flag > 1){
        k = flag + 3
      }else{
        k = 3
      }
    }
    
    # separate conduit and matrix
    if(k >= N-1){
      x$Qc <- x[,3]
    }else{
      
      x$Qc[1:k] <- x[1:k, 3]
      k = k+1
      
      for(i in k:N){
        # if(abs((x[i-1, 3] - x[i, 3])/x[i-1, 3]) <= dQ){
        if(abs(x[i, 4]) < abs(x[i-1, 4]) || abs(x[i, 4]) <= (1+dQ) * abs(x[i-1, 4])){
          x$Qm[i] <- x[i, 3]
        }else{
          x$Qc[i] <- x[i, 3]
        }
      }
    }
    
    # clean-up: define start of continuous conduit and matrix recession
    for(i in 1:nrow(x)){
      # if(!is.na(x$Qm[i] && x$Qm[i+1] && x$Qm[i+2])){
      if(!is.na(x$Qm[i] && x$Qm[i+1])){
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
    plot.name <- "Recession segment extraction by Brutseart"
    remPlot(Qts=Qts, segList=segList, plot.name=plot.name, save.plot=save.plot)
    
  }
  
  return(segList)
}
 
  