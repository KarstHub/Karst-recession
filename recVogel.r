 # VOGEL RECESSION EXTRACTION METHOD
  # decreasing spring discharge Q on a 3 days (default) moving average
  # exclusion: if Qi - Qi+1 / Qi+1 > 30% (default)
  # filter: removes first 30% data points
  
  recVogel <- function(x, MA, dQ, lambda, len, plot){
    #x = hydrogr; MA = 3; dQ = 0.3; lambda = 0.3; len = 14; plot = T

    # sort time series data, smoothing with moving day average
    x$Date <- x[,1]
    x$Q <- x[,2]
    x$Qma <- filter(x[,2], rep(1/MA, MA), sides = 2, circular = F)
    
    # define variables
    recSeg <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(x))), c("Date", "Q"))	#dataframe for each extraxted recession segment
    segList <- list(recSeg)
    RP <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(x))), c("Date", "Q", "Qc", "Qm"))	#dataframe with conduit and matrix segement separated
    RPlist <- list(RP)	#list of all RP segments
    a = 0
    b = 0
    # extract all recession segments of hydrograph
    # overiew: select all recession points after smoothing with moving days
    xReces <- subset(x, subset = diff(x$Qma) <= 0)
    
    for(i in 1:nrow(x)){
      
      if(isTRUE(x$Qma[i] >= x$Qma[i+1])){ # select values with decreasing moving average
        recSeg$Date[i] <- format(x[i, 1], format = "%Y-%m-%d")
        recSeg$Q[i] <- x$Q[i]
        
      }else{
        
        dates <- as.Date(x[i, 1], format = "%Y-%m-%d")
        yr <- as.numeric(format(dates, format = "%Y"))
        mnt <- format(dates, format = "%m")
        x$date <- as.Date(x$Date, format = "%Y-%m-%d")
        
        # calculate discharge volume for hydrological year
        if(mnt < 11){
          begin <- as.Date(paste0(yr-1,"-11-01"))
          end <- as.Date(paste0(yr,"-10-31"))
          b2e <- seq.Date(begin, end, "day")
          hyd_yr <- subset(x, date %in% b2e)
          recSeg$Vt <- (sum(hyd_yr[,2])/nrow(hyd_yr))*365*86400
        }else{
          begin <- as.Date(paste0(yr,"-11-01"))
          end <- as.Date(paste0(yr+1,"-10-31"))
          b2e <- seq.Date(begin, end, "day")
          hyd_yr <- subset(x, date %in% b2e)
          recSeg$Vt <- (sum(hyd_yr[,2])/nrow(hyd_yr))*365*86400
        }
        
        a <- a+1
        recSeg <- recSeg[recSeg$Q != 0, ] # select only non-zero discharge value
        recSeg <- recSeg[!is.na(recSeg[,2]), ]     # remove NA values
        segList[[a]] <- recSeg
        recSeg <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(x))), c("Date", "Q"))
        
        }
		  }
		
		for(j in 1:length(segList)){
		  N <- nrow(segList[[j]]) 
		  lambdaDays <- round(lambda * N) #% of days to be removed as number of days
		  
		  if(isTRUE(lambdaDays > 1)){
		    RP <- setNames(data.frame(matrix(ncol = 5, nrow = N)), c("Date", "Q", "Qc", "Qm", "Vt"))
		    RP$Date[1:N] <- segList[[j]][,1]
		    RP$Q[1:N] <- segList[[j]][,2]
		    RP$Qc[1:lambdaDays] <- segList[[j]][1:lambdaDays, 2]
		    RP$Vt <- segList[[j]][,3]
		    
		    for(k in (lambdaDays+1):nrow(segList[[j]])-1){
		      if(isTRUE((segList[[j]][k, 2] - segList[[j]][k+1, 2])/segList[[j]][k, 2] <= dQ)){
		        RP$Qm[k] <- segList[[j]][k, 2]
		      }else{
		        RP$Qc[k] <- segList[[j]][k, 2]
		      }
		        
		    }
		    
		    b <- b+1
		    m <- which(!is.na(RP$Qm))[1] #get postion matrix drainage starts, its assume as mixed drainage position
		    RP$Qc[m] <- RP$Qm[m]
		    
		    # fit a straight linear to matrix recession to get intercept Qro, required for mangin model
		    RP$Index <- seq_along(RP$Date)
		    if(!is.na(m)){
		      z <- summary(lm(Q~Index, data = RP[m:N, ]))
		      RP$Qro <- z$coefficients[1,1]
		    }else{
		      RP$Qro <- NA
		    }
		    
		    RP$ti <- m
		    RP$Date <- as.Date(RP$Date)
		    RPlist[[b]] <- RP
		      
		  }else{
		    next
		    
		  }
		}

    # retain only recession segments longer then len (5 days default) days
    #cumRP <- lapply(RPlist, function(x){
    #  if(nrow(x) >= len && !is.na(x[1,6])) # minimum length
    #    return(x)
    #  else 
    #    next
    #})
    cumRP <- list(); b = 0
    for(i in 1:length(RPlist)){
      if(nrow(RPlist[[i]]) >= len && !is.na(RPlist[[i]][1,6])){
        b = b+1
        cumRP[[b]] <- RPlist[[i]]
      }else{
        next
      }
    }
    
    if(plot == T){
	   plot.name <- "Recession segment extraction by Vogel"
	   recPlot(x=x, cumRP=cumRP, plot.name=plot.name)
	
    }
    return(cumRP)
 
 }