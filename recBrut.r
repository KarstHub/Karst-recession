# BRUTSEART RECESSION EXTRACTION METHOD
 # decreasing spring discharge Q values
 # exlude 3 days after peak discharge and 4 days in case of major events (Qf)
 # 30% quantile was used for default Qf
 # exlude recession point with difference more than dQ (30% default)
  ##-------------------------------------------------------------------------------------------------------------------------------------------
  recBrut <- function(DATA, Qf, dQ, len, plot){
    x <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(DATA)-1)), c("Date", "Q"))
    x[,1] <- DATA[2:(nrow(DATA)),1]
    x[,2] <- (DATA[1:(nrow(DATA)-1),2] + DATA[2:nrow(DATA),2])/2
    
    # define variables
    recSeg <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(x))), c("Date", "Q"))	#dataframe for each extraxted recession segment
    segList <- list(recSeg)
    RP <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(x))), c("Date", "Q", "Qc", "Qm"))  #dataframe with conduit and matrix segement separated
    RPlist <- list(RP)  #list of all RP segments
    a = 0
    b = 0
    
    for(i in 1:nrow(x)){
      if(isTRUE(as.numeric(x$Q[i+1] - x$Q[i]) < 0)){  # select only value of dQ/dt < 0
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
      
      if(isTRUE(N > 4)){
        
        RP <- setNames(data.frame(matrix(ncol = 5, nrow = N)), c("Date", "Q", "Qc", "Qm", "Vt"))
        RP$Date[1:N] <- segList[[j]][,1]
        RP$Q[1:N] <- segList[[j]][,2]
        RP$Vt <- segList[[j]][,3]
        
        
        if(isTRUE(RP$Q[1] <= Qf)){  # check if Q exceed percentile for major event
          jstart <- 3
          
        }else{
          jstart <- 4
        }
        
        RP$Qc[1:jstart] <- RP$Q[1:jstart]
        
        for(k in (jstart+1):nrow(RP)){  # remove last 2 data points
          if(isTRUE((RP$Q[k] - RP$Q[k-1])/RP$Q[k-1] >= -dQ)){
            RP$Qm[k] <- RP$Q[k]
            
          }else{
            RP$Qc[k] <- RP$Q[k]
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
      plot.name <- "Recession segment extraction by Brutseart"
      recPlot(x=x, cumRP=cumRP, plot.name=plot.name)
    }
    return(cumRP)	
  }