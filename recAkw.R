# AKSOY AND WITTENBERG RECESSION EXTRACTION METHOD
# Aksoy and Wittenberg (2011) Nonlinear baseflow recession analysis in watershed with inttermittent streamflow
# Select receeding segment of hydrograph
# Select part of segment with CV =< 0.1

##-------------------------------------------------------------------------------------------------------------------------------------------
recAkw <- function(DATA, CoV, len, plot){
  #DATA = hydrogr; CoV = 0.01; len = 10; plot = T
  if(len<5){
    stop("Specified recession days less than 5")
  }
  
  # refine input data
  x <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(DATA)-1)), c("Date", "Q"))
  x[,1] <- DATA[2:(nrow(DATA)),1]
  x[,2] <- (DATA[1:(nrow(DATA)-1),2] + DATA[2:nrow(DATA),2])/2
  
  # define dataframes
  recSeg <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(x))), c("Date", "Q"))	#dataframe for each extraxted recession segment
  segList <- list(recSeg)
  RP <- setNames(data.frame(matrix(ncol = 4, nrow = nrow(x))), c("Date", "Q", "Qc", "Qm"))  #dataframe with conduit and matrix segement separated
  RPlist <- list(RP)  #list of all RP segments
  a = 0
  n = 0
  
  # select all recession segments
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
  
  # Activate the following lines to apply moving average 
  # to smooth timeseries before recession extraction
  #-----------------------------------------------------------
  #sort time series data, smoothing with moving day average
  # x = hydrogr; MA = 3; dQ = 0.3; lambda = 0.3; len = 30; plot = T
  # x$Date <- x[,1]
  # x$Q <- x[,2]
  # x$Qma <- filter(x[,2], rep(1/MA, MA), sides = 2, circular = T)
  # 
  # for(i in 1:(nrow(x)-1)){
  # 
  #   if(isTRUE(x$Qma[i] >= x$Qma[i+1]*1)){ # select values with decreasing moving average
  #     recSeg$Date[i] <- format(x[i, 1], format = "%Y-%m-%d")
  #     recSeg$Q[i] <- x$Q[i]
  #     #RP$Date[i] <- format(x$Date[i], format = "%Y-%m-%d")
  #     #RP$Q[i] <- x$Q[i]
  # 
  #   }else{
  #     a <- a+1
  #     recSeg <- recSeg[recSeg$Q != 0, ] # select only non-zero discharge value
  #     recSeg <- recSeg[!is.na(recSeg[,2]), ]     # remove NA values
  #     segList[[a]] <- recSeg
  #     recSeg <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(x))), c("Date", "Q"))
  # 
  #   }
  # }
  #---------------------------------------------------------
  
  coefvar <- function(x,y,D){
    N <- nrow(D)
    Q0 <- as.numeric(D[1,2])
    Q <- with(data=D, Q0*(1 + (1-y)*Q0^(1-y)/(x*y)*Index)^(1/(y-1)))
    coef <- (N^2/(N-1) * sum((D[,2]-Q)^2) / (sum(D[,2]))^2)^0.5
  }
  
  # fcal <- function(data, par){
  #   Q0 = Q0; #b = b
  #   alp <- par[1] #baseflow
  #   #b <- par[2]
  #   Qs <- with(data, Q0*(1 + ((1-b)*Q0^(1-b)/(alp*b))*Index)^(1/(b-1))) 
  #   rss <- with(data, lsfit(x = Q, y = Qs))$coefficients[2]
  #   #rss <- with(data, sum((Q - Qs)^2))
  #   return(rss)
  # }
  
  
  # separate segments to quick and slow flow components
  for(j in 1:length(segList)){
    #j=414
    Nr <- nrow(segList[[j]])
    b <- 1
    
    if(Nr>5){
      R <- segList[[j]]
      
      for(k in 1:nrow(R)){
        rp <- R[k:Nr,]
        N <- nrow(rp)
        CV <- 10
        
        if(N<5){
          break
        }
        
        # rp$Index <- seq_along(rp$Date)-1
        rp$Index <- seq(0, N-1, 1)
        Q0 <- as.numeric(rp[1,2])
        a <- sum(rp[1:(N-1),2] + rp[2:N,2])/(2*sum(rp[1:(N-1),2]^b - rp[2:N,2]^b))
        CV <- coefvar(x=a, y=b, D=rp)
        
        # z <- summary(nls(Q ~ Q0*exp(-alpha*Index), data = rp, start = list(alpha=0.001)))
        # alpha <- coef(z)["alpha","Estimate"]
        # f_opt = optim(par = c(0.1,0.9), fn = fcal, data = rp, method = "BFGS", 
        #               control = c(trace = T, maxit = 30000))
        # alp = f_opt$par[1]
        # b = f_opt$par[2]
        # bflow <- with(rp, Q0*exp(-alpha*Index))
        # bflow <- with(rp, Q0*(1 + ((1-b)*Q0^(1-b)/(alp*b))*Index)^(1/(b-1)))
        # CV <- sqrt((N^2 * sum((rp[,2]-bflow)^2)) / ((N-1) * (sum(rp[,2]))^2))
        # 
        
        # select segment when CV <= CoV and break loop
        if(CV <= CoV && N >= 5){
          RP <- setNames(data.frame(matrix(ncol = 5, nrow = Nr)), c("Date", "Q", "Qc", "Qm", "Vt"))
          RP$Date <- R$Date
          RP$Q <- R$Q
          RP$Qm[k:Nr] <- R$Q[k:Nr]
          RP$Qc[1:k] <- R$Q[1:k]
          RP$Vt <- R$Vt
          
          # fit linear model to matrix recession component to get intercept Qro, required for mangin model
          RP$Index <- seq_along(RP$Date)-1
          z <- summary(lm(Q~Index, data = RP[k:N, ]))
          RP$Qro <- z$coefficients[1,1]
          RP$ti <- k-1
          RP$Date <- as.Date(RP$Date)
          
          # append segment to list of recessions
          n <- n+1
          RPlist[[n]] <- RP
          break
        }else{
          next
        }
      }
      
    }else{
      next
    }
    
  }
  
  # further clean up of recession list, chunk can be removed!
  cumRP <- list(); b = 0
  for(i in 1:length(RPlist)){
    if(nrow(RPlist[[i]]) >= len && !is.na(RPlist[[i]][1,6])){
      b = b+1
      cumRP[[b]] <- RPlist[[i]]
    }else{
      next
    }
  }
  
  # plot recession segments
  if(plot == T){
    plot.name <- "Recession segment extraction by Aksoy"
    recPlot(x=x, cumRP=cumRP, plot.name=plot.name)
  }
  return(cumRP)	
}
