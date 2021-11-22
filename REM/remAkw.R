# Author: Tunde Olarinoye, University of Freiburg
# Modified Aksoy and Wittenberg recession extraction method for karst spring recession
# Reference: Aksoy and Wittenberg (2011) Nonlinear baseflow recession analysis in 
# watershed with inttermittent streamflow

# Methodology
# Select receeding segment of hydrograph and filter out segment with CV =< 0.1

# Qts=Qts; CoV = 0.1; len = 10; plot = T

# define function
remAkw <- function(Qts, CoV, len, plot, save.plot){
  
  # Input parameters
  # Qts = A cell array file in the format [date, Q, Qavg, dQdt] for each spring hydrograph; Date format = "yyyy-mm-dd"
  # CoV = maximum coefficient of variation (CV) allowable for spring discharge [0:1]
  # len = Minimum length of recession event to retain for output
  # plot = logical input to plot extracted recession segments
  
  # Output dataframe
  # segList = List of dataframes of extracted recession segments ["Date","Q","Qc","Qm","ti","index","Qro"]
  #-------------
  
  # recession segment must be longer than 5 days
  if(len<5){
    stop("Specified recession days less than 5")
  }
  
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
  
  for(i in 1:(N-2)){
    
    if(!is.na(nQts[i,5]) || (is.na(nQts[i,5]) && !is.na(nQts[i+1,5]))
       # || (is.na(nQts[i,5]) && is.na(nQts[i+1,5]) && !is.na(nQts[i+2,5]))   
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
  
  coefvar <- function(a,b,D){
    nD <- nrow(D)
    Q0 <- as.numeric(D[1,3])
    Q <- with(data=D, Q0*(1 + ((1-b)*Q0^(1-b)/(a*b))*Index)^(1/(b-1)))
    coef <- (nD^2/(nD-1) * sum((D[,3]-Q)^2) / (sum(D[,3]))^2)^0.5
  }
  
  # separate conduit and matrix segment
  segList <- lapply(segList, function(x){
    
    # x <- segList[[21]]
    # use restriction rules to filter conduit and matrix
    x$Qc <- NA; x$Qm <- NA
    N <- nrow(x)
    b <- 1
    
    # slowflow (matrix) separation
    for(i in 1:nrow(x)){
      
      if(N-i > 5){
        
        rp <- x[i:N,]
        Nrp <- nrow(rp)
        
        rp$Index <- seq(0, Nrp-1, 1)
        Q0 <- as.numeric(rp[1,3])
        a <- sum(rp[1:(Nrp-1),3] + rp[2:Nrp,3])/(2*sum(rp[1:(Nrp-1),3]^b - rp[2:Nrp,3]^b))
        CV <- coefvar(a=a, b=b, D=rp)
        
        if(CV <= CoV){
          x$Qc[1:i] <- x[1:i, 3]
          x$Qm[(i+1):N] <- x[(i+1):N, 3]
          break
          
        }else{
          next
        }
        
      }else{
        x$Qc <- x[ ,3]
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

# fcal <- function(data, par){
#   Q0 = Q0; #b = b
#   alp <- par[1] #baseflow
#   #b <- par[2]
#   Qs <- with(data, Q0*(1 + ((1-b)*Q0^(1-b)/(alp*b))*Index)^(1/(b-1))) 
#   rss <- with(data, lsfit(x = Q, y = Qs))$coefficients[2]
#   #rss <- with(data, sum((Q - Qs)^2))
#   return(rss)
# }

# z <- summary(nls(Q ~ Q0*exp(-alpha*Index), data = rp, start = list(alpha=0.001)))
# alpha <- coef(z)["alpha","Estimate"]
# f_opt = optim(par = c(0.1,0.9), fn = fcal, data = rp, method = "BFGS", 
#               control = c(trace = T, maxit = 30000))
# alp = f_opt$par[1]
# b = f_opt$par[2]
# bflow <- with(rp, Q0*exp(-alpha*Index))
# bflow <- with(rp, Q0*(1 + ((1-b)*Q0^(1-b)/(alp*b))*Index)^(1/(b-1)))
# CV <- sqrt((N^2 * sum((rp[,2]-bflow)^2)) / ((N-1) * (sum(rp[,2]))^2))


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
          
          
       
      
        
        
  

