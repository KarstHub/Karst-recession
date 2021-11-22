##---------------------------------------------------------------------------------------------------
# Mangin model for karst spring hydrograph recession analysis
# Model use two different types of functions: ψt representing contribution from unsaturated zone
# and ϕt representing contribution from the saturated zone (Mangin 1975). Spring dsicharge can be
# modelled by: Qt = ϕt + ψt; ϕt is the baseflow component modelled with Maillet exponential model
# and ψt is the flood recession component modelled with a non-exponential model.
##---------------------------------------------------------------------------------------------------
# The function fit Mangin model to extracted recession hydrographs to calibraterecession parameters. 
# This is a one-step optimisation procedure. The three recession parameters of matrix and conduit
# models are optimised together. Note that N which is 1/ti is not fixed here
# Define function, parameters and variables
# Define function, parameters and variables
source("./Recession analysis/Recession model/KGE.R")
#source("./Recession analysis/Recession model/setMin.R")
mangin4 <- function(RP, Q0 = NULL, plotSim = TRUE, plotName = NULL, orig = 0){
  #Q0 = NULL; #Qro = NULL; q0 = NULL
  #RP <- tryme
  # RP is the extracted recession points in data frame ["Date","Q","Qc","Qm","Qro","ti"]
  # Q0 corresponds to peak discharge at time t = 0
  # Qro corresponds to baseflow component of discharge at inflexion point (change from fast to slow flow)
  # q0 is the diffrence between peak discharge Q0 and Qro
  
  # define output variables
  coeff <- data.frame("alpha"="", "beta"="")
  
  if(is.null(Q0)){
    Q0 = RP[1,2]
  }

  RP$t_index <- seq(0, nrow(RP)-1, 1)
  ti <- which(!is.na(RP$Qm))[1] 
  NR <- nrow(RP)
  
  setMin2 <- function(expr){
    pmax(0, eval(expr))
  }
  
  # optimise recession model
  Nup <- 1/(NR-2)
  qfb = function(data, par){
    Q0 = Q0
    alpha <- par[1]
    N <- par[2]
    E <- par[3]
    Qr <- par[4]
    Qs <- with(data, Qr * exp(-(alpha)*t_index) + 
                 setMin2((Q0 - Qr) * ((1 - (N*t_index))/(1 + (E*t_index)))))
    rss <- with(data, sum((Q - Qs)^2))
    return(rss)
  }
  q <- 0.3*Q0
  qfbFit = optim(par = c(0.1,0.5,0.1,q), fn = qfb, data = RP, method = "L-BFGS-B", 
                 lower = c(0,Nup,0,0), upper = c(Inf,1,10,Q0), control = c(trace = T, maxit = 30000))
  alpha = qfbFit$par[1]
  N = qfbFit$par[2]; tm <- round(1/N)
  E = qfbFit$par[3]
  Qr = qfbFit$par[4]
  rss = qfbFit$value
  q0 = Q0 - Qr
  
  # simulated Q
  RP$simQ = with(data = RP, Qr * exp(-(alpha)*t_index) + setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
  simQb <- with(data = RP, Qr * exp(-(alpha)*t_index))
  simQf <- with(data = RP, setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
  kge <- KGE(RP$Q,RP$simQ)
  
  
  #---------------------------------------------------------
  # if non-modified extracted REM points are considered
  if(orig == 1){
    qb_df <- RP[ti:NR, c("Q","t_index")]
    qb_df[,2] <- qb_df$t_index - qb_df$t_index[1]
    Q0 <- qb_df$Q[1]
    
    qfb = function(data, par){
      Q0 = Q0
      alpha <- par[1]
      N <- par[2]
      E <- par[3]
      Qr <- par[4]
      Qs <- with(data, Qr * exp(-(alpha)*t_index) + 
                   setMin2((Q0 - Qr) * ((1 - (N*t_index))/(1 + (E*t_index)))))
      rss <- with(data, sum((Q - Qs)^2))
      return(rss)
    }
    q <- 0.3*Q0
    qfbFit = optim(par = c(0.1,0.5,0.1,q), fn = qfb, data = qb_df, method = "L-BFGS-B", 
                   lower = c(0,0,0,0), upper = c(Inf,1,10,Q0), control = c(trace = T, maxit = 30000))
    alpha = qfbFit$par[1]
    N = qfbFit$par[2]; tm <- round(1/N)
    E = qfbFit$par[3]
    Qr = qfbFit$par[4]
    rss = qfbFit$value
  }
  #----------------------------------------------------------------
  
  # output df
  RP$alpha[1] <- alpha
  RP$E[1] <- E
  RP$N[1] <- N; 
  RP$tm[1] <- tm
  RP$KGE[1] <- kge
  RP$Qr[1] <- Qr
  RP$N_len[1] <- NR-1
  
  # calculate dynamic volume discarge
  Vdyn <- 86400*(Qr/alpha)
  RP$Vdyn[1] <- Vdyn
  
  # plot measured and predicted baseflow recession
  if(isTRUE(plotSim)){
    dev.new()
    par(xaxs='i', yaxs='i')
    xlim <- c(min(RP$t_index), max(RP$t_index))
    ylim <- c(min(RP$Q), max(RP$Q))
    ylim <- c(0, max(RP$Q))
    plot(Q ~ t_index, data = RP, xlim = xlim, ylim = ylim, pch = 18)
    
    par(new=T, xaxs='i', yaxs='i')
    ylim <- c(min(RP$simQ), max(RP$simQ))
    ylim <- c(0, max(RP$Q))
    plot(simQ ~ t_index, data = RP, type = "l", lwd = 3, xlim = xlim, ylim = ylim, col = "black", axes = F, xlab = "", ylab = "")
    lines(RP$t_index, simQf, xlim = xlim, ylim = ylim, col = "green", axes = F, xlab = "", ylab = "")
    lines(RP$t_index, simQb, xlim = xlim, ylim = ylim, col = "blue", axes = F, xlab = "", ylab = "")
    legend("topright", bty = "n", pch = c(18,NA,NA,NA), lty = c(NA,1,1,1), col = c("black","black","green","blue"), 
           legend = c("Observed discharge","Simulated discharge","Conduit drainage","Matrix drainage"))
    abline(v = ti, lty = "solid", lwd = 2, col ="red")
    abline(v = tm, lty = "dashed", lwd = 2, col = "red")
    
  }
  
  # save plot if plotName argument is provided, plotName should include dir and output file format
  # plotName must include path and file extension of the image to be saved
  if(!is.null(plotName)){
    pdf(paste0("D:/Project_2/Results/Fitted recession/",plotName), height = 5, width = 7)
    par(xaxs='i', yaxs='i')
    xlim <- c(min(RP$t_index), max(RP$t_index))
    ylim <- c(min(RP$Q), max(RP$Q))
    ylim <- c(0, max(RP$Q))
    plot(Q ~ t_index, data = RP, xlim = xlim, ylim = ylim, pch = 18, xlab = "t[days]", ylab = "Q[m3/s]", main = plotName)
    
    par(new=T, xaxs='i', yaxs='i')
    ylim <- c(min(RP$simQ), max(RP$simQ))
    ylim <- c(0, max(RP$Q))
    plot(simQ ~ t_index, data = RP, type = "l", lwd = 3, xlim = xlim, ylim = ylim, col = "black", axes = F, xlab = "", ylab = "")
    lines(RP$t_index, simQf, xlim = xlim, ylim = ylim, col = "green", axes = F, xlab = "", ylab = "")
    lines(RP$t_index, simQb, xlim = xlim, ylim = ylim, col = "blue", axes = F, xlab = "", ylab = "")
    legend("topright", bty = "n", pch = c(18,NA,NA,NA), lty = c(NA,1,1,1), col = c("black","black","green","blue"), 
           legend = c("Observed discharge","Simulated discharge","Conduit drainage","Matrix drainage"))
    abline(v = ti, lty = "solid", lwd = 2, col ="red")
    abline(v = tm, lty = "dashed", lwd = 2, col = "red")
    dev.off()
  }
  
  return(RP)
}
