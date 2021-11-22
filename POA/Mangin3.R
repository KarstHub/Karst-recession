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
mangin3 <- function(RP, Q0 = NULL, Qro = NULL, q0 = NULL, plotSim = TRUE, plotName = NULL){
  #Q0 = NULL; Qro = NULL; q0 = NULL
  #RP <- XX
  # RP is the extracted recession points in data frame ["Date","Q","Qc","Qm","Qro","ti"]
  # Q0 corresponds to peak discharge at time t = 0
  # Qro corresponds to baseflow component of discharge at inflexion point (change from fast to slow flow)
  # q0 is the diffrence between peak discharge Q0 and Qro
  
  # define output variables
  coeff <- data.frame("alpha"="", "beta"="")
  
  if(is.null(Q0)){
    Q0 = RP[1,2]
  }
  
  if(is.null(Qro)){
    Qr = RP$Qro[1]
  }
  
  if(Qr > Q0){
    Qr = Q0
  }
  
  q0 = Q0 - Qr
  RP$t_index <- seq(0, nrow(RP)-1, 1)
  ti <- which(!is.na(RP$Qm))[1] #end of conduit and start of matrix (baseflow) drainge
  NR <- nrow(RP)
  
  setMin2 <- function(expr){
    pmax(0, eval(expr))
  }
  
  # optimise recession model
  Nup <- 1/(NR-2)
  qfb <- function(data, par){
    Qr = Qr; q0 = q0
    alpha <- par[1]
    N <- par[2]
    E <- par[3]
    Qmat <- with(data, Qr * exp(-(alpha)*t_index))
    Qcon <- with(data, setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
    rss <- with(data, sum((Q - (Qmat + Qcon))^2))
    return(rss)
  }
  qfbFit = optim(par = c(0.1,0.5,0.1), fn = qfb, data = RP, method = "L-BFGS-B", 
                    lower = c(0,Nup,0), upper = c(Inf,1,10), control = c(trace = T, maxit = 30000))
  alpha = qfbFit$par[1]
  N = qfbFit$par[2]; tm <- round(1/N)
  E = qfbFit$par[3]
  rss = qfbFit$value
  
  # simulated Q
  RP$simQ = with(data = RP, Qr * exp(-(alpha)*t_index) + setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
  #RP$simQ[tm:NR] <- with(data = RP[tm:NR, ], Qr * exp(-(alpha)*t_index))
  #RP$simQ[1:tm] <- with(data = RP[1:tm, ], Qr * exp(-(alpha)*t_index) + (q0 * (1 - (N*t_index))/(1 + (E*t_index))))
  kge <- KGE(RP$Q,RP$simQ)
  
  # output df
  RP$alpha[1] <- alpha
  RP$E[1] <- E
  RP$N[1] <- N; #ti <- round(1/N)
  RP$tm[1] <- tm
  RP$KGE[1] <- kge
  RP$Qr[1] <- Qr
  
  # plot measured and predicted baseflow recession
  if(isTRUE(plotSim)){
    dev.new()
    xlim <- c(min(RP$t_index), max(RP$t_index))
    ylim <- c(0, max(RP$Q))
    plot(Q ~ t_index, data = RP, xlim = xlim, ylim = ylim, pch = 18)
    
    par(new=T)
    plot(simQ ~ t_index, data = RP, type = "l", xlim = xlim, ylim = ylim, col = "red", axes = F, ylab = "")
    abline(v = ti, lty = "solid", lwd = 2)
    abline(v = tm, lty = "dashed", lwd = 2, col = "green")
    
  }
  
  # save plot if plotName argument is provided, plotName should include dir and output file format
  # plotName must include path and file extension of the image to be saved
  xlim <- c(min(RP$t_index), max(RP$t_index))
  ylim <- c(min(RP$Q), max(RP$Q))
  
  if(!is.null(plotName)){
    png(plotName)
    plot(Q ~ t_index, data = RP, xlim = xlim, ylim = ylim, pch = 18)
    par(new=T)
    plot(simQ ~ t_index, data = RP, type = "l", xlim = xlim, ylim = ylim, col = "red", axes = F, ylab = "")
    abline(v = ti, lty = "dashed", lwd = 2)
    dev.off()
  }
  
  return(RP)
}
