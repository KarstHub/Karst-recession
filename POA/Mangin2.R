##---------------------------------------------------------------------------------------------------
# Mangin model for karst spring hydrograph recession analysis
# Model use two different types of functions: ψt representing contribution from unsaturated zone
# and ϕt representing contribution from the saturated zone (Mangin 1975). Spring dsicharge can be
# modelled by: Qt = ϕt + ψt; ϕt is the baseflow component modelled with Maillet exponential model
# and ψt is the flood recession component modelled with a non-exponential model.
##---------------------------------------------------------------------------------------------------
# The function fit Mangin model to extracted recession hydrographs to calibraterecession parameters. 
# This is a two-step optimisation procedure and one parameter is optimised in first step of matrix model 
# optimisation while two parameters are optimised in conduit model. Note that N which is 1/ti is not fixed here
# Define function, parameters and variables
source("./Recession analysis/Recession model/KGE.R")
#source("./Recession analysis/Recession model/setMin.R")
mangin2 <- function(RP, Q0 = NULL, Qro = NULL, q0 = NULL, plotSim = TRUE, plotName = NULL, orig = 0){
  #Q0 = NULL; Qro = NULL; q0 = NULL
  #RP <- res2[[18]]
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
  
  # define and optimise the matrix drainage model
  # Baseflow (matrix) recession, contribution from conduit assummed to be zero; ψt = 0
  qb_df <- RP[ti:NR, c("Q","t_index")]
  qb <- function(data, par){
    Qr = Qr
    alpha <- par
    Qmat <- with(data, Qr * exp(-(alpha)*t_index))
    with(data, sum((Q - Qmat)^2))
  }
  qbFit <- optim(par = 0.1, fn = qb, data = qb_df, method = "L-BFGS-B", 
                 lower = 0, upper = Inf, control = c(trace = T, maxit = 20000))
  alpha = qbFit$par
  
  # Quick-flow (conduit) recession
  #qFlow <- RP[1:ti, c("Q","t_index")]
  Nup <- 1/(NR-2)
  qf <- function(data, par){
    Qr = Qr; alpha = alpha; q0 = q0
    N <- par[1]
    E <- par[2]
    Qmat <- with(data, Qr * exp(-(alpha)*t_index))
    Qcon <- with(data, setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
    rss <- with(data, sum((Q - (Qmat + Qcon))^2))
    return(rss)
  }
  qfFit <- optim(par = c(0.5,0.1), fn = qf, data = RP, method = "L-BFGS-B", 
                 lower = c(Nup,0), upper = c(1,10), control = c(trace = T, maxit = 20000))      
  N = qfFit$par[1]; tm <- round(1/N)
  E = qfFit$par[2] 
  
  # simulated Q
  RP$simQ[tm:NR] <- with(data = RP[tm:NR, ], Qr * exp(-(alpha)*t_index))
  RP$simQ[1:tm] <- with(data = RP[1:tm, ], Qr * exp(-(alpha)*t_index) + (q0 * (1 - (N*t_index))/(1 + (E*t_index))))
  simQb <- with(data = RP, Qr * exp(-(alpha)*t_index))
  simQf <- with(data = RP, setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
  kge <- KGE(RP$Q,RP$simQ)
  
  #---------------------------------------------------------
  # if non-modified extracted REM points are considered.
  if(orig == 1){
    
    qb_df <- RP[ti:NR, c("Q","t_index")]
    qb_df[,2] <- qb_df$t_index - qb_df$t_index[1]
    Qr <- qb_df$Q[1]
    q0 <- Qr - tail(qb_df$Q,1)
    
    qb <- function(data, par){
      Qr = Qr
      alpha <- par
      rss <- with(data, sum((Q - (Qr * exp(-(alpha)*t_index)))^2))
      return(rss)
    }
    qbFit <- optim(par = 0.01, fn = qb, data = qb_df, method = "L-BFGS-B", 
                   lower = 0, upper = 0.5, control = c(trace = T, maxit = 20000))
    alpha = qbFit$par
    
    # Quick-flow (conduit) recession
    #qFlow <- RP[1:ti, c("Q","t_index")]
    Nup <- 1/nrow(qb_df)
    qf <- function(data, par){
      Qr = Qr; alpha = alpha; q0 = q0
      N <- par[1]
      E <- par[2]
      Qmat <- with(data, Qr * exp(-(alpha)*t_index))
      Qcon <- with(data, setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
      rss <- with(data, sum((Q - (Qmat + Qcon))^2))
      return(rss)
    }
    qfFit <- optim(par = c(0.5,0.1), fn = qf, data = qb_df, method = "L-BFGS-B", 
                   lower = c(0,0), upper = c(1,10), control = c(trace = T, maxit = 20000))      
    N = qfFit$par[1]; tm <- round(1/N)
    E = qfFit$par[2] 
  }
  #----------------------------------------------------------------
  
  # output df
  RP$alpha[1] <- alpha
  RP$E[1] <- E
  RP$N[1] <- N; tm <- round(1/N)
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
    plot(Q ~ t_index, data = RP, xlim = xlim, ylim = ylim, pch = 18, xlab = "t[days]", ylab = "Q[m3/s]")
    
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
    pdf(paste0("D:/Project_2/Results/Fitted recession/",plotName), height = 5, width = 8)
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
