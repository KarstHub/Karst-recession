##---------------------------------------------------------------------------------------------------
# Mangin model for karst spring hydrograph recession analysis
# Model use two different types of functions: ψt representing contribution from unsaturated zone
# and ϕt representing contribution from the saturated zone (Mangin 1975). Spring dsicharge can be
# modelled by: Qt = ϕt + ψt; ϕt is the baseflow component modelled with Maillet exponential model
# and ψt is the flood recession component modelled with a non-exponential model.
##---------------------------------------------------------------------------------------------------
# The function fit Mangin model to extracted recession hydrographs to calibrate
# recession parameters. This is a two-step optimisation procedure and one parameter
# is optimised at each step of matrix and conduit model optimisation
# Define function, parameters and variables
source("./Recession analysis/Recession model/KGE.R")
#source("./Recession analysis/Recession model/setMin.R")
mangin1 <- function(RP, Q0 = NULL, Qro = NULL, q0 = NULL, plotSim = TRUE, plotName = NULL){
  #Q0 = NULL; Qro = NULL; q0 = NULL; plotSim = TRUE
  #RP <- hyd
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
  tm <- ti
  NR <- nrow(RP)
  
  setMin2 <- function(expr){
    pmax(0, eval(expr))}
  
  # its recommended to fit the linear reservior model first for the Mangin approach
  # Baseflow (matrix) recession, contribution from conduit assummed to be zero; ψt = 0
  qb_df <- RP[ti:NR, c("Q","t_index")]
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
  qf_df <- RP[1:ti, c("Q","t_index")]
  N <- 1/ti 
  qf <- function(data, par){
    Qr = Qr; alpha = alpha; q0 = q0; N = N
    E <- par
    rss <- with(data, sum((Q - (Qr * exp(-(alpha)*t_index) + (q0 * (1 - (N*t_index))/(1 + (E*t_index)))))^2))
    return(rss)
  }
  qfFit <- optim(par = 0.1, fn = qf, data = qf_df, method = "L-BFGS-B", lower = 0, upper = 5,
                    control = c(trace = T, maxit = 20000))      
  E = qfFit$par
  
  # simulated Q
  RP$simQ[ti:NR] <- with(data = qb_df, Qr * exp(-(alpha)*t_index))
  simQb <- with(data = RP, Qr * exp(-(alpha)*t_index))
  simQf <- with(data = RP, setMin2(q0 * ((1 - (N*t_index))/(1 + (E*t_index)))))
  RP$simQ[1:ti] <- with(data = qf_df, Qr * exp(-(alpha)*t_index) + (q0 * (1 - (N*t_index))/(1 + (E*t_index))))
  kge <- KGE(RP$Q,RP$simQ)
  
  # output df
  RP$alpha[1] <- alpha
  RP$E[1] <- E
  RP$N[1] <- N; tm <- round(1/N)
  RP$tm[1] <- tm
  RP$KGE[1] <- kge
  RP$Qr[1] <- Qr
  
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
    plot(Q ~ t_index, data = RP, xlim = xlim, ylim = ylim, pch = 18, xlab = "t[days]", ylab = "Q[m3/s]", main = plotName)
    
    par(new=T, xaxs='i', yaxs='i')
    ylim <- c(min(RP$simQ), max(RP$simQ))
    ylim <- c(0, max(RP$Q))
    plot(simQ ~ t_index, data = RP, type = "l", lwd = 2, xlim = xlim, ylim = ylim, col = "black", axes = F, xlab = "", ylab = "")
    lines(RP$t_index, simQf, xlim = xlim, ylim = ylim, col = "green", axes = F, xlab = "", ylab = "")
    lines(RP$t_index, simQb, xlim = xlim, ylim = ylim, col = "blue", axes = F, xlab = "", ylab = "")
    legend("topright", bty = "n", pch = c(18,NA,NA,NA), lty = c(NA,1,1,1), col = c("black","black","green","blue"), 
           legend = c("Observed discharge","Simulated discharge","Conduit drainage","Matrix drainage"))
    abline(v = ti, lty = "solid", lwd = 2, col ="red")
    abline(v = tm, lty = "dashed", lwd = 2, col = "red")
    #dev.off()
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
