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
mangin1 <- function(RP, Q0 = NULL, Qro = NULL, q0 = NULL, plotSim = TRUE, plotName = NULL, orig = 0){
  
  #EXT <- match.arg(EXT)
  Q0 = NULL; Qro = NULL; q0 = NULL; plotSim = TRUE; plotName = NULL;EXT = "original"
  RP <- summerEv[["Qachqouch@Akw"]][[1]]
  # RP is the extracted recession points in data frame ["Date","Q","Qavg","dQdt","Qc","Qm","ti","index","Qro","vt"]
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
  
  # slow-flow (matrix) recession
  qb_df <- RP[which(!is.na(RP$Qm)), ]
  
  # Quick-flow (conduit) recession
  qf_df <- RP[which(!is.na(RP$Qc)), ]
  
  if(Qr > Q0){
    Qr = qf_df$Qm[nrow(qf_df)]
  }
  q0 = Q0 - Qr
  
  # RP$t_index <- seq(0, nrow(RP)-1, 1)
  ti <- RP$ti[1] - 1 #end of conduit and start of matrix (baseflow) drainge
  tm <- ti
  NR <- nrow(RP)
  
  setMin2 <- function(expr){
    pmax(0, eval(expr))}
  
  # its recommended to fit the linear reservior model first for the Mangin approach
  # Baseflow (matrix) recession, contribution from conduit assummed to be zero; ψt = 0
  
  # check and verify approach ---
  # qb_df <- qb_df[which(qb_df$dQdt<0), ]
  # a <- summary(lm(log10(-dQdt) ~ log10(Qavg), data = qb_df))$coefficients[1]
  # a <- exp(-a)
  # plot(qb_df$Qavg, -qb_df$dQdt, log="xy")
  #-----------------------------------------
 
  # estimate slow-flow parameter, alpha
  qb <- function(data, par){
    Qr = Qr
    alpha <- par
    rss <- with(data, sum((Q - (Qr * exp(-(alpha)*index)))^2))
    return(rss)
  }
  qbFit <- optim(par = 0.0001, fn = qb, data = qb_df, method = "L-BFGS-B", 
                 lower = 0, upper = 0.5, control = c(trace = T, maxit = 20000))
  alpha = qbFit$par
  
  # estimate quick-flow recession parameter, E 
  N <- 1/ti 
  qf <- function(data, par){
    Qr = Qr; alpha = alpha; q0 = q0; N = N
    E <- par
    rss <- with(data, sum((Q - (Qr * exp(-(alpha)*index) + (q0 * (1 - (N*index))/(1 + (E*index)))))^2))
    return(rss)
  }
  qfFit <- optim(par = 0.01, fn = qf, data = qf_df, method = "L-BFGS-B", lower = 0, upper = 10,
                 control = c(trace = T, maxit = 200000))      
  E = qfFit$par
 
  # simulated Q
  RP$simQ <- NA
  simQb <- with(data = RP, Qr * exp(-(alpha)*index))
  simQf <- with(data = RP, setMin2(q0 * ((1 - (N*index))/(1 + (E*index)))))
  
  RP$simQ[(ti+1):NR] <- simQb[(ti+1):NR]
  RP$simQ[1:ti] <- with(data = qf_df, Qr * exp(-(alpha)*index) + (q0 * (1 - (N*index))/(1 + (E*index))))
  
  kge <- KGE(RP$Q,RP$simQ)
  
  #---------------------------------------------------------
  # if non-modified extracted REM points are considered.
  # if(orig == 1){
  #   qb_df <- RP[ti:NR, c("Q","t_index")]
  #   qb_df[,2] <- qb_df$t_index - qb_df$t_index[1]
  #   Qr <- qb_df$Q[1]
  #   q0 <- Qr - tail(qb_df$Q,1)
  #   
  #   qb <- function(data, par){
  #     Qr = Qr
  #     alpha <- par
  #     rss <- with(data, sum((Q - (Qr * exp(-(alpha)*t_index)))^2))
  #     return(rss)
  #   }
  #   qbFit <- optim(par = 0.01, fn = qb, data = qb_df, method = "L-BFGS-B", 
  #                  lower = 0, upper = 0.5, control = c(trace = T, maxit = 20000))
  #   alpha = qbFit$par
  #   
  #   # fit baseflow segment and quickflow
  #   N <- 1/tail(qb_df$t_index,1) 
  #   qf <- function(data, par){
  #     Qr = Qr; alpha = alpha; q0 = q0; N = N
  #     E <- par
  #     rss <- with(data, sum((Q - (Qr * exp(-(alpha)*t_index) + (q0 * (1 - (N*t_index))/(1 + (E*t_index)))))^2))
  #     return(rss)
  #   }
  #   qfFit <- optim(par = 0.1, fn = qf, data = qb_df, method = "L-BFGS-B", lower = 0, upper = 5,
  #                  control = c(trace = T, maxit = 20000))      
  #   E = qfFit$par
  # }
  #----------------------------------------------------------------
  
  # output df
  RP$alpha[1] <- alpha
  RP$E[1] <- E
  RP$N[1] <- N
  RP$Qr[1] <- Qr
  RP$N_len[1] <- NR-1
  RP$ti[1] <- ti
  RP$KGE[1] <- kge
  
  # calculate dynamic volume discarge
  Vdyn <- 86400*(Qr/alpha)
  RP$Vdyn[1] <- Vdyn
  
  # plot measured and predicted baseflow recession
  if(isTRUE(plotSim)){
    dev.new()
    par(xaxs='i', yaxs='i')
    xlim <- c(min(RP$index), max(RP$index))
    ylim <- c(min(RP$Q), max(RP$Q))
    ylim <- c(0, max(RP$Q))
    plot(Q ~ index, data = RP, xlim = xlim, ylim = ylim, pch = 18, xlab = "t[days]", ylab = "Q[m3/s]", main = plotName)
    
    par(new=T, xaxs='i', yaxs='i')
    ylim <- c(min(RP$simQ), max(RP$simQ))
    ylim <- c(0, max(RP$Q))
    plot(simQ ~ index, data = RP, type = "l", lwd = 2, xlim = xlim, ylim = ylim, col = "black", axes = F, xlab = "", ylab = "")
    lines(RP$index, simQf, xlim = xlim, ylim = ylim, col = "green", xlab = "", ylab = "")
    lines(RP$index, simQb, xlim = xlim, ylim = ylim, col = "blue", xlab = "", ylab = "")
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
      xlim <- c(min(RP$index), max(RP$index))
      ylim <- c(min(RP$Q), max(RP$Q))
      ylim <- c(0, max(RP$Q))
      plot(Q ~ index, data = RP, xlim = xlim, ylim = ylim, pch = 18, xlab = "t[days]", ylab = "Q[m3/s]", main = plotName)
      
      par(new=T, xaxs='i', yaxs='i')
      ylim <- c(min(RP$simQ), max(RP$simQ))
      ylim <- c(0, max(RP$Q))
      plot(simQ ~ index, data = RP, type = "l", lwd = 3, xlim = xlim, ylim = ylim, col = "black", axes = F, xlab = "", ylab = "")
      lines(RP$index, simQf, xlim = xlim, ylim = ylim, col = "green", axes = F, xlab = "", ylab = "")
      lines(RP$index, simQb, xlim = xlim, ylim = ylim, col = "blue", axes = F, xlab = "", ylab = "")
      legend("topright", bty = "n", pch = c(18,NA,NA,NA), lty = c(NA,1,1,1), col = c("black","black","green","blue"), 
             legend = c("Observed discharge","Simulated discharge","Conduit drainage","Matrix drainage"))
      abline(v = ti, lty = "solid", lwd = 2, col ="red")
      abline(v = tm, lty = "dashed", lwd = 2, col = "red")
      dev.off()
  }
  
  return(RP)
}
