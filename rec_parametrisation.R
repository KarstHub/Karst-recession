#### ----------- PARAMETER FITTING FOR RECESSION ANALYSIS METHODS -----------------###
# calculate linear matrix recession coefficients from martrix recession segments of the hydrograph
# Input: list of recession segment in daily time setps (output of rec.extract function)
# Bsp x [[1]]
# Date      , Q,        Qc ,            Qm ,        Vt ,     Index,  Qro,       ti
# 2017-03-16 0.8387450 0.8387450        NA      20652133     1      0.3741067   74
# ...x[[2]].
# 3 different regressions methods can be choosen for parametrisation of the linear recession model
# 1) the linear approch fit a linear regression line by ordinary least squares to the recession plot (-dQdt gegen Q) 
# 2) the secound method uses a lower enverlope regression to account for unavoidable errors
# 3) a linear approch using the root mean square error for optimisation [work in progress]
# Output:  the intersect of the regresssion line at x = 0 gives the recession rate.
#          the recession rate (rr) acn be transformt to the recession constant or coefficient (rc) by  rc = 1/ exp(rr)
#          resulting in the linear recession coefficient in days 
# 
# library(quantreg)
# library(pspline)
# 
x <- singel_model_run(3932, t.plot= T) %>%
    dplyr::filter(Date >= start_date + years(warmup)) %>%
    rename(Q_out_model = Q) %>%
    rename(Q_Conduit_model= Q_C) %>%
    rename(Q_Matrix_model = Q_M)


x <- rec.extract(x,  model = "Vogel", plot = T)
# #ganez recession
# yy <- rbind(x[[60]][2]$Q)
# xx <- rbind(x[[60]][1]$Date)
# plot(yy ~ xx)
#matrix recession
df <- cbind(x[[21]][4], x[[21]][1])
df$dQdt <- c(NA, diff(df$Qm, 1))
df <- df %>% filter( dQdt != is.na(Qm))
plot(df$Qm ~df$Date)
plot(log10(df$Qm) ~log10(df$Date))
# 
 plot(-df$dQdt ~ df$Qm)
# #recession plot
 plot(log(-df$dQdt) ~ log(df$Qm))
 lm1 <- lm(log(-df$dQdt) ~ log(df$Qm))
 abline(lm1, col= "red")
rc <- 1/exp(lm1$coefficients[1])

rec_coeff_para <- function(x, plot=TRUE, method=c("linear", "lower_env", "RMS"))
    {    
    x <- rbindlist(x) %>%   
        rename(Q_rec = Q) %>%
        rename(Q_Matrix_rec = Qm) %>%
        rename(Q_Conduit_rec = Qc) %>%
        select(-c(Vt, Index, Qro, ti)) %>% 
        mutate("dQdt_matrix" = lag(Q_Matrix_rec, n = 1) - Q_Matrix_rec ) %>% # dQ/dt: Ableitung der matrix recession über zwei 1 Zeitschritt
        mutate("dQdt_matrix2" = lag(Q_Matrix_rec, n = 1) - lead(Q_Matrix_rec, n=1)/2) # dQ/dt: Ableitung der matrix recession über zwei 2 zeitschritte
    
    if (method=="lower_env"){
        rec_mod <- rq(log(-dQdt_matrix) ~log(Q_Matrix_rec), data = x, tau = 0.05)
        rec_coeff <- 1/exp(rec_mod$coefficients[1]) 
        
        if (plot== TRUE){
            #plot(-x$dQdt_matrix ~ x$Q_Matrix_rec)
            
            par(mfrow=c(1,1), xpd=F, mar=c(4, 4, 2, 2), cex=1.5, mgp=c(3,1,0))
            plot(
                log(-x$dQdt_matrix) ~ log(x$Q_Matrix_rec),
                xlab = "log(Q)",
                ylab = "log(-dQ/dt )",
                ylim = range(log(-x$dQdt_matrix),na.rm=T),
                xlim = range(log(x$Q_Matrix_rec),na.rm=T),
                main= "")
            abline(rec_mod, col="green", lwd=2)
            legend("bottomright", legend=paste(method, "regression line"),
                   bty = "n", col="green", lwd=2, xpd=T)

        }
    }
    
    if (method=="linear"){
        rec_mod <- lm(log(-dQdt_matrix) ~log(Q_Matrix_rec), data = x)
        rec_coeff <- 1/exp(rec_mod$coefficients[1]) 
        
        if (plot== TRUE){
            #plot(-x$dQdt_matrix ~ x$Q_Matrix_rec)
            
            
            par(mfrow=c(1,1), xpd=F, mar=c(4, 4, 2, 2), cex=1.5, mgp=c(3,1,0))
            plot(
                log(-x$dQdt_matrix) ~ log(x$Q_Matrix_rec),
                xlab = "log(Q)",
                ylab = "log(-dQ/dt )",
                ylim = range(log(-x$dQdt_matrix),na.rm=T),
                xlim = range(log(x$Q_Matrix_rec),na.rm=T),
                main= "")
            abline(rec_mod, col="red", lwd=2)
            legend("bottomright", legend=paste(method, "regression line"),
                   bty = "n", col="red", lwd=2, xpd=T)
            
        }
    }
    if (method=="RMS"){
        rec_mod <- lm(log(-dQdt_matrix) ~log(Q_Matrix_rec), data = x) ## whiche mdoel to use?
        rec_coeff <- 1/exp(rec_mod$coefficients[1]) 
        
        
        if (plot== TRUE){
            #plot(log(x$dQdt_matrix) ~ x$Date)
            
            
            par(mfrow=c(1,1), xpd=F, mar=c(4, 4, 2, 2), cex=1.5, mgp=c(3,1,0))
            plot(
                log(-x$dQdt_matrix) ~ log(x$Q_Matrix_rec),
                xlab = "log(Q)",
                ylab = "log(-dQ/dt )",
                ylim = range(log(-x$dQdt_matrix),na.rm=T),
                xlim = range(log(x$Q_Matrix_rec),na.rm=T),
                main= "Is root-mean square error the obejtive function of the linear model?")
            abline(rec_mod, col="blue", lwd=2)
            legend("bottomright", legend=paste(method, "regression line"),
                   bty = "n", col="blue", lwd=2, xpd=T)
    }
    
    
    
    
    
    }
    return(rec_coeff)
    
}
# rec_coff_para(y, method="linear", plot=T)
# rec_coff_para(y, method="quant5", plot=T)
