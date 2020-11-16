Sys.setenv(language="en")
require(hydroGOF)
##------------------------------------------------------------------------------------------------------------------------------------------
# Recession extraction and analysis for karst spring hydrograph
# Extraction of karst spring hydrograph recessions with tradtional extraction methods
# The traditional extraction methods are developed to extract baseflow, in this case they extract matrix recessions
# By kicking out the restrictions in traditonal extraction methods, conduit component can also be extracted
# Recession extraction methods: Vogel and Kroll (1992)
#                               Brutsaert and Nieber (1997)
#                               Aksoy and Wittenberg (2011)

# Mangin and Malliet karst spring recession model for conduit and matrix recession coefficients
# Mangin recession model: insert citation here
# Malliet recession model: insert citation here
source("./Recession analysis/Recession model/recVogel.R")
source("./Recession analysis/Recession model/recBrut.R")
source("./Recession analysis/Recession model/recAkw.R")
source("./Recession analysis/Recession model/recClassical.R")
source("./Recession analysis/Recession model/recPlot.R")
##==============================================
##    Hydrograph recession extraction model   ==
##==============================================
rec.extract <- function(x, model = c("Vogel", "Brut", "Akw", "Classical"), par1 = NULL, par2 = NULL, 
                        CoV = 0.10, lambda = 0.3, fit = 0.85, len = 5, plot = FALSE)
  { 
  # include par1 and par2 instead par1 and PAR
  # extract recession segments of hydrograph with different methods
  
  # INPUT PARAMETERS
  #   x --> A cell array file in the format [date, Q] for each spring hydrograph
  #   model == 1 Vogel and Kroll
  #   model == 2 Brutsaert and Neiber
  #   model == 3 Classical (combine vogel and linear curve fitting)
  #   model == 4 Aksoy and Wittenberg
  #   par1 --> if model == Vogel, is the moving average days for smoothing spring hydrograph
  #            if model == Brut, is define as the streamflow percentile for major events [0:100]
  #            if model == Classical, moving average days as in model 1
  #   Par2 --> if model == Vogel, is the maximum percentage difference between consecutive discharge value [0:100]
  #           if model == Brut, maximum percentage difference (-dQ/dt) for unspurious baseflow condition [0:100]
  #           if model == 4, maximum coefficient of variation (CV) allowable for spring discharge [0:100]
  #   fit = for Classical model, is the R squared value [0:100] of linear model fitted to extracted recession points
  #	  lambda = Fraction of Q influenced by spurious flow to be removed, 0.3 default and applicable only to model 1
  
  # OUTPUT
  #   cumRP --> dataframes [Date, Q, Qc, Qm] array of extracted recession segments of hydrograph for model 1,2 and 4
  #             dataframes [Date, Q, Qc, Qm, Qro, ti, Index] for model 3
  #             Qc = conduit and Qm = matrix recession, Qro = intercept of linear model, ti = start time of conduit drainage
  # define input data columns name and format date
  
  # x=Q; model=1; par1=NULL; par2=NULL; len=5; plot=T
  
  Date <- format(x[,1], format = "%Y-%m-%d")
  Q <- x[, 2]
  
  x <- data.frame(Date, Q) # create new dataframe
  model <- match.arg(model)
##--------------------------------------------------------------------------------------------------------------------------------------------
  # set default variables
  if(model == "Vogel"){ # select vogel extraction method
    MA = 3  # default moving average days for smoothing hydrograph
    if(!is.null(par1)){
      MA = par1
    }
    dQ = 0.3 # percentage difference (30%) between consecutive spring discharge (default)
    if(!is.null(par2)){
      dQ = par2/100
    }
    recession <- recVogel(x=x, MA=MA, dQ=dQ, lambda=lambda, len=len, plot=plot)
    return(recession)
  }
  
  if(model == "Brut"){ # select Brutsaert extraction method
      Qf = quantile(x$Q, probs = 0.3, na.rm = T) # maximum percentile for spring discharge
      if(!is.null(par1)){
        Qf = quantile(x$Q, probs = par1/100, na.rm = T)
      }
      dQ = 0.3
      if(!is.null(par2)){
        dQ = par2/100
      }
      recession <- recBrut(DATA=x, Qf=Qf, dQ=dQ, len=len, plot=plot)
      return(recession)
  }
    
  if(model == "Classical"){ # select new model
      MA = 3  # moving average days
      if(!is.null(par1)){
        MA = par1
      }
      recession <- recClassical(x=x, MA=MA, fit=fit, len=len, plot=plot)
      return(recession)
  }
    
  if(model == "Akw"){ # select Aksoy extraction method
      v = 0.1  # maximum coeeficient of variation 20% (default)
      if(!is.null(CoV)){
        v = CoV
      }
      recession <- recAkw(DATA=x, CoV=v, len=len, plot=plot)
      return(recession)
  }

}
