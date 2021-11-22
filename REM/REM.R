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
source("./Recession analysis/Recession model/REM/remVogel.R")
source("./Recession analysis/Recession model/REM/remBrut.R")
source("./Recession analysis/Recession model/REM/remAkw.R")
source("./Recession analysis/Recession model/REM/remPlot.R")
##==============================================
##    Hydrograph recession extraction model   ==
##==============================================
rem.extract <- function(x, model = c("Vogel", "Brut", "Akw"), MA = NULL, dQ = NULL, Qf = NULL,
                        CoV = 0.10, lambda = 0.3, fit = 0.85, len = 5, plot = FALSE, save.plot=NULL)
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
  #           if model == Akw, maximum coefficient of variation (CV) allowable for spring discharge [0:1]
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
  
  model <- match.arg(model)
##--------------------------------------------------------------------------------------------------------------------------------------------
  # set default variables
  if(model == "Vogel"){ # select vogel extraction method
    if(!is.null(MA)){
      movg = MA
    }else{
	  movg = 3
	}
   
    if(!is.null(dQ)){
      Qdiff = dQ/100
    }else{
	  Qdiff = 0.3
	}
    recession <- remVogel(Qts=x, MA=movg, dQ=Qdiff, lambda=lambda, len=len, plot=plot, save.plot=save.plot)
    return(recession)
  }
  
  if(model == "Brut"){ # select Brutsaert extraction method
      
      if(!is.null(Qf)){
        Qmajor = quantile(x[,2], probs = (100-Qf)/100, na.rm = T)
      }else{
	    Qmajor = quantile(x[,2], probs = 0.7, na.rm = T) # upper 30% percentile for major events
	  }
      
      if(!is.null(dQ)){
        Qdiff = dQ/100
      }else{
	    Qdiff <- 0.3
	  }
      recession <- remBrut(Qts=x, Qf=Qmajor, dQ=Qdiff, len=len, plot=plot, save.plot=save.plot)
      return(recession)
  }
    
  if(model == "Akw"){ # select Aksoy extraction method
      v = 0.1  # maximum coeeficient of variation 10% (default)
      if(!is.null(CoV)){
        v = CoV
      }
      recession <- remAkw(Qts=x, CoV=v, len=len, plot=plot, save.plot=save.plot)
      return(recession)
  }

}
