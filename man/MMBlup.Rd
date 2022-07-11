\name{MMBlup}
\alias{MMBlup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  BLUP from MM4LMM results
  %%  ~~function to do ... ~~
}
\description{
  This function computes the BLUP for each random vector included in the \code{\link{MMEst}} output. Note that this function can be called only if the argument \code{X} of \code{\link{MMEst}} was set to \code{NULL}
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
    %%, as in the default vcov() function for lme4,
  %%, as in predict() function for ASReml
}
\usage{
  MMBlup(Y,Cofactor = NULL, X = NULL, fmla = NULL,ZList=NULL,VarList,ResMM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
    The vector of response values used in the function \code{\link{MMEst}}.
  }
  \item{Cofactor}{
    The incidence matrix corresponding to fixed effects common to all models to be adjusted used in the function \code{\link{MMEst}}. If \code{NULL}, a vector full of 1 is used.
  }
  \item{X}{
   Must be \code{NULL}. 
    %%     ~~Describe \code{X} here~~
  }
  \item{fmla}{
    The formula object specifying the fixed effect part of all models separated by + operators used in the function \code{\link{MMEst}} (default is \code{NULL}).
  }
  \item{ZList}{
    The list of incidence matrices associated with random and residual effects used in the function \code{\link{MMEst}} (default is \code{NULL}).
    
  }
  \item{VarList}{
    The list of covariance matrices associated with random and residual effects used in the function \code{\link{MMEst}}.
  }

  \item{ResMM}{
    A list as displayed by the \code{\link{MMEst}} function for a Variance Component Analysis.
  }
  
}

\value{
  The function returns a list where each element corresponds to the Best Linear Unbiased Prediction of a random component of the model.
  
  %%	The output of the function is a matrix. This matrix is the variance covariance matrix of the variance estimates of the model fitted in \code{ResMM}. To extract the standard errors of variance estimates, one has to extract the diagonal of this matrix.
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
\author{GQMS CoreFunctions Team
  %%  ~~who you are~~
}
%%\note{
  %%  ~~further notes~~
    %%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  %%\seealso{
    %% ~~objects to See Also as \code{\link{help}}, ~~~
      %%}
\examples{
  require('MM4LMM')
  data(VarianceComponentExample)
  DataHybrid <- VarianceComponentExample$Data
  KinF <- VarianceComponentExample$KinshipF
  KinD <- VarianceComponentExample$KinshipD
  
  ##Build incidence matrix for each random effect
  Zf <- t(sapply(as.character(DataHybrid$CodeFlint), function(x)
    as.numeric(rownames(KinF)==x)))
  Zd <- t(sapply(as.character(DataHybrid$CodeDent), function(x)
    as.numeric(rownames(KinD)==x)))
  
  ##Build the VarList and ZList objects
  VL = list(Flint=KinF , Dent=KinD , Error = diag(1,nrow(DataHybrid)))
  ZL <- list(Flint=Zf , Dent=Zd , Error = diag(1,nrow(DataHybrid)))
  
  ##Perform inference
  #A first way to call MMEst
  ResultVA <- MMEst(Y=DataHybrid$Trait  , Cofactor = matrix(DataHybrid$Trial)
                    , ZList = ZL  ,  VarList = VL)
  BlupVA <- MMBlup(Y=DataHybrid$Trait  , Cofactor = matrix(DataHybrid$Trial)
                    , ZList = ZL  ,  VarList = VL , ResMM=ResultVA)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ blup estimator }
\keyword{ random effect }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
