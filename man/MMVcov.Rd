\name{MMVcov}
\alias{MMVcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Covariance Matrix for variance estimators.
%%  ~~function to do ... ~~
}
\description{
This function computes the covariance matrix of variance estimators using either the inverse of the Expected Hessian Matrix or the inverse of the Average Information Matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
%%, as in the default vcov() function for lme4,
%%, as in predict() function for ASReml
}
\usage{
MMVcov(ResMM , Y , Cofactor = NULL , formula = NULL,
    ZList = NULL , VarList , information="Expected")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ResMM}{
	A list as displayed by the \code{\link{MMEst}} function for a Variance Component Analysis (only the first element of the list will be analyzed).
}
  \item{Y}{
	The vector of response values used in the function \code{\link{MMEst}}.
}
  \item{Cofactor}{
  The incidence matrix corresponding to fixed effects common to all models to be adjusted used in the function \code{\link{MMEst}}. If \code{NULL}, a vector full of 1 is used.
}
  \item{formula}{
  The formula object specifying the fixed effect part of all models separated by + operators used in the function \code{\link{MMEst}} (default is \code{NULL}).
}
  \item{ZList}{
	The list of incidence matrices associated with random and residual effects used in the function \code{\link{MMEst}} (default is \code{NULL}).

}
  \item{VarList}{
		The list of covariance matrices associated with random and residual effects used in the function \code{\link{MMEst}}.
}

  \item{information}{
	A string specifying the method used to approximate the covariance matrix. It can be either "Expected" (default) to use the Expected Hessian Matrix or "AI" to use the Average Information Matrix. The AI matrix is always computed using Reml estimates whereas the expected hessian matrix could also be used for ML estimates.
}

}
\details{
	If \code{information} is not specified then the algorithm computes the covariance matrix using the Expected matrix based on the inference method (Reml or ML) used in the first item of \code{ResMM}. If \code{information} is equal to "AI" then it computes the AI matrix to calculate the covariance matrix. Only the first item of \code{ResMM} is analyzed.

	
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  The output of the function is a list:
  \item{vcov }{The covariance matrix of the variance estimators}
  \item{SE }{The standard errors of the variance estimators (the square root of the covariance matrix diagonal)}

%%	The output of the function is a matrix. This matrix is the variance covariance matrix of the variance estimates of the model fitted in \code{ResMM}. To extract the standard errors of variance estimates, one has to extract the diagonal of this matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{F. Laporte and T. Mary-Huard
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
  Expected_vcov <- MMVcov(ResMM=ResultVA,Y=DataHybrid$Trait,
                  Cofactor = matrix(DataHybrid$Trial),
                  , ZList = ZL  ,  VarList = VL)
  AI_vcov <- MMVcov(ResMM=ResultVA,Y=DataHybrid$Trait,
                  Cofactor = matrix(DataHybrid$Trial),
                  , ZList = ZL  ,  VarList = VL , information = "AI")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ covariance matrix }
\keyword{ standard errors }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
