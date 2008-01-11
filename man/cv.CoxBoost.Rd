\name{cv.CoxBoost}
\alias{cv.CoxBoost}
\title{Determines the optimal number of boosting steps by cross-validation}
\description{
Performs a K-fold cross-validation for \code{\link{CoxBoost}} in search for the optimal number of boosting steps.  
}
\usage{
cv.CoxBoost(time,status,x,maxstepno=100,K=10,type=c("verweij","naive"),
            trace=FALSE,\dots) 
}
\arguments{
\item{time}{vector of length \code{n} specifying the observed times.}
\item{status}{censoring indicator, i.e., vector of length \code{n} with entries \code{0} for censored observations and \code{1} for uncensored observations.}
\item{x}{\code{n * p} matrix of covariates.}
\item{maxstepno}{maximum number of boosting steps to evaluate, i.e, the returned ``optimal'' number of boosting steps will be in the range \code{[1,maxstepno]}.}
\item{K}{number of folds to be used for cross-validation. If \code{K} is larger or equal to the number of non-zero elements in \code{status}, leave-one-out cross-validation is performed.}
\item{type}{way of calculating the partial likelihood contribution of the observation in the hold-out folds: \code{"verweij"} uses the more appropriate method described in Verweij and van Houwelingen (1996), \code{"naive"} uses the approach where the observations that are not in the hold-out folds are ignored (often found in other R packages).}
\item{trace}{logical value indicating whether progress in estimation should be indicated by printing the number of the cross-validation fold and the index  of the covariate updated.}
\item{\dots}{miscellaneous parameters for the calls to \code{\link{CoxBoost}}}
}
\value{
List with the following components:
\item{mean.logplik}{vector of length \code{maxstepno} with the mean partial log-likelihood for each boosting step}
\item{se.logplik}{vector with standard error estimates for the mean partial log-likelihood criterion for each boosting step.}
\item{optimal.step}{optimal boosting step number, i.e., with minimum mean partial log-likelihood.}
}
\author{
Harald Binder \email{binderh@fdm.uni-freiburg.de}
}
\references{
Verweij, P. J. M. and van Houwelingen, H. C. (1996). Cross-validation in survival analysis. Statistics in Medicine, 12(24):2305-2314.
}
\seealso{
\code{\link{CoxBoost}}, \code{\link{optimCoxBoostPenalty}}
}
\examples{
\dontrun{
#   Generate some survival data with 10 informative covariates 
n <- 200; p <- 100
beta <- c(rep(1,10),rep(0,p-10))
x <- matrix(rnorm(n*p),n,p)
real.time <- -(log(runif(n)))/(10*exp(drop(x \%*\% beta)))
cens.time <- rexp(n,rate=1/10)
status <- ifelse(real.time <= cens.time,1,0)
obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)


#  10-fold cross-validation

cv.res <- cv.CoxBoost(time=obs.time,status=status,x=x,maxstepno=500,
                      K=10,type="verweij",penalty=100) 

#   examine mean partial log-likelihood in the course of the boosting steps
plot(cv.res$mean.logplik)   

#   Fit with optimal number of boosting steps

cbfit <- CoxBoost(time=obs.time,status=status,x=x,stepno=cv.res$optimal.step,
                  penalty=100) 
summary(cbfit)

}

}
\keyword{models} \keyword{regression} \keyword{survial}
