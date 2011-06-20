\name{CDVineMLE}        
\alias{CDVineMLE}

\title{Maximum likelihood estimation of C- and D-vine copula models}

\description{
This function calculates the MLE of C- or D-vine copula model parameters using sequential estimates as initial values (if not provided).
}

\usage{
CDVineMLE(data, family, start=NULL, start2=NULL, type, maxit=200,
          max.df=30, max.BB=list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1)))
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{family}{A d*(d-1)/2 integer vector of C-/D-vine pair-copula families with values \cr
		\code{0} = independence copula \cr
	        \code{1} = Gaussian copula \cr
	        \code{2} = Student t copula (t-copula) \cr
	        \code{3} = Clayton copula \cr
	        \code{4} = Gumbel copula \cr
	        \code{5} = Frank copula \cr
	        \code{6} = Joe copula \cr 
		\code{7} = BB1 copula \cr
		\code{8} = BB6 copula \cr
		\code{9} = BB7 copula \cr
		\code{10} = BB8 copula \cr
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr 
		\code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
		\code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
		\code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
		\code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{27} = rotated BB1 copula (90 degrees) \cr
		\code{28} = rotated BB6 copula (90 degrees) \cr
		\code{29} = rotated BB7 copula (90 degrees) \cr
		\code{30} = rotated BB8 copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees) \cr
		\code{37} = rotated BB1 copula (270 degrees) \cr
		\code{38} = rotated BB6 copula (270 degrees) \cr
		\code{39} = rotated BB7 copula (270 degrees) \cr
		\code{40} = rotated BB8 copula (270 degrees)
		}
  \item{start}{A d*(d-1)/2 numeric vector of starting values for C-/D-vine pair-copula parameters
    (optional; otherwise they are calculated via \code{\link{CDVineSeqEst}}; default: \code{start = NULL}).}
  \item{start2}{A d*(d-1)/2 numeric vector of starting values for second C-/D-vine pair-copula parameters
    (optional; otherwise they are calculated via \code{\link{CDVineSeqEst}}; default: \code{start2 = NULL}).}
  \item{type}{Type of the vine model:\cr
    \code{1} or \code{"CVine"} = C-vine\cr
		\code{2} or \code{"DVine"} = D-vine}
  \item{maxit}{The maximum number of iteration steps (optional; default: \code{maxit = 200}).}
  \item{max.df}{Numeric; upper bound for the estimation of the degrees of freedom parameter of the t-copula
    (default: \code{max.df = 30}; for more details see \code{\link{BiCopEst}}).}
  \item{max.BB}{List; upper bounds for the estimation of the two parameters (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr
    (default: \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).}
}

\value{
  \item{par}{Estimated (first) C-/D-vine pair-copula parameters.}
  \item{par2}{Estimated second C-/D-vine pair-copula parameters for families with two parameters (t, BB1,BB6, BB7, BB8). All other entries are zero.}
  \item{loglik}{Optimized log-likelihood value corresponding to the estimated pair-copula parameters.}
  \item{convergence}{An integer code indicating either successful convergence (\code{convergence = 0})
    or an error (cp. \code{\link{optim}}; the CDVineMLE-function uses the "L-BFGS-B" method):\cr 
    \code{1} = the iteration limit \code{maxit} has been reached \cr
    \code{51} = a warning from the "L-BFGS-B" method; see component \code{message} for further details \cr
    \code{52} = an error from the "L-BFGS-B" method; see component \code{message} for further details}
  \item{message}{A character string giving any additional information returned by \code{\link{optim}}, or \code{NULL}.} 
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.
}

\author{Carlos Almeida, Ulf Schepsmeier}

\seealso{\code{\link{CDVineLogLik}}, \code{\link{CDVineSeqEst}}}

\examples{
## Example 1: 4-dimensional D-vine model with Gaussian pair-copulas
data(worldindices)
Data = as.matrix(worldindices)[,1:4]
fam = rep(1,6)

# maximum likelihood estimation
CDVineMLE(Data,family=fam,type=2,maxit=100)


## Example 2: 4-dimensional D-vine model with mixed pair-copulas
fam2 = c(5,1,3,14,3,2)

# sequential estimation
m = CDVineSeqEst(Data,family=fam2,type=2)
m

# calculate the log-likelihood
LogLik0 = CDVineLogLik(Data,fam2,m$par,m$par2,type=2)
LogLik0$loglik

# maximum likelihood estimation
CDVineMLE(Data,family=fam2,type=2,maxit=5)  # 5 iterations
CDVineMLE(Data,family=fam2,type=2)  # default: 200 iterations
}

