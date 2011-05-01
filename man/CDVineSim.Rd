\name{CDVineSim}               
\alias{CDVineSim}
                                    
\title{Simulation from C- and D-vine copula models}

\description{
This function simulates from given C- and D-vine copula models.
}

\usage{
CDVineSim(N, family, par, par2=rep(0,length(family)), type)
}

\arguments{
  \item{N}{Number of d-dimensional observations simulated.}
  \item{family}{A d*(d-1)/2 integer vector of C-/D-vine pair-copula families with values \cr
		\code{0} = independence copula \cr
    \code{1} = Gaussian copula \cr
    \code{2} = Student t copula (t-copula) \cr
    \code{3} = Clayton copula \cr
    \code{4} = Gumbel copula \cr
    \code{5} = Frank copula \cr
    \code{6} = Joe copula \cr 
		\code{7} = BB1 copula \cr
		\code{9} = BB7 copula \cr
		\code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
		\code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
		\code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
		\code{23} = rotated Clayton copula (90 degrees) \cr
		\code{24} = rotated Gumbel copula (90 degrees) \cr
		\code{26} = rotated Joe copula (90 degrees) \cr
		\code{33} = rotated Clayton copula (270 degrees) \cr
		\code{34} = rotated Gumbel copula (270 degrees) \cr
		\code{36} = rotated Joe copula (270 degrees)}
  \item{par}{A d*(d-1)/2 vector of pair-copula parameters.}
  \item{par2}{A d*(d-1)/2 vector of second parameters for pair-copula families with two parameters (t, BB1, BB7; default: \code{par2 = rep(0,length(family))}).}
  \item{type}{Type of the vine model:\cr
    \code{1} or \code{"CVine"} = C-vine\cr
		\code{2} or \code{"DVine"} = D-vine}
}

\value{
An \code{N} x d matrix of data simulated from the given C- or D-vine copula model.
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.
}

\author{Carlos Almeida, Ulf Schepsmeier, Eike Brechmann, Jakob Stoeber }

\seealso{\code{\link{CDVineLogLik}}, \code{\link{BiCopSim}}}

\examples{
## Example 1: simulate from a 6-dimensional C-vine model 
## with Gaussian pair-copulas
d = 6
dd = d*(d-1)/2
fam1 = rep(1,dd)
par1 = c(0.2,0.69,0.73,0.22,-0.09,0.51,0.32,0.01,0.82,0.01,
         -0.2,-0.32,-0.19,-0.17,-0.06)
N = 100
U1 = CDVineSim(N,fam1,par1,type=1)
head(U1)

# calculate the log-likelihood
logLik = CDVineLogLik(U1,fam1,par1,type=1)
logLik$loglik


## Example 2: simulate from a 6-dimensional C-vine model 
## with Student t pair-copulas each with three degrees of freedom
fam2 = rep(2,dd)
par2 = rep(0.5,dd)
nu2 = rep(3,dd)
U2 = CDVineSim(N,fam2,par2,nu2,type=1)

# calculate the log-likelihood
logLik2 = CDVineLogLik(U2,fam2,par2,nu2,type=1)
logLik2$loglik


## Example 3: simulate from a 6-dimensional D-vine model
## with Student t pair-copulas each with five degrees of freedom
fam3 = rep(2,dd)
par3 = rep(0.4,dd)
nu3 = rep(5,dd)
U3 = CDVineSim(N,fam3,par3,nu3,type=2)

# calculate the log-likelihood
logLik3 = CDVineLogLik(U3,fam3,par3,nu3,type=2)
logLik3$loglik
}

