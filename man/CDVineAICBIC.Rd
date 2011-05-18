\name{CDVineAIC/BIC}        
\alias{CDVineAIC}
\alias{CDVineBIC}                

\title{AIC and BIC of C- and D-vine copula models}

\description{
These functions calculate the Akaike and Bayesian Information criteria of d-dimensional C- and D-vine copula models for a given copula data set.
}

\usage{
CDVineAIC(data, family, par, 
          par2=rep(0,dim(data)[2]*(dim(data)[2]-1)/2), type)
CDVineBIC(data, family, par, 
          par2=rep(0,dim(data)[2]*(dim(data)[2]-1)/2), type)          
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
  \item{par2}{A d*(d-1)/2 vector of second parameters for two parameter pair-copula families
  (default: \code{par2 = rep(0,dim(data)[2]*(dim(data)[2]-1)/2)}).}
  \item{type}{Type of the vine model:\cr
    \code{1} or \code{"CVine"} = C-vine\cr
		\code{2} or \code{"DVine"} = D-vine}
}

\value{
  \item{AIC, BIC}{The computed AIC or BIC value, respectively.}
  \item{pair.AIC, pair.BIC}{An array of individual contributions to the AIC or BIC value
    for each pair-copula, respectively. Note: \code{AIC = sum(pair.AIC)} and similarly \code{BIC = sum(pair.BIC)}.}
}

\details{
If \eqn{k} denotes the number of parameters of a C-vine copula model with log-likelihood \eqn{l_{CVine}} and parameter set \eqn{\boldsymbol{\theta}}{\theta},
then the Akaike Information Criterion (AIC) by Akaike (1973) is defined as
\deqn{
AIC := -2 l_{CVine}\left(\boldsymbol{\theta}|\boldsymbol{u}\right) + 2 k,
}{
AIC := -2 l_{CVine}(\theta|u) + 2 k,
}
for observations \eqn{\boldsymbol{u}=(\boldsymbol{u}_1^\prime,...,\boldsymbol{u}_N^\prime)^\prime}{u=(u'_1,...u'_N)}.

Similarly, the Bayesian Information Criterion (BIC) by Schwarz (1978) is given by
\deqn{
BIC := -2 l_{CVine}\left(\boldsymbol{\theta}|\boldsymbol{u}\right) + \log(N) k.
}{
BIC := -2 l_{CVine}(\theta|u) + log(N) k.
}

The AIC and BIC expressions for D-vine copula models are defined accordingly.
}

\references{ 
Akaike, H. (1973).
Information theory and an extension of the maximum likelihood principle.
In B. N. Petrov and F. Csaki (Eds.), Proceedings of the Second International
Symposium on Information Theory Budapest, Akademiai Kiado, pp. 267-281.

Schwarz, G. E. (1978).
Estimating the dimension of a model.
Annals of Statistics 6 (2), 461-464.
}

\author{Eike Brechmann}

\seealso{\code{\link{CDVineLogLik}}, \code{\link{CDVineVuongTest}}, \code{\link{CDVineClarkeTest}}}

\examples{
## Example 1: 3-dimensional D-vine model with Gaussian pair-copulas
data(worldindices)
Data = as.matrix(worldindices)[,1:3]
fam1 = c(1,1,1) 
par1 = c(0.2,0.3,0.4)

# calculate AIC and BIC
CDVineAIC(Data,fam1,par1,type=2)
CDVineBIC(Data,fam1,par1,type=2)


## Example 2: 6-dimensional C-vine model with Student t pair-copulas
## with 5 degrees of freedom
data(worldindices)
Data = as.matrix(worldindices)
dd = dim(Data)[2]*(dim(Data)[2]-1)/2
fam2 = rep(2,dd)
par2 = rep(0.5,dd)
nu2 = rep(5,dd)

# calculate AIC and BIC
CDVineAIC(Data,fam2,par2,nu2,type=1)
CDVineBIC(Data,fam2,par2,nu2,type=1)


## Example 3: 4-dimensional C-vine model with mixed pair-copulas
fam3 = c(5,1,3,14,3,2)
par3 = c(0.9,0.3,0.2,1.1,0.2,0.7)
nu3 = c(0,0,0,0,0,7)

# calculate AIC and BIC
CDVineAIC(Data[,1:4],fam3,par3,nu3,type=1)
CDVineBIC(Data[,1:4],fam3,par3,nu3,type=1)
}


