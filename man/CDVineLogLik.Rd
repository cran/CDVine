\name{CDVineLogLik}        
\alias{CDVineLogLik}                

\title{Log-likelihood of C- and D-vine copula models}

\description{
This function calculates the log-likelihood of d-dimensional C- and D-vine copula models for a given copula data set.
}

\usage{
CDVineLogLik(data, family, par, 
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

\details{
Let \eqn{\boldsymbol{u}=(\boldsymbol{u}_1^\prime,...,\boldsymbol{u}_N^\prime)^\prime}{u=(u'_1,...,u'_N)'} be \eqn{d}-dimensional observations
with \eqn{\boldsymbol{u}_i=(u_{i,1},...,u_{i,d})^\prime\in[0,1]^d,\ i=1,...,N}{u_i=(u_{i,1},...,u_{i,d})' in [0,1]^d}.
Then the log-likelihood of a C-vine copula is given by
\deqn{
\texttt{loglik}:=l_{CVine}\left(\boldsymbol{\theta}|\boldsymbol{u}\right)=
\sum_{i=1}^N \sum_{j=1}^{d-1} \sum_{k=1}^{d-j}
\ln\left[c_{j,j+k|1,...,j-1}\right],
}{
loglik:=l_{CVine}(\theta|u)=
\sum_{i=1}^N \sum_{j=1}^{d-1} \sum_{k=1}^{d-j}
ln[c_{j,j+k|1,...,j-1}],
}
where
\deqn{
c_{j,j+k|1,...,j-1}:=c_{j,j+k|1:(j-1)}(F(u_{i,j}|u_{i,1},...,u_{i,j-1}),F(u_{i,j+k}|u_{i,1},...,u_{i,j-1})|\boldsymbol{\theta}_{j,j+k|1,...,j-1})
}{
c_{j,j+k|1,...,j-1}:=c_{j,j+k|1:(j-1)}(F(u_{i,j}|u_{i,1},...,u_{i,j-1}),F(u_{i,j+k}|u_{i,1},...,u_{i,j-1})|\theta_{j,j+k|1,...,j-1})
}
denote pair-copulas with parameter(s) \eqn{\boldsymbol{\theta}_{j,j+k|1,...,j-1}}{\theta_{j,j+k|1,...,j-1}}.

Similarly, the log-likelihood of a \eqn{d}-dimensional D-vine copula is
\deqn{
\texttt{loglik}:=l_{DVine}\left(\boldsymbol{\theta}|\boldsymbol{u}\right)=
\sum_{i=1}^N \sum_{j=1}^{d-1} \sum_{k=1}^{d-j}
\ln\left[c_{k,k+j|k+1,...,k+j-1}\right],
}{
loglik:=l_{DVine}(\theta|u)=
\sum_{i=1}^N \sum_{j=1}^{d-1} \sum_{k=1}^{d-j}
ln[c_{k,k+j|k+1,...,k+j-1}],
}
again with pair-copula densities denoted by
\deqn{
c_{k,k+j|k+1,...,k+j-1}:=
}\deqn{
c_{k,k+j|k+1,...,k+j-1}(F(u_{i,k}|u_{i,k+1},...,u_{i,k+j-1}),F(u_{i,k+j}|u_{i,k+1},...,u_{i,k+j-1})|\boldsymbol{\theta}_{k,k+j|k+1,...,k+j-1}).
}{
c_{k,k+j|k+1,...,k+j-1}(F(u_{i,k}|u_{i,k+1},...,u_{i,k+j-1}),F(u_{i,k+j}|u_{i,k+1},...,u_{i,k+j-1})|\theta_{k,k+j|k+1,...,k+j-1}).
}

Conditional distribution functions in both expressions are obtained recursively using the relationship
\deqn{
h(u|\boldsymbol{v},\boldsymbol{\theta}) := F(u|\boldsymbol{v}) =
\frac{\partial C_{uv_j|\boldsymbol{v}_{-j}}(F(u|\boldsymbol{v}_{-j}),F(v_j|\boldsymbol{v}_{-j}))}{\partial F(v_j|\boldsymbol{v}_{-j})},
}{
h(u|v,\theta) := F(u|v) =
\partial C_{uv_j|v_{-j}}(F(u|v_{-j}),F(v_j|v_{-j})) / \partial F(v_j|v_{-j}),
}
where \eqn{C_{uv_j|\boldsymbol{v}_{-j}}}{C_{uv_j|v_{-j}}} is a bivariate copula distribution function with parameter(s) \eqn{\boldsymbol{\theta}}{\theta}
and \eqn{\boldsymbol{v}_{-j}}{v_{-j}} denotes a vector with the \eqn{j}-th component \eqn{v_j} removed.  
The notation of h-functions is introduced for convenience. For more details see Aas et al. (2009).

Note that both log-likelihoods can also be written as \eqn{\texttt{loglik}=\sum_{k=1}^{d(d-1)/2}ll_{k}}{loglik=\sum_{k=1}^{d(d-1)/2}ll_{k}},
where \eqn{ll_{k}} are the individual contributions to the log-likelihood of each pair-copula.
}

\value{
  \item{loglik}{The calculated log-likelihood value of the C- or D-vine copula model.}
  \item{ll}{An array of individual contributions to the log-likelihood
    for each pair-copula. Note: \code{loglik = sum(ll)}.}
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.
}

\author{Carlos Almeida, Ulf Schepsmeier}

\seealso{\code{\link{CDVineAIC}}, \code{\link{CDVineBIC}}, \code{\link{BiCopHfunc}}}

\examples{
## Example 1: 3-dimensional D-vine model with Gaussian pair-copulas
data(worldindices)
Data = as.matrix(worldindices)[,1:3]
fam1 = c(1,1,1) 
par1 = c(0.2,0.3,0.4)

# calculate the log-likelihood
logLik1 = CDVineLogLik(Data,fam1,par1,type=2)

# check the above formula
sum(logLik1$ll)
logLik1$loglik


## Example 2: 6-dimensional C-vine model with Student t pair-copulas
## with 5 degrees of freedom
data(worldindices)
Data = as.matrix(worldindices)
dd = dim(Data)[2]*(dim(Data)[2]-1)/2
fam2 = rep(2,dd)
par2 = rep(0.5,dd)
nu2 = rep(5,dd)

# calculate the log-likelihood
logLik2 = CDVineLogLik(Data,fam2,par2,nu2,type=1)
logLik2$loglik


## Example 3: 4-dimensional C-vine model with mixed pair-copulas
fam3 = c(5,1,3,14,3,2)
par3 = c(0.9,0.3,0.2,1.1,0.2,0.7)
nu3 = c(0,0,0,0,0,7)

# calculate the log-likelihood
logLik3 = CDVineLogLik(Data[,1:4],fam3,par3,nu3,type=2)
logLik3$loglik
}


