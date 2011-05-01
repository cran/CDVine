\name{BiCopHfunc}           
\alias{BiCopHfunc}

\title{Conditional distribution function (h-function) of a bivariate copula}

\description{
This function evaluates the conditional distribution function (h-function) of a given parametric bivariate copula.
}

\usage{
BiCopHfunc(u1, u2, family, par, par2=0)
}

\arguments{
  \item{u1,u2}{Numeric vectors of equal length with values in [0,1].}
  \item{family}{An integer defining the bivariate copula family: \cr
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
  \item{par}{Copula parameter.}
  \item{par2}{Second parameter for bivariate copulas with two parameters (t, BB1, BB7; default: \code{par2 = 0}).}
}

\details{
The h-function is defined as the conditional distribution function of a bivariate copula, i.e.,
\deqn{
h(u|v,\boldsymbol{\theta}) := F(u|v) =
\frac{\partial C(u,v)}{\partial v},
}{
h(u|v,\theta) := F(u|v) =
\partial C(u,v) / \partial v,
}
where \eqn{C} is a bivariate copula distribution function with parameter(s) \eqn{\boldsymbol{\theta}}{\theta}.  
For more details see Aas et al. (2009).
}


\value{
\item{hfunc1}{Numeric vector of the conditional distribution function (h-function) evaluated at \code{u2} given \code{u1}, i.e., \eqn{h(\code{u2}|\code{u1},\boldsymbol{\theta})}{h(u2|u1,\theta)}.}
\item{hfunc2}{Numeric vector of the conditional distribution function (h-function) evaluated at \code{u1} given \code{u2}, i.e., \eqn{h(\code{u1}|\code{u2},\boldsymbol{\theta})}{h(u1|u2,\theta)}.}
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.
}


\author{Ulf Schepsmeier}

\seealso{\code{\link{BiCopPDF}}, \code{\link{BiCopCDF}}, \code{\link{CDVineLogLik}}, \code{\link{CDVineSeqEst}}}

\examples{
## Example 1: 4-dimensional C-vine model with mixed pair-copulas
data(worldindices)
Data = as.matrix(worldindices)[,1:4]
d = dim(Data)[2]
fam = c(5,1,3,14,3,2)

# sequential estimation
seqpar1 = CDVineSeqEst(Data,fam,type=1,method="itau")

# calculate the inputs of the second tree using h-functions
h1 = BiCopHfunc(Data[,1],Data[,2],fam[1],seqpar1$par[1])
h2 = BiCopHfunc(Data[,1],Data[,3],fam[2],seqpar1$par[2])
h3 = BiCopHfunc(Data[,1],Data[,4],fam[3],seqpar1$par[3])

# compare estimated parameters
BiCopEst(h1$hfunc1,h2$hfunc1,fam[4],method="itau")
seqpar1$par[4]

BiCopEst(h1$hfunc1,h3$hfunc1,fam[5],method="itau")
seqpar1$par[5]


## Example 2: 4-dimensional D-vine model with mixed pair-copulas

# sequential estimation
seqpar2 = CDVineSeqEst(Data,fam,type=2,method="itau")

# calculate the inputs of the second tree using h-functions
h1 = BiCopHfunc(Data[,1],Data[,2],fam[1],seqpar2$par[1])
h2 = BiCopHfunc(Data[,2],Data[,3],fam[2],seqpar2$par[2])
h3 = BiCopHfunc(Data[,3],Data[,4],fam[3],seqpar2$par[3])

# compare estimated parameters
BiCopEst(h1$hfunc2,h2$hfunc1,fam[4],method="itau")
seqpar2$par[4]

BiCopEst(h2$hfunc2,h3$hfunc1,fam[5],method="itau")
seqpar2$par[5]
}
