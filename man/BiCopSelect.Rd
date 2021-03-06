\name{BiCopSelect}      
\alias{BiCopSelect}                     

\title{Selection and maximum likelihood estimation of bivariate copula families}

\description{
This function selects an appropriate bivariate copula family for given bivariate copula data
using one of a range of methods.
The corresponding parameter estimates are obtained by maximum likelihood estimation.
}

\usage{
BiCopSelect(u1, u2, familyset=NA, selectioncrit="AIC",
            indeptest=FALSE, level=0.05)
}

\arguments{
  \item{u1,u2}{Data vectors of equal length with values in [0,1].}
  \item{familyset}{Vector of bivariate copula families to select from
    (the independence copula MUST NOT be specified in this vector, otherwise it will be selected).
    The vector has to include at least one bivariate copula family that allows for positive and one that allows for negative dependence.
    If \code{familyset = NA} (default), selection among all possible families is performed.
    Coding of bivariate copula families: \cr
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
  \item{selectioncrit}{Character indicating the criterion for bivariate copula selection. Possible choices: \code{selectioncrit = "AIC"} (default) or \code{"BIC"}.}
  \item{indeptest}{Logical; whether a hypothesis test for the independence of \code{u1} and \code{u2} is performed before bivariate copula selection
    (default: \code{indeptest = FALSE}; cp. \code{\link{BiCopIndTest}}).
    The independence copula is chosen if the null hypothesis of independence cannot be rejected.}
  \item{level}{Numeric; significance level of the independence test (default: \code{level = 0.05}).}
}

\value{
  \item{family}{The selected bivariate copula family.}
  \item{par, par2}{The estimated bivariate copula parameter(s).}
  \item{p.value.indeptest}{P-value of the independence test if performed.}
}

\details{
Copulas can be selected according to the Akaike and Bayesian Information Criteria (AIC and BIC, respectively).
First all available copulas are fitted using maximum likelihood estimation.
Then the criteria are computed for all available copula families
(e.g., if \code{u1} and \code{u2} are negatively dependent, Clayton, Gumbel, Joe, BB1, BB6, BB7 and BB8 and their survival copulas are not considered)
and the family with the minimum value is chosen.
For observations \eqn{u_{i,j},\ i=1,...,N,\ j=1,2,}{u_{i,j}, i=1,...,N,\ j=1,2,} the AIC of a bivariate copula family \eqn{c} with parameter(s) \eqn{\boldsymbol{\theta}} is defined as
\deqn{
AIC := -2 \sum_{i=1}^N \ln[c(u_{i,1},u_{i,2}|\boldsymbol{\theta})] + 2k,
}{
AIC := -2 \sum_{i=1}^N ln[c(u_{i,1},u_{i,2}|\theta)] + 2k,
}
where \eqn{k=1} for one parameter copulas and \eqn{k=2} for the two parameter t-, BB1, BB6, BB7 and BB8 copulas.
Similarly, the BIC is given by
\deqn{
BIC := -2 \sum_{i=1}^N \ln[c(u_{i,1},u_{i,2}|\boldsymbol{\theta})] + \ln(N)k.
}{
BIC := -2 \sum_{i=1}^N ln[c(u_{i,1},u_{i,2}|\theta)] + ln(N)k.
}
Evidently, if the BIC is chosen, the penalty for two parameter families is stronger than when using the AIC.

Additionally a test for independence can be performed beforehand.
}

\note{
When the bivariate t-copula is considered and the degrees of freedom are estimated to be larger than 30,
then the bivariate Gaussian copula is taken into account instead.
Similarly, when BB1 (Clayton-Gumbel), BB6 (Joe-Gumbel), BB7 (Joe-Clayton) or BB8 (Joe-Frank) copulas are considered
and the parameters are estimated to be very close to one of their boundary cases,
the respective one parameter copula is taken into account instead. 
}

\author{Eike Brechmann, Jeffrey Dissmann}

\seealso{\code{\link{CDVineCopSelect}}, \code{\link{BiCopIndTest}}}

\references{
Akaike, H. (1973).
Information theory and an extension of the maximum likelihood principle.
In B. N. Petrov and F. Csaki (Eds.), Proceedings of the Second International
Symposium on Information Theory Budapest, Akademiai Kiado, pp. 267-281.

Brechmann, E. C. (2010).
Truncated and simplified regular vines and their applications.
Diploma thesis, Technische Universitaet Muenchen.\cr
\url{http://mediatum.ub.tum.de/doc/1079285/1079285.pdf}.

Manner, H. (2007).
Estimation and model selection of copulas with an application to exchange rates.
METEOR research memorandum 07/056, Maastricht University.

Schwarz, G. E. (1978).
Estimating the dimension of a model.
Annals of Statistics 6 (2), 461-464.
}

\examples{
## Example 1: Gaussian copula with large dependence parameter
par1 = 0.7
fam1 = 1
dat1 = BiCopSim(500,fam1,par1)

# select the bivariate copula family and estimate the parameter(s)
cop1 = BiCopSelect(dat1[,1],dat1[,2],familyset=c(1:10),indeptest=FALSE,level=0.05)
cop1$family
cop1$par
cop1$par2


## Example 2: Gaussian copula with small dependence parameter
par2 = 0.01
fam2 = 1
dat2 = BiCopSim(500,fam2,par2)

# select the bivariate copula family and estimate the parameter(s)
cop2 = BiCopSelect(dat2[,1],dat2[,2],familyset=c(1:10),indeptest=TRUE,level=0.05)
cop2$family
cop2$par
cop2$par2

\dontrun{
## Example 3: empirical data
data(worldindices)
cop3 = BiCopSelect(worldindices[,1],worldindices[,4],familyset=c(1:10,13,14,16,23,24,26))
cop3$family
cop3$par
cop3$par2
}
}
