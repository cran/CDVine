\name{BiCopPar2Tau}      
\alias{BiCopPar2Tau} 

\title{Kendall's tau value of a bivariate copula}

\description{
This function computes the theoretical Kendall's tau value of a bivariate copula for given parameter values. 
}

\usage{
BiCopPar2Tau(family, par, par2=0)
}

\arguments{
  \item{family}{An integer defining the bivariate copula family:\cr
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
  \item{par2}{Second parameter for the two parameter BB1 and BB7 copulas (default: \code{par2 = 0}).
    Note that the degrees of freedom parameter of the t-copula does not need to be set,
    because the theoretical Kendall's tau value of the t-copula is independent of this choice.}
}                                                                                  

\value{
Theoretical value of Kendall's tau corresponding to the bivariate copula family and parameter(s)
(\eqn{\theta} for one parameter families and the first parameter of the t-copula,
\eqn{\theta} and \eqn{\delta} for the two parameter BB1 and BB2 copulas). 
\tabular{ll}{
No. \tab Kendall's tau \cr
\code{1, 2} \tab \eqn{\frac{2}{\pi}\arcsin(\theta)}{2 / \pi arcsin(\theta)} \cr
\code{3, 13} \tab \eqn{\frac{\theta}{\theta+2}}{\theta / (\theta+2)} \cr
\code{4, 14} \tab \eqn{1-\frac{1}{\theta}}{1-1/\theta} \cr
\code{5} \tab \eqn{1-\frac{4}{\theta}+4\frac{D_1(\theta)}{\theta}}{1-4/\theta + 4 D_1(\theta)/\theta} \cr
	 \tab with \eqn{D_1(\theta)=\int_0^\theta \frac{c/\theta}{\exp(x)-1}dx}{D_1(\theta)=\int_0^\theta (c/\theta)/(exp(x)-1)dx} (Debey function) \cr
\code{6, 16} \tab \eqn{\frac{2\theta-4+2\gamma +2\log2+\Psi(\frac{1}{\theta})+\Psi(\frac{2+\theta}{2\theta})}{\theta-2}}{(2\theta-4+2\gamma + 2log2+\Psi(1/\theta)+\Psi((2+\theta)/(2\theta))) / (\theta-2)} \cr
	 \tab with \eqn{\gamma=\lim_{n\to\infty}(\sum_{i=1}^n \frac{1}{i} -\log n)\approx 0.57}{\gamma=lim_{n->\infty}(\sum_{i=1}^n 1/i - log n) ~= 0.57} (Euler's const.) \cr
	 \tab and \eqn{\Psi(x)=\frac{d}{dx}\log(\Gamma(x))}{\Psi(x)=d/dx log(\Gamma(x))} (Digamma function) \cr
\code{7} \tab \eqn{1-\frac{2}{\delta(\theta+2)}}{1-2/(\delta(\theta+2))} \cr
\code{9} \tab \eqn{1-\frac{2}{\delta(2-\theta)}+\frac{4}{\theta^2 \delta}B(\frac{2-\theta}{\theta},\delta+2)}{1-2/(\delta(2-\theta)) + 4/(\theta^2 \delta) B((2-\theta)/\theta , \delta+2)} \cr
	 \tab with \eqn{B(x,y)=\int_0^1 t^{x+1}(t-1)^{y-1}dt} (Beta function) \cr
\code{23, 33} \tab \eqn{\frac{\theta}{-\theta+2}}{\theta/(-\theta+2)} \cr
\code{24, 34} \tab \eqn{-1-\frac{1}{\theta}}{-1-1/\theta} \cr
\code{26, 36} \tab \eqn{\frac{2\theta+4-2*\gamma-2*\log2-\Psi(\frac{1}{-\theta})-\Psi(\frac{2-\theta}{-2\theta})}{\theta+2}}{(2\theta + 4 - 2*\gamma - 2*log2-\Psi(1/-\theta) - \Psi((2-\theta)/(-2\theta))) / (\theta+2)} \cr
}
}

\author{Ulf Schepsmeier}

\references{
Joe, H. (1997).
Multivariate Models and Dependence Concepts.
Chapman and Hall, London.

Czado, C., U. Schepsmeier, and A. Min (2011).
Maximum likelihood estimation of mixed C-vines with application to exchange rates.
Submitted for publication.
\url{http://www-m4.ma.tum.de/Papers/Schepsmeier/Paper.pdf}.
}

\seealso{\code{\link{CDVinePar2Tau}}, \code{\link{BiCopTau2Par}}}

\examples{
## Example 1: Gaussian copula
tt1 = BiCopPar2Tau(1,0.7)

# transform back
BiCopTau2Par(1,tt1)


## Example 2: Clayton copula
BiCopPar2Tau(3,1.3)
}

