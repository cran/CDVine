\name{BiCopSim}        
\alias{BiCopSim}

\title{Simulation from a bivariate copula}

\description{
This function simulates from a given parametric bivariate copula.
}

\usage{
BiCopSim(N, family, par, par2=0)
}

\arguments{
  \item{N}{Number of bivariate observations simulated.}
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

\value{
An \code{N} x 2 matrix of data simulated from the bivariate copula.
}

\author{Ulf Schepsmeier}

\seealso{\code{\link{BiCopCDF}}, \code{\link{BiCopPDF}}, \code{\link{CDVineSim}}}

\examples{
# simulate from a bivariate t-copula
simdata = BiCopSim(300,2,-0.7,par2=4)
}

