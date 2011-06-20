\name{CDVinePar2Tau}        
\alias{CDVinePar2Tau}

\title{Kendall's tau values of a vine copula model}

\description{
This function computes the values of Kendall's tau corresponding to the parameters of a C- or D-vine copula model.
}

\usage{
CDVinePar2Tau(family, par, par2=rep(0,length(family)))
}

\arguments{
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
  \item{par}{A d*(d-1)/2 vector of pair-copula parameters.}
  \item{par2}{A d*(d-1)/2 vector of second parameters for pair-copula families with two parameters (t, BB1, BB6, BB7, BB8; 
  default: \code{par2 = rep(0,length(family))}).}
}

\value{
A d*(d-1)/2 vector of theoretical Kendall's tau values corresponding to the given pair-copula families and parameters.
}

\author{Eike Brechmann}

\seealso{\code{\link{BiCopPar2Tau}}}

\examples{
# specify the vine
fam3 = c(5,1,3,14,3,2)
par3 = c(0.9,0.3,0.2,1.1,0.2,0.7)
nu3 = c(0,0,0,0,0,7)

# compute the corresponding Kendall's tau values 
tau = CDVinePar2Tau(fam3,par3,nu3)
}

