\name{CDVineCopSelect}       
\alias{CDVineCopSelect}

\title{Sequential copula selection and estimation of C- and D-vine copula models}

\description{
This function fits either a C- or a D-vine copula model to a d-dimensional copula data set.
Appropriate pair-copula families are selected using \code{\link{BiCopSelect}} and estimated sequentially.
}

\usage{
CDVineCopSelect(data, familyset=NA, type, selectioncrit="AIC",
                indeptest=FALSE, level=0.05)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{familyset}{An integer vector of pair-copula families to select from (the independence copula MUST NOT be specified in this vector unless one wants to fit an independence vine!).
    The vector has to include at least one pair-copula family that allows for positive and one that allows for negative dependence.
    If \code{familyset = NA} (default), selection among all possible families is performed.
    The coding of pair-copula families is shown below.}
  \item{type}{Type of the vine model:\cr
		\code{1} or \code{"CVine"} = C-vine\cr
		\code{2} or \code{"DVine"} = D-vine}  
  \item{selectioncrit}{Character indicating the criterion for pair-copula selection. Possible choices: \code{selectioncrit = "AIC"} (default) or \code{"BIC"} (see \code{\link{BiCopSelect}}).}
  \item{indeptest}{Logical; whether a hypothesis test for the independence of \code{u1} and \code{u2} is performed before bivariate copula selection
    (default: \code{indeptest = FALSE}; cp. \code{\link{BiCopIndTest}}).
    The independence copula is chosen for a (conditional) pair if the null hypothesis of independence cannot be rejected.}
  \item{level}{Numeric; significance level of the independence test (default: \code{level = 0.05}).}
}


\value{
  \item{family}{A d*(d-1)/2 vector of pair-copula families with values\cr
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
  \item{par2}{A d*(d-1)/2 vector of second pair-copula parameters for the t-, BB1, BB6, BB7 and BB8 copulas.}
}

\author{Eike Brechmann}

\seealso{\code{\link{BiCopSelect}}, \code{\link{CDVineSeqEst}}}

\examples{
# simulate from a 4-dimensional D-vine copula with mixed pair-copulas
d = 4
dd = d*(d-1)/2
family = c(1,2,3,4,7,3)
par = c(0.5,0.4,2,1.5,1.2,1.5)
par2 = c(0,5,0,0,2,0)
type = 2
simdata = CDVineSim(1000,family,par,par2,type)

# determine appropriate pair-copula families and parameters
# of a D-vine structure
CDVineCopSelect(simdata,type=2,familyset=c(1:10,13,14,23,24))
}

