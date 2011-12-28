\name{CDVineSeqEst}       
\alias{CDVineSeqEst}

\title{Sequential estimation of C- and D-vine copula models}

\description{
This function sequentially estimates the pair-copula parameters of d-dimensional C- or D-vine copula models.
}

\usage{
CDVineSeqEst(data, family, type, method="mle", se=FALSE, max.df=30,
             max.BB=list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1)),
             progress=FALSE)
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
  \item{type}{Type of the vine model:\cr
    \code{1} or \code{"CVine"} = C-vine\cr
		\code{2} or \code{"DVine"} = D-vine}
  \item{method}{Character indicating the estimation method:
    either pairwise maximum likelihood estimation (\code{method = "mle"}; default) or inversion of Kendall's tau (\code{method = "itau"}; see \code{\link{BiCopEst}}).
		For \code{method = "itau"} only one parameter pair-copula families can be used 
		(\code{family = 1}, \code{3}, \code{4}, \code{5}, \code{6}, \code{13}, \code{14}, \code{16}, \code{23}, \code{24}, \code{26}, \code{33}, \code{34} or \code{36}).}
  \item{se}{Logical; whether standard errors are estimated (default: \code{se=FALSE}).}
  \item{max.df}{Numeric; upper bound for the estimation of the degrees of freedom parameter of the t-copula
    (default: \code{max.df = 30}; for more details see \code{\link{BiCopEst}}).}
  \item{max.BB}{List; upper bounds for the estimation of the two parameters (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr
    (default: \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).}
  \item{progress}{Logical; whether the pairwise estimation progress is printed (default: \code{progress = FALSE}).}    
}

\details{
The pair-copula parameter estimation is performed tree-wise,
i.e., for each C-/D-vine tree the results from the previous tree(s) are used
to calculate the new copula parameters using \code{\link{BiCopEst}}.
}


\value{
  \item{par}{Estimated (first) C-/D-vine pair-copula parameters.}
  \item{par2}{Estimated second C-/D-vine pair-copula parameters for families with two parameters (t, BB1, BB6, BB7, BB8). All other entries are zero.}
  \item{se}{Estimated standard errors of the (first) pair-copula parameter estimates \cr
  (if \code{se = TRUE}).}
  \item{se2}{Estimated standard errors of the second pair-copula parameter estimates \cr
  (if \code{se = TRUE}).}
}    

\references{
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.

Czado, C., U. Schepsmeier, and A. Min (2011).
Maximum likelihood estimation of mixed C-vines with application to exchange rates.
Statistical Modelling, to appear.
}

\author{Carlos Almeida, Ulf Schepsmeier}

\seealso{\code{\link{BiCopEst}}, \code{\link{BiCopHfunc}}, \code{\link{CDVineLogLik}}, \code{\link{CDVineMLE}}}

\examples{
## Example 1: 4-dimensional D-vine model with Gaussian pair-copulas
data(worldindices)
Data = as.matrix(worldindices)[,1:4]
d = dim(Data)[2]
fam = rep(1,d*(d-1)/2)

# sequential estimation 
CDVineSeqEst(Data,fam,type=2,method="itau")$par
CDVineSeqEst(Data,fam,type=2,method="mle")$par


## Example 2: 4-dimensional D-vine model with mixed pair-copulas
fam2 = c(5,1,3,14,3,2)

# sequential estimation
CDVineSeqEst(Data,fam2,type=2,method="mle",se=TRUE,progress=TRUE)
}


