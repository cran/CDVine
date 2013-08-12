\name{CDVineClarkeTest}     
\alias{CDVineClarkeTest}  

\title{Clarke test comparing two vine copula models}

\description{
This function performs a Clarke test between two d-dimensional C- or D-vine copula models, respectively.
}

\usage{
CDVineClarkeTest(data, Model1.order=1:dim(data)[2],
                 Model2.order=1:dim(data)[2], Model1.family,
                 Model2.family, Model1.par, Model2.par,
                 Model1.par2=rep(0,dim(data)[2]*(dim(data)[2]-1)/2),
                 Model2.par2=rep(0,dim(data)[2]*(dim(data)[2]-1)/2),
                 Model1.type, Model2.type)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins).}
  \item{Model1.order, Model2.order}{Two numeric vectors giving the order of the variables in the first D-vine trees
    or of the C-vine root nodes in models 1 and 2
    (default: \code{Model1.order} and \code{Model2.order = 1:dim(data)[2]}, i.e., standard order).}
  \item{Model1.family, Model2.family}{Two d*(d-1)/2 numeric vectors of the pair-copula families of models 1 and 2, respectively, with values\cr		
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
  \item{Model1.par, Model2.par}{Two d*(d-1)/2 numeric vectors of the (first) copula parameters of models 1 and 2, respectively.}
  \item{Model1.par2, Model2.par2}{Two d*(d-1)/2 numeric vectors of the second copula parameters of models 1 and 2, respectively;
    necessary for t, BB1, BB6, BB7 and BB8 copulas. If no such families are included in \code{Model1.family}/\code{Model2.family},
    these arguments do not need to be specified
    (default: \code{Model1.par2} and \code{Model2.par2 = rep(0,dim(data)[2]*(dim(data)[2]-1)/2)}).}
  \item{Model1.type, Model2.type}{Type of the respective vine model:\cr
    \code{1} or \code{"CVine"} = C-vine\cr
		\code{2} or \code{"DVine"} = D-vine}
}

\details{
The test proposed by Clarke (2007) allows to compare non-nested models.
For this let \eqn{c_1} and \eqn{c_2} be two competing vine copulas in terms of their densities
and with estimated parameter sets \eqn{\hat{\boldsymbol{\theta}}_1}{\theta_1} and \eqn{\hat{\boldsymbol{\theta}}_2}{\theta_2}.
The null hypothesis of statistical indistinguishability of the two models is 
\deqn{
H_0: P(m_i > 0) = 0.5\ \forall i=1,..,N,
}{
H_0: P(m_i > 0) = 0.5 forall i=1,..,N,
}
where \eqn{m_i:=\log\left[\frac{c_1(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_1)}{c_2(\boldsymbol{u}_i|\hat{\boldsymbol{\theta}}_2)}\right]}{m_i:=log[ c_1(u_i|\theta_1) / c_2(u_i|\theta_2) ]}
for observations \eqn{\boldsymbol{u}_i,\ i=1,...,N}{u_i, i=1,...,N}.

Since under statistical equivalence of the two models the log likelihood ratios of the single observations
are uniformly distributed around zero and in expectation \eqn{50\%} of the log likelihood ratios greater than zero,
the tets statistic
\deqn{
\texttt{statistic} := B = \sum_{i=1}^N \mathbf{1}_{(0,\infty)}(m_i),
}{
statistic := B = \sum_{i=1}^N 1_{(0,\infty)}(m_i),
}
where \eqn{\mathbf{1}}{1} is the indicator function, is distributed Binomial with parameters \eqn{N} and \eqn{p=0.5},
and critical values can easily be obtained.
Model 1 is interpreted as statistically equivalent to model 2 if \eqn{B} is not significantly different from the expected value \eqn{Np = \frac{N}{2}}{np=N/2}.

Like AIC and BIC, the Clarke test statistic may be corrected for the number of parameters used in the models.
There are two possible corrections; 
the Akaike and the Schwarz corrections, which correspond to the penalty terms in the AIC and the BIC, respectively. 
}

\value{
  \item{statistic, statistic.Akaike, statistic.Schwarz}{Test statistics without correction, with Akaike correction and with Schwarz correction.}
  \item{p.value, p.value.Akaike, p.value.Schwarz}{P-values of tests without correction, with Akaike correction and with Schwarz correction.}
}

\seealso{\code{\link{CDVineVuongTest}}, \code{\link{CDVineAIC}}, \code{\link{CDVineBIC}}}

\author{Jeffrey Dissmann, Ulf Schepsmeier, Eike Brechmann}

\references{
Clarke, K. A. (2007).
A Simple Distribution-Free Test for Nonnested Model Selection.
Political Analysis, 15, 347-363.
}

\examples{
\dontrun{
# load data set
data(worldindices)
d = dim(worldindices)[2]

# select the C-vine families and parameters
cvine = CDVineCopSelect(worldindices,c(1:6),type="CVine")

# select the D-vine families and parameters
dvine = CDVineCopSelect(worldindices,c(1:6),type="DVine")

# compare the two models based on the data
clarke = CDVineClarkeTest(worldindices,1:d,1:d,cvine$family,dvine$family,
                          cvine$par,dvine$par,cvine$par2,dvine$par2,
                          Model1.type=1,Model2.type=2)
clarke$statistic
clarke$statistic.Schwarz
clarke$p.value
clarke$p.value.Schwarz
}
}

