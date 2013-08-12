\name{CDVineTreePlot}    
\alias{CDVineTreePlot}

\title{Plot function for C- or D-vine trees}

\description{
This function plots one or all trees of a given C- and D-vine copula model. 
}

\usage{
CDVineTreePlot(data=NULL, family,
               par=rep(0,length(family)), par2=rep(0,length(family)),
               names=NULL, type, method="mle", max.df=30,
               max.BB=list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1)),
               tree="ALL", edge.labels=c("family"), P=NULL, ...)
}

\arguments{
  \item{data}{An N x d data matrix (with uniform margins); default: \code{data = NULL}.}
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
  \item{par}{A d*(d-1)/2 vector of pair-copula parameters\cr
    (optional; default: \code{par = rep(0,length(family))}.}
  \item{par2}{A d*(d-1)/2 vector of second parameters for pair-copula families with two parameters
    (optional; default: \code{par2 = rep(0,length(family))}.}
  \item{names}{A vector of names for the d variables.
    If \code{names = NULL} (default), column names of \code{data} (if available) or simply numbers are used.}
  \item{type}{Type of the vine model:\cr
    \code{1} or \code{"CVine"} = C-vine\cr
		\code{2} or \code{"DVine"} = D-vine}
  \item{method}{Character indicating the estimation method:
    either maximum likelihood estimation (\code{method = "mle"}; default) or inversion of Kendall's tau (\code{method = "itau"}).} 
  \item{max.df}{Numeric; upper bound for the estimation of the degrees of freedom parameter of the t-copula
    (default: \code{max.df = 30}; for more details see \code{\link{BiCopEst}}).}
  \item{max.BB}{List; upper bounds for the estimation of the two parameters (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr
    (default: \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).}
  \item{tree}{Number of the tree to be plotted or \code{tree = "ALL"} (default) to plot all trees.}
  \item{edge.labels}{Vector of edge labels. Possible choices:\cr
    \code{FALSE}: no edge labels\cr
    \code{"family"}: pair-copula families (default)\cr
    \code{"par"}: pair-copula parameters\cr
    \code{"par2"}: second pair-copula parameters\cr
    \code{"theotau"}: theoretical Kendall's tau values corresponding to pair-copula families and parameters (see \code{\link{BiCopPar2Tau}})\cr
    \code{"emptau"}: empirical Kendall's tau values (only if data is provided!)}
  \item{P}{A list of matrices with two columns for the x-y-coordinates of the nodes in the plot(s) (optional; default: \code{P = NULL}).}
  \item{...}{further graphical parameters}
}

\note{
The function computes the positions of the nodes automatically with the Fruchterman-Reingold algorithm (see \code{\link{plot.igraph}} for a detailed description).
If one would like to set the positions manually, one has to specify a list of matrices \code{P} in the argument list.
A good starting point may be to run the function \code{\link{CDVineTreePlot}} and manipulate the returning matrix P.

The user can set the copula parameters \code{par} and \code{par2}.
If not set and data is provided, the parameters of the C- or D-vine copula model are estimated sequentially using
\code{\link{CDVineSeqEst}}/\code{\link{BiCopEst}}.
Then the edge width is chosen according to the empirical Kendall's tau values. Otherwise theoretical values are used.
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.
}

\author{Ulf Schepsmeier}

\seealso{\code{\link{CDVineMLE}}}

\examples{
# simulate from a 6-dimensional C-vine model with Gaussian pair-copulas
d = 6
dd = d*(d-1)/2
par1 = c(0.2,0.69,0.73,0.22,-0.09,0.51,0.32,0.01,0.82,0.01,-0.2,
         -0.32,-0.19,-0.17,-0.06);
fam1 = rep(1,dd)
N = 100
U = CDVineSim(N,fam1,par1,type=1)

# plot the first tree with pair-copula families and
# empirical Kendall's tau values as edge labels
CDVineTreePlot(U,fam1,type=1,tree=1,edge.labels=c("family","emptau"))

# plot all trees without edge labels and without sequential estimation
CDVineTreePlot(data=NULL,family=fam1,par=par1,type=1, edge.labels=FALSE)
}
