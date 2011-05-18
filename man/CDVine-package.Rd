\name{CDVine-package}
\alias{CDVine-package}        
\alias{CDVine}
\docType{package}

\title{Statistical inference of C- and D-vine copulas}

\description{
This package provides functions for statistical inference of canonical vine (C-vine) and D-vine copulas.
It contains tools for bivariate exploratory data analysis and for bivariate as well as vine copula selection.
Models can be estimated either sequentially or by joint maximum likelihood estimation.
Sampling algorithms and plotting methods are also included.
Data is assumed to lie in the unit hypercube (so-called copula data).
}

\section{Bivariate copula families}{
In this package several bivariate copula families are included for bivariate analysis as well as for multivariate analysis using vine copulas.
It provides functionality of elliptical (Gaussian and Student t) as well as Archimedean (Clayton, Gumbel, Frank, Joe, BB1 and BB7) copulas
to cover a large bandwidth of possible dependence structures.
For the one parameter Archimedean copula families rotated versions are included to cover negative dependence too.
The two parameter BB1 and BB7 copulas are however numerically instable for large parameters.
From experience we recommend that the respective first parameters are chosen smaller than 5 and the second parameters smaller than 6.
In general, the user should be careful with extreme parameter choices.

The following table shows the parameter ranges of bivariate copula families with parameters \code{par} and \code{par2}:
\tabular{lll}{
Copula family \tab \code{par} \tab \code{par2} \cr
Gaussian \tab \eqn{(-1,1)} \tab - \cr
Student t \tab \eqn{(-1,1)} \tab \eqn{(1,\infty)} \cr
(Survival) Clayton \tab \eqn{(0,\infty)} \tab - \cr
(Survival) Gumbel \tab \eqn{[1,\infty)} \tab - \cr
Frank \tab \eqn{R\backslash\{0\}} \tab - \cr
(Survival) Joe \tab \eqn{(1,\infty)} \tab - \cr
Rotated Clayton (90 and 270 degrees) \tab \eqn{(-\infty,0)} \tab - \cr
Rotated Gumbel (90 and 270 degrees)\tab \eqn{(-\infty,-1]} \tab - \cr
Rotated Joe (90 and 270 degrees) \tab \eqn{(-\infty,-1)} \tab - \cr
Clayton-Gumbel (BB1) \tab \eqn{(0,\infty)} \tab \eqn{[1,\infty)} \cr
Joe-Clayton (BB7) \tab \eqn{[1,\infty)} \tab \eqn{(0,\infty)} \cr
} 
}

\section{C- and D-vine copula models}{
When specifying C- and D-vine copula models, one has to select an order of the variables.
For a D-vine the order of the variables in the first tree has to be chosen and for a C-vine the root nodes for each tree need to be determined.
Functions for inference of C- and D-vine copula models in this package assume that the order of the variables in the data set under investigation
exactly corresponds to this C- or D-vine order.
E.g., in a C-vine the first column of a data set is the first root node, the second column the second root node, etc.
According to this order arguments have to be provided to functions for C- and D-vine copula inference.
After choosing the type of the vine model,
the copula families (\code{family}) and parameters (\code{par} and \code{par2}) have to be specified as vectors of length \eqn{d(d-1)/2},
where \eqn{d} is the number of variables.
In a C-vine, the entries of this vector correspond to the following pairs and associated pair-copula terms\cr\cr
\eqn{(1,2), (1,3), (1,4), ..., (1,d),}
\eqn{(2,3|1), (2,4|1), ..., (2,d|1),}
\eqn{(3,4|1,2), (3,5|1,2), ..., (3,d|1,2),}
\eqn{...,}
\eqn{(d-1,d|1,...,d-2).}\cr\cr
Similarly, the pairs of a D-vine are denoted in the following order:\cr\cr
\eqn{(1,2), (2,3), (3,4), ..., (d-1,d),}
\eqn{(1,3|2), (2,4|3), ..., (d-2,d|d-1),}
\eqn{(1,4|2,3), (2,5|3,4), ..., (d-3,d|d-2,d-1),}
\eqn{...,}
\eqn{(1,d|2,...,d-1).} 
}

\section{Empirical Kendall's tau}{
We like to thank Shing (Eric) Fu, Feng Zhu, Guang (Jack) Yang, and Harry Joe for providing their implementation
of the method by Knight (1966) for efficiently computing the empirical Kendall's tau.
}

\details{
\tabular{ll}{
Package: \tab CDVine\cr
Type: \tab Package\cr
Version: \tab 1.0-1\cr
Date: \tab 2011-05-18\cr
License: \tab GPL (>=2)\cr
Depends: \tab R (>= 2.8.1), MASS, mvtnorm, igraph \cr
LazyLoad: \tab yes
}
}

\author{Ulf Schepsmeier, Eike Christian Brechmann <CDVine@ma.tum.de>
}

\references{ 
Aas, K., C. Czado, A. Frigessi, and H. Bakken (2009).
Pair-copula constructions of multiple dependence.
Insurance: Mathematics and Economics 44 (2), 182-198.

Bedford, T. and R. M. Cooke (2001).
Probability density decomposition for conditionally dependent random variables modeled by vines.
Annals of Mathematics and Artificial intelligence 32, 245-268.

Bedford, T. and R. M. Cooke (2002).
Vines - a new graphical model for dependent random variables.
Annals of Statistics 30, 1031-1068.

Czado, C., U. Schepsmeier, and A. Min (2011).
Maximum likelihood estimation of mixed C-vines with application to exchange rates.
Submitted for publication.
\url{http://www-m4.ma.tum.de/Papers/Schepsmeier/Paper.pdf}.

Joe, H. (1996).
Families of m-variate distributions with given margins and m(m-1)/2 bivariate dependence parameters.
In L. Rueschendorf, B. Schweizer, and M. D. Taylor (Eds.), Distributions with fixed marginals and related topics, pp. 120-141.
Hayward: Institute of Mathematical Statistics.

Joe, H. (1997).
Multivariate Models and Dependence Concepts.
Chapman and Hall, London.

Knight, W. R. (1966).
A computer method for calculating Kendall's tau with ungrouped data.
Journal of the American Statistical Association 61 (314), 436-439.

Kurowicka, D. and R. M. Cooke (2006).
Uncertainty Analysis with High Dimensional Dependence Modelling.
Chichester: John Wiley.

Kurowicka, D. and H. Joe (Eds.) (2011).
DEPENDENCE MODELING: Vine Copula Handbook.
Singapore: World Scientific Publishing Co.
}
