\name{worldindices}      
\docType{data}
\alias{worldindices}

\title{Major World Indices}

\description{
This data set contains transformed standardized residuals of daily log returns of major world stock indices in 2009 and 2010. 
The considered indices are the leading stock exchanges of the six largest economies in the world:
the US American S&P 500, the Japanese Nikkei 225, the Chinese SSE Composite Index, the German DAX, the French CAC 40 and the British FTSE 100 Index.
Each time series is filtered using an ARMA(1,1)-GARCH(1,1) model with Student t innovations.
}

\format{
A data frame with 396 observations on 6 variables. Column names correspond to ticker symbols of the indices.
}

\source{Yahoo! Finance}

\examples{
# load the data set
data(worldindices)

# compute the empirical Kendall's tau matrix
cor(worldindices,method="kendall")
}
