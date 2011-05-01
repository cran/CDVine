\name{BiCopName}       
\alias{BiCopName}

\title{Bivariate copula family names}

\description{
This function transforms the bivariate copula family number into its character expression and vice versa.
}

\usage{
BiCopName(family, short=TRUE)
}

\arguments{
\item{family}{Bivariate copula family, either its number or its character expression (see table below).
\tabular{rll}{
No. \tab Short name \tab Long name \cr
\code{0} \tab \code{"I"} \tab \code{"Independent"} \cr
\code{1} \tab \code{"N"} \tab \code{"Gaussian"} \cr
\code{2} \tab \code{"t"} \tab \code{"t"} \cr
\code{3} \tab \code{"C"} \tab \code{"Clayton"} \cr
\code{4} \tab \code{"G"} \tab \code{"Gumbel"} \cr
\code{5} \tab \code{"F"} \tab \code{"Frank"} \cr
\code{6} \tab \code{"J"} \tab \code{"Joe"} \cr
\code{7} \tab \code{"BB1"} \tab \code{"Clayton-Gumbel"} \cr
\code{9} \tab \code{"BB7"} \tab \code{"Joe-Clayton"} \cr
\code{13} \tab \code{"SC"} \tab \code{"Survival Clayton"} \cr
\code{14} \tab \code{"SG"} \tab \code{"Survival Gumbel"} \cr
\code{16} \tab \code{"SJ"} \tab \code{"Survival Joe"} \cr
\code{23} \tab \code{"C90"} \tab \code{"Rotated Clayton 90 degrees"} \cr
\code{24} \tab \code{"G90"} \tab \code{"Rotated Gumbel 90 degrees"} \cr
\code{26} \tab \code{"J90"} \tab \code{"Rotated Joe 90 degrees"} \cr
\code{33} \tab \code{"C270"} \tab \code{"Rotated Clayton 270 degrees"} \cr
\code{34} \tab \code{"G270"} \tab \code{"Rotated Gumbel 270 degrees"} \cr
\code{36} \tab \code{"J270"} \tab \code{"Rotated Joe 270 degrees"}
}
}
\item{short}{Logical; if the number of a bivariate copula family is used and \code{short = TRUE} (default),
  a short version of the corresponding character expression is returned, otherwise the long version.}
}

\value{
The transformed bivariate copula family (see table above).
}

\author{Ulf Schepsmeier}

\seealso{\code{\link{CDVineTreePlot}}}

\examples{
# family as number
family = 1
BiCopName(family,short=TRUE)	 # short version
BiCopName(family,short=FALSE)	# long version

# family as character expression (short version)
family = "C"
BiCopName(family)	# as number

# long version
family = "Clayton"
BiCopName(family)	# as number
}
