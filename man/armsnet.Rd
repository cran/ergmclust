\name{armsnet}
\alias{armsnet}
\docType{data}
\title{
Arms Trade Network Data in 2003.
}
\description{
The directed network on all transfers of major conventional weapons internationally. We define the edges as \eqn{y_{ij}=1}, if the volume of international transfers of arms, measured by Trend Indicator Value (TIV) from country i to country j exceeds 1 million dollars, and \eqn{y_{ij}=0} otherwise.
}
\usage{data(armsnet)}
\format{
The format is a 69 \eqn{\times} 69 network adjacency matrix.
}
\source{
\url{https://www.sipri.org/databases/armstransfers}
}
\references{
Akerman, A., & Seim, A. L. (2014) \emph{The global arms trade network 1950–2007}, Journal of Comparative Economics, Vol. 42(3), 535-551

\url{https://www.sciencedirect.com/journal/journal-of-comparative-economics/vol/42/issue/3}
}
\examples{data(armsnet)}