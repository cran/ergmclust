\name{tradenet}
\alias{tradenet}
\docType{data}
\title{
International Trade Network Data in 1991.
}
\description{
The undirected network on all trade relations internationally among 58 countries. We define the edges as \eqn{y_{ij}=1}, if there is a bilateral trade between country \eqn{i} and \eqn{j}, and \eqn{y_{ij}=0} otherwise.
}
\usage{data(tradenet)}
\format{
The format is a 58 \eqn{\times} 58 network adjacency matrix.
}
\source{
\url{https://projecteuclid.org/euclid.aoas/1310562208#supplemental}
}
\references{
Westveld, A. H. and Hoff, P. D. (2011) \emph{A mixed effects model for longitudinal relational and network data, with applications to international trade and conflict}, The Annals of Applied Statistics 5(2A), 843–872 

\url{https://projecteuclid.org/euclid.aoas/1310562208}
}
\examples{data(tradenet)}