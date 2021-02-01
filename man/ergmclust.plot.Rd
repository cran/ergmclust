\name{ergmclust.plot}
\alias{ergmclust.plot}
\title{
Visualization For Model-Based Clustering of Large Networks.
}
\description{
Visualization of the network data with the clusters node colors representing different clusters in the Exponential-Family Random Graph Models (ERGMs) clustered network.
}
\usage{
ergmclust.plot(adjmat, K, directed = FALSE, thresh = 1e-06, 
iter.max = 200, coef.init = NULL, node.labels = NULL)
}
\arguments{
  \item{adjmat}{
  An object of class matrix of dimension (N x N) containing the adjacency matrix, where N is the number of nodes in the network.
}
  \item{K}{
  Number of clusters in the mixed membership Exponential-Family Random Graph Models (ERGMs).
}
  \item{directed}{
  If \code{TRUE}, the network is supposed to be directed (and therefore adjmat is must be asymmetric in general). By default, this is set as \code{FALSE}.
}
  \item{thresh}{
  Optional user-supplied convergence threshold for relative error in the objective in Variational Expectation-Maximization (VEM) algorithm. The default value is set as 1e-06. 
}
  \item{iter.max}{
  The maximum number of iterations after which the algorithm is terminated. The default value is set as 200.
}
  \item{coef.init}{
  Optional user-supplied network canonical parameter vector (K-dimensional); default is \code{NULL}, and ergmclust chooses a random perturbation around zero vector.
}
  \item{node.labels}{
  Optional user-supplied network node names character vector (N-dimensional); default is \code{NULL}.
}
}
\details{
ergmclust.plot provides the visualization tool for network data clustered through mixed membership Exponential-Family Random Graph Models (ERGMs). The optional argument \code{node.labels} could help track the cluster membership of specific nodes.
}
\value{
Returns a plot of network object with colored nodes corresponding to \code{K} clusters.
}
\references{
Vu D. Q., Hunter, D. R., and Schweinberger, M. (2013) \emph{Model-based Clustering of Large Networks}, The Annals of Applied Statistics, Vol. 7(2), 1010-1039

\url{https://projecteuclid.org/euclid.aoas/1372338477}
}
\author{
Authors: Amal Agarwal [aut, cre], Kevin Lee [aut], Lingzhou Xue [aut, cre], Anna Yinqi Zhang [cre]

Maintainer: Amal Agarwal <amalag.19@gmail.com>
}
\examples{
## undirected network:
data(tradenet)
## Plotting clustered network
ergmclust.plot(adjmat = tradenet, K = 2, directed = FALSE, 
thresh = 1e-06)

## directed network:
data(armsnet)
## Plotting clustered network
ergmclust.plot(adjmat = armsnet, K = 2, directed = TRUE, 
thresh = 1e-06)
}