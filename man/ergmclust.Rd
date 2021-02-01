\name{ergmclust}
\alias{ergmclust}
\title{
Model-Based Clustering of Large Networks Through ERGMs.
}
\description{
Model-based clustering and cluster-specific parameter estimation through the mixed membership Exponential-Family Random Graph Models (ERGMs) using Variational Expectation-Maximization  algorithm.
}
\usage{
ergmclust(adjmat, K, directed = FALSE, thresh = 1e-06, iter.max = 200, coef.init = NULL)
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
  Optional user-supplied network canonical parameter vector (K-dimensional). Default is \code{NULL} and ergmclust chooses a random perturbation around zero vector.
}
}
\details{
ergmclust is an R implementation for the model-based clustering through the mixed membership Exponential-Family Random Graph Models (ERGMs) with undirected and directed network data. It uses the Variational Expectation-Maximization algorithm to solve the approximate maximum likelihood estimation. 
}
\value{
Returns a list of \code{ergmclust} object. Each object of class \code{ergmclust} is a \code{list} with the following components:
  \item{coefficients}{An object of class vector of size (K x 1) containing the canonical network parameters in Exponential-Family Random Graph Models (ERGMs).}
  \item{probability}{An object of class matrix of size (N x K) containing the mixed membership probabilities of the model for N nodes distributed in K clusters.}
  \item{clust.labels}{An object of class vector of size (N x 1) containing the cluster membership labels in \{1, ..., K\} for N nodes.}
  \item{ICL}{Integrated Classification Likelihood (ICL) score calculated from completed data log-likelihood and penalty terms.}
}
\references{
Agarwal, A. and Xue, L. (2019) \emph{Model-Based Clustering of Nonparametric Weighted Networks With Application to Water Pollution Analysis}, Technometrics, to appear

\url{https://amstat.tandfonline.com/doi/abs/10.1080/00401706.2019.1623076}

Blei, D. M. , Kucukelbir, A., and McAuliffe, J. D. (2017), \emph{Variational Inference: A Review for Statisticians}, Journal of the American Statistical Association, Vol. 112(518), 859-877

\url{https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1285773}

Lee, K. H., Xue, L, and Hunter, D. R. (2017) \emph{Model-Based Clustering of Time-Evolving Networks through Temporal Exponential-Family Random Graph Models}, Journal of Multivariate Analysis, to appear

\url{https://arxiv.org/abs/1712.07325}

Vu D. Q., Hunter, D. R., and Schweinberger, M. (2013) \emph{Model-based Clustering of Large Networks}, The Annals of Applied Statistics, Vol. 7(2), 1010-1039

\url{https://projecteuclid.org/euclid.aoas/1372338477}
}
\author{
Authors: Amal Agarwal [aut, cre], David R. Hunter [aut], Kevin Lee [aut], Lingzhou Xue [aut, cre], Anna Yinqi Zhang [cre]

Maintainer:  Amal Agarwal <amalag.19@gmail.com>
}
\examples{
## undirected network:
data(tradenet)
## clustering and estimation for K = 2 groups
ergmclust(adjmat = tradenet, K = 2, directed = FALSE, 
thresh = 1e-06, iter.max = 120, coef.init = NULL)

## directed network:
data(armsnet)
## clustering and estimation for K = 2 groups
ergmclust(adjmat = armsnet, K = 2, directed = TRUE, 
thresh = 1e-06, iter.max = 120, coef.init = NULL)
}