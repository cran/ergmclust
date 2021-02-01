\name{ergmclust.ICL}
\alias{ergmclust.ICL}
\title{
Model Selection Based On Integrated Classification Likelihood.
}
\description{
Model-based clustering and cluster-specific parameter estimation through the mixed membership Exponential-Family Random Graph Models (ERGMs) for the different number of clusters. Model selection is based on maximum value of Integrated Classification Likelihood (ICL).
}
\usage{
ergmclust.ICL(adjmat, Kmax = 5, directed = FALSE,
    thresh = 1e-06, iter.max = 200, coef.init = NULL)
}
\arguments{
  \item{adjmat}{
  An object of class matrix of dimension (N x N) containing the adjacency matrix, where N is the number of nodes in the network.
}
  \item{Kmax}{
  Maximum number of clusters.
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
  Optional user-supplied network canonical parameter vector (K-dimensional). Default is \code{NULL}, and ergmclust chooses a random perturbation around zero vector.
}
}
\details{
ergmclust.ICL is an R implementation for the model selection for an appropriate number of clusters in the mixed membership Exponential-Family Random Graph Models (ERGMs). The Integrated Classification Likelihood (ICL) was proposed by Biernacki et al. (2000) and Daudin, et. al. (2008) to assess the model-based clustering. 
}
\value{
Returns a list of \code{ergmclust} object. Each object of class \code{ergmclust} is a \code{list} with the following components:
  \item{Kselect}{Optimum number of clusters chosen after model selection through Integrated Classification Likelihood (ICL).}
  \item{coefficients}{An object of class vector of size (Kselect x 1) containing the canonical network parameters of the model.}
  \item{probability}{An object of class matrix of size (N x Kselect) containing the mixed membership probabilities of the model for N nodes distributed in Kselect clusters.}
  \item{clust.labels}{An object of class vector of size (N x 1) containing the cluster membership labels in \{1, ..., Kselect\} for N nodes.}
  \item{ICL}{Integrated Classification Likelihood (ICL) score calculated from completed data log-likelihood and penalty terms.}
}
\references{
Biernacki, C., Celeux, G., and Govaert, G. (2000) \emph{Assessing a mixture model for clustering with the integrated completed likelihood}, IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 22(7), 719-725

\url{https://ieeexplore.ieee.org/document/865189}

Daudin, J. J., Picard, F., and Robin, S. (2008) \emph{A Mixture Model for Random Graphs}, Statistics and Computing, Vol. 18(2), 173–183

\url{https://link.springer.com/article/10.1007/s11222-007-9046-7}
}
\author{
Authors: Amal Agarwal [aut, cre], Kevin Lee [aut], Lingzhou Xue [aut, cre], Anna Yinqi Zhang [cre]

Maintainer: Lingzhou Xue <lzxue@psu.edu>
}
\examples{
## undirected network:
data(tradenet)
## Model selection for Kmax = 3
ergmclust.ICL(adjmat = tradenet, Kmax = 3, directed = FALSE, 
thresh = 1e-06, iter.max = 120, coef.init = NULL)

## directed network:
data(armsnet)
## Model selection for Kmax = 3
ergmclust.ICL(adjmat = armsnet, Kmax = 3, directed = TRUE,
thresh = 1e-06, iter.max = 60, coef.init = NULL)
}