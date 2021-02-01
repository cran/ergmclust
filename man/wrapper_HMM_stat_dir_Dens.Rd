\name{wrapper_HMM_stat_dir_Dens}
\alias{wrapper_HMM_stat_dir_Dens}
\title{
An estimation framework for static binary directed networks.
}
\description{
The inner wrapper function wrapper_HMM_stat_dir_Dens() is an estimation framework for binary static directed network based on variational Expectation-Maximization algorithm.
}
\usage{
wrapper_HMM_stat_dir_Dens(adjmat, K, thresh, iter.max, theta_init)
}
\arguments{
  \item{adjmat}{
  Array of the network adjacency matrices.
}
  \item{K}{
  Number of clusers for which the estimation must be performed.
}
  \item{thresh}{
  Convergence threshold for estimation, the default being 1e-6.
}
  \item{iter.max}{
  The maximum number of iterations after which the algorithm is terminated. The default value is set as 200.
  }
  \item{theta_init}{
  The initial value of the density, stability, transitivity parameters for which clustering is desired, the default being 1e-1.
}
}
\value{
Returns a list of \code{ergmclust} object. Each object of class \code{ergmclust} is a \code{list} with the following components:
  \item{Converged Parameters}{The first elements of the list in sequential order are gamma, alpha, pi, tau, theta for dynamic networks. For static networks, there is no pi or tau.}
  \item{Estimated Cluster IDs}{The second element of the list is a vector (matrix) of estimated cluster memberships for a static (dynamic) case. This is absent for all K=1 cases.}
  \item{ICL Values}{The third element of the list is the Integrated Classification Likelihood value that can be used for model selection, i.e. selecting the appropriate number of clusters.}
}
\keyword{internal}