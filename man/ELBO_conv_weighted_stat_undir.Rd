\name{ELBO_conv_weighted_stat_undir}
\alias{ELBO_conv_weighted_stat_undir}
\title{
C++ function for updating the ELBO convergence function in a weighted static binary undirected network.
}
\description{
The C++ function ELBO_conv_weighted_stat_undir() updates the ELBO convergence function when the number of clusters is greater than or equal to 2 in a weighted static binary undirected network.
}
\usage{
ELBO_conv_weighted_stat_undir(gamma, pi, theta, block_dens_mat, adjmat, N, K)
}
\arguments{
  \item{gamma}{
  Vector of current \code{gamma} (variational parameters) values.
  }
  \item{pi}{
  Vector of current \code{pi} (mixture proportions) values.
  }
  \item{theta}{
  Vector of current \code{theta} (network canonical parameters) values.
  }
  \item{block_dens_mat}{
  Array of a network block density matrices.
  }
  \item{adjmat}{
  Array of a network adjacency matrices.
  }
  \item{N}{
  Number of nodes.
  }
  \item{K}{
  Number of clusters.
  }
}
\value{
Returns the value of ELBO for calculating relative error from the previous iterate that is subsequently used in the stopping criteria of the algorithm.
}
\keyword{internal}