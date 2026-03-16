\name{grad_theta_weighted_stat_undir_K1}
\alias{grad_theta_weighted_stat_undir_K1}
\title{
C++ function for updating the gradient function in a weighted static binary undirected network.
}
\description{
The C++ function grad_HMM_stat_undir() updates the gradient function when the number of clusters is 1 in a weighted static binary undirected network.
}
\usage{
grad_theta_weighted_stat_undir_K1(theta, adjmat, N)
}
\arguments{
  \item{theta}{
  Vector of current \code{theta} (network canonical parameters) values.
  }
  \item{adjmat}{
  Array of a network adjacency matrices.
  }
  \item{N}{
  Number of nodes.
  }
}
\value{
Returns the gradient vector of ELBO with respect to theta parameters for Newton-Raphson update in the M-step.
}
\keyword{internal}