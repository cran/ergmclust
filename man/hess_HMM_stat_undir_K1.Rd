\name{hess_HMM_stat_undir_K1}
\alias{hess_HMM_stat_undir_K1}
\title{
C++ function for updating the Hessian matrix in a static binary undirected network.
}
\description{
The C++ function hess_HMM_stat_undir() updates the Hessian matrix when the number of cluster is 1 in a static binary undirected network.
}
\usage{
hess_HMM_stat_undir_K1(theta, N)
}
\arguments{
  \item{theta}{
  Vector of current \code{theta} (network canonical parameters) values.
}
  \item{N}{
  Number of nodes.
}
}
\value{
Returns the Hessian matrix of ELBO with respect to theta parameters for Newton-Raphson update in M-step.
}
\keyword{internal}