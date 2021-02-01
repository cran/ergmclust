\name{gamma_update_HMM_stat_dir}
\alias{gamma_update_HMM_stat_dir}
\title{
C++ function for updating the variable \code{gamma} (variational parameters) in a static binary directed network.
}
\description{
The C++ function gamma_update_HMM_stat_dir() updates the value of variable \code{gamma} (variational parameters) in a static binary directed network.
}
\usage{
gamma_update_HMM_stat_dir(gamma, pi, theta, network, N, K)
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
  \item{network}{
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
Returns all quadratic and linear term coefficients for the quadratic programs corresponding to all nodes.
}
\keyword{internal}