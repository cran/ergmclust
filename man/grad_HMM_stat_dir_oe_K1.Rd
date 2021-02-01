\name{grad_HMM_stat_dir_oe_K1}
\alias{grad_HMM_stat_dir_oe_K1}
\title{
C++ function for updating the gradient function in a static binary directed network.
}
\description{
The C++ function grad_HMM_stat_dir_oe_K1() updates the gradient function when the number of cluster is 1 in a static binary directed network. (The above description may need further improvement.)
}
\usage{
grad_HMM_stat_dir_oe_K1(theta, network, N)
}
\arguments{
  \item{theta}{
  Vector of current \code{theta} (network canonical parameters) values.
}
  \item{network}{
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