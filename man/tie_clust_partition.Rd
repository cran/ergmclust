\name{tie_clust_partition}
\alias{tie_clust_partition}
\title{
C++ function to acquire the y_ij for all pairs of clusters.
}
\description{
The C++ function retrieves y_ij for all pairs of clusters in a weighted static binary undirected network.
}
\usage{
tie_clust_partition(clust_est, adjmat, wtmat, N, K)
}
\arguments{
  \item{clust_est}{
  A vector of estimated cluster membership.
  }
  \item{adjmat}{
  Array of a network adjacency matrices.
  }
  \item{wtmat}{
  Array of a network weight matrices.
  }
  \item{N}{
  Number of nodes.
  }
  \item{K}{
  Number of clusters.
  }
}
\value{
Returns y_ij (0, 1) for all pairs of clusters
}
\keyword{internal}