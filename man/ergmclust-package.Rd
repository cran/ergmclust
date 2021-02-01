\name{ergmclust-package}
\alias{ergmclust-package}
\docType{package}
\title{
ERGM-based network clustering
}
\description{
Clustering and estimation of parameters in ERGMs for static undirected and directed networks with inference based on VEM algorithm.
}
\details{
The ergmclust package is an R implementation that serves as an estimation framework for static binary networks, in both undirected and directed cases. Its main functions include ergmclust for clustering and parameter estimation, ergmclust.ICL for model selection, and ergmclust.plot for visualizing the clustered network. The package is based on VEM algorithm (Vu et. al., 2013) and works well with both simulated and real-world data.
}
\author{
Authors: Amal Agarwal [aut, cre], Kevin Lee [aut], Lingzhou Xue [aut, cre], Anna Yinqi Zhang [cre]

Maintainer: Amal Agarwal <amalag.19@gmail.com>
}
\references{
Agarwal, A. and Xue, L. (2019) \emph{Model-Based Clustering of Nonparametric Weighted Networks With Application to Water Pollution Analysis}, Technometrics, to appear

\url{https://amstat.tandfonline.com/doi/abs/10.1080/00401706.2019.1623076}

Biernacki, C., Celeux, G., and Govaert, G. (2000) \emph{Assessing a mixture model for clustering with the integrated completed likelihood}, IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 22(7), 719-725

\url{https://ieeexplore.ieee.org/document/865189}

Blei, D. M. , Kucukelbir, A., and McAuliffe, J. D. (2017), \emph{Variational Inference: A Review for Statisticians}, Journal of the American Statistical Association, Vol. 112(518), 859-877

\url{https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1285773}

Daudin, J. J., Picard, F., and Robin, S. (2008) \emph{A Mixture Model for Random Graphs}, Statistics and Computing, Vol. 18(2), 173–183

\url{https://link.springer.com/article/10.1007/s11222-007-9046-7}

Lee, K. H., Xue, L, and Hunter, D. R. (2017) \emph{Model-Based Clustering of Time-Evolving Networks through Temporal Exponential-Family Random Graph Models}, Journal of Multivariate Analysis, to appear

\url{https://arxiv.org/abs/1712.07325}

Vu D. Q., Hunter, D. R., and Schweinberger, M. (2013) \emph{Model-based Clustering of Large Networks}, The Annals of Applied Statistics, Vol. 7(2), 1010-1039

\url{https://projecteuclid.org/euclid.aoas/1372338477}
}