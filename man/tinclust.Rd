% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tinclust.R
\name{tinclust}
\alias{tinclust}
\title{Model based clustering with MTIN distribution}
\usage{
tinclust(
  x,
  G = 1,
  max_iter = 100,
  tol = 10^-1,
  init_method = c("mclust", "kmedoids", "kmeans", "heirarchical"),
  ...,
  verbose = FALSE,
  plot_likelihood = TRUE
)
}
\arguments{
\item{x}{A data matrix.}

\item{G}{An integer representing the number of clusters.}

\item{max_iter}{An integer representing the maximum number of iterations.}

\item{tol}{Float representing tolerance for Aitken's acceleration criterion.}

\item{init_method}{Model used for initialization: \code{'mclust'}, \code{"kmedoids"}, \code{"kmeans"}, \code{"heirarchical"}}

\item{...}{Additional arguments passed to \code{dmtin}}

\item{verbose}{Print log-likelihood at every iteration.}

\item{plot_likelihood}{Plot log-likelihoods after convergence.}
}
\value{
A list with the following elements:
\item{Mu}{List of cluster means.}
\item{Sigma}{List of cluster scale matrices.}
\item{Theta}{List of cluster inflation parameters.}
\item{Pi}{Mixing proportions.}
\item{Z}{Expected value of Z.}
\item{W_y}{Expected value of W.}
\item{cluster}{Membership labels.}
\item{L}{List of log-likelihoods at each iteration.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{KIC}{Kullback information criterion.}??
}
\description{
Model based clustering with MTIN distribution
}
\examples{
\donttest{
data <- palmerpenguins::penguins
keep <- !rowSums(is.na(data[,3:6])) == 4
data <- data[keep,]

res <- tinclust(x = data[,3:6], G = 3)
plot(res$X[,3], res$X[,1], col = res$cluster)
}



}
