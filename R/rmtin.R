#' Generate random numbers from the MTIN Distribution.
#'
#' @param n Number of points to be generated.
#' @param mu Mean vector.
#' @param Sigma Scale matrix.
#' @param theta Inflation parameter.
#' @param norm.dens Method used to calculate the density of the normal distribution: \code{"dmnorm"} (see \code{\link[mnormt]{dmnorm}}) or \code{"dmvnorm"} (see \code{\link[mvtnorm]{Mvnorm}}).
#'
#' @return A list with the following elements:
#' \item{X}{Data matrix.}
#' \item{w}{Vector of weights.}
#' @examples
#' d <- 3
#' rmtin(10,mu=rep(0,d),Sigma=diag(d),theta=0.4)
#'
#' @export
rmtin <- function(n, mu = rep(0,d), Sigma, theta = 0.01, norm.dens = c("dmnorm", "dmvnorm", "Rfast")){
     norm.dens <- match.arg(norm.dens)

  if(missing(Sigma))
    stop("Sigma is missing")
  if(theta <= 0 | theta >= 1)
    stop("theta must be in the interval (0,1)")

  if(is.matrix(Sigma))
    d <- ncol(Sigma)
  if(!is.matrix(Sigma))
    d <- 1

  X <- matrix(0,n,d)
  w <- stats::runif(n=n,min=1-theta,1)

  if(norm.dens == 'Rfast') {
       for (i in 1:n) {
            X[i,] <- Rfast::rmvnorm(n = 1, mu = mu, sigma = Sigma)
       }
  }
  if(norm.dens == "dmnorm")
    for(i in 1:n)
      X[i,] <- mnormt::rmnorm(n = 1, mean = mu, varcov = Sigma/w[i])
  if(norm.dens == "dmvnorm")
    for(i in 1:n)
      X[i,] <- mvtnorm::rmvnorm(n = 1, mean = mu, sigma = Sigma/w[i])

  return(list(X=X,w=w))

}
