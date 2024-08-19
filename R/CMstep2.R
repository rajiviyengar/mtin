#' Conditional weighted ML estimation of theta given mu and Sigma.
#'
#' @param X Data matrix.
#' @param weights Vector of weights.
#' @param mu Mean vector.
#' @param Sigma Scalar matrix.
#' @param naMat Matrix of missing value indicators for data matrix.
#' @param obsInd Indices of rows containing complete observations.
#' @param misInd Indices of rows containing missing observations.
#' @param ... Additional arguments passed to \code{dmtin}
#'
#' @return
#' \item{theta}{Estimated inflation parameter.}
#' @examples
#' d <- 2
#' mu <- rep(0,d)
#' Sigma <- diag(d)
#' set.seed(5)
#' X <- rmtin(n=100,mu,Sigma,theta=0.4)$X
#' CMstep2(X=X,mu=mu,Sigma=Sigma)
#'
#' @export
CMstep2 <- function(X,weights=NULL,mu,Sigma, naMat = NULL, obsInd = NULL, misInd = NULL, ...){


  if(is.vector(X))
    X <- matrix(X,ncol=1)
  if(is.data.frame(X))
    X <- as.matrix(X)
  if(any(is.na(X)))
    stop('No NAs allowed.')

  n <- nrow(X)

  if(is.null(weights))
    weights <- rep(1,n)

  # objective function

  if (is.null(naMat)) {
       f <- function(par,weights,X,mu,Sigma,formula){

            # ----- #
            # theta #
            # ----- #

            theta <- par

            # -------------------- #
            # pseudo log-likelihood #
            # -------------------- #

            pll <- sum(weights*log(dmtin(x=X,mu=mu,Sigma=Sigma,theta=theta,...)))

            return(pll)
       }
  }
  else {
       f <- function(par,weights,X,mu,Sigma,...){

            # ----- #
            # theta #
            # ----- #

            theta <- par

            # -------------------- #
            # pseudo loglikelihood #
            # -------------------- #

            pll <- sum(weights[obsInd] * log(dmtin(x = X[obsInd,], mu = mu, Sigma = Sigma, theta = theta)))

            for (i in misInd) {
                 o <- !naMat[i,]

                 pll <- pll + weights[i] * log(dmtin(x = X[i,o], mu = mu[o], Sigma = Sigma[o,o], theta = theta, ...))
            }

            return(pll)

       }
  }

  suppressWarnings({
       res <- stats::optimize(f=f, interval=c(0,1), weights=weights, X=X, mu=mu, Sigma=Sigma, maximum=TRUE, ...)
  })



  theta <- res$maximum

  return(theta)

}
