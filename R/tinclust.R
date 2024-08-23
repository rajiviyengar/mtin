#' Model based clustering with MTIN distribution
#'
#' @param x A data matrix.
#' @param G An integer representing the number of clusters.
#' @param max_iter An integer representing the maximum number of iterations.
#' @param tol Float representing tolerance for Aitken's acceleration criterion.
#' @param init_method Model used for initialization: \code{'mclust'}, \code{"kmedoids"}, \code{"kmeans"}, \code{"heirarchical"}
#' @param verbose Print log-likelihood at every iteration.
#' @param plot_likelihood Plot log-likelihoods after convergence.
#' @param ... Additional arguments passed to \code{dmtin}
#' @return A list with the following elements:
#' \item{Mu}{List of cluster means.}
#' \item{Sigma}{List of cluster scale matrices.}
#' \item{Theta}{List of cluster inflation parameters.}
#' \item{Pi}{Mixing proportions.}
#' \item{Z}{Expected value of Z.}
#' \item{W_y}{Expected value of W.}
#' \item{cluster}{Membership labels.}
#' \item{L}{List of log-likelihoods at each iteration.}
#' \item{AIC}{Akaike information criterion.}
#' \item{BIC}{Bayesian information criterion.}
#' \item{KIC}{Kullback information criterion.}??
#'
#' @examples
#' \donttest{
#' data <- palmerpenguins::penguins
#' keep <- !rowSums(is.na(data[,3:6])) == 4
#' data <- data[keep,]
#'
#' res <- tinclust(x = data[,3:6], G = 3)
#' plot(res$X[,3], res$X[,1], col = res$cluster)
#' }
#'
#'
#'
#' @export
tinclust <- function(x, G = 1, max_iter = 100, tol = 10^-1, init_method = c("mclust", "kmedoids", "kmeans", "heirarchical"), ..., verbose = FALSE, plot_likelihood = TRUE) {

     {
          ####################
          # Input Formatting #
          ####################
          x <- as.matrix(x)
          M <- is.na(x)
          n <- nrow(x)
          p <- ncol(x)
          mis <- which(as.logical(rowSums(M)))
          obs <- which(!as.logical(rowSums(M)))
          nMisPatterns <- integer(0)
          if (p == 1 && any(is.na(x))) stop("NAs not allowed with univariate data.")
          if (any(rowSums(M) == p)) stop("Remove observations that contain only NAs.")

          if (any(M)) {
               misPatterns <- matrix(unique(matrix(M[mis,], ncol = p)), ncol = p)
               nMisPatterns <- 1:nrow(misPatterns)
               whichMis <- rep(list(0), length(nMisPatterns))

               for (j in nMisPatterns) {
                    whichMis[[j]] <- which(apply(M, 1, function (row) all(row == misPatterns[j,])))
               }

          }

          ####################
          # Model parameters #
          ####################
          Mus <- rep(list(matrix(0, nrow = 1, ncol = p)), G)
          Sigmas <- rep(list(matrix(0, nrow = p, ncol = p)), G)
          Thetas <- rep(0.1, G)
          Pis <- rep(0, G)
          Z <- matrix(0, nrow = n, ncol = G)
          W <- matrix(0, nrow = n, ncol = G)
          li <- matrix(0, nrow = n, ncol = G)
          iter <- 0

          #########################
          # Store Log-Likelihoods #
          #########################
          L <- rep(0, max_iter)
          diffs <- rep(0, max_iter)

          #########################
          # Initialize Parameters #
          #########################
          Z <- init_clusters(x, G, obsInd = obs, init_method = init_method)

          for (g in 1:G) {
               Mus[[g]] <- matrix(Rfast::colmeans(matrix(x[Z[,g] == 1,], ncol = p)), ncol = p)
               Sigmas[[g]] <- Rfast::cova(matrix(x[Z[,g] == 1,], ncol = p))
               Pis[g] <- sum(Z[,g]) / length(obs)
               li[obs,g] <- Pis[g] * dmtin(x[obs,], Mus[[g]], Sigmas[[g]], Thetas[g], ...)
          }

          for (j in nMisPatterns) {
               o <- !misPatterns[j,]
               pObs <- sum(o)

               for (g in 1:G) {
                    li[whichMis[[j]],g] <- Pis[g] * dmtin(x[whichMis[[j]],o], Mus[[g]][o], Sigmas[[g]][o,o], Thetas[g], ...)
               }
          }

          ##################
          # ECME Algorithm #
          ##################
          for (k in 1:max_iter) {
               #############
               # E-Step: W #
               #############
               for (g in 1:G) {
                    delta <- Rfast::mahala(matrix(x[obs,], ncol = p), Mus[[g]], Sigmas[[g]])
                    num <- Rfast::Pmax(2 * ((zipfR::Igamma((p/2)+2, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((p/2)+2, delta/2, lower = FALSE))), rep(10^(-322), length(obs)))
                    den <- Rfast::Pmax(delta * ((zipfR::Igamma((p/2)+1, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((p/2)+1, delta/2, lower = FALSE))), num)
                    W[obs,g] <- num/den
               }

               for (j in nMisPatterns) {
                    o <- !misPatterns[j,]
                    pObs <- sum(o)

                    for (g in 1:G) {
                         delta <- Rfast::mahala(matrix(x[whichMis[[j]],o], ncol = pObs), Mus[[g]][o], Sigmas[[g]][o,o])
                         num <- Rfast::Pmax(2 * ((zipfR::Igamma((p/2)+2, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((p/2)+2, delta/2, lower = FALSE))), rep(10^(-322), length(whichMis[[j]])))
                         den <- Rfast::Pmax(delta * ((zipfR::Igamma((p/2)+1, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((p/2)+1, delta/2, lower = FALSE))), num)
                         W[whichMis[[j]],g] <- num / den
                    }
               }

               #############
               # E-Step: Z #
               #############
               for (g in 1:G) {
                    Z[,g] <- li[,g] / Rfast::rowsums(li)
               }

               ###########
               # CM-Step #
               ###########
               for (g in 1:G) {

                    for (j in nMisPatterns) {
                         o <- !misPatterns[j,]
                         m <- !o
                         correction <- Sigmas[[g]][m,o] %*% solve(Sigmas[[g]][o,o])

                         for (i in whichMis[[j]]) {
                              x[i,m] <- Mus[[g]][m] + correction %*% (x[i,o] - Mus[[g]][o])
                         }
                    }

                    Mus[[g]] <- Rfast::colsums(Z[,g] * W[,g] * x) / sum(Z[,g] * W[,g])
                    sig_temp <- Rfast::Crossprod(Z[obs,g] * W[obs,g] * Rfast::eachrow(matrix(x[obs,], ncol = p), Mus[[g]], '-'), Rfast::eachrow(matrix(x[obs,], ncol = p), Mus[[g]], '-'))

                    for (j in nMisPatterns) {
                         m <- misPatterns[j,]
                         o <- !m
                         pObs <- sum(o)
                         pMis <- sum(m)

                         cross <- matrix(data = rep(0, p^2), nrow = p)
                         cross[o,o] <- Rfast::Crossprod(Z[whichMis[[j]],g] * W[whichMis[[j]],g] * Rfast::eachrow(matrix(x[whichMis[[j]],o], ncol = pObs), Mus[[g]][o], '-'), Rfast::eachrow(matrix(x[whichMis[[j]],o], ncol = pObs), Mus[[g]][o], '-'))
                         cross[o,m] <- Rfast::Crossprod(Z[whichMis[[j]],g] * W[whichMis[[j]],g] * Rfast::eachrow(matrix(x[whichMis[[j]],o], ncol = pObs), Mus[[g]][o], '-'), Rfast::eachrow(matrix(x[whichMis[[j]],m], ncol = pMis), Mus[[g]][m], '-'))
                         cross[m,o] <- t(cross[o,m])
                         cross[m,m] <- sum(Z[whichMis[[j]],g]) * (Sigmas[[g]][m,m] - Sigmas[[g]][m,o] %*% solve(Sigmas[[g]][o,o]) %*% Sigmas[[g]][o,m]) + Rfast::Crossprod(Z[whichMis[[j]],g] * W[whichMis[[j]],g] * Rfast::eachrow(matrix(x[whichMis[[j]],m], ncol = pMis), Mus[[g]][m], '-'), Rfast::eachrow(matrix(x[whichMis[[j]],m], ncol = pMis), Mus[[g]][m], '-'))

                         sig_temp <- sig_temp + cross
                    }

                    Sigmas[[g]] <- sig_temp / sum(Z[,g])
                    Thetas[g] <- CMstep2(x, Z[,g], Mus[[g]], Sigmas[[g]], naMat = M, obsInd = obs, misInd = mis, ...)
               }

               ######################
               # Compute Likelihood #
               ######################
               for (g in 1:G) {
                    li[obs,g] <- sum(Z[,g]) * dmtin(x[obs,], Mus[[g]], Sigmas[[g]], Thetas[g], ...) / n
               }

               for (j in nMisPatterns) {
                    o <- !misPatterns[j,]
                    pObs <- sum(o)

                    for (g in 1:G) {
                         li[whichMis[[j]],g] <- sum(Z[,g]) * dmtin(x[whichMis[[j]],o], Mus[[g]][o], Sigmas[[g]][o,o], Thetas[g], ...) / n
                    }
               }

               L[k] <- sum(log(Rfast::rowsums(li)))
               if (verbose) print(L[k])

               ########################
               # Aitken's Convergence #
               ########################
               if (k > 2) {
                    ak <- (L[k] - L[k-1])/(L[k-1] - L[k-2])
                    linf <- L[k-1] + 1/(1 - ak) * (L[k] - L[k-1])

                    diffs[k] <- linf - L[k-1]
                    if (diffs[k] > 0) {
                         if (diffs[k] < tol) {
                              iter <- k
                              break
                         }
                    }
               }
          }

          ##################
          # Prepare Output #
          ##################
          cluster <- apply(Z, 1, which.max)
          if (iter == 0) {
               iter <- max_iter
          }

          for (i in mis) {
               g <- cluster[i]
               m <- M[i,]
               o <- !m

               x[i,m] <- Mus[[g]][m] + Sigmas[[g]][m,o] %*% solve(Sigmas[[g]][o,o]) %*% (x[i,o] - Mus[[g]][o])
          }

          L <- L[1:iter]

          npar <- list(
               pi    = G - 1,
               mu    = G * p,
               Sigma = G * p * (p + 1) / 2,
               theta = G
          )

          total_params <- Reduce('+', npar)

          AIC <- -2 * L[iter] + 2 * total_params
          BIC <- -2 * L[iter] + total_params * log(n)

          KIC  <- -2 * L[iter] + 3 * (total_params + 1)
          KICc <- -2 * L[iter] + 2 * (total_params + 1) * n/(n-total_params -2) - n * digamma((n-total_params)/2) + n * log(n/2)

          AIC3 <- -2 * L[iter] + 3 * total_params
          CAIC <- -2 * L[iter] + total_params * (1 + log(n))
          AICc <- -2 * L[iter] + 2 * total_params * n/(n - total_params - 1)

          if (iter == max_iter) warning("Convergence criterion not met. Try increasing max_iter.")
          if (!all(L == sort(L))) warning("Likelihood did not increase monotonically. Try using formula = 'indirect'.")
          if (plot_likelihood) {
               plot(1:iter, L, xlab = "Iteration", ylab = "Log-likelihood")
          }

          return(list(Mu = Mus, Sigma = Sigmas, Theta = Thetas, X = x, Pi = Pis, Z = Z, W = W, L = L, cluster = cluster, iterations = iter, npar = total_params, AIC = AIC, BIC = BIC, KIC = KIC, KICc = KICc, AIC3 = AIC3, CAIC = CAIC, AICc = AICc))
     } # Multivariate Case
}
