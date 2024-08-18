#' Fit cluster weighted regression model using the MTIN distribution
#'
#' @param x A data matrix.
#' @param G An integer representing the number of clusters.
#' @param max_iter An integer representing the maximum number of iterations.
#' @param tol Float representing tolerance for Aitken's acceleration criterion.
#' @param formula Method used to calculate density passed to \code{dmtin}: \code{"direct"}, \code{"indirect"}, \code{"series"}
#' @param init_method Model used for initialization: \code{'mclust'}, \code{"kmedoids"}, \code{"kmeans"}, \code{"heirarchical"}
#' @param verbose Print log-likelihood at every iteration.
#' @param plot_likelihood Plot log-likelihoods after convergence.
#'
#' @return A list with the following elements:
#' \item{Beta}{List of regression coefficients.}
#' \item{Mu_x}{List of mean vectors for independent variable.}
#' \item{Sigma_x}{List of scale matrices for independent variable.}
#' \item{Theta_x}{List of tailedness parameters for independent variable.}
#' \item{Sigma_y}{List of scale matrices for dependent variable.}
#' \item{Theta_y}{List of tailedness parameters for dependent variable.}
#' \item{Pi}{Mixing proportions.}
#' \item{Z}{Expected value of Z.}
#' \item{W_y}{Expected value of W for dependent variable.}
#' \item{W_X}{Expected value of W for independent variable.}
#' \item{cluster}{Membership labels.}
#' \item{L}{List of log-likelihoods at each iteration.}
#' \item{Model}{Type of model fitted.}
#' \item{AIC}{Akaike information criterion.}
#' \item{BIC}{Bayesian information criterion.}
#' \item{KIC}{Kullback information criterion.}??
#'
#' @export
tinclust <- function(x, G = 1, max_iter = 100, tol = 10^-1, formula = c("direct", "indirect", "series"), init_method = c("mclust", "kmedoids", "kmeans", "heirarchical"), verbose = FALSE, plot_likelihood = TRUE) {
     formula <- match.arg(formula)

     ####################
     # Input Formatting #
     ####################
     x <- as.matrix(x)
     M <- is.na(x)
     n <- nrow(x)
     p <- ncol(x)
     mis <- which(as.logical(rowSums(M)))
     obs <- which(!as.logical(rowSums(M)))


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
          Mus[[g]] <- colMeans(x[Z[,g] == 1,])
          Sigmas[[g]] <- cov(x[Z[,g] == 1,])
          Pis[g] <- sum(Z[,g]) / length(obs)
          li[obs,g] <- Pis[g] * dmtin(x[obs,], Mus[[g]], Sigmas[[g]], Thetas[g], formula = formula)
     }

     for (i in mis) {
          m <- M[i,]
          o <- !m

          li[i,g] <- Pis[g] * dmtin(x[i,o], Mus[[g]][o], Sigmas[[g]][o,o], Thetas[g], formula = formula)
     }


     ##################
     # ECME Algorithm #
     ##################
     for (k in 1:max_iter) {
          #############
          # E-Step: W #
          #############
          for (g in 1:G) {
               delta <- mahalanobis(as.matrix(x[obs,], ncol = p), Mus[[g]], Sigmas[[g]], tol = 1e-20)
               num <- pmax(2 * ((zipfR::Igamma((p/2)+2, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((p/2)+2, delta/2, lower = FALSE))), rep(10^(-322), length(obs)))
               den <- pmax(delta * ((zipfR::Igamma((p/2)+1, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((p/2)+1, delta/2, lower = FALSE))), num)
               W[obs,g] <- num/den
          }

          for (i in mis) {
               m <- M[i,]
               o <- !M[i,]
               pObs <- sum(o)

               for (g in 1:G) {
                    delta <- mahalanobis(as.matrix(x[i,o], ncol = pObs), Mus[[g]][o], Sigmas_y[[g]][o,o])
                    num <- max(2 * ((zipfR::Igamma((pObs/2)+2, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((pObs/2)+2, delta/2, lower = FALSE))), 10^(-322))
                    den <- max(delta * ((zipfR::Igamma((pObs/2)+1, (1-Thetas[g])*delta/2, lower = FALSE) - zipfR::Igamma((pObs/2)+1, delta/2, lower = FALSE))), num)
                    W[i,g] <- num / den
               }
          }

          #############
          # E-Step: Z #
          #############
          for (g in 1:G) {
               Z[,g] <- li[,g] / rowSums(li)
          }

          ###########
          # CM-Step #
          ###########
          for (g in 1:G) {
               for (i in mis) {
                    m <- M[i,]
                    o <- !m

                    x[i,m] <- Mus[[g]][m] + Sigmas[[g]][m,o] %*% solve(Sigmas[[g]][o,o]) %*% (x[i,o] - Mus[[g]][o])
               }

               Mus[[g]] <- colSums(Z[,g] * W[,g] * x) / sum(Z[,g] * W[,g])
               Sigmas[[g]] <- crossprod(sqrt(Z[,g] * W[,g]) * sweep(x, 2, Mus[[g]]))

               for (i in mis) {
                    m <- M[i,]
                    o <- !m

                    cross <- matrix(data = rep(0, p^2), nrow = p)
                    cross[o,o] <- Z[i,g] * W[i,g] * (x[i,o] - Mus[[g]][o]) %*% t(x[i,o] - Mus[[g]][o])
                    cross[o,m] <- Z[i,g] * W[i,g] * (x[i,o] - Mus[[g]][o]) %*% t(x[i,m] - Mus[[g]][m])
                    cross[m,o] <- t(cross[o,m])
                    cross[m,m] <- Z[i,g] * (Sigmas[[g]][m,m] - Sigmas[[g]][m,o] %*% solve(Sigmas[[g]][o,o]) %*% Sigmas[[g]][o,m]) + Z[i,g] * W[i,g] * (x[i,m] - Mus[[g]][m]) %*% t(x[i,m] - Mus[[g]][m])

                    Sigmas[[g]] <- Sigmas[[g]] + cross
               }

               Sigmas[[g]] <- Sigmas[[g]] / sum(Z[,g])
               Thetas[g] <- CMstep2(x, Z[,g], Mus[[g]], Sigmas[[g]], formula = formula, naMat = M, obsInd = obs, misInd = mis)
          }

          ######################
          # Compute Likelihood #
          ######################
          for (g in 1:G) {
               li[,g] <- sum(Z[,g]) * dmtin(x[obs,], Mus[[g]], Sigmas[[g]], Thetas[g], formula = formula) / n
          }

          for (i in mis) {
               o <- !M[i,]

               for (g in 1:G) {
                    li[i,g] <- sum(Z[,g]) * dmtin(x[i,o], Mus[[g]][o], Sigmas[[g]][o,o], Thetas[g], formula = formula)
               }
          }

          L[k] <- sum(log(rowSums(li)))
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



















}
