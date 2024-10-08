#' Fit cluster weighted regression model using the MTIN distribution
#'
#' @param x A matrix whose columns are the predictor variables.
#' @param y A matrix whose columns are the response variables.
#' @param G An integer representing the number of clusters.
#' @param max_iter An integer representing the maximum number of iterations.
#' @param tol Float representing tolerance for Aitken's acceleration criterion.
#' @param model Type of regression model to fit: \code{"fixed"}, \code{"random"}
#' @param formula Method used to calculate density passed to dmtin: \code{"direct"}, \code{"indirect"}, \code{"series"}
#' @param init_method Model used for initialization: \code{'mclust'}, \code{"kmedoids"}, \code{"kmeans"}, \code{"heirarchical"}
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
#' @examples
#' \donttest{
#' data <- palmerpenguins::penguins
#' keep <- !rowSums(is.na(data[,3:6])) == 4
#' data <- data[keep,]
#'
#' res <- tinreg(x = data$flipper_length_mm, y = data$bill_length_mm, G = 3, model = 'random')
#' plot(res$X, res$Y, col = res$cluster)
#' abline(a = res$Beta[[1]][1,1], b = res$Beta[[1]][2,1], col = 'black')
#' abline(a = res$Beta[[2]][1,1], b = res$Beta[[2]][2,1], col = 'red')
#' abline(a = res$Beta[[3]][1,1], b = res$Beta[[3]][2,1], col = 'green')
#' }
#'
#' @export
tinreg <- function(x, y, G = 1, max_iter = 100, tol = 10^-1, model = c("fixed", "random"), init_method = c("mclust", "kmedoids", "kmeans", "heirarchical"), verbose = FALSE, plot_likelihood = FALSE, ...) {
     model <- match.arg(model)

     if(model == 'fixed') {
          ####################
          # Input formatting #
          ####################
          x <- as.matrix(x)
          if (any(is.na(x))) stop("Missing values in X are not allowed.")
          y <- as.matrix(y)
          x <- cbind(rep(1, nrow(x)), x)
          M_y <- is.na(y)
          obs_y <- which(!as.logical(rowSums(is.na(y))))
          mis_y <- which((as.logical(rowSums(is.na(y)))))
          n <- nrow(x)
          p_y <- ncol(y)
          p_x <- ncol(x)
          if (any(rowSums(M_y) == p_y)) stop("Remove observations containing only NAs")

          ############################
          # Declare Model Parameters #
          ############################
          Betas <- rep(list(matrix(nrow = p_x, ncol = p_y, data = rep(0, p_y*p_x))))
          Sigmas_y <- rep(list(matrix(data = rep(0, p_y^2), nrow = p_y)))
          Pis <- rep(0,G)
          Thetas_y <- rep(0.2, G)
          Z <- matrix(0, nrow = n, ncol = G)
          W_y <- matrix(0, nrow = n, ncol = G)
          y_center <- rep(list(matrix(data = rep(0, p_y*n), nrow = n)), G)
          iter <- 0

          #########################
          # Store Log-Likelihoods #
          #########################
          li <- matrix(0, nrow = n, ncol = G)
          L <- rep(0, max_iter)
          diffs <- rep(0, max_iter)

          #########################
          # Initialize Parameters #
          #########################
          Z <- init_clusters(cbind(x,y), G, obsInd = obs_y, init_method = init_method)

          for (g in 1:G) {
               Betas[[g]] <- Rfast::mat.mult(solve(Rfast::Crossprod(Z[obs_y,g] * as.matrix(x[obs_y,], ncol = p_x), matrix(x[obs_y,], ncol = p_x))), Rfast::Crossprod(Z[obs_y,g] * matrix(x[obs_y,], ncol = p_x), matrix(y[obs_y,], ncol = p_y)))
               y_center[[g]][obs_y,] <- y[obs_y,] - Rfast::mat.mult(matrix(x[obs_y,], ncol = p_x), Betas[[g]])
               Sigmas_y[[g]] <- Rfast::cova(matrix(y_center[[g]][as.logical(Z[, g]), ], ncol = p_y))
               li[obs_y,g] <- sum(Z[,g]) * dmtin(matrix(y_center[[g]][obs_y,], ncol = p_y), rep(0, p_y), Sigmas_y[[g]], Thetas_y[g], ...) / n
          }

          for (i in mis_y) {
               m <- M_y[i,]
               o <- !M_y[i,]
               pObs <- sum(o)

               for (g in 1:G) {
                    li[i,g] <- sum(Z[,g]) * dmtin(y_center[[g]][i,o], rep(0, pObs), Sigmas_y[[g]][o,o], Thetas_y[g], ...) / n
               }
          }


          ##################
          # ECME Algorithm #
          ##################
          for (k in 1:max_iter) {
               #############
               # E-Step: Y #
               #############
               for (i in mis_y) {
                    m <- M_y[i,]
                    o <- !m

                    for (g in 1:G) {
                         Mu <- x[i,] %*% Betas[[g]]
                         y_center[[g]][i,m] <- Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])
                         y_center[[g]][i,o] <- y[i,o] - Mu[o]
                    }

               }

               #############
               # E-Step: W #
               #############
               for (g in 1:G) {
                    delta <- Rfast::mahala(matrix(y_center[[g]][obs_y,], ncol = p_y), rep(0, p_y), Sigmas_y[[g]])
                    num <- Rfast::Pmax(2 * ((zipfR::Igamma((p_y/2)+2, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((p_y/2)+2, delta/2, lower = FALSE))), rep(10^(-322), length(obs_y)))
                    den <- Rfast::Pmax(delta * ((zipfR::Igamma((p_y/2)+1, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((p_y/2)+1, delta/2, lower = FALSE))), num)
                    W_y[obs_y,g] <- num/den
               }

               for (i in mis_y) {
                    m <- M_y[i,]
                    o <- !M_y[i,]
                    pObs <- sum(!M_y[i,])

                    for (g in 1:G) {
                         delta <- Rfast::mahala(matrix(y_center[[g]][i,o], ncol = pObs), rep(0, pObs), Sigmas_y[[g]][o,o])
                         num <- max(2 * ((zipfR::Igamma((pObs/2)+2, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((pObs/2)+2, delta/2, lower = FALSE))), 10^(-322))
                         den <- max(delta * ((zipfR::Igamma((pObs/2)+1, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((pObs/2)+1, delta/2, lower = FALSE))), num)
                         W_y[i,g] <- num / den
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

                    for (i in mis_y) {
                         m <- M_y[i,]
                         o <- !m

                         Mu <- x[i,] %*% Betas[[g]]

                         y[i,m] <- Mu[m] + Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])
                    }

                    Betas[[g]] <- Rfast::mat.mult(solve(Rfast::Crossprod(Z[,g] * W_y[,g] * x, x)), Rfast::Crossprod(x, Z[,g] * W_y[,g] * matrix(y, ncol = p_y)))
                    y_center[[g]] <- y - Rfast::mat.mult(x, Betas[[g]])
                    sig_temp_y <- Rfast::Crossprod(Z[obs_y,g] * W_y[obs_y,g] * matrix(y_center[[g]][obs_y,], ncol = p_y), matrix(y_center[[g]][obs_y,], ncol = p_y))

                    for (i in mis_y) {
                         m <- M_y[i,]
                         o <- !M_y[i,]

                         Mu <- x[i,] %*% Betas[[g]]

                         y[i,m] <- Mu[m] + Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])

                         cross <- matrix(data = rep(0, p_y^2), nrow = p_y)
                         cross[o,o] <- Z[i,g] * W_y[i,g] * (y[i,o] - Mu[o]) %*% t(y[i,o] - Mu[o])
                         cross[o,m] <- Z[i,g] * W_y[i,g] * (y[i,o] - Mu[o]) %*% t(y[i,m] - Mu[m])
                         cross[m,o] <- t(cross[o,m])
                         cross[m,m] <- Z[i,g] * (Sigmas_y[[g]][m,m] - Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% Sigmas_y[[g]][o,m]) + Z[i,g] * W_y[i,g] * (y[i,m] - Mu[m]) %*% t(y[i,m] - Mu[m])

                         sig_temp_y <- sig_temp_y + cross
                    }

                    Sigmas_y[[g]] <- sig_temp_y / sum(Z[,g])
                    Thetas_y[g] <- CMstep2(matrix(y_center[[g]], ncol = p_y), Z[,g], rep(0, p_y), Sigmas_y[[g]], naMat = M_y, obsInd = obs_y, misInd = mis_y, ...)
                    Pis[g] <- sum(Z[,g]) / n
               }

               ######################
               # Compute Likelihood #
               ######################
               for (g in 1:G) {
                    li[obs_y,g] <- sum(Z[,g]) * dmtin(matrix(y_center[[g]][obs_y,], ncol = p_y), rep(0, p_y), Sigmas_y[[g]], Thetas_y[g], ...) / n

               }

               for(i in mis_y) {
                    o <- !M_y[i,]
                    pObs <- sum(o)

                    for (g in 1:G) {
                         li[i,g] <- sum(Z[,g]) * dmtin(y_center[[g]][i,o], rep(0, pObs), Sigmas_y[[g]][o,o], Thetas_y[g], ...) / n
                    }
               }

               L[k] <- sum(Rfast::Log(Rfast::rowsums(li)))

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

          if (!(all(L == sort(L)))) warning("Likelihood did not increase monotonically. Try using 'formula = indirect'.")


          ##################
          # Prepare Output #
          ##################
          for (i in mis_y) {
               g <- which.max(Z[i,])
               m <- M_y[i,]
               o <- !m
               Mu <- x[i,] %*% Betas[[g]]

               y[i,m] <- Mu[m] + Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])
          }

          cluster <- apply(Z, 1, which.max)
          if (iter == 0) iter <- max_iter
          if (plot_likelihood) {
               plot(1:iter, L[1:iter], xlab = "Iteration", ylab = "Log-Likelihood")
          }

          L <- L[1:iter]

          npar <- list(
               pi    = G - 1,
               beta = G * length(Betas[[1]]),
               Sigma_y = G * p_y * (p_y + 1) / 2,
               theta_y = G
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

          return(list(Beta = Betas, Theta_y = Thetas_y, X = x, Y = y, Sigma_y = Sigmas_y, Pi = Pis, Z = Z, W_y = W_y, L = L, cluster = cluster, Model = model, iterations = iter, npar = total_params, AIC = AIC, BIC = BIC, KIC = KIC, KICc = KICc, AIC3 = AIC3, CAIC = CAIC, AICc = AICc))

     }

     if (model == 'random') {

          ####################
          # Input Formatting #
          ####################
          x <- as.matrix(x)
          y <- as.matrix(y)
          x <- cbind(rep(1, nrow(x)), x)
          M_y <- is.na(y)
          if (any(is.na(x))) stop("Missing values in X are not allowed.")
          obs_y <- which(!as.logical(rowSums(is.na(y))))
          mis_y <- which((as.logical(rowSums(is.na(y)))))
          n <- nrow(x)
          p_y <- ncol(y)
          p_x <- ncol(x)
          iter <- 0
          if (any(rowSums(M_y) == p_y)) stop("Remove observations containing only NAs")

          ####################
          # Model Parameters #
          ####################
          Betas <- rep(list(matrix(nrow = p_x, ncol = p_y, data = rep(0, p_y*p_x))))
          Sigmas_y <- rep(list(matrix(data = rep(0, p_y^2), nrow = p_y)))
          Sigmas_x <- rep(list(matrix(data = rep(0, p_x^2), nrow = p_x)))
          Pis <- rep(0,G)
          Thetas_y <- rep(0.2, G)
          Thetas_x <- rep(0.2, G)
          Z <- matrix(0, nrow = n, ncol = G)
          W_y <- matrix(0, nrow = n, ncol = G)
          W_x <- matrix(0, nrow = n, ncol = G)
          y_center <- rep(list(matrix(data = rep(0, p_y*n), nrow = n)), G)
          Mus_x <- rep(list(matrix(nrow = 1, ncol = p_x-1, data = rep(0, p_x-1))))

          #########################
          # Store Log-Likelihoods #
          #########################
          li <- matrix(0, nrow = n, ncol = G)
          L <- rep(0, max_iter)
          diffs <- rep(0, max_iter)

          #########################
          # Initialize Parameters #
          #########################
          Z <- init_clusters(cbind(x,y), G, obsInd = obs_y, init_method = init_method)

          for (g in 1:G) {
               Betas[[g]] <- Rfast::mat.mult(solve(Rfast::Crossprod(Z[obs_y,g] * as.matrix(x[obs_y,], ncol = p_x), matrix(x[obs_y,], ncol = p_x))), Rfast::Crossprod(Z[obs_y,g] * matrix(x[obs_y,], ncol = p_x), matrix(y[obs_y,], ncol = p_y)))
               Mus_x[[g]] <- Rfast::colsums(as.matrix(x[as.logical(Z[, g]),2:p_x], ncol = p_x - 1)) / sum(Z[,g])
               y_center[[g]][obs_y,] <- as.matrix(y[obs_y,], ncol = p_y) - Rfast::mat.mult(x[obs_y,], Betas[[g]])
               Sigmas_y[[g]] <- Rfast::cova(as.matrix(y_center[[g]][as.logical(Z[, g]), ], ncol = p_y))
               Sigmas_x[[g]] <- Rfast::cova(as.matrix(x[as.logical(Z[, g]), 2:p_x], ncol = p_x))

               li[,g] <- sum(Z[,g]) * dmtin(as.matrix(x[,2:p_x], ncol = p_x - 1), Mus_x[[g]], Sigmas_x[[g]], Thetas_x[g], ...) / n
               li[obs_y,g] <- li[obs_y,g] * dmtin(as.matrix(y_center[[g]][obs_y,], ncol = p_y), rep(0, p_y), Sigmas_y[[g]], Thetas_y[g], ...)
          }

          for (i in mis_y) {
               m <- M_y[i,]
               o <- !M_y[i,]
               pObs <- sum(o)

               for (g in 1:G) {
                    li[i,g] <- li[i,g] * dmtin(y_center[[g]][i,o], rep(0, pObs), Sigmas_y[[g]][o,o], Thetas_y[g], ...)
               }
          }

          ##################
          # ECME Algorithm #
          ##################
          for (k in 1:max_iter) {
               #############
               # E-Step: Y #
               #############
               for (i in mis_y) {
                    m <- M_y[i,]
                    o <- !m

                    for (g in 1:G) {
                         Mu <- x[i,] %*% Betas[[g]]
                         y_center[[g]][i,m] <- Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])
                         y_center[[g]][i,o] <- y[i,o] - Mu[o]
                    }
               }

               #############
               # E-Step: W #
               #############
               for (g in 1:G) {
                    delta <- Rfast::mahala(as.matrix(y_center[[g]][obs_y,], ncol = p_y), rep(0, p_y), Sigmas_y[[g]])
                    num <- Rfast::Pmax(2 * ((zipfR::Igamma((p_y/2)+2, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((p_y/2)+2, delta/2, lower = FALSE))), rep(10^(-322), length(obs_y)))
                    den <- Rfast::Pmax(delta * ((zipfR::Igamma((p_y/2)+1, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((p_y/2)+1, delta/2, lower = FALSE))), num)
                    W_y[obs_y,g] <- num/den

                    delta <- Rfast::mahala(as.matrix(x[,2:p_x], ncol = p_x - 1), Mus_x[[g]], Sigmas_x[[g]])
                    num <- Rfast::Pmax(2 * ((zipfR::Igamma(((p_x - 1)/2)+2, (1-Thetas_x[g])*delta/2, lower = FALSE) - zipfR::Igamma(((p_x - 1)/2)+2, delta/2, lower = FALSE))), rep(10^(-322), n))
                    den <- Rfast::Pmax(delta * ((zipfR::Igamma(((p_x - 1)/2)+1, (1-Thetas_x[g])*delta/2, lower = FALSE) - zipfR::Igamma(((p_x - 1)/2)+1, delta/2, lower = FALSE))), num)
                    W_x[,g] <- num/den
               }

               for (i in mis_y) {
                    o <- !M_y[i,]
                    pObs <- sum(o)

                    for (g in 1:G) {
                         delta <- Rfast::mahala(as.matrix(y_center[[g]][i,o], ncol = pObs), rep(0, pObs), Sigmas_y[[g]][o,o])
                         num <- max(2 * ((zipfR::Igamma((pObs/2)+2, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((pObs/2)+2, delta/2, lower = FALSE))), 10^(-322))
                         den <- max(delta * ((zipfR::Igamma((pObs/2)+1, (1-Thetas_y[g])*delta/2, lower = FALSE) - zipfR::Igamma((pObs/2)+1, delta/2, lower = FALSE))), num)
                         W_y[i,g] <- num / den

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

                    for (i in mis_y) {
                         m <- M_y[i,]
                         o <- !m

                         Mu <- x[i,] %*% Betas[[g]]

                         y[i,m] <- Mu[m] + Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])
                    }

                    Mus_x[[g]] <- Rfast::colsums(Z[,g] * W_x[,g] * as.matrix(x[,2:p_x], ncol = p_x - 1)) / sum(Z[,g] * W_x[,g])
                    Betas[[g]] <- Rfast::mat.mult(solve(Rfast::Crossprod(x, Z[,g] * W_y[,g] * x)), Rfast::Crossprod(x, Z[,g] * W_y[,g] * y))
                    y_center[[g]] <- as.matrix(y, ncol = p_y) - Rfast::mat.mult(x, Betas[[g]])
                    Sigmas_x[[g]] <- Rfast::Crossprod(Z[,g] * W_x[,g] * Rfast::eachrow(as.matrix(x[, 2:p_x], ncol = p_x - 1), Mus_x[[g]], '-'), Rfast::eachrow(as.matrix(x[, 2:p_x], ncol = p_x - 1), Mus_x[[g]], '-')) / sum(Z[,g])
                    sig_temp_y <- Rfast::Crossprod(Z[obs_y,g] * W_y[obs_y,g] * matrix(y_center[[g]][obs_y,], ncol = p_y), matrix(y_center[[g]][obs_y,], ncol = p_y))

                    for (i in mis_y) {
                         m <- M_y[i,]
                         o <- !m

                         Mu <- x[i,] %*% Betas[[g]]

                         y[i,m] <- Mu[m] + Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])

                         cross <- matrix(data = rep(0, p_y^2), nrow = p_y)
                         cross[o,o] <- Z[i,g] * W_y[i,g] * (y[i,o] - Mu[o]) %*% t(y[i,o] - Mu[o])
                         cross[o,m] <- Z[i,g] * W_y[i,g] * (y[i,o] - Mu[o]) %*% t(y[i,m] - Mu[m])
                         cross[m,o] <- t(cross[o,m])
                         cross[m,m] <- Z[i,g] * (Sigmas_y[[g]][m,m] - Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% Sigmas_y[[g]][o,m]) + Z[i,g] * W_y[i,g] * (y[i,m] - Mu[m]) %*% t(y[i,m] - Mu[m])

                         sig_temp_y <- sig_temp_y + cross
                    }

                    Sigmas_y[[g]] <- sig_temp_y / sum(Z[,g])
                    Thetas_y[g] <- CMstep2(as.matrix(y_center[[g]], ncol = p_y), Z[,g], rep(0, p_y), Sigmas_y[[g]], naMat = M_y, obsInd = obs_y, misInd = mis_y, ...)
                    Thetas_x[g] <- CMstep2(as.matrix(x[,2:p_x], ncol = p_x - 1), Z[,g], Mus_x[[g]], Sigmas_x[[g]], ...)
                    Pis[g] <- sum(Z[,g]) / n
               }

               ######################
               # Compute Likelihood #
               ######################
               for (g in 1:G) {
                    li[,g] <- sum(Z[,g]) * dmtin(as.matrix(x[,2:p_x], ncol = p_x - 1), Mus_x[[g]], Sigmas_x[[g]], Thetas_x[g], ...) / n
                    li[obs_y,g] <- li[obs_y,g] * dmtin(as.matrix(y_center[[g]][obs_y,], ncol = p_y), rep(0, p_y), Sigmas_y[[g]], Thetas_y[g], ...)

               }

               for(i in mis_y) {
                    o <- !M_y[i,]
                    pObs <- sum(o)

                    for (g in 1:G) {
                         li[i,g] <- li[i,g] * dmtin(y_center[[g]][i,o], rep(0, pObs), Sigmas_y[[g]][o,o], Thetas_y[g], ...)
                    }
               }

               L[k] <- sum(log(Rfast::rowsums(li)))

               if (verbose) print(L[k])

               ########################
               # Aitken's Convergence #
               ########################
               {
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
          }

          ##################
          # Prepare Output #
          ##################
          for (i in mis_y) {
               g <- which.max(Z[i,])
               m <- M_y[i,]
               o <- !m
               Mu <- x[i,] %*% Betas[[g]]

               y[i,m] <- Mu[m] + Sigmas_y[[g]][m,o] %*% solve(Sigmas_y[[g]][o,o]) %*% (y[i,o] - Mu[o])
          }

          if (!all(L == sort(L))) warning("Likelihood did not increase monotonically. Try using formula = 'indirect'.")

          cluster <- apply(Z, 1, which.max)
          if (iter == 0) {
               iter <- max_iter
          }

          if (plot_likelihood) {
               plot(1:iter, L[1:iter], xlab = "Iteration", ylab = "Log-Likelihood")
          }

          L <- L[1:iter]

          npar <- list(
               pi    = G - 1,
               mu_x    = G * p_x,
               Sigma_x = G * p_x * (p_x + 1) / 2,
               theta_x    = G,
               beta = G * length(Betas[[1]]),
               Sigma_y = G * p_y * (p_y + 1) / 2,
               theta_y = G
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

          return(list(Beta = Betas, Mu_x = Mus_x, Theta_x = Thetas_x, Theta_y = Thetas_y, X = x, Y = y, Sigma_x = Sigmas_x, Sigma_y = Sigmas_y, Pi = Pis, Z = Z, W_x = W_x, W_y = W_y, L = L, cluster = cluster, Model = model, iterations = iter, npar = total_params, AIC = AIC, BIC = BIC, KIC = KIC, KICc = KICc, AIC3 = AIC3, CAIC = CAIC, AICc = AICc))

     }
}
