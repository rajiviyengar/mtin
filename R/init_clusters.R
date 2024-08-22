#' Initialize membership labels
#'
#' @keywords internal
#'
#' @param x A matrix whose columns are variables.
#' @param G An integer representing the number of clusters.
#' @param obsInd Indices of rows containing complete observations.
#' @param labels Vector of membership labels to be included if using manual initialization.
#' @param init_method Model used for initialization: \code{'mclust'}, \code{"kmedoids"}, \code{"kmeans"}, \code{"heirarchical"}
#'
#' @return
#' \item{Z}{An (n x G) matrix representing membership labels.}

init_clusters <- function(X, G, obsInd = NULL, init_method = c("mclust", "kmedoids", "kmeans", "heirarchical")) {
     init_method <- match.arg(init_method)

     Z <- matrix(0, ncol = G, nrow = nrow(X)) # Membership matrix to be return

     if (init_method == "mclust") {
          init <- mclust::Mclust(data = X[obsInd,], G = G)

          for (i in 1:length(obsInd)) {
               Z[obsInd[i], init$classification[i]] <- 1
          }

          return(Z)
     }

     if (init_method == "kmedoids") {
          init <- cluster::pam(x = X[obsInd,], k = G)

          for (i in 1:length(obsInd)) {
               Z[obsInd[i], init$clustering[i]] <- 1
          }

          return(Z)
     }

     if (init_method == "kmeans") {
          init <- stats::kmeans(x = X[obsInd,], centers = G)

          for (i in 1:length(obsInd)) {
               Z[obsInd[i], init$cluster[i]] <- 1
          }

          return(Z)
     }

     if (init_method == "heirarchical") {
          d <- stats::dist(X[obsInd,])
          tree <- stats::hclust(d)
          init <- stats::cutree(tree = tree, k = G)

          for (i in 1:length(obsInd)) {
               Z[obsInd[i], init[i]] <- 1
          }

          return(Z)
     }
}
