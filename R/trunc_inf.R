# ----- functions to test the means of clusters -----

#' Exact significance test for hierarchical clustering
#'
#' This tests the null hypothesis of no difference in means between 
#' clusters \code{k1} and \code{k2} at level \code{K} in a hierarchical clustering. 
#'(The \code{K} clusters are numbered as per the results of the \code{cutree} 
#' function in the \code{stats} package.)
#' 
#' In order to account for the fact that the clusters have been estimated from the data, 
#' the p-values are computed conditional on the fact that those clusters were estimated. 
#' This function computes p-values exactly via an analytic characterization of the conditioning set. 
#' 
#' Currently, this function supports squared Euclidean distance as a measure of dissimilarity 
#' between observations, and the following six linkages: single, average, centroid, Ward, 
#' McQuitty (also known as WPGMA), and median (also kown as WPGMC). 
#' 
#' By default, this function assumes that the covariance matrix of the features is isotropic 
#' i.e. \eqn{Cov(X_i) = \sigma^2 I_p}. Setting \code{iso} to \code{FALSE} instead assumes that 
#' \eqn{Cov(X_i) = \Sigma}. If known, \eqn{\sigma} can be passed in using the \code{sigma} argument 
#' or \eqn{\Sigma^{-1}} can be passed in the \code{SigInv} argument; otherwise, an 
#' estimate of \eqn{\sigma} or \eqn{\Sigma} will be used. 
#'
#' @export
#'
#' @param X \eqn{n} by \eqn{p} matrix containing numeric data.
#' @param link String selecting the linkage. Supported options are 
#' \code{"single", "average", "centroid", "ward.D", "median"}, and \code{"mcquitty"}.
#' @param hcl Object of the type \code{hclust} containing the hierarchical clustering of X. 
#' @param K Integer selecting the total number of clusters.
#' @param k1,k2 Integers selecting the clusters to test, as indexed by the results of \code{cutree(hcl, K)}.
#' @param iso Boolean. If \code{TRUE}, isotropic covariance matrix model, otherwise not.
#' @param sig Optional scalar specifying \eqn{\sigma}, relevant if \code{iso} is \code{TRUE}.
#' @param SigInv Optional matrix specifying \eqn{\Sigma^{-1}}, relevant if \code{iso} is \code{FALSE}.
#' @param dist The distances of matrix X
#' @return
#' \item{stat}{the test statistic: the Euclidean distance between the mean of cluster \code{k1} and the mean of cluster \code{k2}  }
#' \item{pval}{the p-value}
#' \item{trunc}{object of the type \code{Intervals} containing the conditioning set}
#'
#' @examples
#' # Simulates a 100 x 2 data set with three clusters
#' set.seed(123)
#' dat <- rbind(c(-1, 0), c(0, sqrt(3)), c(1, 0))[rep(1:3, length=100), ] + 
#' matrix(0.2*rnorm(200), 100, 2)
#'
#' # Average linkage hierarchical clustering
#' hcl <- hclust(dist(dat, method="euclidean")^2, method="average")
#' 
#' # plot dendrograms with the 1st and 2nd clusters (cut at the third split) 
#' # displayed in blue and orange 
#' plot(hcl)
#' rect_hier_clusters(hcl, k=3, which=1:2, border=c("blue", "orange"))
#' 
#' # tests for a difference in means between the blue and orange clusters
#' test_hier_clusters_exact(X=dat, link="average", hcl=hcl, K=3, k1=1, k2=2)
#' 
#' @seealso \code{\link{rect_hier_clusters}} for visualizing clusters \code{k1} and \code{k2} in the dendrogram;
#' 
#' \code{\link{test_complete_hier_clusters_approx}} for approximate p-values for complete linkage hierarchical clustering;
#' 
#' \code{\link{test_clusters_approx}} for approximate p-values for a user-specified clustering function;
#' 
#' @references Lucy L. Gao et al. "Selective inference for hierarchical clustering". 
test_hier_clusters_exact <- function(X, link, hcl, K, k1, k2, iso=TRUE, sig=NULL, SigInv=NULL, dist=NULL) {
    # error checking 
    if(!is.matrix(X)) stop("X should be a matrix")
  
    n <- nrow(X)
    q <- ncol(X)
    
    if(link == "complete") stop("Exact p-value not supported. See 'test_complete_hier_clusters_approx' for an approximate p-value.")
    if(!link %in% c("single", "average", "centroid", "ward.D", "mcquitty", "median")) stop("Linkage should be 'single', 'average', 'centroid', 'ward.D', 'mcquitty', or 'median'")
    if(!is_integer_between_a_b(K, 2, n)) stop("number of clusters (K) should be between 2 and n")
    if(!is_integer_between_a_b(k1, 1, K) | !is_integer_between_a_b(k2, 1, K)) stop(paste("cluster indices should be between 1 and K", sep=""))
    if((iso != TRUE) & (iso != FALSE)) stop("iso should be TRUE or FALSE")
    
    # hierarchical clustering with squared Euclidean distance and specified linkage
    hcl_at_K <- stats::cutree(hcl, K)

    n1 <- sum(hcl_at_K == k1)
    n2 <- sum(hcl_at_K == k2)
    squared_norm_nu <- 1/n1 + 1/n2
    diff_means <- colMeans(X[hcl_at_K == k1, , drop=FALSE]) - colMeans(X[hcl_at_K == k2, , drop=FALSE])
    
    if(is.null(dist)) dist <- stats::dist(X)^2
    
    if(iso) {
      if(is.null(sig)) {
        sig <- sqrt(sum(scale(X, scale=FALSE)^2)/(n*q - q))
      }
      
      # compute test statistic
      stat <- norm_vec(diff_means)

      # compute truncation set
      if(link == "single") S <- compute_S_single(X, hcl, K, k1, k2, dist)
      if(link == "average") S <- compute_S_average(X, hcl, K, k1, k2, dist)
      if(link == "centroid") S <-  compute_S_centroid(X, hcl, K, k1, k2, dist)
      if(link == "ward.D") S <-  compute_S_ward(X, hcl, K, k1, k2, dist)
      if(link == "mcquitty") S <-  compute_S_mcquitty(X, hcl, K, k1, k2, dist)
      if(link == "median") S <-  compute_S_median(X, hcl, K, k1, k2, dist)

      # set distribution of phi
      scale_factor <- squared_norm_nu*sig^2

    } else {
      if(is.null(SigInv)) {
        Sig <- stats::cov(scale(X, scale=F))
        SigInv <- solve(Sig)
      }
      
      # compute test statistic
      stat <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))

      # compute truncation set
      if(link == "single") S <- compute_S_single_gencov(X, hcl, K, k1, k2, stat)
      if(link == "average") S <- compute_S_average_gencov(X, hcl, K, k1, k2, stat, dist)
      if(link == "centroid") S <-  compute_S_centroid_gencov(X, hcl, K, k1, k2, stat, dist)
      if(link == "ward.D") S <-  compute_S_ward_gencov(X, hcl, K, k1, k2, stat, dist)
      if(link == "mcquitty") S <-  compute_S_mcquitty_gencov(X, hcl, K, k1, k2, stat, dist)
      if(link == "median") S <-  compute_S_median_gencov(X, hcl, K, k1, k2, stat, dist)

      # set distribution of phi
      scale_factor <- squared_norm_nu
    }

    # compute p-value using truncated chi-squared distribution
    gestat <- intervals::Intervals(c(stat^2/scale_factor, Inf))
    denom <- S^2/scale_factor
    numer <- suppressWarnings(intervals::interval_intersection(gestat, denom))
    pval <- TChisqRatioApprox(q, numer, denom)

    return(list(stat=stat, pval=pval, trunc=S))
}

#' Monte Carlo significance test for complete linkage hierarchical clustering
#'
#' This tests the null hypothesis of no difference in means between 
#' clusters \code{k1} and \code{k2} at level \code{K} in a complete 
#' linkage hierarchical clustering. (The \code{K} clusters are numbered as per 
#' the results of the \code{cutree} function in the \code{stats} package.)
#' 
#' Important note: Before calling \code{hclust} and this function, make sure to 
#' load the package \code{fastcluster}. This is because the p-value approximation 
#' procedure requires running hierarchical clustering on a large number of simulated 
#' data sets, and the version of \code{hclust} in the \code{fastcluster} package
#' is much faster than the version of \code{hclust} in \code{stats}.  
#' 
#' In order to account for the fact that the clusters have been estimated from the data, 
#' the p-values are computed conditional on the fact that those clusters were estimated. 
#' This function approximates p-values via importance sampling. 
#' 
#' Currently, this function supports squared Euclidean distance as a measure of dissimilarity 
#' between observations. (Note that complete linkage is invariant under monotone transformations 
#' of the measure of dissimilarity between observations, so unsquared Euclidean distance 
#' would produce the same hierarchical clustering.)
#' 
#' By default, this function assumes that the covariance matrix of the features is isotropic 
#' i.e. \eqn{Cov(X_i) = \sigma^2 I_p}. Setting \code{iso} to false instead assumes that 
#' \eqn{Cov(X_i) = \Sigma}. If known, \eqn{\sigma} can be passed in using the \code{sigma} argument 
#' or \eqn{\Sigma^{-1}} can be passed in the \code{SigInv} argument; otherwise, an 
#' estimate of \eqn{\sigma} or \eqn{\Sigma} will be used. 
#'
#' @export
#'
#' @param X \eqn{n} by \eqn{p} matrix containing numeric data.
#' @param hcl An object of the type \code{hclust} containing the hierarchical clustering of X.  
#' @param K Integer selecting the total number of clusters.
#' @param k1,k2 Integers selecting the clusters to test.
#' @param iso Boolean. If \code{TRUE}, isotropic covariance matrix model, otherwise not.
#' @param sig Optional scalar specifying \eqn{\sigma}, relevant if \code{iso} is \code{TRUE}.
#' @param SigInv Optional matrix specifying \eqn{\Sigma^{-1}}, relevant if \code{iso} is \code{FALSE}.
#' @param ndraws Integer selecting the number of importance samples, default of 2000.
#'
#' @return
#' \item{stat}{the test statistic: the Euclidean distance between the mean of cluster \code{k1} and the mean of cluster \code{k2}  }
#' \item{pval}{the approximate p-value}
#' \item{stderr}{estimated standard error of the p-value estimate}
#'
#' @examples
#' # Simulates a 100 x 2 data set with no clusters
#' set.seed(1)
#' dat <- matrix(rnorm(200), 100, 2)
#'
#' # Complete linkage hierarchical clustering
#' library(fastcluster)
#' hcl <- hclust(dist(dat, method="euclidean")^2, method="complete")
#' 
#' # plot dendrograms with the 1st and 2nd clusters (cut at the third level) 
#' # displayed in blue and orange 
#' plot(hcl)
#' rect_hier_clusters(hcl, k=3, which=1:2, border=c("blue", "orange"))
#' 
#' # Monte Carlo test for a difference in means between the blue and orange clusters
#' test_complete_hier_clusters_approx(X=dat, hcl=hcl, K=3, k1=1, k2=2, ndraws=1000)
#' 
#' @seealso \code{\link{rect_hier_clusters}} for visualizing clusters \code{k1} and \code{k2} in the dendrogram;
#' 
#' \code{\link{test_hier_clusters_exact}} for exact p-values for hierarchical clustering with other linkages;
#' 
#' \code{\link{test_clusters_approx}} for approximate p-values for a user-specified clustering function;
#' 
#' @references Lucy L. Gao et al. "Selective inference for hierarchical clustering". 
test_complete_hier_clusters_approx <- function(X, hcl, K, k1, k2, iso=TRUE, sig=NULL, SigInv=NULL, ndraws=2000) {
  # error checking 
  if(!is.matrix(X)) stop("X should be a matrix")
  
  n <- nrow(X)
  q <- ncol(X)
  
  if(!is_integer_between_a_b(K, 2, n)) stop("number of clusters (K) should be between 2 and n")
  if(!is_integer_between_a_b(k1, 1, K) | !is_integer_between_a_b(k2, 1, K)) stop(paste("cluster indices should be between 1 and K", sep=""))
  if((iso != TRUE) & (iso != FALSE)) stop("iso should be TRUE or FALSE")
  if(!("fastcluster" %in% .packages())) stop("The fastcluster package must be loaded before calling hclust and before calling this function!")
  
  
  hcl_at_K <- stats::cutree(hcl, K)

  k1_obs <- which(hcl_at_K == k1)
  k2_obs <- which(hcl_at_K == k2)
  n1 <- length(k1_obs)
  n2 <- length(k2_obs)
  squared_norm_nu <- 1/n1 + 1/n2
  diff_means <- colMeans(X[k1_obs, , drop=FALSE]) - colMeans(X[k2_obs, , drop=FALSE])

  prop_k2 <- n2/(n1+n2)

  if(iso) {
    if(is.null(sig)) {
      sig <- sqrt(sum(scale(X, scale=FALSE)^2)/(n*q - q))
    }

    scale_factor <- squared_norm_nu*sig^2
    # compute test statistic
    stat <- norm_vec(diff_means)
  } else {
    if(is.null(SigInv)) {
      Sig <- stats::cov(scale(X, scale=FALSE))
      SigInv <- solve(Sig)
    }

    scale_factor <- squared_norm_nu

    # compute test statistic
    stat <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
  }

  scale_factor <- sqrt(scale_factor)
  log_survives <- rep(NA, ndraws)
  phi <- stats::rnorm(ndraws)*scale_factor + stat
  
  k1_constant <- prop_k2*diff_means/stat
  k2_constant <- (prop_k2 - 1)*diff_means/stat
  orig_k1 <- t(X[k1_obs, ])
  orig_k2 <- t(X[k2_obs, ])
  
  Xphi <- X
  
  for(j in 1:ndraws) {
    if(phi[j] < 0) next
    
    # Compute perturbed data set
    phi_minus_stat <- phi[j] - stat
    Xphi[k1_obs, ] <- t(orig_k1 + k1_constant*phi_minus_stat)
    Xphi[k2_obs, ] <- t(orig_k2 + k2_constant*phi_minus_stat)
    
    # Recluster the perturbed data set
    hcl_Xphi <- fastcluster::hclust(stats::dist(Xphi)^2, method="complete")
    clusters_Xphi <- stats::cutree(hcl_Xphi, K)
    if(same_cl(hcl_at_K, clusters_Xphi, K)) {
      log_survives[j] <- -(phi[j]/scale_factor)^2/2 + (q-1)*log(phi[j]/scale_factor) - (q/2 - 1)*log(2) - log(gamma(q/2)) - log(scale_factor) -
        stats::dnorm(phi[j], mean=stat, sd=scale_factor, log=TRUE)
    }
  }

  # Trim down to only survives
  phi <- phi[!is.na(log_survives)]
  log_survives <- log_survives[!is.na(log_survives)]

  survives <- length(log_survives)

  # Return nothing if nothing survives
  if(survives == 0) {
    warning("Oops - we didn't generate any samples that preserved the clusters! Try re-running with a larger value of ndraws.")
    return(list(stat=stat, pval=NA, stderr=NA))
  }
  
  #  Approximate p-values
  log_survives_shift <- log_survives - max(log_survives)
  props <- exp(log_survives_shift)/sum(exp(log_survives_shift))
  pval <- sum(props[phi >= stat])
  
  var_pval <- (1 - pval)^2*sum(props[phi >= stat]^2) + pval^2*sum(props[phi < stat]^2)
  
  return(list(stat=stat, pval=pval, stderr=sqrt(var_pval)))
}

#' Monte Carlo significance test for any clustering method
#'
#' This function performs a user-specified clustering method \code{cl_fun} on the rows of a 
#' data matrix to obtain \code{K} clusters, and tests the null hypothesis of no difference in means 
#' between clusters \code{k1} and \code{k2}. 
#' 
#' In order to account for the fact that the clusters have been estimated from the data, 
#' the p-values are computed conditional on the fact that those clusters were estimated. 
#' This function approximates p-values via importance sampling. 
#' 
#' This function assumes that \code{cl_fun} takes a \eqn{n \times p} numeric data matrix as input 
#' and outputs integer assignments to clusters 1 through \code{K}. 
#' 
#' Thank you to August Guang for providing code to speed-up the function by 
#' parallelizing via the \code{future} package.
#'
#' @export
#'
#' @param X \eqn{n} by \eqn{p} matrix containing numeric data.
#' @param k1,k2 Integers selecting the clusters to test.
#' @param iso Boolean. If \code{TRUE}, isotropic covariance matrix model, otherwise not.
#' @param sig Optional scalar specifying \eqn{\sigma}, relevant if \code{iso} is \code{TRUE}.
#' @param SigInv Optional matrix specifying \eqn{\Sigma^{-1}}, relevant if \code{iso} is \code{FALSE}.
#' @param ndraws Integer selecting the number of importance samples, default of 2000.
#' @param cl_fun Function returning assignments to clusters 1 through \code{K}. 
#' @param cl Optionally pass in the results of calling \code{cl_fun} on your data. This is for
#' efficiency and reproducibility (when the clustering function is non-deterministic). 
#'
#' @return
#' \item{stat}{the test statistic: the Euclidean distance between the mean of cluster \code{k1} and the mean of cluster \code{k2}  }
#' \item{pval}{the approximate p-value}
#' \item{stderr}{standard error of the p-value estimate}
#' \item{clusters}{the estimated cluster assignments}
#' 
#' @examples 
#' # Simulates a 100 x 2 data set with three clusters
#' set.seed(123)
#' dat <- rbind(c(-1, 0), c(0, sqrt(3)), c(1, 0))[rep(1:3, length=100), ] + 
#' matrix(0.2*rnorm(200), 100, 2)
#' 
#' # Function to run k-means clustering w/ k = 3 and 50 random starts
#' km_cluster <- function(X) { 
#'  km <- kmeans(X, 3, nstart=50)
#'  return(km$cluster)
#' }
#' 
#' # Cluster data using k-means 
#' clusters <- km_cluster(dat)
#' table(rep(1:3, length=100), clusters)
#'
#' # tests for a difference in means between clusters 1 and 2
#' # We pass in earlier k-means clustering results from earlier
#' results <- test_clusters_approx(dat, k1=1, k2=2, cl_fun=km_cluster, ndraws=500, cl=clusters)
#' results$stat
#' results$pval
#' results$stderr
#' 
#' @references Lucy L. Gao et al. "Selective inference for hierarchical clustering". 
test_clusters_approx <- function(X, k1, k2, iso=TRUE, sig=NULL, SigInv=NULL, ndraws=2000, cl_fun, cl=NULL) {
  if(!is.matrix(X)) stop("X should be a matrix")
  
  n <- nrow(X)
  q <- ncol(X)
  
  if(is.null(cl)) cl <- cl_fun(X)
  K <- length(unique(cl))
  
  if(!is_integer_between_a_b(K, 2, n)) stop("number of clusters (K) should be between 2 and n")
  if(!is_integer_between_a_b(k1, 1, K) | !is_integer_between_a_b(k2, 1, K)) stop(paste("cluster indices should be between 1 and K", sep=""))
  if((iso != TRUE) & (iso != FALSE)) stop("iso should be TRUE or FALSE")
  
  n1 <- sum(cl == k1)
  n2 <- sum(cl == k2)
  squared_norm_nu <- 1/n1 + 1/n2
  diff_means <- colMeans(X[cl == k1, , drop=FALSE]) - colMeans(X[cl == k2, , drop=F])
  
  prop_k2 <- n2/(n1+n2)
  
  
  if(iso) {
    if(is.null(sig)) {
      sig <- sqrt(sum(scale(X, scale=FALSE)^2)/(n*q - q))
    }
    
    scale_factor <- squared_norm_nu*sig^2
    # compute test statistic
    stat <- norm_vec(diff_means)
  } else {
    if(is.null(SigInv)) {
      Sig <- stats::cov(scale(X, scale=FALSE))
      SigInv <- solve(Sig)
    }

    scale_factor <- squared_norm_nu
    
    # compute test statistic
    stat <- sqrt(as.numeric(t(diff_means)%*%SigInv%*%diff_means))
  }
  
  scale_factor <- sqrt(scale_factor)
  phi <- stats::rnorm(ndraws)*scale_factor + stat
  
  
  k1_constant <- prop_k2*diff_means/stat
  k2_constant <- (prop_k2 - 1)*diff_means/stat
  orig_k1 <- t(X[cl == k1, ])
  orig_k2 <- t(X[cl == k2, ])
  
  Xphi <- X
  
  log_survives <- unlist(future.apply::future_lapply(X = 1:ndraws, FUN = function(j) {
    if(phi[j] < 0) return(NA)
    
    # Compute perturbed data set
    Xphi <- X
    Xphi[cl == k1, ] <- t(orig_k1 + (phi[j] - stat)*k1_constant)
    Xphi[cl == k2, ] <- t(orig_k2 + (phi[j] - stat)*k2_constant)
    
    # Recluster the perturbed data set
    cl_Xphi <- cl_fun(Xphi)
    if(preserve_cl(cl, cl_Xphi, k1, k2)) {
      log_survives <- -(phi[j]/scale_factor)^2/2 + (q-1)*log(phi[j]/scale_factor) - (q/2 - 1)*log(2) - log(gamma(q/2)) - log(scale_factor) -
        stats::dnorm(phi[j], mean=stat, sd=scale_factor, log=TRUE)
      return(log_survives)
    }
    
    return(NA)
    
  }, future.seed=TRUE))
  
  # Trim down to only survives
  phi <- phi[!is.na(log_survives)]
  log_survives <- log_survives[!is.na(log_survives)]
  
  survives <- length(log_survives)
  
  # Return nothing if nothing survives
  if(survives == 0) {
    warning("Oops - we didn't generate any samples that preserved the clusters! Try re-running with a larger value of ndraws.")
    return(list(stat=stat, pval=NA, stderr=NA, clusters=cl))
  }
  
  #  Approximate p-values
  log_survives_shift <- log_survives - max(log_survives)
  props <- exp(log_survives_shift)/sum(exp(log_survives_shift))
  pval <- sum(props[phi >= stat])

  var_pval <- (1 - pval)^2*sum(props[phi >= stat]^2) + pval^2*sum(props[phi < stat]^2)
  
  return(list(stat=stat, pval=pval, stderr=sqrt(var_pval), clusters=cl))
}