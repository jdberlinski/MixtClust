#' Robust Clustering for Complete and Incomplete Data
#'
#' @description Clustering using finite mixture of multivariate t distributions
#'   including handling of incomplete data.
#'
#' @param x A matrix with \eqn{n} observations (rows), \eqn{p} columns
#'   (dimensions), and missing entries set to \code{NA}.
#' @param initial.values Either \code{"emEM"} specifiying the emEM
#'   initialization strategy (see \code{emEM.args} for additional arguments),
#'   \code{"kmeans"} specifiying use of kmeans to generate an initial partition,
#'   a vector of integers specifying an initial partition,
#'   or a named list of initial parameter values (see details).
#' @param nclusters Positive integer. The assumed number of clusters if initial
#'   values are not provided.
#' @param max.iter Positive integer. The maximum number of EM iterations
#'   allowed.
#' @param tol Positive scalar. The desired stopping criterion value.
#' @param convergence Either \code{"lop"} specifying use of relative change
#'   in loglikelihood as the convergence criterion, or \code{"aitkens"} specifying
#'   Aitken's acceleration (default).
#' @param sigma.constr Character specifying constrains on dispersion matrices
#'   for each group. See details for further explanation.
#' @param df.constr Logical. Should the degrees of freedom \eqn{\nu_k} be held
#'   constant over \eqn{k = 1,\dots,K} all clusters?
#' @param approx.df Logical. If \code{approx.df = TRUE}, a numerical
#'   approximation of the objective function used to estimate the degrees of
#'   freedom \eqn{\nu_k} for \eqn{k = 1,\dots,K}.
#' @param method How should missing entries be handled? Must be either
#' \code{"fullEM"} to include missing entries in the EM algorithm following Lin 2009,
#'  \code{"marginalization"} to integrating out missing entries, or \code{"deletion"}
#'  to analyze complete cases only.
#' @param verbose Logical. Should progress be periodically reported to the
#'   screen?
#' @param emEM.args A named list of options utilized if \code{initial.values =
#'   "emEM"} (see details).
#' @param scaled Logical variable that indicates if computations for multi-dimensional datasets should be done after scaling the dataset. Note that the resulting parameters are scaled back and so should not theoretically have much effect on the performance, except to potentially offer stability in numerical computations.
#' @param labels Numeric vector of size \eqn{n} which contains the known labels for the data, to be used in semi-supervised
#'    clustering. Labels of unknown observations should be \code{"NA"}.
#'
#' @details
#'
#' Model-based clustering using finite mixtures of t distributions, with
#' handling of incomplete data using either marginalization or the EM algorithm.
#' If supplying initial values, format as a named list with elements: \itemize{
#' \item{"pi"}{ Mixing proportions. A vector of length \eqn{K} that sums to
#' one.} \item{"nu"}{ Degrees of freedom. A vector of length \eqn{K} with
#' entries at least equal to three (thus requiring the existance of the first
#' two moments.)} \item{"mu"}{ Locations. A \eqn{K \times p} matrix, where the
#' \eqn{k}-th row is the location \eqn{\mu_k \in R^p} for cluster \eqn{k}.}
#' \item{"Sigma"}{ Dispersions. A \eqn{p \times p \times K} array, where the
#' \eqn{k}-th slice is the \eqn{p \times p} positive-definite dispersion matrix
#' \eqn{\Sigma_k} for cluster \eqn{k}.} }
#'
#' The arguments for emEM specified in the list \code{emEM.args} are: \itemize{
#' \item{"nstarts"}{ Positive integer. The number of randomly generated initial
#' starting parameter values under consideration.} \item{"em.iter"}{ Positive
#' integer. The number of short EM iterations to be performed on each set of
#' initial starting parameter values.} \item{"nbest"}{ Positive integer. After
#' \code{em.iter} EM iterations are performed in each of the \code{nstarts}
#' initial values, the number of top ranking (according to loglikelihood)
#' parameter values on which to run the long EM either to convergence (specified
#' by \code{tol}) or maximum number of iterations (specified by
#' \code{max.iter}). If \code{nbest} is greater than one, the long EM run
#' achieving the largest loglikelihood will be returned.}
#'
#' Note that because of how the 'marginalization' algorithm is derived, it is
#' recommended to use em.iter = 1,  this results in the Rnd-EM algorithm of
#' Maitra (2009) when using 'marginalization'
#' }
#'
#' Constraints on the dispersion matrices for each group are specified according
#' to the decomposition of the dispersion matrices \deqn{\Sigma_k = \zeta_k \Gamma_k \Lambda_k \Gamma_k^\top.}
#' The nomenclature is a three letter string following that of other popular packages, such as
#' \code{mclust}. The three letters specify constrains on the volume
#' (\eqn{\zeta_k}), shape (\eqn{\Lambda_k}), and orientation (\eqn{\Gamma_k}),
#' respectively. Each letter is either "V" or "E", specifying that the
#' corresponding component is allowed to vary or must remain equal across
#' groups, respectively. Additional specifications include "VII", "EII", "VVI",
#' "VEI", "EVI", and "EEI", where the "I" indicates the corresponding component
#' is an identity matrix of appropriate dimensions.
#'
#' @references Emily Goren & Ranjan Maitra, 2021.
#' "Model-based clustering of partial records," arXiv:2103.16336.
#' https://arxiv.org/abs/2103.16336
#'
#' Ranjan Maitra, 2009. "Initializing optimization partitioning algorithms," ACM/IEEE Transactions on Computational Biology and Bioinformatics, 6:1:144–157.
#'
#' Tsung-I Lin & Hsiu Ho & Pao Shen, 2009. "Computationally
#' efficient learning of multivariate t mixture models with missing
#' information," Computational Statistics, 24(3): 375-392.
#'#'
#' @return A list containing: \itemize{ \item{"estimates"}{ A list of the final
#'   estimates "pi", "nu", "mu", and "Sigma" containing the MLEs for the mixing
#'   proportions, degrees of freedom, locations, and dispersions, respectively.}
#'   \item{"iterations"}{ Number of EM iterations performed (long EM run only;
#'   if emEM was performed, this excludes the short em run iterations specified
#'   in emEM.args$em.iter).} \item{"Zs"}{ A \eqn{n \times K} matrix where the
#'   \eqn{i}-th row contains the posterior probabilities of membership in
#'   cluster \eqn{1, \dots, K} for the \eqn{i}-th observation
#'   (\eqn{i=1,\dots,n}).} \item{"class"}{ A vector of length \eqn{n} with the
#'   predicted class memberships for each observation.} \item{"loglik"}{ The
#'   log likelihood at each (long EM run) iteration.} \item{"loglik"}{The
#'   log likelihood at the last iteration, computed for all cases (including
#'   those with missing values when \code{method = "deletion"})}
#'   \item{"bic"}{ The BIC for
#'   the final fitted model.} \item{"EM.time"}{ Runtime for the long EM run(s).}
#'   \item{"em.time"}{ Runtime for the short em run(s) when \code{initial.values
#'   = "emEM"}.} \item{"total.time"}{ Runtime for the entire function call.}
#'   \item{"call"}{ Supplied function call.} \code{npar}{The number of model parameters.}}
#'
#' @examples
#' set.seed(20180626)
#' # Use iris data.
#' d <- subset(iris, select = -Species)
#' # Create missing data -- MCAR with 10% chance of missingness.
#' missing <- matrix(rbinom(n = ncol(d)*nrow(d), size = 1, prob = 0.1), ncol = ncol(d))
#' x <- d; x[missing == 1] <- NA
#' # Run EM with emEM initialization strategy for candidate clusters K = 2, 3, 4.
#' Ks <- 2:4
#' ans <- lapply(Ks, function(K) {
#'     MixtClust(x, nclusters = K, emEM.args = list(nstarts=K*10, em.iter=5, nbest=1))
#' })
#' # Get BIC for each K.
#' BICs <- sapply(ans, function(f) f$bic)
#' # Plot BIC by K.
#' plot(BICs ~ Ks, pch = 20, xlab = 'Number of Clusters', ylab = 'BIC')
#'
#' @author Emily Goren, \email{emily.goren@gmail.com}
#'
#' @export
#'
MixtClust <- function(x,
                      initial.values = "emEM",
                      nclusters = NULL,
                      max.iter = 1e3,
                      tol = 1e-3,
                      convergence = "aitkens",
                      sigma.constr = "VVV",
                      df.constr = FALSE,
                      approx.df = TRUE,
                      method = "marginalization",
                      verbose = TRUE,
                      scaled = FALSE,
                      labels = rep(NA, nrow(x)),
                      emEM.args = list(nstarts = round(sqrt(nclusters*prod(dim(x)))), em.iter = 1, nbest = 10))
{
  mf <- match.call(expand.dots = FALSE)
  #######################################################
  # Check supplied arguments
  #######################################################
  if (max.iter < 1 | ceiling(max.iter) != max.iter | !is.numeric(max.iter))
    stop("Please supply a positive integer value for max.iter")
  if (max.iter < 2 & convergence == 'aitkens')
    stop("max.iter must be at least 2 in order to implemement Aitken's convergence criterion")
  if (!(convergence %in% c("lop", "aitkens")))
    stop("Supplied convergence not recognized, must be either lop or aitkens")
  if (tol <= 0 | !is.numeric(tol))
    stop("Please supply a positive numeric value for tol")
  if (!is.logical(df.constr))
    stop("Please supply a logical value for df.constr")
  if (!is.logical(approx.df))
    stop("Please supply a logical value for approx.df")
  if (!(method %in% c("fullEM", "marginalization", "deletion")))
    stop("method must be one of fullEM, marginalization, or deletion")
  marginalization <- ifelse(method == "marginalization", TRUE, FALSE)
  if (is.null(x))
    stop('No data was supplied!')

  valid_constr <- c("VVV", "EEE", "VII", "EII", "EEI", "VVI", "EVI",
                    "EVV", "VEE", "VEV", "VEI", "EVE", "EEV", "VVE")
  if (is.character(sigma.constr))
    sigma.constr <- toupper(sigma.constr)
  if (!(sigma.constr %in% valid_constr)) {
    if (is.logical(sigma.constr)) {
      new_constr <- ifelse(sigma.constr, "EEE", "VVV")
      warning(paste0("sigma.constr = ", sigma.constr, " is depreciated. Use one of: ",
          paste(valid_constr, collapse = ", "), " instead. \n", "Converting sigma.constr = ",
          sigma.constr, " to sigma.constr = ", new_constr))
      sigma.constr <- new_constr
    } else {
      stop("sigma.constr must be one of: ", paste(valid_constr, collapse = ", "))
    }
  }

  if (length(labels) != nrow(x))
    stop("length of `labels` does not match number of rows in data")



  #######################################################
  # Set up data
  #######################################################
  x <- as.matrix(x)
  R <- is.na(x)
  # Remove observations with no observed coordinates and coordinates with no observations.
  bad.rows <- rowSums(R) == ncol(x)
  bad.cols <- colSums(R) == nrow(x)

  if (scaled) {
      # if scaled computations are desired, make sure that the data are scaled for each coordinate to have sd 1.
      x <- scale(x, center = FALSE, scale = T)
      v <- attr(x,"scaled:scale")
  }

  CC <- (rowSums(R) == 0)
  if (any(bad.rows) & (verbose))
    message(paste("Removing rows with no observations:", which(bad.rows)))
  if (any(bad.cols) & (verbose))
    message(paste("Removing columns with no observations:", which(bad.cols)))
  if (method == "deletion") { # Also remove non-complete cases
    Y <- X <- x[CC, !bad.cols]
    R <- R[CC, !bad.cols]
  } else {
    Y <- X <- x[!bad.rows, !bad.cols]
    R <- R[!bad.rows, !bad.cols]
  }
  n <- nrow(X); p <- ncol(X)
  # Missingness indicator matrix.
  A <- matrix(as.numeric(!R), ncol = p, nrow = n)
  ps <- rowSums(A)
  # Set missing values to 0 for cpp code. --------------------------------------
  X[R] <- 0
  # Unique patterns of missingness.
  R.unique <- unique(R, MARGIN = 1)
  # Missingness pattern of each observation.
  miss.grp <- apply(R, 1, row_match, R.unique)
  Ru <- 1*!R.unique

  # Set class indicators
  # recode labels to start at 1 and increase sequentially
  labels <- as.numeric(as.factor(labels))
  class_indicators <- matrix(0, nrow = nrow(x), ncol = nclusters)
  for (i in 1:nclusters) class_indicators[labels == i, i] <- 1

  labeled_observation <- !is.na(labels)

  #######################################################
  # Initialize
  #######################################################
  parnames <- c("pi", "nu", "mu", "Sigma")
  if (initial.values[1] == "emEM") {
    emEM <- TRUE
    nstarts <- emEM.args$nstarts
    em.iter <- emEM.args$em.iter
    nbest <- emEM.args$nbest
    if (nbest > nstarts) {
      message("Number of em runs is less than number of EM runs, increasing nstarts to nbest.")
      nstarts <- nbest
    }
    if (is.null(nclusters))
      stop("Please provide the number of clusters")
    if (nclusters < 1 | ceiling(nclusters) != nclusters)
      stop("Please provide a positive integer value for the number of clusters")
    if (verbose)
      message(paste0('Generating initial values and running short em assuming ', nclusters, ' clusters ...'))
    ptm <- proc.time()
    shortEMest <- rep(list(NA), nstarts)
    shortEMll <- rep(NA, nstarts)
    for (st in 1:nstarts) {
      tmp <- run.em(
        nclusters, X,
        miss.grp, A, R, Ru, ps,
        em.iter, sigma.constr, df.constr, marginalization,
        init = "smart-random", labeled_observation, class_indicators)
      shortEMll[st] <- tmp$loglik
      # save estimates if candidates for nbest logliks
      if (st <= nbest) {
        shortEMest[[st]] <- tmp$estimates
        cutoff <- min(shortEMll, na.rm = TRUE)
      } else {
        if (tmp$loglik > cutoff) {
          shortEMest[[st]] <- tmp$estimates
          cutoff <- sort(shortEMll, decreasing = TRUE)[nbest]
        }
      }
    }
    em.time <- proc.time() - ptm
    if (verbose)
      message(paste0('Determining best ', nbest, ' values from ', nstarts, ' short em runs assuming ', nclusters, ' clusters...'))
    best <- order(shortEMll, decreasing = TRUE)[1:nbest]
    initial.values <- shortEMest[best]
  } else if (initial.values[1] == "kmeans") {
    emEM <- FALSE
    initial.values <- run.em(
      nclusters, X,
      miss.grp, A, R, Ru, ps, em.iter,
      sigma.constr, df.constr,
      marginalization, init = "kmeans", labeled_observation, class_indicators)
    em.time <- NA
    initial.values <- list(initial.values$estimates)
  } else if (is(initial.values, "list")) {
    emEM <- FALSE
    # check supplied values
    if (!(all(names(initial.values) %in% parnames & parnames %in% names(initial.values))))
      stop("initial.values must be a list with named elements 'pi', 'nu', 'mu', and 'Sigma'")
    nclustersold <- nclusters
    nclusters <- length(initial.values$pi)
    if (nclustersold != nclusters)
      message(paste0('Length of supplied initial.values$pi is ', nclustersold, '. Setting number of clusters to ', 'K = ', nclusters))
    if (length(initial.values$nu) != nclusters)
      stop(paste0('Length of supplied initial.values$nu must be ', nclusters, '.'))
    if (ncol(initial.values$mu) != p | nrow(initial.values$mu) != nclusters)
      stop(paste0('Dimension of supplied initial.values$mu must be ', nclusters, ' by ', p))
    if (!all(dim(initial.values$Sigma) == c(p,p,nclusters)))
      stop(paste0('Dimension of supplied initial.values$Sigma must be ', p, ' by ', p, ' by ', nclusters))
    em.time <- NA
    initial.values <- list(initial.values)
  } else {
    if (!all.equal(labels[labeled_observation], initial.values[labeled_observation]))
      stop("supplied `labels` do not agree with starting values `initial.values`")
    emEM <- FALSE
    init.ids <- initial.values
    if (length(init.ids) != nrow(X))
      stop("initial.values are detected to be the initial partition, but is not of length n")
    nclustersold <- nclusters
    nclusters <- length(unique(init.ids))
    if (nclustersold != nclusters)
      message(paste0('Length of supplied initial.values$pi is ', nclustersold, '. Setting number of clusters to ', 'K = ', nclusters))
    if (!all(init.ids %in% 1:nclusters))
      stop("initial.values are detected to be the initial partition, but contains values not in 1:nclusters")
    zz <- sapply(1:nclusters, function(k) as.numeric(init.ids == k))
    initial.values <- list(get.init.val(X, R, nclusters, df.constr, sigma.constr, "ids", zz))
    em.time <- NA
  }

  #######################################################
  # Run long EM
  #######################################################
  # Number of parameters = |pi|-1 + |mu| + |nu| + |Sigma|
  sigp <- (p*(p+1))/2
  # npar <- (nclusters-1) + nclusters*p + ifelse(df.constr, 1, nclusters) + ifelse(sigma.constr, sigp, nclusters*sigp)

  # assuming K == nclusters
  # Volume: K parameters for varying, 1 for equal
  # Shape: K*(p - 1) for varying, p - 1 for equal (one is left to match desired volume)
  # Orientation: K*p*(p-1)/2 for varying, p*(p-1)/2 for equal number of parameters in an orthogonal matrix
  npar <- (nclusters - 1) + nclusters*p + ifelse(df.constr, 1, nclusters)
  if (sigma.constr == "VVV") {
    npar <- npar + nclusters * (p * (p + 1)) / 2
  } else if (sigma.constr == "EEE") {
    npar <- npar + (p * (p + 1)) / 2
  } else if (sigma.constr == "VII") {
    npar <- npar + nclusters
  } else if (sigma.constr == "EII") {
    npar <- npar + 1
  } else if (sigma.constr == "EEI") {
    npar <- npar + p
  } else if (sigma.constr == "VVI") {
    npar <- npar + p*nclusters
  } else if (sigma.constr == "EVI") {
    npar <- npar + (p - 1)*nclusters + 1
  } else if (sigma.constr == "EVV") {
    npar <- npar + nclusters * p * (p - 1) / 2 + nclusters * (p - 1) + 1
  } else if (sigma.constr == "VEE") {
    npar <- npar + p * (p + 1) / 2 + nclusters - 1
  } else if (sigma.constr == "VEV") {
    npar <- npar + nclusters * p * (p - 1) / 2 + p + nclusters - 1
  } else if (sigma.constr == "VEI") {
    npar <- npar + p + nclusters - 1
  } else if (sigma.constr == "EVE") {
    npar <- npar + nclusters * (p - 1) + 1 + p * (p - 1) / 2
  } else if (sigma.constr == "EEV") {
    npar <- npar + p + nclusters * p * (p - 1) / 2
  } else if (sigma.constr == "VVE") {
    npar <- npar + nclusters * p + p * (p - 1) / 2
  }

  multres <- rep(list(NA), length(initial.values))
  ptm <- proc.time()
  for (i in 1:length(initial.values)) {
    if (verbose & emEM)
      message(paste0('Running long EM from best em runs: ', i, ' of ', nbest, ' total assuming ', nclusters, ' clusters...'))
    multres[[i]] <- run.EM(
      initial.values[[i]], nclusters, X,
      miss.grp, A, Ru, ps,
      max.iter, tol, convergence,
      sigma.constr, df.constr,
      approx.df, marginalization, npar, labeled_observation, class_indicators)
  }
  EM.time <- proc.time() - ptm

  #######################################################
  # Output.
  # If more than one long EM run, only return best one.
  #######################################################
  if (verbose)
    message(paste0('Finished EM assuming ', nclusters, ' clusters. \n\n'))
  best <- sapply(multres, function(out) tail(out$loglik, n = 1L))
  o <- multres[[which.max(best)]]
  o$class <- rep(NA, ncol(x))
  if (method == "deletion") {
    idx <- (!bad.rows & !CC)
    if (nclusters == 1) {
      o$class[!bad.rows] <- 1
    } else {
      o$class[CC] <- apply(o$Zs, 1, which.max)
      o$class[idx] <- apply(x[idx, !bad.cols], 1, classify, params = o$estimates)
    }
    o$loglikn <- tail(o$loglik, n=1L)
    if (sum(idx) > 0)
      o$loglikn = o$loglikn + loglikelihood(x[idx, !bad.cols], o$estimates)
  } else {
    o$class[!bad.rows] <- apply(o$Zs, 1, which.max)
    o$loglikn <- tail(o$loglik, n=1L)
  }
  if (scaled) {
      ## scale the estimates back
      o$estimates$mu <- t(v * t(o$estimates$mu))
      for (k in 1:nclusters) {
          o$estimates$Sigma[ , , k]  <- diag(v) %*%  o$estimates$Sigma[ , , k] %*% diag(v)
      }
  }
  o$EM.time <- EM.time
  o$em.time <- em.time
  o$total.time <- EM.time + em.time
  o$call <- mf
  o$npar <- npar
  o$sigma.constr <- sigma.constr
  o
}
