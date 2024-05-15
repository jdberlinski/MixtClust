## Row match: return row numbers of matrix that have rows equal to row.
row_match <- function(row, matrix) {
    ans <- apply(matrix, 1, function(r) all.equal(row, r) == TRUE)
    return(which(ans == TRUE))
}

## Perform one EM iteration (of all ECM steps).
EM_iter <- function(oldpars, x, A, Ru, miss.grp, ps, sigma.constr, df.constr, approx.df, marginalization, labeled_obs, class_indicators) {
    k <- length(oldpars$pi)
    ## 1st cycle -- update cluster proportions, means, dfs
    ## E-step followed by CM-step for pi
    ##    print(are.sigmas.valid(oldpars$Sigma))
    w <- up_W(x, oldpars$mu, oldpars$Sigma, oldpars$nu, miss.grp, Ru)
    if (k > 1) {
        z <- up_Z(x, oldpars$mu, oldpars$Sigma, oldpars$nu, oldpars$pi, miss.grp, Ru, labeled_obs, class_indicators)
        pisnew <- as.numeric(up_pi(z))
    } else {
        z <- matrix(1, nrow(x), k)
        pisnew <- 1
    }
    if (marginalization) {
        ## update locations
        musnew <- up_mu(x, z, w, A)
    } else {
        K <- ncol(z)
        M <- length(unique(miss.grp))
        p <- ncol(x)
        ## EM for missing components
        SOiOEOO <- lapply(1:K, function(k) SOiOEOOk(oldpars$Sigma[,,k], Ru))
        xhat <- lapply(1:K, function(k) xhatk(x, oldpars$mu[k,], miss.grp, M, SOiOEOO[[k]]))
                                        # update locations
        musnew <- up_mu_Lin(p, z, w, xhat)
    }
    ## update nu's
    nusnew <- as.numeric(up_nu(z, w, oldpars$nu, ps, df.constr, approx.df))
    bidx <- !is.finite(nusnew); nusnew[bidx] <- oldpars$nu[bidx] # fix any NaN
    ## 2nd cycle -- update dispersions
    if (k > 1)
        z <- up_Z(x, musnew, oldpars$Sigma, nusnew, pisnew, miss.grp, Ru, labeled_obs, class_indicators)
    w <- up_W(x, musnew, oldpars$Sigma, nusnew, miss.grp, Ru)
    if (marginalization) {
        Sigmasnew <- up_Sigma(x, z, w, musnew, A, sigma.constr)
    } else {
        xhat_l <- array(unlist(xhat), dim = c(nrow(xhat[[1]]), ncol(xhat[[1]]), length(xhat)))
        Sigmasnew <- up_Sigma_Lin(z, w, musnew, oldpars$Sigma, xhat_l, miss.grp, sigma.constr, Ru)
    }
    ## Output
    newpars <- list(pisnew, nusnew, musnew, Sigmasnew)
    names(newpars) <- c('pi', 'nu', 'mu', 'Sigma')
    return(newpars)
}

## Function to generate a set of initial parameter values.
get.init.val <- function(X, R, K, df.constr, sigma.constr, init = "smart-random", labeled_obs = NULL, class_indicators = NULL, Z = NULL) {
    ##if (df.constr) nus <- rep(runif(1, 5, 25), K) else nus <- runif(K, 5, 25)
    ## Follow teigen: set dfstart to 50
    ##nus <- runif(K, min = 10, max = 50)
    nus <- rep(50, K)
    if (df.constr) nus <- rep(nus[1], K)
    n <- nrow(X); p <- ncol(X)

    if (init == "smart-random")
        old.inits <- FALSE else
                               old.inits  <- TRUE

    if (any(labeled_obs)) {
      if (old.inits) stop("only emEM is supported in the semi-supervised case")
      y_obs <- X[labeled_obs, ]
      R_obs <- R[labeled_obs, ]
      y_obs[R_obs] <- NA
      y_unobs <- X[!labeled_obs, ]
      R_unobs <- R[!labeled_obs, ]
      y_unobs[R_unobs] <- NA

      obs_class <- apply(class_indicators[labeled_obs, ], 1, which.max)
      M <- max(obs_class)

      Sigmas <- array(0, dim=c(p, p, K))

      # first, initialize the observed classes
      Sigmas_obs <- array(0, dim = c(p, p, M))
      mus_obs <- matrix(0, nrow = M, ncol = p)
      nks_obs <- numeric(M)
      pis_obs <- numeric(M)
      for (k in 1:M) {
        # rare case: single observation in a labeled cluster
        if (sum(obs_class == k) > 1)
          Sigmas_obs[, , k] <- cov(y_obs[obs_class == k, ], use = "pairwise.complete.obs")
        Sigmas_obs[, , k] <- Sigmas_obs[, , k] + 1e-3 * diag(p)
        mus_obs[k, ] <- colMeans(y_obs[obs_class == k, ], na.rm = T)
        nks_obs[k] <- sum(obs_class == k)
        pis_obs[k] <- sum(obs_class == k) / n
      }

      # now, initialize the remaining classes using the cases with unobserved labels
      Sigmas_unobs <- array(0, dim = c(p, p, K - M))
      if ((K - 1) == M) {
        Sigmas_unobs[, , 1] <- cov(y_unobs, use = "pairwise.complete.obs")
        mus_unobs <- matrix(colMeans(y_unobs, na.rm = T), nrow = 1)
        nks_unobs <- nrow(y_unobs)
        pis_unobs <- nrow(y_unobs) / n
      } else if (K != M) {
        minnks <- 0
        while (minnks <= p) {
          ## Let us at least put in at least p + 1 observation pairs in each group
          res <- kmmeans::kmmeans(data = as.data.frame(y_unobs), K = K - M, n.init = 1, kmmns.iter = 0)
          res$partition <- res$partition + 1
          minks <- sapply(1:(K - M), function(x)(min(crossprod(!R_unobs[res$partition==x,]))))
          minnks <- min(minks)
        }

        nks_unobs <- table(res$partition)
        pis_unobs <- nks_unobs / n
        mus_unobs <- res$centers
        for (k in 1:(K - M)) {
          Sigmas_unobs[,,k] <- cov(y_unobs[res$partition==k,], use = "pairwise.complete.obs")
          Sigmas_unobs[,,k]  <- Sigmas_unobs[,,k] + 1e-3 * diag(p)
        }
      } else {
        # if K == M, we don't need to consider the unobserved cases for initialization
        nks_unobs <- NULL
        pis_unobs <- NULL
        mus_unobs <- NULL
      }

      # combine results
      for (k in 1:K) {
        if (k <= M)
          Sigmas[, , k] <- Sigmas_obs[, , k]
        else
          Sigmas[, , k] <- Sigmas_unobs[, , k - M]
      }

      mus <- rbind(mus_obs, mus_unobs)
      nks <- c(nks_unobs, nks_obs)
      pis <- c(pis_unobs, pis_obs)
    }

    else if (!old.inits) {
        ## smart random initialization chooses random points in sequence from the dataset with probability inversely proportional to the distance from (the closest of) the previous points. We use the kmmeans with 0 iterations to achieve this. Note that this function can handle missing data and so does not particularly care over missing cases, and can handle them fine.
        ## because we use pairwise complete observations to calculate the estimates of the Sigmas, it is possible that these are singular. So, instead of checking for stability as we were doing earlier, we simply add a small value to the doagonal (0.001 times the identity matrix) to make the matrix more non-singular).
        y <- X
        y[R] <- NA # recreate the data matrix with NAs for call to kmmeans
        Sigmas <- array(0, dim=c(p, p, K))
        if (K == 1) {
            Sigmas[,,1]  <- cov(y, use = "pairwise.complete.obs")
            pis <- 1
            mus <- matrix(colMeans(y,na.rm=T), nrow = 1)
        } else {
          # TODO: This should be done only on the unknown observations, with the number of clusters initialized being
          # the difference between the number of observed clusters and the total number of clusters
          #   - What if these are equal, and there are less than p-1 observations in one of the observed clusters?
          #   - More speficically, if I have M groups in the known labels, what if I have less than (K - M)(p + 1)
          #   unlabeled observations?
            ##            sigmas.stab <- F
            ##            grand.iter  <- 0
            ##            while ((!sigmas.stab) & (grand.iter < 5)) {
            minnks <- 0
            while (minnks <= p) {
                ## Let us at least put in at least p + 1 observation pairs in each group
                ## cat("minnks = ", minnks, "\n")
                res <- kmmeans::kmmeans(data = as.data.frame(y), K = K, n.init = 1, kmmns.iter = 0)
                res$partition <- res$partition + 1
                minks <- sapply(1:K, function(x)(min(crossprod(!R[res$partition==x,]))))
                minnks <- min(minks)
                ## cat("minnks = ", minnks, "\n")
            }

            nks <- table(res$partition)
            pis <- nks/n
            mus <- res$centers
            for (k in 1:K) {
                Sigmas[,,k] <- cov(y[res$partition==k,], use = "pairwise.complete.obs")
                Sigmas[,,k]  <- Sigmas[,,k] + 1e-3 * diag(p)
            }
        }
    } else     {
        minstable <- p/n
        CC <- rowSums(R) == 0
        ## Uniform initialization of z.
        if (K > 1 & is.null(Z)) {
            nonstable <- TRUE
            while(nonstable) {
                Z <- switch(init,
                            uniform = {
                                matrix(runif(n*K, 0, 1), n, K)
                            },
                            kmeans = {
                                id <- kmeans(X[CC, ], nstart = 1e3, centers = K)
                                sapply(1:K, function(k) as.numeric(id$cluster == k))
                            })
                Z <- Z/rowSums(Z)
                pis <- colMeans(Z)
                if (all(pis > minstable)) nonstable <- FALSE
            }
        } else if (K > 1 & !is.null(Z)) {
            Z <- Z/rowSums(Z)
            pis <- colMeans(Z)
        } else {
            Z <- matrix(1, n, K)
            pis <- 1
        }
        ## Use updated IDs to set remaining parameters.
        mus <- matrix(0, K, p)
        Sigmas <- array(0, dim=c(p, p, K))
        for (k in 1:K) {
            wtk <- switch(init,
                          uniform = {
                              Z[,k] * CC # set wt to zero if not complete case
                          },
                          kmeans = {
                              out <- rep(0, n)
                              out[CC] <- Z[,k]
                              out
                          },
                          ids = {
                              Z[,k] * CC
                          })
            wtcov <- cov.wt(X, wt = wtk, method = "ML")
            mus[k,] <- wtcov$center
            Sigmas[,,k] <- wtcov$cov
        }
    }
    if ((K > 1) & sigma.constr == "EEE") {
        wtdSigmas <- lapply(1:K, function(k) pis[k]*Sigmas[,,k])
        S <- Reduce('+', wtdSigmas)
        for (k in 1:K) Sigmas[,,k] <- S
    }
    out <- list(pi = pis, nu = nus, mu = mus, Sigma = Sigmas)
    out
}

## Run EM.
run.EM <- function(init, nclusters, X, miss.grp, A, Ru, ps, max.iter, tol, convergence, sigma.constr, df.constr, approx.df, marginalization, npar, labeled_obs, class_indicators) {
    ##print("entering run.EM")
    old <- init
    del <- 1e6 # Initialize convergence holder.
    iter <- 0
    initL <- sapply(1:nclusters, function(k) {old$pi[k] * h(X, old$mu[k,], old$Sigma[,,k], old$nu[k], miss.grp, Ru)})
    LLs <- c(sum(log(rowSums(initL))), rep(NA, max.iter)) # Store loglikelihood at each iteration.
    while (del > tol) {
        iter <- iter + 1
        new <- EM_iter(old, X, A, Ru, miss.grp, ps, sigma.constr, df.constr, approx.df, marginalization, labeled_obs, class_indicators)
        newL <- sapply(1:nclusters, function(k) {new$pi[k] * h(X, new$mu[k,], new$Sigma[,,k], new$nu[k], miss.grp, Ru)})
        newLLn <- sum(log(rowSums(newL)))
        ##    cat("ITER =", iter, "newLLn" = newLLn, "\n")
        LLs[iter+1] <- newLLn
        old <- new
                                        # Stop if loglik didn't change or went down
        if (newLLn <= LLs[iter]) break
                                        # Otherwise check specificed convergence critrion
        if (convergence == 'lop') {
            del <- ((LLs[iter] - newLLn) / LLs[iter]) # Relative change in loglikelihoods.
        } else if (convergence == 'aitkens') {
            if (iter > 1 ) {
                a <- (newLLn - LLs[iter]) / (LLs[iter] - LLs[iter - 1])
                LLinf <- LLs[iter] + ((newLLn - LLs[iter]) / (1 - a))
                del <- abs(LLinf - newLLn)
            } else {
                del <- 1e6
            }
        } else {
            stop('Specified convergence criterion not recognized.')
        }
        if (iter == max.iter) break
    }
                                        # Final posterior probabilities of cluster membership
    Zs <- up_Z(X, new$mu, new$Sigma, new$nu, new$pi, miss.grp, Ru, labeled_obs, class_indicators)
    BIC <- 2*LLs[iter+1] - npar*log(nrow(X))
    res <- list(estimates = new,
                iterations = iter,
                Zs = Zs,
                loglik = LLs[2:(iter+1)],
                bic = BIC)
#    print("exiting run.EM")
#    print(res)
    return(res)
}

## Run an initialization with short em -- don't keep track of LL, BIC, etc to save time
run.em <- function(nclusters, X, miss.grp, A, R, Ru, ps, niter, sigma.constr, df.constr, marginalization, init = "uniform", labeled_obs, class_indicators) {
  old <- get.init.val(X, R, nclusters, df.constr, sigma.constr, init, labeled_obs, class_indicators)
  iter <- 0
  while (iter <= niter) {
    iter <- iter + 1
    new <- EM_iter(old, X, A, Ru, miss.grp, ps, sigma.constr, df.constr, approx.df = TRUE, marginalization, labeled_obs, class_indicators)
    old <- new
  }
  # Final loglik
  L <- sapply(1:nclusters, function(k) {new$pi[k] * h(X, new$mu[k,], new$Sigma[,,k], new$nu[k], miss.grp, Ru)})
  res <- list(estimates = new, loglik = sum(log(rowSums(L))))
  return(res)
}
##
## check for positive definiteness of the Sigmas
##

are.sigmas.valid <- function(Sigma, eps = 1e-4)
{
    all(sapply(1:dim(Sigma)[3], function(k) {
        if ((max(diag(Sigma[,,k])) > eps*0.01)  & (min(diag(Sigma[,,k]) > 0)))
            eigen(cov2cor(Sigma[,,k]),symmetric = TRUE, only.values = TRUE)$values[dim(Sigma)[1]] > eps else FALSE
    }))
}

## Update Sigma -- R version used for debugging only
up_SigmaR <- function(x, z, w, mus, A, constr = FALSE) {
    k <- ncol(z)
    p <- ncol(x)
    n <- nrow(x)
    S <- sapply(1:k, function(kk) {
        muk <- mus[kk,]
        x.c <- sweep(x, 2, muk)
        Lk <- lapply(1:n, function(i) {
            Acent <- as.numeric(diag(A[i,]) %*% x.c[i,])
            return(z[i,kk] * w[i,kk] * Acent %o% Acent)
        })
        Lks <- Reduce("+", Lk)
        Rk <- lapply(1:n, function(i) {
            return((z[i,kk]) * A[i,] %o% A[i,])
        })
        Rks <- Reduce("+", Rk)
        return(Lks / Rks)
    }, simplify = 'array')
    if (constr) stop('Dispersion constraints not yet implemented')
    return(S)
}

                                        # Update mu -- R version used for debugging only
up_muR <- function(x, z, w, A) {
    k <- ncol(z)
    p <- ncol(x)
    n <- nrow(x)
    mu <- t(sapply(1:k, function(kk) {
        Wk <- lapply(1:n, function(i) z[i,kk] * w[i, kk] * diag(A[i,]))
        Lks <- solve(Reduce("+", Wk))
        Rk <- lapply(1:n, function(i) Wk[[i]] %*% x[i,])
        Rks <- Reduce("+", Rk)
        return(as.numeric(Lks %*% Rks))
    }))
    return(mu)
}

                                        # Loglik
loglikelihood <- function(x, ests) {
    X <- as.matrix(x)
    R <- is.na(X)
                                        # Remove observations with no observed coordinates and coordinates with no observations.
    bad.rows <- rowSums(R) == ncol(X)
    bad.cols <- colSums(R) == nrow(X)
    X[R] <- 0
    n <- nrow(X); p <- ncol(X)
                                        # Set missing values to 0 for cpp code.
                                        # Unique patterns of missingness.
    R.unique <- unique(R, MARGIN = 1)
                                        # Missingness pattern of each observation.
    miss.grp <- apply(R, 1, row_match, R.unique)
    Ru <- 1*!R.unique
    nclusters <- length(ests$pi)

    L <- sapply(1:nclusters, function(k) {ests$pi[k] * h(X, ests$mu[k,], ests$Sigma[,,k], ests$nu[k], miss.grp, Ru)})
    LL <- sum(log(rowSums(L)))
    return(LL)
}


                                        #h.try <- function(x, mu, Sigma, nu, miss.grp, Ru) {
                                        #    try.h <- try(h(x, mu, Sigma, nu, miss.grp, Ru))
                                        #    print(try.h)
                                        #    print(length(class(try.h))) # appears to be matrix(array)
                                        #    print(length(class(try.h)))
                                        #    if (length(class(try.h)==1)) {
                                        #        if (class(try.h) == "try-error")
                                        #            browser()
                                        #    }  else
                                        #        try.h
                                        #    h(x, mu, Sigma, nu, miss.grp, Ru)
                                        #}

