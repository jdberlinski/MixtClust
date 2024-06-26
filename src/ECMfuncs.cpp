#define ARMA_DONT_PRINT_ERRORS
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// Make a vector into an array (in Armadillo speak, a cube)
// [[Rcpp::export]]
arma::cube to_array(NumericVector sigmas, int K, int p) {
  NumericVector vSigmas(sigmas);
  arma::cube Sigmas(vSigmas.begin(), p, p, K, false);
  return Sigmas;
}


// Make a symmetric singular matrix non-singular by adding epsilon(tol) to the
// zero eigenvalues.
// [[Rcpp::export]]
arma::mat fix_var(arma::mat sigma, double tol = 1e-3) {
  int p = sigma.n_rows;
  arma::mat eigvec(p,p), ans(p,p);
  arma::vec eigval(p), failures(p); failures.zeros();
  arma::eig_sym(eigval, eigvec, sigma);
  double counter = 0.0;
  for (int j=0; j<p; j++) {
    bool fails = ((eigval(j) / eigval(p-1)) < tol);
    if (fails) {
      counter += 1.0;
      failures(j) = 1;
    }
  }
  if (counter > 0) {
    for (int j=0; j<p; j++) {
      if (failures(j) == 1) {
        eigval(j) = tol / (counter/p);
      }
    }
    ans = eigvec * arma::diagmat(eigval) * eigvec.t();
  } else {
    ans = sigma;
  }
  return ans;
}


//' Mahalanobis distance
//'
//' @description Compute the squared Mahalanobis distance (fast implementation).
//'
//' @references Mahalanobis, Prasanta Chandra. "On the Generalised Distance in
//'   Statistics." Proceedings of the National Institute of Sciences of India 2
//'   (1936): 49-55.
//'
//' @param x Numeric. A vector (of length \eqn{p}) or matrix (with \eqn{p}
//'   columns).
//' @param mu Numeric. A vector of length \eqn{p}.
//' @param sigma Numeric. A \eqn{p \times p} non-negative definite matrix.
//' @param ischol Logical. Set to \code{TRUE} if \code{sigma} is provided as a
//'   Cholesky decomposition.
//'
//' @return The squared Mahalanobis distance for all rows in \code{x} and the
//'   mean vector \code{mu} with respect to covariance matrix \code{sigma},
//'   defined as \eqn{(x - \mu)' \Sigma^{-1}(x - \mu)}.
//'
//' @author Emily Goren, \email{emily.goren@gmail.com}
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec mahalanobis(arma::mat x, arma::vec mu, arma::mat sigma, bool ischol = false) {
  // Check inputs.
  if (mu.n_elem != sigma.n_cols) {
    Rcpp::stop("The supplied mean vector and covariance matrix have incompatible dimensions.");
  }
  if (x.n_cols != sigma.n_cols)  {
    Rcpp::stop("The supplied data matrix and covariance matrix have incompatible dimensions.");
  }
  // Cholesky decomp of covariance matrix -- take lower triangle.
  int n = x.n_rows, p = x.n_cols;
  arma::mat A(p,p);
  if (!ischol) {
    arma::mat Atmp(p,p);
    bool success = arma::chol(Atmp, sigma);
    if (!success) {
      Atmp = arma::chol(fix_var(sigma));
    }
    A = arma::trimatl(Atmp.t());
  } else {
    A = arma::trimatl(sigma.t());
  }
  arma::vec D = A.diag();
  // Solve linear system.
  arma::vec ans(n);
  for (int i = 0; i < n; i++) {
    arma::vec tmp(p);
    for (int j = 0; j < p; j++) {
      double s = 0.0;
      for (int k = 0; k < j; k++) {
        s += tmp(k) * A(j, k);
      }
      tmp(j) = ( x(i, j) - mu(j) - s ) / D(j);
    }
    ans.at(i) = sum(square(tmp));
  }
  return ans;
}

//' Multivariate t distribution
//'
//' @description Compute the density of the multivariate t distribution (fast
//'   implementation).
//'
//' @param x Numeric. A vector (of length \eqn{p}) or matrix (with \eqn{p}
//'   columns).
//' @param mu Numeric. A vector of length \eqn{p} representing the mean.
//' @param sigma Numeric. A \eqn{p \times p} non-negative definite matrix (or
//'   its Cholesky decomposition).
//' @param nu Numeric. A positive scalar representing the degrees of freedom.
//' @param logans Logical. If \code{TRUE}, the log density is returned.
//' @param ischol Logical. Set to \code{TRUE} if \code{sigma} is provided as a
//'   Cholesky decomposition.
//'
//' @return The multivariate t density for all rows in \code{x} using degrees of
//'   freedom \code{nu}, mean vector \code{mu}, and covariance matrix
//'   \code{sigma}.
//'
//' @author Emily Goren, \email{emily.goren@gmail.com}
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec dMVT(arma::mat x, arma::vec mu, arma::mat sigma, double nu, bool logans = false, bool ischol = false) {
  // Check inputs.
  if (mu.n_elem != sigma.n_cols) {
    Rcpp::stop("The supplied mean vector and covariance matrix have incompatible dimensions.");
  }
  if (x.n_cols != sigma.n_cols)  {
    Rcpp::stop("The supplied data matrix and covariance matrix have incompatible dimensions.");
  }
  int p = x.n_cols, n = x.n_rows;
  arma::mat A(p,p);
  if (!ischol) {
    bool success = arma::chol(A, sigma);
    if (!success) {
      A = arma::chol(fix_var(sigma));
    }
  } else {
    A = sigma;
  }
  arma::vec ans(n);
  arma::vec maha = mahalanobis(x, mu, A, true);
  double logDet = sum(arma::log(A.diag()));
  if (nu <= 0.0) { // MVN
    ans = -0.5*maha - (p/2.0)*std::log(2.0*M_PI) + logDet;
  } else {
    double c = lgamma((nu + p) / 2.0) - lgamma(nu / 2.0) - (p / 2.0) * std::log(nu*M_PI) - logDet;
    for (int i=0; i<n; i++) {
      ans(i) = c - 0.5*(nu+p) * log1p(maha(i)/nu);
    }
  }
  if (!logans) ans = exp(ans);
  return ans;
}


//' Multivariate t distribution with missing data
//'
//' @description Compute the marginal density of the multivariate t distribution
//'   for the observed coordinates.
//'
//' @param x Numeric. A vector (of length \eqn{p}) or matrix (with \eqn{p}
//'   columns). Missing values for row \eqn{i} correspondind to unique pattern
//'   of missingness \eqn{m} are indicated by the \eqn{m^{th}} element of the
//'   selection index list \code{D}.
//' @param mu Numeric. A vector of length \eqn{p} representing the mean.
//' @param sigma Numeric. A \eqn{p \times p} non-negative definite matrix (or
//'   its Cholesky decomposition).
//' @param nu Numeric. A postitive scalar.
//' @param grp Numeric. A vector of length \eqn{n} with elements indicating the
//'   which row of \code{Ru} to use as the observed coordinates for each of
//'   the \eqn{n} observations.
//' @param Ru Binary matrix. Each row corresponds to a unique pattern of missingness,
//' where \eqn{1} indicates observed and \eqn{0} indicates missing coordinates.
//'
//' @return The multivariate t density for all rows in \code{x} using degrees of
//'   freedom \code{nu}, mean vector \code{mu}, covariance matrix \code{sigma},
//'   and observed coordinates specified by rows in \code{Ru}.
//'
//' @author Emily Goren, \email{emily.goren@gmail.com}
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec h(arma::mat x, arma::vec mu, arma::mat sigma, double nu, arma::vec grp, arma::umat Ru) {
  int n = x.n_rows, M = Ru.n_rows;
  arma::vec ans(n);
  for (int m=0; m<M; m++) {
    int g = m+1;
    // get mu, cholesky decom of sigma for this missingness pattern
    arma::uvec oidx = arma::find(Ru.row(m) == 1);
    int pg = oidx.size();
    arma::vec mug = mu.elem(oidx);
    arma::mat sigmag = sigma.submat(oidx, oidx);
    arma::mat Rg(pg,pg);
    bool success = arma::chol(Rg, sigmag);
    if (!success) {
      Rg = arma::chol(fix_var(sigmag));
    }
    // get obs for this missingness pattern
    arma::uvec gidx = arma::find(grp == g);
    arma::mat xg = x.submat(gidx, oidx);
    ans(gidx) = dMVT(xg, mug, Rg, nu, false, true);
  }
  return ans;
}


// E-step: update Z.
// [[Rcpp::export]]
arma::mat up_Z(arma::mat x, arma::mat mus, NumericVector sigmas, arma::vec nus, arma::vec pis, arma::vec grp, arma::umat Ru,
    arma::vec labeled_obs, arma::mat class_indicators) {
  int K = mus.n_rows, n = x.n_rows, p = x.n_cols;
  arma::cube Sigmas = to_array(sigmas, K, p);
  arma::mat ans(n, K);
  for (int k=0; k<K; k++) {
    ans.col(k) = pis(k) * h(x, mus.row(k).t(), Sigmas.slice(k), nus(k), grp, Ru);
  }
  for (int i=0; i<n; i++) {
    double rowsum = 0.0;
    for (int k=0; k<K; k++) {
      rowsum += ans(i,k);
    }
    // only update the Z's for unlabeled observations
    if (labeled_obs(i))
      ans.row(i) = class_indicators.row(i);
    else
      ans.row(i) = ans.row(i) / rowsum;
  }
  return ans;
}


// E-step: update W.
// [[Rcpp::export]]
arma::mat up_W(arma::mat x, arma::mat mus, NumericVector sigmas, arma::vec nus, arma::vec grp, arma::umat Ru) {
  int K = mus.n_rows, n = x.n_rows, p = x.n_cols, M = Ru.n_rows;
  arma::cube Sigmas = to_array(sigmas, K, p);
  arma::mat ans(n, K);
  for (int m=0; m<M; m++) {
    int g = m+1;
    // get obs for this missingness pattern
    arma::uvec oidx = arma::find(Ru.row(m) == 1);
    arma::uvec gidx = arma::find(grp == g);
    int pg = oidx.size();
    int ng = gidx.size();
    arma::mat xg = x.submat(gidx, oidx);
    for (int k=0; k<K; k++) {
      arma::vec muk = mus.row(k).t();
      arma::mat sigmak = Sigmas.slice(k);
      // get mu, cholesky decom of sigma for this missingness pattern
      arma::vec mukg = muk.elem(oidx);
      arma::mat sigmakg = sigmak.submat(oidx, oidx);
      arma::mat Rkg(pg,pg);
      bool success = arma::chol(Rkg, sigmakg);
      if (!success) {
        Rkg = arma::chol(fix_var(sigmakg));
      }
      arma::vec maha = mahalanobis(xg, mukg, Rkg, true);
      for (int i=0; i<ng; i++) {
        ans(gidx(i), k) = ( nus(k) + pg ) / ( nus(k) + maha(i) );
      }
    }
  }
  return ans;
}


// CM-step 1: update pis.
// [[Rcpp::export]]
arma::vec up_pi(arma::mat z) {
  int n = z.n_rows, K = z.n_cols;
  arma::vec ans(K); ans.zeros();
  for (int k=0; k<K; k++) {
    for (int i=0; i<n; i++) {
      ans(k) += z(i,k);
    }
    ans(k) = ans(k) / n;
  }
  return ans;
}


// CM-step 2: update mus.
// [[Rcpp::export]]
arma::mat up_mu(arma::mat x, arma::mat z, arma::mat w, arma::mat A) {
  int p = x.n_cols, n = x.n_rows, K = z.n_cols;
  arma::mat ans(K, p);
  for (int k=0; k<K; k++) {
    arma::mat L(p,p); L.zeros();
    arma::vec R(p); R.zeros();
    for (int i=0; i<n; i++) {
      arma::mat dA = arma::diagmat(A.row(i));
      L += z(i,k) * w(i,k) * dA;
      R += z(i,k) * w(i,k) * dA * x.row(i).t();
    }
    ans.row(k) = arma::solve(L, R).t();
  }
  return ans;
}


// CM-step 2: update Sigmas.
// [[Rcpp::export]]
arma::cube up_Sigma(arma::mat x, arma::mat z, arma::mat w, arma::mat mus, arma::mat A, String constr) {
  int p = x.n_cols, n = x.n_rows, K = z.n_cols;
  arma::cube Sigmas(p,p,K);
  //TODO: add parameter for maximum iterations for iterative solutions
  // same with tolerance
  int iter_max = 100;
  double tol = 0.01;
  // constraints with closed form solutions
  if (constr == "VVV" || constr == "EEE") {
    for (int k=0; k<K; k++) {
      arma::mat L(p,p); L.zeros();
      arma::mat R(p,p); R.zeros();
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
      Sigmas.slice(k) = L / R;
    }
    if (constr == "EEE") {
      arma::vec pis = up_pi(z);
      arma::mat S(p,p); S.zeros();
      for (int k=0; k<K; k++) {
        S += pis(k) * Sigmas.slice(k);
      }
      for (int k=0; k<K; k++) {
        Sigmas.slice(k) = S;
      }
    }
  }
  else if (constr == "VII") {
    for (int k = 0; k < K; k++) {
      arma::mat L(p,p); L.zeros();
      arma::mat R(p,p); R.zeros();
      arma::vec rowsums(p);
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
      rowsums = arma::sum(R, 1);
      R.zeros();
      R.diag() = rowsums;
      arma::mat Iden(p, p, arma::fill::eye);
      double zeta = arma::trace(L);
      Sigmas.slice(k) = zeta * Iden * R.i();
    }
  }
  else if (constr == "EII") {
    arma::mat L(p,p); L.zeros();
    arma::mat R(p,p); R.zeros();
    arma::vec rowsums(p);
    for (int k = 0; k < K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
    }
    rowsums = arma::sum(R, 1);
    R.zeros();
    R.diag() = rowsums;
    arma::mat Iden(p, p, arma::fill::eye);
    double zeta = arma::trace(L);
    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = zeta * Iden * R.i();
    }
  }
  else if (constr == "EEI") {
    arma::mat L(p,p); L.zeros();
    arma::mat R(p,p); R.zeros();
    for (int k = 0; k < K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
    }

    arma::mat Lambda = diagmat(L);
    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = Lambda / R;
    }
  }
  else if (constr == "VVI") {
    for (int k = 0; k < K; k++) {
      arma::mat L(p,p); L.zeros();
      arma::mat R(p,p); R.zeros();
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
      // R = diagmat(R);
      arma::mat Lambda = diagmat(L);
      Sigmas.slice(k) = Lambda / R;
      // Sigmas.slice(k) = Lambda * R.i();
    }
  }
  else if (constr == "EVI") {
    double zeta = 0.0;
    double detval = 0.0;
    arma::mat R(p,p); R.zeros();
    for (int k = 0; k < K; k++) {
      arma::mat L(p, p); L.zeros();
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
      arma::mat Lambda = arma::diagmat(L);
      detval = pow(arma::det(Lambda), 1.0 / p);
      zeta += detval;
      Lambda = Lambda / detval;

      Sigmas.slice(k) = Lambda;
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = Sigmas.slice(k) * zeta / R;
    }
  }
  else if (constr == "EVV") {
    double zeta = 0.0;
    double detval = 0.0;
    arma::mat R(p,p); R.zeros();
    for (int k = 0; k < K; k++) {
      arma::mat L(p, p); L.zeros();
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
      detval = pow(arma::det(L), 1.0 / p);
      zeta += detval;
      L = L / detval;

      Sigmas.slice(k) = L;
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = zeta * Sigmas.slice(k) / R;
    }
  }
  // constraints with iterative solutions
  else if (constr == "VEE") {
    // first, we initialize the matrix C
    // using the EEE solution
    arma::vec zeta(K); zeta.zeros();
    arma::mat C(p, p); C.zeros();
    arma::cube L(p,p, K); L.zeros();
    arma::cube R(p,p, K); R.zeros();
    arma::vec rowsum(p);
    for (int k=0; k<K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L.slice(k) += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R.slice(k) += z(i,k) * A.row(i).t() * A.row(i);
      }
      rowsum = sum(R.slice(k), 1);
      R.slice(k).zeros();
      R.slice(k).diag() = rowsum;
    }
    /* arma::mat Rsum; */
    /* Rsum = arma::sum(R, 2); */
    /* C = arma::sum(L, 2) * Rsum.i(); */
    C.eye();
    // start iterative procedure
    arma::mat W(p, p); W.zeros(); // accumulator for C
    arma::mat C_old(p, p);
    arma::vec zeta_old(K);
    double c_diff = 0;
    double c_denom = 0;
    double zeta_diff = 0;
    double zeta_denom = 0;
    for (int iter = 0; iter < iter_max; iter++) {
      // first, find optimal zeta's with respect to C
      W.zeros();
      C_old = C;
      zeta_old = zeta;
      for (int k = 0; k < K; k++) {
        zeta(k) = arma::trace(L.slice(k) * C.i());
        W += L.slice(k) / zeta(k) * R.slice(k);
      }
      C = W / pow(arma::det(W), 1.0 / p);

      // check termination condition
      if (iter > 1) {
        c_diff = arma::accu(arma::square(C_old - C));
        zeta_diff = arma::accu(arma::square(zeta_old - zeta));

        c_denom = arma::accu(arma::square(C_old));
        zeta_denom = arma::accu(arma::square(zeta_old));

        if (abs(zeta_diff / zeta_denom - 1) < tol && abs(c_diff / c_denom - 1) < tol) {
          break;
        }
      }
    }
    for (int k = 0; k < K; k++) {
      // Sigmas.slice(k) = zeta(k) * C * R.slice(k).i();
      Sigmas.slice(k) = arma::symmatu(zeta(k) * C * R.slice(k).i());
    }
  }
  else if (constr == "VEV") {
    // initialize solution using EEE
    Sigmas = up_Sigma(x, z, w, mus, A, constr = "EEE");

    arma::cube L(p,p, K); L.zeros();
    arma::cube R(p,p, K); R.zeros();

    arma::mat eigvec(p, p);
    arma::vec eigval(p);
    arma::vec rowsum(p);

    arma::cube Gamma(p, p, K);
    arma::cube Omega(p, p, K);

    for (int k = 0; k < K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L.slice(k) += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R.slice(k) += z(i,k) * A.row(i).t() * A.row(i);
      }
      rowsum = arma::sum(R.slice(k), 1);
      R.slice(k).zeros();
      R.slice(k).diag() = rowsum;
      arma::eig_sym(eigval, eigvec, L.slice(k)); // to ensure computational symchange2 here
      Gamma.slice(k) = eigvec;
      Omega.slice(k) = arma::diagmat(eigval);
    }

    arma::vec zeta(K); zeta.zeros();
    arma::mat Lambda(p, p); Lambda.eye();

    arma::mat Lambda_old(p, p);
    arma::vec zeta_old(K);

    double lambda_diff, lambda_denom;
    double zeta_diff, zeta_denom;

    for (int iter = 0; iter < iter_max; iter++) {
      Lambda_old = Lambda;
      Lambda.zeros();
      for (int k = 0; k < K; k++) {
        zeta(k) = arma::trace(L.slice(k) * Gamma.slice(k) * Lambda_old.i() * Gamma.slice(k).t());
        Lambda += Omega.slice(k) * R.slice(k) / zeta(k);
      }

      Lambda = Lambda / pow(arma::det(Lambda), 1.0 / p);

      if (iter > 1) {
        zeta_diff = arma::accu(arma::square(zeta_old - zeta));
        lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));

        lambda_denom = arma::accu(arma::square(Lambda_old));
        zeta_denom = arma::accu(arma::square(zeta_old));

        if (abs(zeta_diff / zeta_denom - 1) < tol && abs(lambda_diff / lambda_denom - 1) < tol) {
          break;
        }
      }
    }

    for (int k = 0; k < K; k++) {
      // rowsum = R.slice(k).diag();
      // R.slice(k).zeros();
      // R.slice(k).diag() = rowsum;
      // force symmetry to ensure computational errors don't happen
      Sigmas.slice(k) = arma::symmatu(zeta(k) * Gamma.slice(k) * Lambda * Gamma.slice(k).t() * R.slice(k).i());
    }
  }
  else if (constr == "VEI") {
    // values to be updated
    // start with identity matrix for Lambda
    arma::mat Lambda = arma::eye(p, p);
    arma::vec zeta(K);

    // values to perform calculations with
    arma::mat W_acc(p, p);

    // values to check termination conditions
    double zeta_diff, zeta_denom;
    double lambda_diff, lambda_denom;
    arma::mat Lambda_old(p, p);
    arma::vec zeta_old(K);

    // get values for L and R
    arma::cube L(p, p, K); L.zeros();
    arma::cube R(p, p, K); R.zeros();
    arma::vec rowsum(p);
    for (int k = 0; k < K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L.slice(k) += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R.slice(k) += z(i,k) * A.row(i).t() * A.row(i);
      }
      rowsum = arma::sum(R.slice(k), 1);
      R.slice(k).zeros();
      R.slice(k).diag() = rowsum;
    }

    for (int iter = 0; iter < iter_max; iter++) {
      W_acc.zeros();
      Lambda_old = Lambda;
      zeta_old = zeta;
      // update values, starting with zeta vector

      // important to keep track of the idea that in non-missing cases, zeta
      // should be divided by (pn_k), but we are using the matrix R, so it
      // any time a matrix is multiplied by zeta, it should be divided
      // (element-wise) by R
      for (int k = 0; k < K; k++) {
        zeta(k) = arma::trace(L.slice(k) * Lambda.i());
        W_acc += L.slice(k) * R.slice(k) / zeta(k);
      }
      Lambda = arma::diagmat(W_acc) / pow(arma::det(arma::diagmat(W_acc)), 1.0 / p);

      // check termination condition
      if (iter > 1) {
        zeta_diff = arma::accu(arma::square(zeta_old - zeta));
        lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));

        lambda_denom = arma::accu(arma::square(Lambda_old));
        zeta_denom = arma::accu(arma::square(zeta_old));

        if (abs(zeta_diff / zeta_denom - 1) < tol && abs(lambda_diff / lambda_denom - 1) < tol) {
          break;
        }
      }
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = Lambda * zeta(k) * R.slice(k).i();
      // Sigmas.slice(k) = Lambda * zeta(k) / R.slice(k);
    }
  }
  else if (constr == "EVE") {
    // get values for L and R
    arma::cube L(p, p, K); L.zeros();
    arma::mat R(p, p); R.zeros();

    // values to store eigenvalues to be worked with later
    arma::vec eigval(p);
    arma::mat eigvec(p, p);

    // vector to store largest eigenvalue for each K
    arma::vec L_eig(K);

    for (int k = 0; k < K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L.slice(k) += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
      arma::eig_sym(eigval, eigvec, arma::symmatu(L.slice(k))); // values are in ascending order
      L_eig(k) = eigval(p - 1);
    }

    // initialize lamdba matrices as identity
    arma::cube Lambda(p, p, K);
    for (int k = 0; k < K; k++) {
      Lambda.slice(k).eye();
    }

    arma::mat Gamma(p, p); Gamma.zeros();

    // matrix for MM iterations (and its svd)
    arma::mat F(p, p);
    arma::mat P;
    arma::mat C;
    arma::vec b;

    // vector to store largest eigenvalues of lambda matrices
    arma::vec Lambda_eig(K);

    // working vectors for termination conditions, etc
    double detval = 0.0;
    arma::mat Gamma_old(p, p);
    arma::mat Gamma_old2(p, p);
    arma::cube Lambda_old(p, p, K);
    double Gamma_diff, Gamma_denom;
    double Lambda_diff, Lambda_denom;

    // main iteration loop
    for (int iter = 0; iter < iter_max; iter++) {
      F.zeros();
      Gamma_old = Gamma;
      Lambda_old = Lambda;
      // there needs to be an inner loop here to solve for Gamma
      for (int in = 0; in < iter_max; in++) {
        Gamma_old2 = Gamma;
        if (iter % 2 == 0) {
          for (int k = 0; k < K; k++) {
            F += Lambda.slice(k).i() * Gamma.t() * L.slice(k) - L_eig(k) * Lambda.slice(k).i() * Gamma.t();
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }
        else {
          // we need to find the largest eigenvalue of each Lambda(k);
          for (int k = 0; k < K; k++) {
            arma::eig_sym(eigval, eigvec, arma::symmatu(Lambda.slice(k).i()));
            Lambda_eig(k) = eigval(p - 1);
            F += L.slice(k) * Gamma * Lambda.slice(k).i() - Lambda_eig(k) * L.slice(k) * Gamma;
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }

        // check termination
        if (in > 1) {
          Gamma_diff = arma::accu(arma::square(Gamma_old2 - Gamma));
          Gamma_denom = arma::accu(arma::square(Gamma_old2));

          if (abs(Gamma_diff / Gamma_denom - 1) < tol) {
            break;
          }
        }
      }
      // end of inner loop

      for (int k = 0; k < K; k++) {
        Lambda.slice(k) = arma::diagmat(Gamma.t() * L.slice(k) * Gamma);
        detval = arma::det(Lambda.slice(k));
        Lambda.slice(k) = Lambda.slice(k) / pow(detval, 1.0 / p);
      }

      // check termination
      if (iter > 1) {
        Gamma_diff = arma::accu(arma::square(Gamma_old - Gamma));
        Gamma_denom = arma::accu(arma::square(Gamma_old));

        Lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));
        Lambda_denom = arma::accu(arma::square(Lambda_old));

        if (abs(Gamma_diff / Gamma_denom - 1) < tol && abs(Lambda_diff / Lambda_denom - 1) < tol) {
          break;
        }
      }
    }

    // compute optimal zeta and return resulting matrices
    double zeta = 0.0;
    for (int k = 0; k < K; k++) {
      zeta += arma::trace(Gamma * Lambda.slice(k) * Gamma.t() * L.slice(k));
    }

    arma::vec rowsum(p);
    for (int k = 0; k < K; k++) {
      rowsum = arma::sum(R, 1);
      R.zeros();
      R.diag() = rowsum;
      // Sigmas.slice(k) = zeta * Gamma * Lambda.slice(k) * Gamma.t() * R.i();
      Sigmas.slice(k) = arma::symmatu(zeta * Gamma * Lambda.slice(k) * Gamma.t() * R.i());
    }
  }
  else if (constr == "VVE") {
    // get values for L and R
    arma::cube L(p, p, K); L.zeros();
    arma::cube R(p, p, K); R.zeros();

    // values to store eigenvalues to be worked with later
    arma::vec eigval(p);
    arma::mat eigvec(p, p);

    // vector to store largest eigenvalue for each K
    arma::vec L_eig(K);

    for (int k = 0; k < K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L.slice(k) += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R.slice(k) += z(i,k) * A.row(i).t() * A.row(i);
      }
      arma::eig_sym(eigval, eigvec, arma::symmatu(L.slice(k))); // values are in ascending order
      L_eig(k) = eigval(p - 1);
    }

    // initialize lamdba matrices as identity
    arma::cube Lambda(p, p, K);
    for (int k = 0; k < K; k++) {
      Lambda.slice(k).eye();
    }

    arma::mat Gamma(p, p); Gamma.zeros();

    // matrix for MM iterations (and its svd)
    arma::mat F(p, p);
    arma::mat P;
    arma::mat C;
    arma::vec b;

    // vector to store largest eigenvalues of lambda matrices
    arma::vec Lambda_eig(K);

    // working vectors for termination conditions, etc
    double detval = 0.0;
    arma::mat Gamma_old(p, p);
    arma::mat Gamma_old2(p, p);
    arma::cube Lambda_old(p, p, K);
    double Gamma_diff, Gamma_denom;
    double Lambda_diff, Lambda_denom;

    // main iteration loop
    for (int iter = 0; iter < iter_max; iter++) {
      F.zeros();
      Gamma_old = Gamma;
      Lambda_old = Lambda;
      // there needs to be an inner loop here to solve for Gamma
      for (int in = 0; in < iter_max; in++) {
        Gamma_old2 = Gamma;
        if (iter % 2 == 0) {
          for (int k = 0; k < K; k++) {
            F += Lambda.slice(k).i() * Gamma.t() * L.slice(k) - L_eig(k) * Lambda.slice(k).i() * Gamma.t();
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }
        else {
          // we need to find the largest eigenvalue of each Lambda(k);
          for (int k = 0; k < K; k++) {
            arma::eig_sym(eigval, eigvec, arma::symmatu(Lambda.slice(k).i()));
            Lambda_eig(k) = eigval(p - 1);
            F += L.slice(k) * Gamma * Lambda.slice(k).i() - Lambda_eig(k) * L.slice(k) * Gamma;
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }

        // check termination
        if (in > 1) {
          Gamma_diff = arma::accu(arma::square(Gamma_old2 - Gamma));
          Gamma_denom = arma::accu(arma::square(Gamma_old2));

          if (abs(Gamma_diff / Gamma_denom - 1) < tol) {
            break;
          }
        }
      }
      // end of inner loop

      for (int k = 0; k < K; k++) {
        Lambda.slice(k) = arma::diagmat(Gamma.t() * L.slice(k) * Gamma);
        detval = arma::det(Lambda.slice(k));
        Lambda.slice(k) = Lambda.slice(k) / pow(detval, 1.0 / p);
      }

      // check termination
      if (iter > 1) {
        Gamma_diff = arma::accu(arma::square(Gamma_old - Gamma));
        Gamma_denom = arma::accu(arma::square(Gamma_old));

        Lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));
        Lambda_denom = arma::accu(arma::square(Lambda_old));

        if (abs(Gamma_diff / Gamma_denom - 1) < tol && abs(Lambda_diff / Lambda_denom - 1) < tol) {
          break;
        }
      }
    }

    // compute optimal zeta and return resulting matrices
    arma::vec zeta(K);
    for (int k = 0; k < K; k++) {
      zeta(k) = arma::trace(Gamma * Lambda.slice(k) * Gamma.t() * L.slice(k));
    }

    arma::vec rowsum(p);
    for (int k = 0; k < K; k++) {
      rowsum = arma::sum(R.slice(k), 1);
      R.slice(k).zeros();
      R.slice(k).diag() = rowsum;
      // Sigmas.slice(k) = zeta(k) * Gamma * Lambda.slice(k) * Gamma.t() * R.slice(k).i();
      Sigmas.slice(k) = arma::symmatu(zeta(k) * Gamma * Lambda.slice(k) * Gamma.t() * R.slice(k).i());
    }
  }
  else if (constr == "EEV") {
    // based on the notes, there is no need to calculate the zetas
    arma::cube L(p,p, K); L.zeros();
    arma::mat R(p,p); R.zeros();

    arma::cube Gamma(p, p, K);
    arma::mat Lambda(p, p); Lambda.zeros();

    // work arrays
    arma::mat eigvec(p, p);
    arma::vec eigval(p);
    arma::vec diagval(p);

    for (int k = 0; k < K; k++) {
      for (int i=0; i<n; i++) {
        arma::mat Ai = arma::diagmat(A.row(i));
        arma::vec u = x.row(i).t() - mus.row(k).t();
        L.slice(k) += z(i,k) * w(i,k) * Ai * u * u.t() * Ai;
        R += z(i,k) * A.row(i).t() * A.row(i);
      }
      arma::eig_sym(eigval, eigvec, L.slice(k));
      Gamma.slice(k) = eigvec;
      Lambda += arma::diagmat(eigval);
    }

    for (int k = 0; k < K; k++) {
      // diagval = R.diag();
      // R.zeros();
      // R.diag() = diagval;
      // Sigmas.slice(k) = Gamma.slice(k) * Lambda * Gamma.slice(k).t() * R.i();
      // avoid accumulated computational error
      Sigmas.slice(k) = arma::symmatu(Gamma.slice(k) * Lambda * Gamma.slice(k).t() / R);
    }

  }
  // check if covariance is close to singular, if it is add a bit to the diagonal
  // the only effect this should have is for cases when there is a degenerate solution,
  // which (i hope) is thrown out anyway
  // TODO: maybe delet
  double r;
  for (int k = 0; k < K; k++) {
    r = arma::rcond(Sigmas.slice(k));
    if (r < 1e-16)
      Sigmas.slice(k) += 1e-3 * arma::eye(p, p);
  }
  return Sigmas;
}


// CM-step 3: update nus -- helper functions for updating cluster-specifc
// degrees of freedom.
double approx_nu(NumericVector z, NumericVector w, double nu, NumericVector ps, int n) {
  double out, tmp = 0.0, nk = 0.0;
  for (int i=0; i<n; i++) {
    tmp += z[i] * (log(w[i]) - w[i] + R::digamma((nu + ps[i])/2.0) - log((nu + ps[i])/2.0));
    nk += z[i];
    }
  double cc = - 1.0 - tmp/nk;
  double num = - exp(cc) + 2.0*exp(cc) * (exp(R::digamma(nu/2.0)) - ((nu/2.0) - 0.5));
  double den = 1.0 - exp(cc);
  out = num / den;
  return out;
}
struct nu_pars {
  NumericVector z_k;
  NumericVector w_k;
  double nu_k;
  NumericVector P;
  int n;
};
double objective_nu(double df, void *params) {
  struct nu_pars *pars;
  pars = (struct nu_pars *)params;
  int n = pars->n;
  NumericVector ps(n);
  ps = pars->P;
  double nu_k = pars->nu_k;
  NumericVector z_k(n);
  z_k = pars->z_k;
  NumericVector w_k(n);
  w_k = pars->w_k;
  double zsum = 0, asum = 0;
  for (int i=0; i<n; i++) {
    asum += z_k[i] * (log(w_k[i])-w_k[i]+R::digamma((nu_k+ps[i])/2.0)-log((nu_k+ps[i])/2.0));
    zsum += z_k[i];
  }
  double ans = 1.0 - R::digamma(df/2.0) + log(df/2.0) + (asum/zsum);
  return ans;
}
double rootsolver_nu(NumericVector z_k, NumericVector w_k, double nu_k, NumericVector ps, int n, int iter_max = 1e6, double tol = 1e-3, double min_df = 1e-3, double max_df = 1e4) {
  int status, iter = 0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *solver;
  double ans;
  gsl_function F;
  struct nu_pars params = {z_k, w_k, nu_k, ps, n};
  F.function = &objective_nu;
  F.params = &params;
  T = gsl_root_fsolver_brent;
  solver = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set (solver, &F, min_df, max_df);
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (solver);
    ans = gsl_root_fsolver_root (solver);
    min_df = gsl_root_fsolver_x_lower (solver);
    max_df = gsl_root_fsolver_x_upper (solver);
    status = gsl_root_test_interval (min_df, max_df, 0, tol);
  }
  while (status == GSL_CONTINUE && iter < iter_max);
  gsl_root_fsolver_free (solver);
  return ans;
}


// CM-step 3: update nus -- helper functions for updating constrained (same for
// all clusters) degrees of freedom.
double approx_nu_constr(NumericMatrix z, NumericMatrix w, double nu, NumericVector ps, int n, int K) {
  double out;
  double tmp = 0.0;
  for (int k=0; k<K; k++) {
    double tmpk = 0.0;
    for (int i=0; i<n; i++) {
      tmpk += z(i,k) * (log(w(i,k)) - w(i,k) + R::digamma((nu + ps[i])/2.0) - log((nu + ps[i])/2.0));
    }
    tmp += tmpk;
  }
  double cc = -1.0 - tmp/n;
  double num = -exp(cc) + 2.0*exp(cc) * (exp(R::digamma(nu/2.0)) - ((nu/2.0) - 0.5));
  double den = 1.0 - exp(cc);
  out = num / den;
  return out;
}

struct nu_constr_pars {
  NumericMatrix z;
  NumericMatrix w;
  double nu;
  NumericVector P;
  int K;
  int n;
};
double objective_nu_constr(double df, void *params) {
  struct nu_constr_pars *pars;
  pars = (struct nu_constr_pars *)params;
  int n = pars->n;
  int K = pars->K;
  NumericVector ps(n);
  ps = pars->P;
  double nu = pars->nu;
  NumericMatrix z(n,K);
  z = pars->z;
  NumericMatrix w(n,K);
  w = pars->w;
  double asum = 0;
  for (int k=0; k<K; k++) {
    for (int i=0; i<n; i++) {
      asum += z(i,k) * (log(w(i,k))-w(i,k)+R::digamma((nu+ps[i])/2.0)-log((nu+ps[i])/2.0));
    }
  }
  double ans = 1.0 - R::digamma(df/2.0) + log(df/2.0) + (asum/n);
  return ans;
}
double rootsolver_nu_constr(NumericMatrix z, NumericMatrix w, double nu, NumericVector ps, int n, int K, int iter_max = 1e6, double tol = 1e-3, double min_df = 1e-3, double max_df = 1e4) {
  int status, iter = 0;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *solver;
  double ans;
  gsl_function F;
  struct nu_constr_pars params = {z, w, nu, ps, K, n};
  F.function = &objective_nu_constr;
  F.params = &params;
  T = gsl_root_fsolver_brent;
  solver = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (solver, &F, min_df, max_df);
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (solver);
    ans = gsl_root_fsolver_root (solver);
    min_df = gsl_root_fsolver_x_lower (solver);
    max_df = gsl_root_fsolver_x_upper (solver);
    status = gsl_root_test_interval (min_df, max_df, 0, tol);
  }
  while (status == GSL_CONTINUE && iter < iter_max);
  gsl_root_fsolver_free (solver);
  return ans;
}


// CM-step 3: update nus.
// [[Rcpp::export]]
NumericVector up_nu(NumericMatrix z, NumericMatrix w, NumericVector nus, NumericVector ps, bool constr = false, bool approx = false) {
  int n = z.nrow(), K = z.ncol();
  NumericVector ans(K);
  if (constr) {
    double tmp;
    if (approx) {
      tmp = approx_nu_constr(z, w, nus(0), ps, n, K);
    } else {
      tmp = rootsolver_nu_constr(z, w, nus(0), ps, n, K);
    }
    if (tmp < 3.0) {
      tmp = 3.0;
    }
    if (tmp > 200) {
      tmp = 200;
    }
    for (int k=0; k<K; k++) {
      ans(k) = tmp;
    }
  } else if (!constr) {
    for (int k=0; k<K; k++) {
      if (approx) {
        ans(k) = approx_nu(z(_,k), w(_,k), nus(k), ps, n);
      } else {
        ans(k) = rootsolver_nu(z(_,k), w(_,k), nus(k), ps, n);
      }
      if (ans(k) < 3.0) {
        ans(k) = 3.0;
      }
      if (ans(k) > 200) {
        ans(k) = 200;
      }
    }
  } else {
    Rcpp::stop("Degree of freedom constraint option must be boolean.");
  }
  return ans;
}


// Code to implement full EM (also including y.NA) per Lin's method
// For a given k, compute ES'inv(SES')S,
// where S is the selection matrix based on the observed indices in D
// and E is the dispersion of cluster k.
// [[Rcpp::export]]
arma::cube SOiOEOOk(arma::mat sigma, arma::umat Ru) {
  int p = sigma.n_cols, M = Ru.n_rows;
  arma::cube ans(p,p,M);
  for (int m=0; m<M; m++) {
    arma::uvec oidx = arma::find(Ru.row(m) == 1);
    int pg = oidx.size();
    arma::mat Rg(pg,pg);
    arma::mat sigmag = sigma.submat(oidx, oidx);
    bool success = arma::inv(Rg, sigmag);
    if (!success) {
      Rg = arma::inv(fix_var(sigmag));
    }
      arma::mat R(p,p); R.zeros();
    R(oidx, oidx) = Rg;
    ans.slice(m) = sigma * R;
  }
  return ans;
}
// For a given k, compute xhat.
// [[Rcpp::export]]
arma::mat xhatk(arma::mat x, arma::vec mu, arma::vec grp, int M, NumericVector SOiOEOOk) {
  int p = x.n_cols, n = x.n_rows;
  arma::cube R = to_array(SOiOEOOk, M, p);
  arma::mat xhat(n,p);
  for (int i=0; i<n; i++) {
    int g = grp(i) - 1;
    arma::vec u = x.row(i).t() - mu;
    arma::vec xhati = mu + R.slice(g) * u;
    xhat.row(i) = xhati.t();
  }
  return xhat;
}
// Update mus per Lins method.
// [[Rcpp::export]]
arma::mat up_mu_Lin(int p, arma::mat z, arma::mat w, ListOf<NumericMatrix> xhat) {
  int n = z.n_rows, K = z.n_cols;
  arma::mat ans(K, p);
  for (int k=0; k<K; k++) {
    arma::vec num(p); num.zeros();
    double den = 0, wt = 0;
    NumericMatrix tmp = xhat[k];
    arma::mat xhatk(tmp.begin(), n, p, false);
    for (int i=0; i<n; i++) {
      wt = z(i,k) * w(i,k);
      num += wt*xhatk.row(i).t();
      den += wt;
    }
    ans.row(k) = num.t() / den;
  }
  return ans;
}

// Update Sigmas per Lins method (one k at a time)
// [[Rcpp::export]]
arma::mat up_Sigmak_Lin(int M, arma::vec zk, arma::vec wk, arma::vec mu, arma::mat sigma,
			arma::mat xhatk, arma::vec grp, NumericVector SOiOEOOk) {
  int n = xhatk.n_rows, p = xhatk.n_cols;
  arma::mat I = arma::eye(p,p);
  arma::mat num(p,p); num.zeros();
  double den = 0;
  arma::cube R = to_array(SOiOEOOk, M, p);
  for (int i=0; i<n; i++) {
    int g = grp(i) - 1;
    arma::vec u = xhatk.row(i).t() - mu;
    num += zk(i) * (wk(i) * u * u.t() + (I - R.slice(g)) * sigma);
    den += zk(i);
  }
  arma::mat ans = num / den;
  return ans;
}

// Update Sigmas per Lins method with constraints
// [[Rcpp::export]]
arma::cube up_Sigma_Lin(arma::mat z, arma::mat w, arma::mat mu, arma::cube sigma,
			arma::cube xhat, arma::vec grp, String constr, arma::umat Ru) {
  // TODO: add parameter for maximum iterations for iterative solutions
  // same with tolerance
  int iter_max = 100;
  double tol = 0.01;

  int n = xhat.slice(1).n_rows, p = xhat.slice(1).n_cols, K = z.n_cols;
  arma::cube Sigmas(p, p, K); Sigmas.zeros();
  arma::mat I = arma::eye(p,p);
  // L here is the Omega matrix in Lin's paper
  arma::cube L(p, p, K); L.zeros();
  arma::mat L_tot(p, p); L.zeros();
  arma::vec n_k(K); n_k.zeros();
  for (int k = 0; k < K; k++) {
    arma::cube R = SOiOEOOk(sigma.slice(k), Ru);
    for (int i = 0; i < n; i++) {
      int g = grp(i) - 1;
      arma::vec u = xhat.slice(k).row(i).t() - mu.row(k).t(); // this probably needs to be changed
      L.slice(k) += z(i, k) * (w(i, k) * u * u.t() + (I - R.slice(g)) * sigma.slice(k));
      n_k(k) += z(i, k);
    }
    L_tot += L.slice(k);
    if (constr == "VVV") {
      Sigmas.slice(k) = L.slice(k) / n_k(k);
    }
  }

  if (constr == "EEE") {
    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = L_tot / n;
    }
  }
  else if (constr == "VII") {
    double zeta;
    for (int k = 0; k < K; k++) {
      zeta = arma::trace(L.slice(k));

      Sigmas.slice(k) = zeta * I / (p * n_k(k));
    }
  }
  else if (constr == "EII") {
    double zeta = arma::trace(L_tot);
    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = zeta * I / (p * n);
    }
  }
  else if (constr == "EEI") {
    arma::mat Lambda = arma::diagmat(L_tot);
    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = Lambda / n;
    }
  }
  else if (constr == "VVI") {
    arma::mat Lambda;
    for (int k = 0; k < K; k++) {
      Lambda = arma::diagmat(L.slice(k));
      Sigmas.slice(k) = Lambda / n_k(k);
    }
  }
  else if (constr == "EVI") {
    // the method here differs from what Lin has in his paper.
    // I believe this is the correct solution (per Maitra's notes), see function up_Sigma
    // (also the method in Lin's paper produces an error)
    arma::mat Lambda;
    double zeta = 0.0;
    double detval;
    for (int k = 0; k < K; k++) {
      Lambda = arma::diagmat(L.slice(k));
      detval = pow(arma::det(Lambda), 1.0 / p);
      Sigmas.slice(k) = Lambda / detval;
      zeta += detval;
    }
    zeta /= n;
    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) *= zeta;
    }
  }
  else if (constr == "EVV") {
    double zeta = 0.0;
    double detval;
    for (int k = 0; k < K; k++) {
      detval = pow(arma::det(L.slice(k)), 1.0 / p);
      Sigmas.slice(k) = L.slice(k) / detval;
      zeta += detval;
    }
    zeta /= n;
    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) *= zeta;
    }
  }
  else if (constr == "VEE") {
    arma::mat C = arma::eye(p, p);
    arma::mat C_old(p, p);

    arma::vec zeta(K);
    arma::vec zeta_old(K);

    double c_diff, c_denom, zeta_diff, zeta_denom;

    for (int iter = 0; iter < iter_max; iter++) {
      C_old = C;
      zeta_old = zeta;
      C.zeros();
      for (int k = 0; k < K; k++) {
        zeta(k) = arma::trace(L.slice(k) * C_old.i()) / (p * n_k(k));
        C += L.slice(k) / zeta(k);
      }
      C = C / pow(arma::det(C), 1.0 / p);

      if (iter > 1) {
        c_diff = arma::accu(arma::square(C_old - C));
        zeta_diff = arma::accu(arma::square(zeta_old - zeta));

        c_denom = arma::accu(arma::square(C_old));
        zeta_denom = arma::accu(arma::square(zeta_old));

        if (abs(zeta_diff / zeta_denom - 1) < tol && abs(c_diff / c_denom - 1) < tol) {
          break;
        }
      }
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = zeta(k) * C;
    }
  }
  else if (constr == "VEV") {
    // initialize solution using EEE
    Sigmas = up_Sigma_Lin(z, w, mu, sigma, xhat, grp, "EEE", Ru);

    arma::mat eigvec(p, p);
    arma::vec eigval(p);

    arma::cube Gamma(p, p, K);
    // Omega's contain the eigenvalues of the L matrix
    arma::cube Omega(p, p, K);

    for (int k = 0; k < K; k++) {
      arma::eig_sym(eigval, eigvec, L.slice(k));
      Gamma.slice(k) = eigvec;
      Omega.slice(k) = arma::diagmat(eigval);
    }

    arma::vec zeta(K); zeta.zeros();
    arma::mat Lambda(p, p); Lambda.eye();

    arma::mat Lambda_old(p, p);
    arma::vec zeta_old(K);

    double lambda_diff, lambda_denom;
    double zeta_diff, zeta_denom;

    for (int iter = 0; iter < iter_max; iter++) {
      Lambda_old = Lambda;
      Lambda.zeros();
      for (int k = 0; k < K; k++) {
        zeta(k) = arma::trace(L.slice(k) * Gamma.slice(k) * Lambda_old.i() * Gamma.slice(k).t()) / (p * n_k(k));
        Lambda += Omega.slice(k) / zeta(k);
      }

      Lambda = Lambda / pow(arma::det(Lambda), 1.0 / p);

      if (iter > 1) {
        zeta_diff = arma::accu(arma::square(zeta_old - zeta));
        lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));

        lambda_denom = arma::accu(arma::square(Lambda_old));
        zeta_denom = arma::accu(arma::square(zeta_old));

        if (abs(zeta_diff / zeta_denom - 1) < tol && abs(lambda_diff / lambda_denom - 1) < tol) {
          break;
        }
      }
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = zeta(k) * Gamma.slice(k) * Lambda * Gamma.slice(k).t();
    }

  }
  else if (constr == "VEI") {
    arma::mat Lambda = arma::eye(p, p);
    arma::vec zeta(K);

    arma::mat W_acc(p, p);

    double zeta_diff, zeta_denom;
    double lambda_diff, lambda_denom;
    arma::mat Lambda_old(p, p);
    arma::vec zeta_old(K);

    for (int iter = 0; iter < iter_max; iter++) {
      W_acc.zeros();
      Lambda_old = Lambda;
      zeta_old = zeta;

      for (int k = 0; k < K; k++) {
        zeta(k) = arma::trace(L.slice(k) * Lambda.i()) / (p * n_k(k));
        W_acc += L.slice(k) / zeta(k);
      }
      Lambda = arma::diagmat(W_acc) / pow(arma::det(arma::diagmat(W_acc)), 1.0 / p);

      if (iter > 1) {
        zeta_diff = arma::accu(arma::square(zeta_old - zeta));
        lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));

        lambda_denom = arma::accu(arma::square(Lambda_old));
        zeta_denom = arma::accu(arma::square(zeta_old));

        if (abs(zeta_diff / zeta_denom - 1) < tol && abs(lambda_diff / lambda_denom - 1) < tol) {
          break;
        }
      }
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = zeta(k) * Lambda;
    }
  }
  else if (constr == "EVE") {
    arma::vec eigval(p);
    arma::mat eigvec(p, p);

    // vector to store largest eigenvalue for each K
    arma::vec L_eig(K);
    for (int k = 0; k < K; k++) {
      arma::eig_sym(eigval, eigvec, L.slice(k)); // values are in ascending order
      L_eig(k) = eigval(p - 1);
    }

    // initialize lamdba matrices as identity
    arma::cube Lambda(p, p, K);
    for (int k = 0; k < K; k++) {
      Lambda.slice(k).eye();
    }

    arma::mat Gamma(p, p); Gamma.zeros();

    // matrix for MM iterations (and its svd)
    arma::mat F(p, p);
    arma::mat P;
    arma::mat C;
    arma::vec b;

    // vector to store largest eigenvalues of lambda matrices
    arma::vec Lambda_eig(K);

    // working vectors for termination conditions, etc
    double detval = 0.0;
    arma::mat Gamma_old(p, p);
    arma::mat Gamma_old2(p, p);
    arma::cube Lambda_old(p, p, K);
    double Gamma_diff, Gamma_denom;
    double Lambda_diff, Lambda_denom;

    // main iteration loop
    for (int iter = 0; iter < iter_max; iter++) {
      F.zeros();
      Gamma_old = Gamma;
      Lambda_old = Lambda;
      // there needs to be an inner loop here to solve for Gamma
      for (int in = 0; in < iter_max; in++) {
        Gamma_old2 = Gamma;
        if (iter % 2 == 0) {
          for (int k = 0; k < K; k++) {
            F += Lambda.slice(k).i() * Gamma.t() * L.slice(k) - L_eig(k) * Lambda.slice(k).i() * Gamma.t();
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }
        else {
          // we need to find the largest eigenvalue of each Lambda(k);
          for (int k = 0; k < K; k++) {
            arma::eig_sym(eigval, eigvec, Lambda.slice(k));
            Lambda_eig(k) = eigval(p - 1);
            F += L.slice(k) * Gamma * Lambda.slice(k).i() - Lambda_eig(k) * L.slice(k) * Gamma;
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }

        // check termination
        if (in > 1) {
          Gamma_diff = arma::accu(arma::square(Gamma_old2 - Gamma));
          Gamma_denom = arma::accu(arma::square(Gamma_old2));

          if (abs(Gamma_diff / Gamma_denom - 1) < tol) {
            break;
          }
        }
      }
      // end of inner loop

      for (int k = 0; k < K; k++) {
        Lambda.slice(k) = arma::diagmat(Gamma.t() * L.slice(k) * Gamma);
        detval = arma::det(Lambda.slice(k));
        Lambda.slice(k) = Lambda.slice(k) / pow(detval, 1.0 / p);
      }

      // check termination
      if (iter > 1) {
        Gamma_diff = arma::accu(arma::square(Gamma_old - Gamma));
        Gamma_denom = arma::accu(arma::square(Gamma_old));

        Lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));
        Lambda_denom = arma::accu(arma::square(Lambda_old));

        if (abs(Gamma_diff / Gamma_denom - 1) < tol && abs(Lambda_diff / Lambda_denom - 1) < tol) {
          break;
        }
      }
    }
    // end outer loop
    //
    // compute optimal zeta and return resulting matrices
    double zeta = 0.0;
    for (int k = 0; k < K; k++) {
      zeta += arma::trace(Gamma * Lambda.slice(k) * Gamma.t() * L.slice(k));
    }
    zeta /= (n * p);

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = arma::symmatu(zeta * Gamma * Lambda.slice(k) * Gamma.t());
    }
  }
  else if (constr == "VVE") {
    arma::vec eigval(p);
    arma::mat eigvec(p, p);

    // vector to store largest eigenvalue for each K
    arma::vec L_eig(K);
    for (int k = 0; k < K; k++) {
      arma::eig_sym(eigval, eigvec, L.slice(k)); // values are in ascending order
      L_eig(k) = eigval(p - 1);
    }

    // initialize lamdba matrices as identity
    arma::cube Lambda(p, p, K);
    for (int k = 0; k < K; k++) {
      Lambda.slice(k).eye();
    }

    arma::mat Gamma(p, p); Gamma.zeros();

    // matrix for MM iterations (and its svd)
    arma::mat F(p, p);
    arma::mat P;
    arma::mat C;
    arma::vec b;

    // vector to store largest eigenvalues of lambda matrices
    arma::vec Lambda_eig(K);

    // working vectors for termination conditions, etc
    double detval = 0.0;
    arma::mat Gamma_old(p, p);
    arma::mat Gamma_old2(p, p);
    arma::cube Lambda_old(p, p, K);
    double Gamma_diff, Gamma_denom;
    double Lambda_diff, Lambda_denom;

    // main iteration loop
    for (int iter = 0; iter < iter_max; iter++) {
      F.zeros();
      Gamma_old = Gamma;
      Lambda_old = Lambda;
      // there needs to be an inner loop here to solve for Gamma
      for (int in = 0; in < iter_max; in++) {
        Gamma_old2 = Gamma;
        if (iter % 2 == 0) {
          for (int k = 0; k < K; k++) {
            F += Lambda.slice(k).i() * Gamma.t() * L.slice(k) - L_eig(k) * Lambda.slice(k).i() * Gamma.t();
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }
        else {
          // we need to find the largest eigenvalue of each Lambda(k);
          for (int k = 0; k < K; k++) {
            arma::eig_sym(eigval, eigvec, Lambda.slice(k));
            Lambda_eig(k) = eigval(p - 1);
            F += L.slice(k) * Gamma * Lambda.slice(k).i() - Lambda_eig(k) * L.slice(k) * Gamma;
          }
          arma::svd(P, b, C, F);
          Gamma = C * P.t();
        }

        // check termination
        if (in > 1) {
          Gamma_diff = arma::accu(arma::square(Gamma_old2 - Gamma));
          Gamma_denom = arma::accu(arma::square(Gamma_old2));

          if (abs(Gamma_diff / Gamma_denom - 1) < tol) {
            break;
          }
        }
      }
      // end of inner loop

      for (int k = 0; k < K; k++) {
        Lambda.slice(k) = arma::diagmat(Gamma.t() * L.slice(k) * Gamma);
        detval = arma::det(Lambda.slice(k));
        Lambda.slice(k) = Lambda.slice(k) / pow(detval, 1.0 / p);
      }

      // check termination
      if (iter > 1) {
        Gamma_diff = arma::accu(arma::square(Gamma_old - Gamma));
        Gamma_denom = arma::accu(arma::square(Gamma_old));

        Lambda_diff = arma::accu(arma::square(Lambda_old - Lambda));
        Lambda_denom = arma::accu(arma::square(Lambda_old));

        if (abs(Gamma_diff / Gamma_denom - 1) < tol && abs(Lambda_diff / Lambda_denom - 1) < tol) {
          break;
        }
      }
    }
    // end outer loop

    // compute optimal zeta and return resulting matrices
    arma::vec zeta(K);
    for (int k = 0; k < K; k++) {
      zeta(k) = arma::trace(Gamma * Lambda.slice(k) * Gamma.t() * L.slice(k)) / (p * n_k(k));
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = zeta(k) * Gamma * Lambda.slice(k) * Gamma.t();
    }
  }
  else if (constr == "EEV") {
    // based on the notes, there is no need to calculate the zetas
    arma::cube Gamma(p, p, K);
    arma::mat Lambda(p, p); Lambda.zeros();

    // work arrays
    arma::mat eigvec(p, p);
    arma::vec eigval(p);

    for (int k = 0; k < K; k++) {
      arma::eig_sym(eigval, eigvec, L.slice(k));
      Gamma.slice(k) = eigvec;
      // TODO: should this be reversed...?
      // Lin's paper says that this diagonal matrix should have the eigenvalues
      // in _descending_ order, but armadillo returns them in ascending order...
      // this works and gives decent results for the marginalization method,
      // however
      Lambda += arma::diagmat(eigval);
    }

    for (int k = 0; k < K; k++) {
      Sigmas.slice(k) = Gamma.slice(k) * Lambda * Gamma.slice(k).t() / n;
    }

  }

  return Sigmas;
}

// Compute the Q2 function
// [[Rcpp::export]]
double Q2(arma::mat x, arma::mat z, arma::mat w, NumericVector sigmas, arma::mat mus, arma::vec grp, arma::umat Ru) {
  int K = mus.n_rows, n = x.n_rows, p = x.n_cols, M = Ru.n_rows;
  arma::cube Sigmas = to_array(sigmas, K, p);
  arma::mat ans(n, K);
  for (int m=0; m<M; m++) {
    int g = m+1;
    // get obs for this missingness pattern
    arma::uvec oidx = arma::find(Ru.row(m) == 1);
    arma::uvec gidx = arma::find(grp == g);
    int pg = oidx.size();
    int ng = gidx.size();
    arma::mat xg = x.submat(gidx, oidx);
    for (int k=0; k<K; k++) {
      arma::vec muk = mus.row(k).t();
      arma::mat sigmak = Sigmas.slice(k);
      // get mu, cholesky decom of sigma for this missingness pattern
      arma::vec mukg = muk.elem(oidx);
      arma::mat sigmakg = sigmak.submat(oidx, oidx);
      arma::mat Rkg(pg,pg);
      bool success = arma::chol(Rkg, sigmakg);
      if (!success) {
        Rkg = arma::chol(fix_var(sigmakg));
      }
      double logDet = 2.0 * sum(arma::log(Rkg.diag()));
      arma::vec maha = mahalanobis(xg, mukg, Rkg, true);
      for (int i=0; i<ng; i++) {
        int idx = gidx(i);
        ans(idx, k) = z(idx,k) / 2.0 * ( - logDet - w(idx,k) * maha(i) );
      }
    }
  }
  double out = 0;
  for (int k=0; k<K; k++) {
    for (int i=0; i<n; i++) {
      out += ans(i,k);
    }
  }
  return out;
}

// modification of k(m)means++ to allow for existing means chosen from the labeled observations
//           x : the data with unobserved labels
//       means : vector of existing means (may be empty)
// known_means : number of rows in `means` which come from labeled observations
// [[Rcpp::export]]
arma::mat get_init_centers(arma::mat x, arma::mat means, size_t known_means) {
  int n, p, nclusters, divisor;
  size_t i, j, k;
  double diff, total_weight, u, random_tracker, new_weight;

  n = x.n_rows;
  p = x.n_cols;
  nclusters = means.n_rows;

  arma::vec weight(n, arma::fill::zeros);
  // calculate distance between each observation and the current means
  for (i = 0; i < n; i++) {
    for (k = 0; k < known_means; k++) {
      new_weight = 0.0;
      divisor = 0;
      for (j = 0; j < p; j++) {
        if (isfinite(x(i, j)) && isfinite(means(k, j))) {
          divisor++;
          diff = x(i, j) - means(k, j);
          // weight(i) += diff * diff;
          new_weight += diff * diff;
        }
      }
      if (divisor > 1)
        new_weight /= (double) divisor;
        // weight(i) /= (double) divisor;
      if (k == 0 || new_weight < weight(i))
        weight(i) = new_weight;
    }
  }

  // choose remaining means
  for (k = known_means; k < nclusters; k++) {
    total_weight = 0.0;
    for (i = 0; i < n; i++)
      total_weight += weight(i);

    // select observation with probability proportional to total distance
    u = R::runif(0, 1);
    random_tracker = 0.0;
    for (i = 0; i < n; i++) {
      random_tracker += weight(i) / total_weight;
      if (random_tracker >= u)
        break;
    }

    // set new cluster mean to chosen observation
    means.row(k) = x.row(i);

    // recalculate weights if necessary
    if (k + 1 < nclusters) {
      for (i = 0; i < n; i++) {
        new_weight = 0.0;
        divisor = 0;
        for (j = 0; j < p; j++) {
          if (isfinite(x(i, j)) && isfinite(means(k, j))) {
            divisor++;
            diff = x(i, j) - means(k, j);
            new_weight += diff * diff;
          }
        }
        if (divisor > 1)
          new_weight /= (double) divisor;
        if (new_weight < weight(i))
          weight(i) = new_weight;
      }
    }
  }

  return(means);
}

// [[Rcpp::export]]
arma::vec assign_initial_partitions(arma::mat x, arma::mat means) {
  int n, p, nclusters, divisor;
  size_t i, j, k;
  double diff, dist, smallest_dist;

  n = x.n_rows;
  p = x.n_cols;
  nclusters = means.n_rows;
  arma::vec cluster(n, arma::fill::zeros);

  // assign each point to its closest center
  for (i = 0; i < n; i++) {
    smallest_dist = 0.0;
    for (k = 0; k < nclusters; k++) {
      divisor = 0;
      dist = 0.0;
      for (j = 0; j < p; j++) {
        if (isfinite(x(i, j)) && isfinite(means(k, j))) {
          divisor++;
          diff = x(i, j) - means(k, j);
          dist += diff * diff;
        }
      }
      if (divisor > 1)
        dist /= (double) divisor;
      if (k == 0 || dist < smallest_dist) {
        smallest_dist = dist;
        cluster(i) = k;
      }
    }
  }

  return(cluster);
}
