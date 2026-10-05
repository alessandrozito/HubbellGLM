#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// w_j = j^(1-sigma). sigma is fixed for a whole fit and only alpha changes
// between calls, so recomputing pow() inside the alpha loop (as the current
// polyseries_mean_single does) is pure repeated work.
static inline std::vector<double> make_weights(int max_size, double sigma) {
  std::vector<double> w(std::max(max_size, 1));
  const double p = 1.0 - sigma;
  w[0] = 0.0;
  if (p == 1.0) {                       // sigma == 0: j^1, skip pow entirely
    for (int j = 1; j < max_size; j++) w[j] = (double) j;
  } else {
    for (int j = 1; j < max_size; j++) w[j] = std::pow((double) j, p);
  }
  return w;
}

// mean and variance share the same loop; the IRLS step needs both at the same
// eta, so computing them together halves the passes over j.
// [[Rcpp::export]]
List polyseries_meanvar_fast(NumericVector size, NumericVector alpha,
                             double sigma, int nthreads = 1) {
  const int n = size.size();
  if (n != alpha.size()) stop("'size' and 'alpha' must have the same length.");
  int max_size = 1;
  for (int i = 0; i < n; i++) max_size = std::max(max_size, (int) size[i]);
  std::vector<double> w = make_weights(max_size, sigma);

  NumericVector mu(n), va(n);
  double *pmu = REAL(mu), *pva = REAL(va), *pw = w.data();
  const double *psz = REAL(size), *pal = REAL(alpha);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    const double a = pal[i];
    const int s = (int) psz[i];
    double m = 1.0, v = 0.0;
    for (int j = 1; j < s; j++) {
      const double t = a / (a + pw[j]);
      m += t;
      v += t * (1.0 - t);
    }
    pmu[i] = m;
    pva[i] = v;
  }
  return List::create(_["mean"] = mu, _["var"] = va);
}

// Mean only, for the calls that do not need the variance.
// [[Rcpp::export]]
NumericVector polyseries_mean_fast(NumericVector size, NumericVector alpha,
                                   double sigma, int nthreads = 1) {
  const int n = size.size();
  if (n != alpha.size()) stop("'size' and 'alpha' must have the same length.");
  int max_size = 1;
  for (int i = 0; i < n; i++) max_size = std::max(max_size, (int) size[i]);
  std::vector<double> w = make_weights(max_size, sigma);

  NumericVector mu(n);
  double *pmu = REAL(mu), *pw = w.data();
  const double *psz = REAL(size), *pal = REAL(alpha);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    const double a = pal[i];
    const int s = (int) psz[i];
    double m = 1.0;
    for (int j = 1; j < s; j++) m += a / (a + pw[j]);
    pmu[i] = m;
  }
  return mu;
}

// Invert mu(alpha) with safeguarded Newton.
//
// The derivative is free: with t_j = alpha/(alpha + w_j),
//   dmu/dalpha = sum_j w_j/(alpha + w_j)^2 = Var(alpha)/alpha,
// and Var is already accumulated in the same loop as mu. Newton therefore
// costs one pass per iteration and converges in a handful of steps, where
// uniroot's bisection needs ~40 - and uniroot ran from R, once per observation.
// [[Rcpp::export]]
NumericVector inv_polyseries_fast(NumericVector mu_target, NumericVector size,
                                  double sigma, double tol = 1e-10,
                                  int maxit = 100, int nthreads = 1) {
  const int n = mu_target.size();
  int max_size = 1;
  for (int i = 0; i < n; i++) max_size = std::max(max_size, (int) size[i]);
  std::vector<double> w = make_weights(max_size, sigma);

  NumericVector out(n);
  double *pout = REAL(out), *pw = w.data();
  const double *pmt = REAL(mu_target), *psz = REAL(size);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
#endif
  for (int i = 0; i < n; i++) {
    const int s = (int) psz[i];
    double mt = pmt[i];
    // Same clamping as the R implementation.
    if (mt <= 1.0 + 1e-4)          mt = 1.0 + 1e-4;
    else if (mt >= s - 1e-4)       mt = s - 1e-4;

    double lo = (mt - 1.0) / (s - mt);
    double hi = std::pow((double)(s - 1), 1.0 - sigma) * (mt - 1.0) / (s - mt) + 1e-4;
    if (!(hi > lo)) hi = lo + 1e-4;

    double a = 0.5 * (lo + hi);
    for (int it = 0; it < maxit; it++) {
      double m = 1.0, v = 0.0;
      for (int j = 1; j < s; j++) {
        const double t = a / (a + pw[j]);
        m += t;
        v += t * (1.0 - t);
      }
      const double f = m - mt;
      if (f > 0.0) hi = a; else lo = a;
      if (std::fabs(f) < tol * std::max(1.0, mt)) break;

      const double deriv = v / a;                  // dmu/dalpha
      double a_new = (deriv > 0.0) ? a - f / deriv : 0.5 * (lo + hi);
      // Fall back to bisection whenever Newton leaves the bracket.
      if (!(a_new > lo && a_new < hi) || !std::isfinite(a_new))
        a_new = 0.5 * (lo + hi);
      if (std::fabs(a_new - a) < 1e-15 * std::max(1.0, a)) { a = a_new; break; }
      a = a_new;
    }
    pout[i] = a;
  }
  return out;
}

// Invert the Dirichlet-process mean  mu(alpha) = alpha*(psi(alpha+n) - psi(alpha)).
//
// This is the one that actually dominates a fit: hubbell()'s variance() and
// dev.resids() both call it on every IRLS iteration, through a Vectorize()d
// uniroot, i.e. one R-level root-find per observation per iteration.
// The derivative is already known analytically (deriv_dirichlet_process), so
// Newton applies directly.
// [[Rcpp::export]]
NumericVector inv_mean_dp_fast(NumericVector mu_target, NumericVector size,
                               double tol = 1e-10, int maxit = 100,
                               int nthreads = 1) {
  const int n = mu_target.size();
  NumericVector out(n);
  double *pout = REAL(out);
  const double *pmt = REAL(mu_target), *psz = REAL(size);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    const double s = psz[i];
    double mt = pmt[i];
    if (mt <= 1.0 + 1e-4)    mt = 1.0 + 1e-4;
    else if (mt >= s - 1e-4) mt = s - 1e-4;

    double lo = (mt - 1.0) / (s - mt);
    double hi = (s - 1.0) * (mt - 1.0) / (s - mt) + 1e-4;
    if (!(hi > lo)) hi = lo + 1e-4;

    double a = 0.5 * (lo + hi);
    for (int it = 0; it < maxit; it++) {
      const double dg = R::digamma(a + s) - R::digamma(a);
      const double f  = a * dg - mt;
      if (f > 0.0) hi = a; else lo = a;
      if (std::fabs(f) < tol * std::max(1.0, mt)) break;

      const double d = dg + a * (R::trigamma(a + s) - R::trigamma(a));
      double a_new = (d > 0.0) ? a - f / d : 0.5 * (lo + hi);
      if (!(a_new > lo && a_new < hi) || !std::isfinite(a_new))
        a_new = 0.5 * (lo + hi);
      if (std::fabs(a_new - a) < 1e-15 * std::max(1.0, a)) { a = a_new; break; }
      a = a_new;
    }
    pout[i] = a;
  }
  return out;
}
