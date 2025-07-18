/*
  Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @file   gof.c
 * @author Ivan Kojadinovic, Marius Hofert
 * @date   2008, 2019
 *
 * @brief  - Two-sample goodness-of-fit test statistic of Remillard, Scaillet (2009)
 *         - Goodness-of-fit tests for copulas and EV copulas
 *           based on the parametric bootstrap and on a multiplier
 *           bootstrap
 *
 */

#include "nacopula.h"

#include <R.h>
#include <Rmath.h>

#include "An.h"
#include "gof.h"
#include "empcop.h"
#include "indepTests.h" // for progress bar


/**
 * @title Gaussian mixture kernel
 * @param x see gofMMDtest_c()
 * @param y see gofMMDtest_c()
 * @param d see gofMMDtest_c()
 * @param bandwidth2 see gofMMDtest_c()
 * @param K see gofMMDtest_c()
 * @return Gaussian mixture kernel
 * @author Marius Hofert and Ivan Kojadinovic
 */
double GaussMixKernel(double *x, double *y, int d, double *bandwidth2, int K) {
    int i;
    double diff = 0.0;
    double res = 0.0;
    for (i = 0; i < d; i++)
	diff += R_pow(x[i] - y[i], 2);
    diff /= (double)d; /* average norm */
    for (i = 0; i < K; i++)
	res += exp(diff / (-2.0 * bandwidth2[i]));
    return(res / (double)K); /* mean over all Gaussian kernels */
}

/**
 * @title MMD aggregated two-sample test of Schrab et al. (2024), a multiplier
 *        unbiased-MMD^2 based two-sample bootstrap test
 * @param x (d, n)-matrix of observations (typically from a copula)
 * @param y (d, n)-matrix of observations (typically from a copula)
 * @param n sample size (of both x and y)
 * @param d dimension
 * @param N number of bootstrap replications
 * @param bandwidth2 squared bandwidths of the Gaussian mixture kernel
 * @param K length(bandwidth2)
 * @param MMD2 unbiased MMD^2 (the one which removes the diagonal in
 *        each of the three double-sums, so the 2n-diagonal of K and both
 *        n-diagonals of the (2,1) and (1,2) submatrices of the (2n, 2n) kernel
 *        matrix)
 * @param MMD2H0 bootstrapped test statistics
 * @author Marius Hofert and Ivan Kojadinovic
 */
void gofMMDtest_c(double *x, double *y, int *n, int *d, int *N, double *bandwidth2,
		  int *K, double *MMD2, double *MMD2H0)
{
    int i, l, b;
    double *h = R_Calloc((*n) * (*n), double);
    double *eps = R_Calloc((*n), double);
    double s = 0.0;
    double nn1 = (*n) * (*n - 1);

    /* Compute h and the realized test statistic using the symmetry of the kernel */
    for (i = 0; i < *n; i++)
	for (l = 0; l < i; l++) {
	    h[i + (*n) * l] =
		  GaussMixKernel(&x[(*d) * i], &x[(*d) * l], *d, bandwidth2, *K)
		+ GaussMixKernel(&y[(*d) * i], &y[(*d) * l], *d, bandwidth2, *K)
		- GaussMixKernel(&x[(*d) * i], &y[(*d) * l], *d, bandwidth2, *K)
		- GaussMixKernel(&x[(*d) * l], &y[(*d) * i], *d, bandwidth2, *K);
	    h[l + (*n) * i] = h[i + (*n) * l];
	    s += 2.0 * h[i + (*n) * l];
	}
    *MMD2 = s / nn1;

    GetRNGstate();

    /* Generate N realizations under the null */
    /* *p.value = 0.0; => we do that in R */
    for (b = 0; b < *N; b++) {
	/* Generate n iid random variates from U({-1,1}) */
	for (i = 0; i < *n; i++)
	    eps[i] = (unif_rand() < 0.5) ? -1.0 : 1.0 ;

	/* Realization number b: see Gretton et al. (2024, (11)) */
	/* Below we exploit the symmetry of h */
	MMD2H0[b] = 0.0;
	for (i = 0; i < *n; i++)
	    for (l = 0; l < i; l++)
		if (i != l)
		    MMD2H0[b] += 2.0 * eps[i] * eps[l] * h[i + (*n) * l];
	MMD2H0[b] /= nn1;

	/* /\* Update p-value *\/ => we do that in R */
	/* *p.value += (MMD2H0[b] >= *MMD2); */
    }
    /* *p.value /= (double)(*N); => we do that in R */

    PutRNGstate();

    R_Free(h);
    R_Free(eps);
}

/**
 * @title Two-sample goodness-of-fit test statistic of Remillard, Scaillet (2009)
 * @param u1_ (n1, d) matrix of copula observations
 * @param u2_ (n2, d) matrix of copula observations
 * @return Numeric(1) giving the two-sample goodness-of-fit test statistic
 * @author Marius Hofert
 */
SEXP gofT2stat_c(SEXP u1_, SEXP u2_) {
    /* Setup */
    int *dims1 = INTEGER(getAttrib(u1_, R_DimSymbol));
    int n1 = dims1[0];
    int d = dims1[1];
    int *dims2 = INTEGER(getAttrib(u2_, R_DimSymbol));
    int n2 = dims2[0];
    double* u1 = REAL(u1_);
    double* u2 = REAL(u2_);
    int k1, k2, j; /* k1, k2 run over all rows, j runs over all cols */
    double p, res1, res2, res3; /* for sub-results */
    SEXP res = PROTECT(allocVector(REALSXP, 1)); /* final result object ('numeric(1)') */

    /* Main (note: u1[k, j] = u1[n1 * j + k] (u1 is stored in column-major order)) */

    /* Part 1: u1 with u1 */
    res1 = 0.0;
    for(k1 = 0; k1 < n1; k1++) {
	for(k2 = 0; k2 < n1; k2++) {
	    p = 1.0;
	    for(j = 0; j < d; j++) {
		p *= 1.0 - MAX(u1[n1 * j + k1], u1[n1 * j + k2]);
	    }
	    res1 += p;
        }
    }

    /* Part 2: u1 with u2 */
    res2 = 0.0;
    for(k1 = 0; k1 < n1; k1++) {
	for(k2 = 0; k2 < n2; k2++) {
	    p = 1.0;
	    for(j = 0; j < d; j++) {
		p *= 1.0 - MAX(u1[n1 * j + k1], u2[n2 * j + k2]);
	    }
	    res2 += p;
        }
    }

    /* Part 3: u2 with u2 */
    res3 = 0.0;
    for(k1 = 0; k1 < n2; k1++) {
	for(k2 = 0; k2 < n2; k2++) {
	    p = 1.0;
	    for(j = 0; j < d; j++) {
		p *= 1.0 - MAX(u2[n2 * j + k1], u2[n2 * j + k2]);
	    }
	    res3 += p;
        }
    }

    /* Return */
    REAL(res)[0] = ((res1 / (n1 * n1) - 2 * res2 / (n1 * n2) + res3 / (n2 * n2)) / (1.0/n1 + 1.0/n2));
    UNPROTECT(1);
    return res;
}

/**
 * @title Cramer-von Mises test statistic
 * @param n sample size
 * @param p dimension
 * @param U pseudo-observations
 * @param Ctheta values of the fitted copula at U
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises(int *n, int *p, double *U, double *Ctheta,
		     double *stat) {
  double s = 0.;
  for(int k = 0; k < *n; k++) {
    double diff = multCn(U, *n, *p, U, *n, k, 0.0) - Ctheta[k];
    s += diff * diff;
  }
  *stat = s;
}

/**
 * @title Cramer-von Mises test statistic (grid version)
 * @param p dimension
 * @param U pseudo-observations
 * @param n sample size
 * @param V grid
 * @param m grid size
 * @param Ctheta values of the fitted copula at V
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_grid(int *p, double *U, int *n, double *V, int *m,
			  double *Ctheta, double *stat) {
  double s = 0.;
  for (int k = 0; k < *m; k++) {
    double diff = multCn(U, *n, *p, V, *m, k, 0.0) - Ctheta[k];
    s +=  diff * diff;
  }
  *stat = s * (*n) / (*m);
}

/**
 * @title Multiplier bootstrap for the GoF testing: R's  gofCopula(*, simulation="mult")
 * @param p dimension
 * @param U pseudo-observations (n x p)
 * @param n sample size
 * @param G grid (g x p)
 * @param g grid size
 * @param b bandwidth for partial derivative estimation
 * @param influ influence matrix (g x n)
 * @param denom "denominator" n-vector for Rn approach
 * @param N number of multiplier replications
 * @param s0 N replications of the test statistic
 * @author Ivan Kojadinovic
 */
void multiplier(int *p, double *U, int *n, double *G, int *g, double *b,
		double *influ, double *denom, int *N, double *s0, int *verbose)
{
  int i, j, k, l, ind;
  double invsqrtn = 1.0/sqrt(*n);
  double *influ_mat = R_Calloc((*n) * (*g), double),
    *v1  = R_Calloc(*p, double), *v2 = R_Calloc(*p, double),
    *der = R_Calloc(*p, double);

  /* influence matrix */
  for (j = 0; j < *g; j++) /* loop over the grid points */
    {
      /* derivatives wrt args */
      for (k = 0; k < *p; k++)
        {
	  v1[k] = G[j + k * (*g)];
	  v2[k] = v1[k];
        }
      for (k = 0; k < *p; k++)
        {
	  v1[k] += *b;
	  v2[k] -= *b;
	  der[k] = der_multCn(U, *n, *p, v1, v2, 2.0 * (*b));
	  v1[k] -= *b;
	  v2[k] += *b;
	}

      for (i = 0; i < *n; i++) /* loop over the data */
        {
	  influ_mat[i + j * (*n)] = 0.0;
	  ind = 1;
	  for (k = 0; k < *p; k++)
	    {
	      ind *= (U[i + k * (*n)] <= G[j + k * (*g)]);
	      influ_mat[i + j * (*n)] -= der[k] * (U[i + k * (*n)] <= G[j + k * (*g)]);
	    }
	  influ_mat[i + j * (*n)] += ind; /* - influ[j + i * (*g)]; */
	  influ    [j + i * (*g)] *= invsqrtn;
	  influ_mat[i + j * (*n)] *= invsqrtn;
	}
    }
  R_Free(v1);
  R_Free(v2);
  R_Free(der);

  double *random = R_Calloc(*n, double);

  GetRNGstate();
  /* generate N approximate realizations */
  for (l=0; l < *N; l++)
    {
      /* generate n variates */
      double mean = 0.0;
      for (i=0; i<*n; i++)
	{
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= *n;

      /* realization number l */
      s0[l] = 0.0;
      for (j=0; j<*g; j++)
	{
	  double process = 0.0;
	  for (i=0; i<*n; i++)
	      process += ((random[i] - mean) * influ_mat[i + j * (*n)]
			  - random[i] * influ[j + i * (*g)]) / denom[j];
	  s0[l] += process * process;
	}
      s0[l] /= *g;

      /* update progress bar */
      if (*verbose)
	  progressBar(l, *N, 70);
    }

  PutRNGstate();

  R_Free(influ_mat);
  R_Free(random);
}

/// Goodness-of-fit tests for extreme-value copulas

/**
 * @title Cramer-von Mises test statistic based on the Pickands estimator
 *        used in the GOF for EV copulas
 *        stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_Pickands(int n, int m, double *S,
			      double *T, double *Atheta,
			      double *stat) {
  double t, Ac, Au, dc, du,
    invA0 = biv_invAP(n, S, T, 0.0);

  stat[0] = 0.0; stat[1] = 0.0;
  for (int i = 0; i < m; i++) {
      t = (double)i/(double)(m);
      Au = biv_invAP(n, S, T, t);
      Ac = Au - invA0 + 1.0; // correction
      du = 1 / Au - Atheta[i];
      dc = 1 / Ac - Atheta[i];
      stat[0] += dc * dc;
      stat[1] += du * du;
    }
  stat[0] = stat[0] * (double)(n) / (double)(m);
  stat[1] = stat[1] * (double)(n) / (double)(m);
}

/**
 * @title Cramer-von Mises test statistic based on the CFG estimator
 *        used in the GOF for EV copulas
 *        stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_CFG(int n, int m, double *S,
			 double *T, double *Atheta,
			 double *stat) {
  double t, Au, Ac, dc, du,
    logA0 = biv_logACFG(n, S, T, 0.0);

  stat[0] = 0.0; stat[1] = 0.0;
  for (int i = 0; i < m; i++) {
      t = (double) i / (double) (m);
      Au = biv_logACFG(n, S, T, t);
      Ac = Au - logA0; // correction
      dc = exp(Ac) - Atheta[i];
      du = exp(Au) - Atheta[i];
      stat[0] += dc * dc;
      stat[1] += du * du;
    }
  stat[0] = stat[0] * (double)(n)/(double)(m);
  stat[1] = stat[1] * (double)(n)/(double)(m);
}

/**
 * @title Wrapper for the Cramer-von Mises test statistics
 *        used in the GOF for EV copulas
 *        stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @param CFG if > 0 then CFG, else Pickands
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_Afun(int *n, int *m, double *S,
			  double *T, double *Atheta,
			  double *stat, int *CFG)
{
  if (*CFG)
    cramer_vonMises_CFG(*n, *m, S, T, Atheta, stat);
  else
    cramer_vonMises_Pickands(*n, *m, S, T, Atheta, stat);
}
