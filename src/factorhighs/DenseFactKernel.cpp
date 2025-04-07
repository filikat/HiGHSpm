#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "ReturnValues.h"
#include "util/HighsRandom.h"

// Dense Factorization kernel

std::pair<int, double> maxInCol(int j, int n, int m, double* A, int lda) {
  // Given the symemtric matrix A, of size nxn, accessed with leading dimension
  // lda, in upper triangular format, ignoring rows 0:j-1, find the maximum in
  // row/col m.

  double maxval = -1.0;
  int r = -1;
  // scan column m, rows j:m-1
  for (int row = j; row < m; ++row) {
    double val = std::abs(A[row + lda * m]);
    if (val > maxval) {
      maxval = val;
      r = row;
    }
  }
  // scan row m, columns m+1:n
  for (int col = m + 1; col < n; ++col) {
    double val = std::abs(A[m + lda * col]);
    if (val > maxval) {
      maxval = val;
      r = col;
    }
  }
  return {r, maxval};
}

void staticReg(double& pivot, int sign, double& regul) {
  // apply static regularization

  double old_pivot = pivot;
  if (sign > 0)
    pivot += kDualStaticRegularization;
  else
    pivot -= kPrimalStaticRegularization;
  regul = pivot - old_pivot;
}

bool blockBunchKaufman(int j, int n, double* A, int lda, int* swaps, int* sign,
                       double thresh, double* regul, int sn, int bl,
                       double max_in_R) {
  // Perform Bunch-Kaufman pivoting within a block of the supernode (see Schenk,
  // Gartner, ETNA 2006).
  // It works only for upper triangular A.
  // Return true if 2x2 pivot should be used, false otherwise.
  // Swap of columns may be performed.
  // Regularization of pivot may be performed.

  Clock clock;
  bool flag_2x2 = false;

  /*if (n == 1) {
    // only pivot in the block or supernode
    if (sign[j] * A[j + lda * j] < 0)
      DataCollector::get()->setWrongSign(A[j + lda * j]);

    double old_pivot = A[j + lda * j];
    staticReg(A[j + lda * j], sign[j], regul[j]);
    if (std::abs(A[j + lda * j]) < kAlphaBK * max_in_R * 1e-6) {
      A[j + lda * j] = sign[j] * kAlphaBK * max_in_R * 1e-6;
      DataCollector::get()->countRegPiv();
      // #ifdef PRINT_REGULARIZATION
      printf("%4d %2d %2d: pivot %8.1e set to %8.1e\n", sn, bl, j, old_pivot,
             A[j + lda * j]);
      // #endif
    }
    regul[j] = A[j + lda * j] - old_pivot;
    DataCollector::get()->setMaxReg(std::abs(regul[j]));
  }*/

#if 1
  // Find largest diagonal entry in the residual part of the block
  int ind_max_diag = -1;
  double max_diag = -1.0;
  for (int i = j; i < n; ++i) {
    double val = std::abs(A[i + lda * i]);
    if (val > max_diag) {
      max_diag = val;
      ind_max_diag = i;
    }
  }
  assert(ind_max_diag >= j);

  // put column with max pivot as first column
  swapCols('U', n, A, lda, j, ind_max_diag, swaps, sign);
#endif

  // Max in column j of diagonal block
  auto res = maxInCol(j, n, j, A, lda);
  double gamma_j = res.second;
  int r = res.first;
  double Ajj = sign[j] > 0 ? A[j + lda * j] + kDualStaticRegularization
                           : A[j + lda * j] - kPrimalStaticRegularization;

  if (std::max(std::abs(Ajj), gamma_j) <= thresh || sign[j] * Ajj < 0 ||
      j == n - 1) {
    // Must accept current pivot
    double old_pivot = A[j + lda * j];
    staticReg(A[j + lda * j], sign[j], regul[j]);

    if (sign[j] * A[j + lda * j] < 0) {
      DataCollector::get()->setWrongSign(A[j + lda * j]);
      //A[j + lda * j] = sign[j] * thresh;
    }

    if (std::max(std::abs(Ajj), gamma_j) < thresh) {
      // perturbe pivot
      A[j + lda * j] = sign[j] * thresh;
      DataCollector::get()->countRegPiv();
#ifdef PRINT_REGULARIZATION
      printf("%4d %2d %2d: pivot %8.1e set to %8.1e\n", sn, bl, j, old_pivot,
             A[j + lda * j]);
#endif
    }
    regul[j] = A[j + lda * j] - old_pivot;
    DataCollector::get()->setMaxReg(std::abs(regul[j]));

  } else {
    // Max in column r of diagonal block
    assert(r >= 0);
    res = maxInCol(j, n, r, A, lda);
    double gamma_r = res.second;
    double Arr = sign[r] > 0 ? A[r + lda * r] + kDualStaticRegularization
                             : A[r + lda * r] - kPrimalStaticRegularization;

    if ((std::abs(Ajj) >= kAlphaBK * gamma_j ||
         std::abs(Ajj) * gamma_r >= kAlphaBK * gamma_j * gamma_j)) {
      // Accept current pivot
      staticReg(A[j + lda * j], sign[j], regul[j]);

      if (sign[j] * A[j + lda * j] < 0)
        DataCollector::get()->setWrongSign(A[j + lda * j]);

    } else if (std::abs(Arr) >= kAlphaBK * gamma_r) {
      // Use pivot r

      swapCols('U', n, A, lda, j, r, swaps, sign);
      staticReg(A[j + lda * j], sign[j], regul[j]);

      if (sign[j] * A[j + lda * j] < 0)
        DataCollector::get()->setWrongSign(A[j + lda * j]);

    } else {
      // Use 2x2 pivot (j,r)

      swapCols('U', n, A, lda, j + 1, r, swaps, sign);
      flag_2x2 = true;
      staticReg(A[j + lda * j], sign[j], regul[j]);
      staticReg(A[j + 1 + lda * (j + 1)], sign[j + 1], regul[j + 1]);

      if (sign[j] * A[j + lda * j] < 0)
        DataCollector::get()->setWrongSign(A[j + lda * j]);
      if (sign[j + 1] * A[j + 1 + lda * (j + 1)] < 0)
        DataCollector::get()->setWrongSign(A[j + 1 + lda * (j + 1)]);
    }
  }

  DataCollector::get()->sumTime(kTimeDenseFact_pivoting, clock.stop());
  return flag_2x2;
}

double regularizePivot(double pivot, double thresh, const int* sign,
                       const double* A, int lda, int j, int n, char uplo,
                       int sn, int bl) {
  // add static regularization
  if (sign[j] == 1)
    pivot += kDualStaticRegularization;
  else
    pivot -= kPrimalStaticRegularization;

  double s = (double)sign[j];
  double old_pivot = pivot;

  double spivot = s * pivot;
  double K = 1e12;

  bool adjust = false;
  bool modified_pivot = false;

  if (spivot <= thresh && spivot >= -thresh) {
    // small pivot, lift to thresh
    pivot = s * thresh;
    adjust = true;
    modified_pivot = true;
#ifdef PRINT_REGULARIZATION
    printf("%2d, %2d, %2d: small pivot %e, with sign %d, set to %e\n", sn, bl,
           j, old_pivot, sign[j], pivot);
#endif

  } else if (spivot < -thresh && spivot >= -thresh * K) {
    // wrong sign, lift more
    pivot = s * thresh * 10;
    adjust = true;
    modified_pivot = true;
#ifdef PRINT_REGULARIZATION
    printf("%2d, %2d, %2d: wrong pivot %e, with sign %d, set to %e\n", sn, bl,
           j, old_pivot, sign[j], pivot);
#endif

  } else if (spivot < -thresh * K) {
    // pivot is completely lost
    pivot = s * 1e100;
    modified_pivot = true;
#ifdef PRINT_REGULARIZATION
    printf("%2d, %2d, %2d: disaster pivot %e, with sign %d, set to %e\n", sn,
           bl, j, old_pivot, sign[j], pivot);
#endif
  }

  if (adjust) {
    // compute the minimum pivot required to keep the diagonal of the
    // current block acceptable:
    // b is column below pivot, d is diagonal of block, p is pivot
    // we want d_k - b_k^2 / p \ge thresh
    // i.e.
    // p \ge b_k^2 / (d_k - thresh)
    //
    double required_pivot = pivot;
    for (int k = j + 1; k < n; ++k) {
      double bk = uplo == 'L' ? A[k + j * lda] : A[j + k * lda];
      double dk = A[k + k * lda];

      // if pivot and dk have different sign, skip
      double sk = sign[k];
      if (s * sk < 0) continue;

      double temp = (dk - sk * thresh);
      temp = (bk * bk) / temp;

      if (s > 0) {
        required_pivot = std::max(required_pivot, temp);
        required_pivot = std::min(required_pivot, 1e100);
      } else {
        required_pivot = std::min(required_pivot, temp);
        required_pivot = std::max(required_pivot, -1e100);
      }
    }

    if (required_pivot != pivot) {
      modified_pivot = true;

      if (s > 0)
        pivot = std::max(pivot, required_pivot);
      else
        pivot = std::min(pivot, required_pivot);

#ifdef PRINT_REGULARIZATION
      printf("\t%2d, %2d, %2d: adjust %e to %e\n", sn, bl, j, old_pivot, pivot);
#endif
    }
  }

  if (modified_pivot) DataCollector::get()->countRegPiv();

  return pivot;
}

int denseFactK(char uplo, int n, double* A, int lda, int* pivot_sign,
               double thresh, double* regul, int* swaps, double* pivot_2x2,
               int sn, int bl, double max_in_R) {
  // ===========================================================================
  // Factorization kernel
  // Matrix A is in format F
  // BLAS calls: dscal, dcopy, dger
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n) {
    printf("\ndenseFactK: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

#ifdef FINE_TIMING
  Clock clock;
#endif

  // ===========================================================================
  // LOWER TRIANGULAR
  // ===========================================================================
  // No pivoting performed.
  // This is kept only for reference, to use the formats F and FP.
  if (uplo == 'L') {
    // allocate space for copy of col
    std::vector<double> temp(n - 1);

    for (int j = 0; j < n; ++j) {
      // diagonal element
      double Ajj = A[j + lda * j];

      if (std::isnan(Ajj)) {
        printf("\ndenseFactK: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // add regularization
      double old_pivot = Ajj;
      Ajj =
          regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, sn, bl);
      regul[j] = Ajj - old_pivot;
      DataCollector::get()->setMaxReg(std::abs(regul[j]));

      // save diagonal element
      A[j + lda * j] = Ajj;

      const int M = n - j - 1;
      if (M > 0) {
        // make copy of column
        callAndTime_dcopy(M, &A[j + 1 + j * lda], 1, temp.data(), 1);

        // scale column j
        callAndTime_dscal(M, 1.0 / Ajj, &A[j + 1 + j * lda], 1);

        // update rest of the matrix
        callAndTime_dger(M, M, -1.0, temp.data(), 1, &A[j + 1 + j * lda], 1,
                         &A[j + 1 + (j + 1) * lda], lda);
      }
    }
  }

  // ===========================================================================
  // UPPER TRIANGULAR
  // ===========================================================================
  else {
    if (!swaps || !pivot_2x2) {
      printf("\ndenseFactK: invalid input\n");
      return kRetInvalidInput;
    }

    // initialize order of pivots
    for (int i = 0; i < n; ++i) swaps[i] = i;

    // allocate space for copy of col(s)
    std::vector<double> temp(n - 1);
    std::vector<double> temp2(n - 1);

    int step = 1;

    for (int j = 0; j < n; j += step) {
      bool flag_2x2 = false;

#ifdef PIVOTING
      flag_2x2 = blockBunchKaufman(j, n, A, lda, swaps, pivot_sign, thresh,
                                   regul, sn, bl, max_in_R);
#endif

      // cannot do 2x2 pivoting on last column
      assert(j < n - 1 || flag_2x2 == false);

      if (!flag_2x2) {
        // 1x1 pivots
        step = 1;

        // diagonal element
        double Ajj = A[j + lda * j];

        if (std::isnan(Ajj)) {
          printf("\ndenseFactK: invalid pivot %e\n", Ajj);
          return kRetInvalidPivot;
        }

#ifndef PIVOTING
        // add regularization
        double old_pivot = Ajj;
        Ajj = regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, sn,
                              bl);
        regul[j] = Ajj - old_pivot;
        DataCollector::get()->setMaxReg(std::abs(regul[j]));
#endif

        // save reciprocal of pivot
        A[j + lda * j] = 1.0 / Ajj;

        const int M = n - j - 1;
        if (M > 0) {
          // make copy of row
          callAndTime_dcopy(M, &A[j + (j + 1) * lda], lda, temp.data(), 1);

          // scale row j
          callAndTime_dscal(M, 1.0 / Ajj, &A[j + (j + 1) * lda], lda);

          // update rest of the matrix
          callAndTime_dger(M, M, -1.0, temp.data(), 1, &A[j + (j + 1) * lda],
                           lda, &A[j + 1 + (j + 1) * lda], lda);
        }
      } else {
        // 2x2 pivots
        step = 2;

        // diagonal pivot elements
        const double d1 = A[j + lda * j];
        const double d2 = A[j + 1 + lda * (j + 1)];

        if (std::isnan(d1) || std::isnan(d2)) {
          printf("\ndenseFactK: invalid pivot %e %e\n", d1, d2);
          return kRetInvalidPivot;
        }

        // off-diagonal pivot element
        const double offd = A[j + lda * (j + 1)];
        A[j + lda * (j + 1)] = 0.0;

        // compute coefficients of 2x2 inverse
        const double denom = d1 * d2 - offd * offd;
        const double i_d1 = d2 / denom;
        const double i_d2 = d1 / denom;
        const double i_off = -offd / denom;

        // save them in place of pivots
        A[j + lda * j] = i_d1;
        A[j + 1 + lda * (j + 1)] = i_d2;
        pivot_2x2[j] = i_off;

        const int M = n - j - 2;
        if (M > 0) {
          double* r1 = &A[j + (j + 2) * lda];
          double* r2 = &A[j + 1 + (j + 2) * lda];

          // make a copy of first row
          callAndTime_dcopy(M, r1, lda, temp.data(), 1);

          // make a copy of second row
          callAndTime_dcopy(M, r2, lda, temp2.data(), 1);

          // solve rows j,j+1
          callAndTime_dscal(M, i_d1, r1, lda);
          callAndTime_daxpy(M, i_off, temp2.data(), 1, r1, lda);
          callAndTime_dscal(M, i_d2, r2, lda);
          callAndTime_daxpy(M, i_off, temp.data(), 1, r2, lda);

          // update rest of the matrix
          callAndTime_dger(M, M, -1.0, temp.data(), 1, r1, lda,
                           &A[j + 2 + (j + 2) * lda], lda);
          callAndTime_dger(M, M, -1.0, temp2.data(), 1, r2, lda,
                           &A[j + 2 + (j + 2) * lda], lda);
        }

        DataCollector::get()->count2x2();
      }
    }
  }

#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeDenseFact_kernel, clock.stop());
#endif

  return kRetOk;
}
