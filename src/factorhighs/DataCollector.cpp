#include "DataCollector.h"

namespace highspm {

// instance of DataCollector
DataCollector* DataCollector::ptr_ = nullptr;

DataCollector::DataCollector() {
#ifdef DATA_COLLECTION
  times_.resize(kTimeSize);
  blas_calls_.resize(kTimeBlasEnd - kTimeBlasStart + 1);
#endif
}

DataCollector* DataCollector::get() { return ptr_; }

void DataCollector::start() {
  if (!ptr_) ptr_ = new DataCollector();
}
void DataCollector::destruct() {
  delete ptr_;
  ptr_ = nullptr;
}
void DataCollector::append() {
  // add an empty IterData object to the record
  iter_data_record_.push_back(IterData());
}
IterData& DataCollector::back() {
  // access most recent record of data
  return iter_data_record_.back();
}

void DataCollector::sumTime(TimeItems i, double t) {
#ifdef DATA_COLLECTION
  // Keep track of times and blas calls.
  std::lock_guard<std::mutex> lock(times_mutex_);
  times_[i] += t;
#ifdef BLAS_TIMING
  if (i >= kTimeBlasStart && i <= kTimeBlasEnd)
    ++blas_calls_[i - kTimeBlasStart];
#endif
#endif
}

void DataCollector::countSolves() {
#ifdef DATA_COLLECTION
  std::lock_guard<std::mutex> lock(times_mutex_);
  ++total_solves_;
  ++back().num_solves;
#endif
}

void DataCollector::setExtremeEntries(double minD, double maxD, double minoffD,
                                      double maxoffD) {
#ifdef DATA_COLLECTION
  // Store max and min entries of D and L.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  back().minD = std::min(back().minD, minD);
  back().maxD = std::max(back().maxD, maxD);
  back().minL = std::min(back().minL, minoffD);
  back().maxL = std::max(back().maxL, maxoffD);
#endif
}

void DataCollector::countRegPiv() {
#ifdef DATA_COLLECTION
  // Increase the number of dynamically regularized pivots.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_reg_piv;
#endif
}

void DataCollector::countSwap() {
#ifdef DATA_COLLECTION
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_swap;
#endif
}

void DataCollector::count2x2() {
#ifdef DATA_COLLECTION
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_2x2;
#endif
}

void DataCollector::setWrongSign(double p) {
#ifdef DATA_COLLECTION
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_wrong_sign;
  back().max_wrong_sign = std::max(back().max_wrong_sign, std::abs(p));
#endif
}

void DataCollector::setMaxReg(double new_reg) {
#ifdef DATA_COLLECTION
  // Keep track of maximum regularization used.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  back().max_reg = std::max(back().max_reg, new_reg);
#endif
}

void DataCollector::printTimes() const {
#ifdef COARSE_TIMING
  printf("----------------------------------------------------\n");
  printf("Analyse time            \t%8.4f\n", times_[kTimeAnalyse]);

#ifdef FINE_TIMING
  printf("\tMetis:                  %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseMetis],
         times_[kTimeAnalyseMetis] / times_[kTimeAnalyse] * 100);
  printf("\tTree:                   %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseTree],
         times_[kTimeAnalyseTree] / times_[kTimeAnalyse] * 100);
  printf("\tCounts:                 %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseCount],
         times_[kTimeAnalyseCount] / times_[kTimeAnalyse] * 100);
  printf("\tSupernodes:             %8.4f (%4.1f%%)\n", times_[kTimeAnalyseSn],
         times_[kTimeAnalyseSn] / times_[kTimeAnalyse] * 100);
  printf("\tReorder:                %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseReorder],
         times_[kTimeAnalyseReorder] / times_[kTimeAnalyse] * 100);
  printf("\tSn sparsity pattern:    %8.4f (%4.1f%%)\n",
         times_[kTimeAnalysePattern],
         times_[kTimeAnalysePattern] / times_[kTimeAnalyse] * 100);
  printf("\tRelative indices:       %8.4f (%4.1f%%)\n",
         times_[kTimeAnalyseRelInd],
         times_[kTimeAnalyseRelInd] / times_[kTimeAnalyse] * 100);
#endif

  printf("----------------------------------------------------\n");
  printf("Factorise time          \t%8.4f\n", times_[kTimeFactorise]);

#ifdef FINE_TIMING
  printf("\tPrepare:                %8.4f (%4.1f%%)\n",
         times_[kTimeFactorisePrepare],
         times_[kTimeFactorisePrepare] / times_[kTimeFactorise] * 100);
  printf("\tAssembly original:      %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseAssembleOriginal],
         times_[kTimeFactoriseAssembleOriginal] / times_[kTimeFactorise] * 100);
  printf("\tAssemble children in F: %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseAssembleChildrenFrontal],
         times_[kTimeFactoriseAssembleChildrenFrontal] /
             times_[kTimeFactorise] * 100);
  printf("\tAssemble children in C: %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseAssembleChildrenClique],
         times_[kTimeFactoriseAssembleChildrenClique] / times_[kTimeFactorise] *
             100);
  printf("\tDense factorisation:    %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseDenseFact],
         times_[kTimeFactoriseDenseFact] / times_[kTimeFactorise] * 100);
  printf("\t\tmain:           %8.4f\n", times_[kTimeDenseFact_main]);
  printf("\t\tSchur:          %8.4f\n", times_[kTimeDenseFact_schur]);
  printf("\t\tkernel:         %8.4f\n", times_[kTimeDenseFact_kernel]);
  printf("\t\tconvert:        %8.4f\n", times_[kTimeDenseFact_convert]);
  printf("\t\tpivoting:       %8.4f\n", times_[kTimeDenseFact_pivoting]);
  printf("\tTerminate:              %8.4f (%4.1f%%)\n",
         times_[kTimeFactoriseTerminate],
         times_[kTimeFactoriseTerminate] / times_[kTimeFactorise] * 100);
#endif

  printf("----------------------------------------------------\n");
  printf("Solve time              \t%8.4f (%d calls)\n", times_[kTimeSolve],
         total_solves_);
  printf("----------------------------------------------------\n");

#ifdef BLAS_TIMING

  double total_blas_time =
      times_[kTimeBlas_copy] + times_[kTimeBlas_axpy] + times_[kTimeBlas_scal] +
      times_[kTimeBlas_swap] + times_[kTimeBlas_gemv] + times_[kTimeBlas_trsv] +
      times_[kTimeBlas_tpsv] + times_[kTimeBlas_ger] + times_[kTimeBlas_trsm] +
      times_[kTimeBlas_syrk] + times_[kTimeBlas_gemm];

  printf("BLAS time               \t%8.4f\n", total_blas_time);
  printf("\tcopy:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_copy], times_[kTimeBlas_copy] / total_blas_time * 100,
         blas_calls_[kTimeBlas_copy - kTimeBlasStart]);
  printf("\taxpy:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_axpy], times_[kTimeBlas_axpy] / total_blas_time * 100,
         blas_calls_[kTimeBlas_axpy - kTimeBlasStart]);
  printf("\tscal:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_scal], times_[kTimeBlas_scal] / total_blas_time * 100,
         blas_calls_[kTimeBlas_scal - kTimeBlasStart]);
  printf("\tswap:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_swap], times_[kTimeBlas_swap] / total_blas_time * 100,
         blas_calls_[kTimeBlas_swap - kTimeBlasStart]);
  printf("\tgemv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_gemv], times_[kTimeBlas_gemv] / total_blas_time * 100,
         blas_calls_[kTimeBlas_gemv - kTimeBlasStart]);
  printf("\ttrsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_trsv], times_[kTimeBlas_trsv] / total_blas_time * 100,
         blas_calls_[kTimeBlas_trsv - kTimeBlasStart]);
  printf("\ttpsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_tpsv], times_[kTimeBlas_tpsv] / total_blas_time * 100,
         blas_calls_[kTimeBlas_tpsv - kTimeBlasStart]);
  printf("\tger:            \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_ger], times_[kTimeBlas_ger] / total_blas_time * 100,
         blas_calls_[kTimeBlas_ger - kTimeBlasStart]);
  printf("\ttrsm:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_trsm], times_[kTimeBlas_trsm] / total_blas_time * 100,
         blas_calls_[kTimeBlas_trsm - kTimeBlasStart]);
  printf("\tsyrk:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_syrk], times_[kTimeBlas_syrk] / total_blas_time * 100,
         blas_calls_[kTimeBlas_syrk - kTimeBlasStart]);
  printf("\tgemm:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times_[kTimeBlas_gemm], times_[kTimeBlas_gemm] / total_blas_time * 100,
         blas_calls_[kTimeBlas_gemm - kTimeBlasStart]);
  printf("----------------------------------------------------\n");
#endif
#endif
}
void printMemory(double mem) {
  if (mem < 1024)
    printf("%.1f B\n", mem);
  else if (mem < 1024 * 1024)
    printf("%.1f KB\n", mem / 1024);
  else if (mem < 1024 * 1024 * 1024)
    printf("%.1f MB\n", mem / 1024 / 1024);
  else
    printf("%.1f GB\n", mem / 1024 / 1024 / 1024);
}
void DataCollector::printSymbolic(bool verbose) const {
  printf("\nStatistic of Factor L\n");
  printf("size            : %.2e\n", (double)n_);
  printf("nnz             : %.2e\n", nz_);
  printf("fill-in         : %.2f\n", fillin_);
  printf("serial memory   : ");
  printMemory(serial_storage_);

  printf("dense  ops      : %.1e\n", dense_ops_);
  printf("sparse ops      : %.1e\n", sparse_ops_);
  printf("critical ops    : %.1e\n", critical_ops_);
  printf("max tree speedup: %.2f\n", dense_ops_ / critical_ops_);

  if (verbose) {
    printf("artificial nz   : %.1e (%.1f%%)\n", artificial_nz_,
           artificial_nz_ / nz_ * 100);
    printf("artificial ops  : %.1e (%.1f%%)\n", artificial_ops_,
           artificial_ops_ / dense_ops_ * 100);
    printf("largest front   : %5d\n", largest_front_);
    printf("largest sn      : %5d\n", largest_sn_);
    printf("supernodes      : %5d\n", sn_);
    printf("sn size <=   1  : %5d\n", sn_size_1_);
    printf("sn size <=  10  : %5d\n", sn_size_10_);
    printf("sn size <= 100  : %5d\n", sn_size_100_);
    printf("sn avg size     : %5.1f\n", (double)n_ / sn_);
  }
  printf("\n");
}

void DataCollector::printIter() const {
#ifdef DATA_COLLECTION
  printf(
      "\niter |    min D     max D     min L     max L  |"
      "    reg   swap    2x2     ws | "
      "  max_reg  solv   omega     back_err nw/cw      max_ws |\n");
  for (Int i = 0; i < iter_data_record_.size(); ++i) {
    const IterData& iter = iter_data_record_[i];
    printf(
        "%3d  |"
        " %9.1e %9.1e %9.1e %9.1e |"
        " %6d %6d %6d %6d |"
        " %9.1e %4d %9.1e %9.1e %9.1e",
        i, iter.minD, iter.maxD, iter.minL, iter.maxL, iter.n_reg_piv,
        iter.n_swap, iter.n_2x2, iter.n_wrong_sign, iter.max_reg,
        iter.num_solves, iter.omega, iter.nw_back_err, iter.cw_back_err);

    if (iter.max_wrong_sign > 0.0)
      printf(" %9.1e |\n", iter.max_wrong_sign);
    else
      printf(" %9s |\n", "-");
  }

  printf(
      "\niter |    norm1    maxdiag |"
      "    min T     max T  |"
      "     x_j * z_j / mu    small/large |"
      " corr   sigma af/co |\n");
  for (Int i = 0; i < iter_data_record_.size(); ++i) {
    const IterData& iter = iter_data_record_[i];
    printf(
        "%3d  | %9.1e %9.1e |"
        " %9.1e %9.1e |"
        " %9.1e %9.1e %6d %5d  |"
        " %3d %7.2f %6.2f |\n",
        i, iter.M_norm1, iter.M_maxdiag, iter.min_theta, iter.max_theta,
        iter.min_prod, iter.max_prod, iter.num_small_prod, iter.num_large_prod,
        iter.correctors, iter.sigma_aff, iter.sigma);
  }

  return;

  printf(
      "\niter |"
      "        xl       |       dxl       |"
      "        xu       |       dxu       |"
      "        zl       |       dzl       |"
      "        zu       |       dzu       |"
      "     alpha p     |     alpha d     |"
      "\n");

  for (Int i = 0; i < iter_data_record_.size(); ++i) {
    const IterData& iter = iter_data_record_[i];
    printf(
        "%3d  |"
        " %.1e %.1e | %.1e %.1e |"
        " %.1e %.1e | %.1e %.1e |"
        " %.1e %.1e | %.1e %.1e |"
        " %.1e %.1e | %.1e %.1e |"
        " %.1e %.1e | %.1e %.1e |"
        "\n",
        i, iter.min_xl, iter.max_xl, iter.min_dxl, iter.max_dxl, iter.min_xu,
        iter.max_xu, iter.min_dxu, iter.max_dxu, iter.min_zl, iter.max_zl,
        iter.min_dzl, iter.max_dzl, iter.min_zu, iter.max_zu, iter.min_dzu,
        iter.max_dzu, iter.p_limit_x, iter.p_limit_dx, iter.d_limit_z,
        iter.d_limit_dz);
  }
#endif
}

}  // namespace highspm