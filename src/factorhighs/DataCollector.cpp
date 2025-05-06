#include "DataCollector.h"

namespace highspm {

// instance of DataCollector
DataCollector* DataCollector::ptr_ = nullptr;

DataCollector::DataCollector() {
#ifdef DEBUG
  counter_data_.times.resize(kTimeSize);
  counter_data_.blas_calls.resize(kTimeBlasEnd - kTimeBlasStart + 1);
#endif
}

DataCollector* DataCollector::get() { return ptr_; }

void DataCollector::initialise() {
#ifdef DEBUG
  if (!ptr_) ptr_ = new DataCollector();
#endif
}
void DataCollector::terminate() {
#ifdef DEBUG
  delete ptr_;
  ptr_ = nullptr;
#endif
}
void DataCollector::append() {
#ifdef DEBUG
  // add an empty IterData object to the record
  iter_data_record_.push_back(IterData());
#endif
}
IterData& DataCollector::back() {
  // access most recent record of data
  return iter_data_record_.back();
}

void DataCollector::sumTime(TimeItems i, double t) {
#ifdef DEBUG
  // Keep track of times and blas calls.
  std::lock_guard<std::mutex> lock(times_mutex_);
  counter_data_.times[i] += t;
#ifdef BLAS_TIMING
  if (i >= kTimeBlasStart && i <= kTimeBlasEnd)
    ++counter_data_.blas_calls[i - kTimeBlasStart];
#endif
#endif
}
void DataCollector::countSolves() {
#ifdef DEBUG
  std::lock_guard<std::mutex> lock(times_mutex_);
  ++counter_data_.solves;
  ++back().num_solves;
#endif
}
void DataCollector::setExtremeEntries(double minD, double maxD, double minoffD,
                                      double maxoffD) {
#ifdef DEBUG
  // Store max and min entries of D and L.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  back().minD = std::min(back().minD, minD);
  back().maxD = std::max(back().maxD, maxD);
  back().minL = std::min(back().minL, minoffD);
  back().maxL = std::max(back().maxL, maxoffD);
#endif
}
void DataCollector::countRegPiv() {
#ifdef DEBUG
  // Increase the number of dynamically regularised pivots.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_reg_piv;
#endif
}
void DataCollector::countSwap() {
#ifdef DEBUG
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_swap;
#endif
}
void DataCollector::count2x2() {
#ifdef DEBUG
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_2x2;
#endif
}
void DataCollector::setWrongSign(double p) {
#ifdef DEBUG
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_wrong_sign;
  back().max_wrong_sign = std::max(back().max_wrong_sign, std::abs(p));
#endif
}
void DataCollector::setMaxReg(double new_reg) {
#ifdef DEBUG
  // Keep track of maximum regularisation used.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  back().max_reg = std::max(back().max_reg, new_reg);
#endif
}

void DataCollector::setOmega(double omega) {
#ifdef DEBUG
  back().omega = std::max(back().omega, omega);
#endif
}
void DataCollector::setNorms(double norm1, double maxdiag) {
#ifdef DEBUG
  back().M_norm1 = norm1;
  back().M_maxdiag = maxdiag;
#endif
}
void DataCollector::setSigma(double sigma, bool affinescaling) {
#ifdef DEBUG
  if (affinescaling)
    back().sigma_aff = sigma;
  else
    back().sigma = sigma;
#endif
}
void DataCollector::setCorrectors(Int correctors) {
#ifdef DEBUG
  back().correctors = correctors;
#endif
}
void DataCollector::setBackError(double nw, double cw, Int large_components) {
#ifdef DEBUG
  back().nw_back_err = std::max(back().nw_back_err, nw);
  back().cw_back_err = std::max(back().cw_back_err, cw);
  back().large_components_cw =
      std::max(back().large_components_cw, large_components);
#endif
}
void DataCollector::setExtremeTheta(const std::vector<double>& scaling) {
#ifdef DEBUG
  back().min_theta = std::numeric_limits<double>::infinity();
  back().max_theta = 0.0;
  for (double d : scaling) {
    if (d != 0.0) {
      back().min_theta = std::min(back().min_theta, 1.0 / d);
      back().max_theta = std::max(back().max_theta, 1.0 / d);
    }
  }
#endif
}
void DataCollector::setProducts(double min_prod, double max_prod, Int num_small,
                                Int num_large) {
#ifdef DEBUG
  back().min_prod = min_prod;
  back().max_prod = max_prod;
  back().num_small_prod = num_small;
  back().num_large_prod = num_large;
#endif
}

void DataCollector::printTimes() const {
#ifdef COARSE_TIMING

  const std::vector<double>& times = counter_data_.times;

  printf("----------------------------------------------------\n");
  printf("Analyse time            \t%8.4f\n", times[kTimeAnalyse]);

#ifdef FINE_TIMING
  printf("\tMetis:                  %8.4f (%4.1f%%)\n",
         times[kTimeAnalyseMetis],
         times[kTimeAnalyseMetis] / times[kTimeAnalyse] * 100);
  printf("\tTree:                   %8.4f (%4.1f%%)\n", times[kTimeAnalyseTree],
         times[kTimeAnalyseTree] / times[kTimeAnalyse] * 100);
  printf("\tCounts:                 %8.4f (%4.1f%%)\n",
         times[kTimeAnalyseCount],
         times[kTimeAnalyseCount] / times[kTimeAnalyse] * 100);
  printf("\tSupernodes:             %8.4f (%4.1f%%)\n", times[kTimeAnalyseSn],
         times[kTimeAnalyseSn] / times[kTimeAnalyse] * 100);
  printf("\tReorder:                %8.4f (%4.1f%%)\n",
         times[kTimeAnalyseReorder],
         times[kTimeAnalyseReorder] / times[kTimeAnalyse] * 100);
  printf("\tSn sparsity pattern:    %8.4f (%4.1f%%)\n",
         times[kTimeAnalysePattern],
         times[kTimeAnalysePattern] / times[kTimeAnalyse] * 100);
  printf("\tRelative indices:       %8.4f (%4.1f%%)\n",
         times[kTimeAnalyseRelInd],
         times[kTimeAnalyseRelInd] / times[kTimeAnalyse] * 100);
#endif

  printf("----------------------------------------------------\n");
  printf("Factorise time          \t%8.4f\n", times[kTimeFactorise]);

#ifdef FINE_TIMING
  printf("\tPrepare:                %8.4f (%4.1f%%)\n",
         times[kTimeFactorisePrepare],
         times[kTimeFactorisePrepare] / times[kTimeFactorise] * 100);
  printf("\tAssembly original:      %8.4f (%4.1f%%)\n",
         times[kTimeFactoriseAssembleOriginal],
         times[kTimeFactoriseAssembleOriginal] / times[kTimeFactorise] * 100);
  printf("\tAssemble children in F: %8.4f (%4.1f%%)\n",
         times[kTimeFactoriseAssembleChildrenFrontal],
         times[kTimeFactoriseAssembleChildrenFrontal] / times[kTimeFactorise] *
             100);
  printf("\tAssemble children in C: %8.4f (%4.1f%%)\n",
         times[kTimeFactoriseAssembleChildrenClique],
         times[kTimeFactoriseAssembleChildrenClique] / times[kTimeFactorise] *
             100);
  printf("\tDense factorisation:    %8.4f (%4.1f%%)\n",
         times[kTimeFactoriseDenseFact],
         times[kTimeFactoriseDenseFact] / times[kTimeFactorise] * 100);
  printf("\t\tmain:           %8.4f\n", times[kTimeDenseFact_main]);
  printf("\t\tSchur:          %8.4f\n", times[kTimeDenseFact_schur]);
  printf("\t\tkernel:         %8.4f\n", times[kTimeDenseFact_kernel]);
  printf("\t\tconvert:        %8.4f\n", times[kTimeDenseFact_convert]);
  printf("\t\tpivoting:       %8.4f\n", times[kTimeDenseFact_pivoting]);
  printf("\tTerminate:              %8.4f (%4.1f%%)\n",
         times[kTimeFactoriseTerminate],
         times[kTimeFactoriseTerminate] / times[kTimeFactorise] * 100);
#endif

  printf("----------------------------------------------------\n");
  printf("Solve time              \t%8.4f (%d calls)\n", times[kTimeSolve],
         counter_data_.solves);

#ifdef FINE_TIMING
  printf("\tPrepare:                %8.4f (%4.1f%%)\n",
         times[kTimeSolvePrepare],
         times[kTimeSolvePrepare] / times[kTimeSolve] * 100);
  printf("\tSolve:                  %8.4f (%4.1f%%)\n", times[kTimeSolveSolve],
         times[kTimeSolveSolve] / times[kTimeSolve] * 100);
  printf("\tResidual:               %8.4f (%4.1f%%)\n",
         times[kTimeSolveResidual],
         times[kTimeSolveResidual] / times[kTimeSolve] * 100);
  printf("\tOmega:                  %8.4f (%4.1f%%)\n", times[kTimeSolveOmega],
         times[kTimeSolveOmega] / times[kTimeSolve] * 100);
#endif
  printf("----------------------------------------------------\n");

#ifdef BLAS_TIMING

  const std::vector<Int>& blas_calls = counter_data_.blas_calls;

  double total_blas_time =
      times[kTimeBlas_copy] + times[kTimeBlas_axpy] + times[kTimeBlas_scal] +
      times[kTimeBlas_swap] + times[kTimeBlas_gemv] + times[kTimeBlas_trsv] +
      times[kTimeBlas_tpsv] + times[kTimeBlas_ger] + times[kTimeBlas_trsm] +
      times[kTimeBlas_syrk] + times[kTimeBlas_gemm];

  printf("BLAS time               \t%8.4f\n", total_blas_time);
  printf("\tcopy:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_copy], times[kTimeBlas_copy] / total_blas_time * 100,
         blas_calls[kTimeBlas_copy - kTimeBlasStart]);
  printf("\taxpy:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_axpy], times[kTimeBlas_axpy] / total_blas_time * 100,
         blas_calls[kTimeBlas_axpy - kTimeBlasStart]);
  printf("\tscal:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_scal], times[kTimeBlas_scal] / total_blas_time * 100,
         blas_calls[kTimeBlas_scal - kTimeBlasStart]);
  printf("\tswap:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_swap], times[kTimeBlas_swap] / total_blas_time * 100,
         blas_calls[kTimeBlas_swap - kTimeBlasStart]);
  printf("\tgemv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_gemv], times[kTimeBlas_gemv] / total_blas_time * 100,
         blas_calls[kTimeBlas_gemv - kTimeBlasStart]);
  printf("\ttrsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_trsv], times[kTimeBlas_trsv] / total_blas_time * 100,
         blas_calls[kTimeBlas_trsv - kTimeBlasStart]);
  printf("\ttpsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_tpsv], times[kTimeBlas_tpsv] / total_blas_time * 100,
         blas_calls[kTimeBlas_tpsv - kTimeBlasStart]);
  printf("\tger:            \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_ger], times[kTimeBlas_ger] / total_blas_time * 100,
         blas_calls[kTimeBlas_ger - kTimeBlasStart]);
  printf("\ttrsm:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_trsm], times[kTimeBlas_trsm] / total_blas_time * 100,
         blas_calls[kTimeBlas_trsm - kTimeBlasStart]);
  printf("\tsyrk:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_syrk], times[kTimeBlas_syrk] / total_blas_time * 100,
         blas_calls[kTimeBlas_syrk - kTimeBlasStart]);
  printf("\tgemm:           \t%8.4f (%4.1f%%) in %10d calls\n",
         times[kTimeBlas_gemm], times[kTimeBlas_gemm] / total_blas_time * 100,
         blas_calls[kTimeBlas_gemm - kTimeBlasStart]);
  printf("----------------------------------------------------\n");
#endif
#endif
}

void DataCollector::printIter() const {
#ifdef DEBUG
  printf(
      "\niter |    min D     max D     min L     max L  |"
      "    reg   swap    2x2     ws | "
      "  max_reg  solv   omega     nw_be     cw_be     cw_large  max_ws |\n");
  for (Int i = 0; i < iter_data_record_.size(); ++i) {
    const IterData& iter = iter_data_record_[i];
    printf(
        "%3d  |"
        " %9.1e %9.1e %9.1e %9.1e |"
        " %6d %6d %6d %6d |"
        " %9.1e %4d %9.1e %9.1e %9.1e %9d",
        i, iter.minD, iter.maxD, iter.minL, iter.maxL, iter.n_reg_piv,
        iter.n_swap, iter.n_2x2, iter.n_wrong_sign, iter.max_reg,
        iter.num_solves, iter.omega, iter.nw_back_err, iter.cw_back_err,
        iter.large_components_cw);

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

#endif
}

}  // namespace highspm