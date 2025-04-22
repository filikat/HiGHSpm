#include "DataCollector.h"

namespace highspm {

// instance of DataCollector
DataCollector* DataCollector::ptr_ = nullptr;

DataCollector::DataCollector() {
#ifdef DATA_COLLECTION
  counter_data_.times.resize(kTimeSize);
  counter_data_.blas_calls.resize(kTimeBlasEnd - kTimeBlasStart + 1);
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

void FactorData::clear() {
  n = 0;
  nz = 0.0;
  sn = 0;
  fillin = 0.0;
  sparse_ops = 0.0;
  dense_ops = 0.0;
  critical_ops = 0.0;
  artificial_nz = 0.0;
  artificial_ops = 0.0;
  serial_storage = 0.0;
  largest_front = 0.0;
  largest_sn = 0.0;
  sn_size_1 = 0.0;
  sn_size_10 = 0.0;
  sn_size_100 = 0.0;
}
void CounterData::clear() {
  times.clear();
  blas_calls.clear();
  solves = 0;
}

void DataCollector::saveAndClear() {
  // Save current factorization and counter data in temporary storage and clear
  // main data. This is useful if analyse phase of both NE and AS is attempted,
  // otherwise the data would get corrupted.

  saved_factor_data_ = factor_data_;
  saved_counter_data_ = counter_data_;

  factor_data_.clear();
  counter_data_.clear();
}

void DataCollector::loadSaved() {
  // Sets the factorization data equal to the one stored in temporary storage.
  factor_data_ = std::move(saved_factor_data_);
  counter_data_ = std::move(saved_counter_data_);
}

void DataCollector::clearSaved() {
  saved_factor_data_.clear();
  saved_counter_data_.clear();
}

void DataCollector::sumTime(TimeItems i, double t) {
#ifdef DATA_COLLECTION
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
#ifdef DATA_COLLECTION
  std::lock_guard<std::mutex> lock(times_mutex_);
  ++counter_data_.solves;
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
  const FactorData& fd = factor_data_;

  printf("\nStatistic of Factor L\n");
  printf("size            : %.2e\n", (double)fd.n);
  printf("nnz             : %.2e\n", fd.nz);
  printf("fill-in         : %.2f\n", fd.fillin);
  printf("serial memory   : ");
  printMemory(fd.serial_storage);

  printf("dense  ops      : %.1e\n", fd.dense_ops);
  printf("sparse ops      : %.1e\n", fd.sparse_ops);
  printf("critical ops    : %.1e\n", fd.critical_ops);
  printf("max tree speedup: %.2f\n", fd.dense_ops / fd.critical_ops);

  if (verbose) {
    printf("artificial nz   : %.1e (%.1f%%)\n", fd.artificial_nz,
           fd.artificial_nz / fd.nz * 100);
    printf("artificial ops  : %.1e (%.1f%%)\n", fd.artificial_ops,
           fd.artificial_ops / fd.dense_ops * 100);
    printf("largest front   : %5d\n", fd.largest_front);
    printf("largest sn      : %5d\n", fd.largest_sn);
    printf("supernodes      : %5d\n", fd.sn);
    printf("sn size <=   1  : %5d\n", fd.sn_size_1);
    printf("sn size <=  10  : %5d\n", fd.sn_size_10);
    printf("sn size <= 100  : %5d\n", fd.sn_size_100);
    printf("sn avg size     : %5.1f\n", (double)fd.n / fd.sn);
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

#endif
}

}  // namespace highspm