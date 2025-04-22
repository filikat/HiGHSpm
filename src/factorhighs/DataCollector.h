#ifndef FACTORHIGHS_DATA_COLLECTOR_H
#define FACTORHIGHS_DATA_COLLECTOR_H

#include <atomic>
#include <mutex>
#include <vector>

#include "Timing.h"
#include "auxiliary/IntConfig.h"

namespace highspm {

struct IterData {
  // data of a given ipm iteration

#ifdef DATA_COLLECTION
  // factorization data
  double minD = std::numeric_limits<double>::max();
  double maxD = 0.0;
  double minL = std::numeric_limits<double>::max();
  double maxL = 0.0;
  double max_reg = 0.0;
  Int n_reg_piv = 0;
  Int n_swap = 0;
  Int n_2x2 = 0;
  Int n_wrong_sign = 0;
  double max_wrong_sign = 0.0;
#endif

  // generic ipm data
  double p_obj;
  double d_obj;
  double p_inf;
  double d_inf;
  double mu;
  double pd_gap;
  double p_alpha;
  double d_alpha;

  // advanced ipm data
  double min_theta;
  double max_theta;
  double min_prod;
  double max_prod;
  Int num_small_prod;
  Int num_large_prod;
  double sigma_aff;
  double sigma;
  Int correctors;
  double omega{};
  double nw_back_err{};
  double cw_back_err{};
  double M_norm1;
  double M_maxdiag;
  Int num_solves = 0;

  // extremes
  double min_xl = std::numeric_limits<double>::max();
  double min_xu = std::numeric_limits<double>::max();
  double min_zl = std::numeric_limits<double>::max();
  double min_zu = std::numeric_limits<double>::max();
  double max_xl = 0.0;
  double max_xu = 0.0;
  double max_zl = 0.0;
  double max_zu = 0.0;
  double min_dxl = std::numeric_limits<double>::max();
  double min_dxu = std::numeric_limits<double>::max();
  double min_dzl = std::numeric_limits<double>::max();
  double min_dzu = std::numeric_limits<double>::max();
  double max_dxl = 0.0;
  double max_dxu = 0.0;
  double max_dzl = 0.0;
  double max_dzu = 0.0;
  double p_limit_x;
  double p_limit_dx;
  double d_limit_z;
  double d_limit_dz;
};

struct FactorData {
  // Symbolic factorization statistics
  Int n{};
  double nz{};
  Int sn{};
  double fillin{};
  double dense_ops{};
  double sparse_ops{};
  double critical_ops{};
  double artificial_nz{};
  double artificial_ops{};
  double serial_storage{};
  Int largest_front{};
  Int largest_sn{};
  Int sn_size_1{};
  Int sn_size_10{};
  Int sn_size_100{};

  void clear();
};

struct CounterData {
  // Record of times and BLAS calls
  std::vector<double> times{};
  std::vector<Int> blas_calls{};
  Int solves{};

  void clear();
};

// DataCollector is a singleton object.
// Only one copy of it can exist and it does not have a public constructor or
// destructor.
// Use DataCollector::start() to allocate the DataCollector.
// Use DataCollector::destruct() to deallocate the DataCollector.
// Use DataCollector::get()->... to access any non-static member function.

class DataCollector {
  // Symbolic factorization data
  FactorData factor_data_;

  // Record of times and BLAS calls
  CounterData counter_data_;

  // record of data of ipm iterations
  std::vector<IterData> iter_data_record_{};

  // Mutexes for concurrent access
  std::mutex times_mutex_;
  std::mutex iter_data_mutex_;

  // Instance of DataCollector
  static DataCollector* ptr_;

  // Data that was set aside
  FactorData saved_factor_data_;
  CounterData saved_counter_data_;

  // Private ctor and dtor
  DataCollector();
  ~DataCollector() = default;

 public:
  // Access to the object
  static DataCollector* get();
  static void start();
  static void destruct();

  // Manage record of data of iterations
  void append();
  IterData& back();

  // Manage factorization data
  void saveAndClear();
  void loadSaved();
  void clearSaved();
  FactorData& factorData() { return factor_data_; }

  // Functions with lock, they can be accessed simultaneously
  void sumTime(TimeItems i, double t);
  void countSolves();
  void countRegPiv();
  void countSwap();
  void count2x2();
  void setWrongSign(double p);
  void setMaxReg(double new_reg);
  void setExtremeEntries(double minD, double maxD, double minoffD,
                         double maxoffD);

  // Const functions
  void printTimes() const;
  void printSymbolic(bool verbose = false) const;
  void printIter() const;
};

}  // namespace highspm

#endif