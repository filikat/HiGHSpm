#ifndef EXPERIMENTDATA_H_
#define EXPERIMENTDATA_H_
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Highs.h"

class ExperimentData {
public:
  std::string decomposer;
  std::string model_name;
  int model_num_col;
  int model_num_row;
  bool newton_solve;
  int system_size;
  //bool is_A_positive_definite;
  int system_nnz;
  int nnz_L;
  double solution_error;
  double residual_error;
  double fill_in_factor;

  //time
  double time_taken;
  double form_time;
  double analysis_time;
  double factorization_time;
  double solve_time;


  void reset(){
    decomposer = "na";
    model_num_col = -1;
    model_num_row = -1;
    newton_solve = false;
    system_size = -1;
    system_nnz = -1;
    nnz_L = -1;
    solution_error = -1;
    residual_error = -1;
    fill_in_factor = -1;
    time_taken = -1;
    analysis_time = -1;
    factorization_time = -1;
    solve_time = -1;
    
  }
  void fillIn_LL();
  void fillIn_LDL();

};

double getWallTime();

std::ostream& operator<<(std::ostream& os, const ExperimentData& data);
void writeDataToCSV(const std::vector<ExperimentData>& data, const std::string& filename);
double residualError(const HighsSparseMatrix& AAT,
		     const std::vector<double>& b,
		     const std::vector<double>& x);
double residualErrorAThetaAT(const HighsSparseMatrix& A,
			     const std::vector<double>& theta,
			     const std::vector<double>& b,
			     const std::vector<double>& x);
#endif
