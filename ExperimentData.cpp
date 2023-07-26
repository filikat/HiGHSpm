#include "Direct.h"
#include "ExperimentData.h"
#include <iomanip>
double getWallTime() {
    using namespace std::chrono;
    using wall_clock = std::chrono::high_resolution_clock;
    return duration_cast<duration<double> >(wall_clock::now().time_since_epoch()).count();
}

std::ostream& operator<<(std::ostream& os, const ExperimentData& data) {
  const int text_width = 20;
  const int num_width = 12;
  const int pct_width = 8;
  double float_dim = double(data.model_size);
  const double aat_density = data.model_size ? 1e2 * double(data.nnz_AAT) / (float_dim * float_dim) : -1;
  const double l_density = data.model_size ? 1e2 * double(data.nnz_L) / (float_dim * double(data.model_size+1) * 0.5) : -1;
  const double sum_time = data.analysis_time + data.factorization_time + data.solve_time;
  const double accounted_time = data.time_taken > 0 ? 1e2 * sum_time / data.time_taken : -1;
  os 
    << std::left << std::setw(text_width) << "decomposer: "
    << std::right << std::setw(num_width) << data.decomposer << "\n" 
    << std::left << std::setw(text_width) << "model name:" 
    << std::right << std::setw(num_width) << data.model_name << "\n" 
    << std::left << std::setw(text_width) << "model size: " 
    << std::right << std::setw(num_width) << data.model_size << "\n" 
    << std::left << std::setw(text_width) << "AAT nnz: " 
    << std::right << std::setw(num_width) << data.nnz_AAT << " ("
    << std::right << std::setw(pct_width) << aat_density << "%)\n" 
    << std::left << std::setw(text_width) << "L nnz: " 
    << std::right << std::setw(num_width) << data.nnz_L << " ("
    << std::right << std::setw(pct_width) << l_density << "%)\n" 
    << std::left << std::setw(text_width) << "solution error: " 
    << std::right << std::setw(num_width) << data.solution_error << "\n" 
    << std::left << std::setw(text_width) << "residual error: " 
    << std::right << std::setw(num_width) << data.residual_error << "\n" 
    << std::left << std::setw(text_width) << "fill-in: " 
    << std::right << std::setw(num_width) << data.fill_in_factor << "\n" 
    << std::left << std::setw(text_width) << "analyse time: " 
    << std::right << std::setw(num_width) << std::setprecision(6) << data.analysis_time << "\n" 
    << std::left << std::setw(text_width) << "factorization time: " 
    << std::right << std::setw(num_width) << std::setprecision(6) << data.factorization_time << "\n" 
    << std::left << std::setw(text_width) << "solve time: " 
    << std::right << std::setw(num_width) << std::setprecision(6) << data.solve_time << "\n" 
    << std::left << std::setw(text_width) << "sum time: " 
    << std::right << std::setw(num_width) << std::setprecision(6) << sum_time << " ("
    << std::right << std::setw(pct_width) << accounted_time << "%)\n" 
    << std::left << std::setw(text_width) << "time taken: " 
    << std::right << std::setw(num_width) << std::setprecision(6) << data.time_taken << "\n";
  return os;
}

void writeDataToCSV(const std::vector<ExperimentData>& data, const std::string& filename)
{
    std::ofstream outputFile;
    outputFile.open(filename);
    
    // Write header
    outputFile << "Decomposer,Model Name,Model Size,NNZ AAT,NNZ L,Solution Error,Residual Error,Fill in Factor,Time Taken, Analyse time, factorization time, solve time\n";
    
    // Write data
    for(const auto& experimentData : data)
    {
        outputFile << experimentData.decomposer << ",";
        outputFile << experimentData.model_name << ",";
        outputFile << experimentData.model_size << ",";
        outputFile << experimentData.nnz_AAT << ",";
        outputFile << experimentData.nnz_L << ",";
        outputFile << experimentData.solution_error << ",";
        outputFile << experimentData.residual_error << ",";
        outputFile << experimentData.fill_in_factor << ",";
        outputFile << experimentData.time_taken << ",";
        outputFile << experimentData.analysis_time << ",";
        outputFile << experimentData.factorization_time << ",";
        outputFile << experimentData.solve_time << "\n";
    }
    
    outputFile.close();
}

double residualError(const HighsSparseMatrix& A,
		     const std::vector<double>& b,
		     const std::vector<double>& x){
  std::vector<double> residual = b;
  A.alphaProductPlusY(-1, x, residual);
  double residual_error = 0;
  for (int ix = 0; ix < b.size(); ix++)
    residual_error = std::max(std::fabs(residual[ix]), residual_error);
  return residual_error;
}

double residualErrorAThetaAT(const HighsSparseMatrix& A,
			     const std::vector<double>& theta,
			     const std::vector<double>& b,
			     const std::vector<double>& x){
  std::vector<double> AThetaATx;
  productAThetaAT(A, theta, x, AThetaATx);
  double residual_error = 0;
  for (int ix = 0; ix < b.size(); ix++)
    residual_error = std::max(std::fabs(AThetaATx[ix]-b[ix]), residual_error);
  return residual_error;
}

double fillIn_LL(int nnz_AAT, int nnz_L, int MatrixSize){
    return (double) (2*nnz_L - MatrixSize)/ nnz_AAT;
}

double fillIn_LDL(int nnz_AAT, int nnz_L, int MatrixSize){
    return (double) (2*nnz_L + MatrixSize)/ nnz_AAT;
}
