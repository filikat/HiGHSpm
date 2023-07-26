#include "Direct.h"
#include "ExperimentData.h"

double getWallTime() {
    using namespace std::chrono;
    using wall_clock = std::chrono::high_resolution_clock;
    return duration_cast<duration<double> >(wall_clock::now().time_since_epoch()).count();
}

std::ostream& operator<<(std::ostream& os, const ExperimentData& data) {
    os << "decomposer: " << data.decomposer <<
      "\n model name:" << data.model_name <<
      "\n model size: " << data.model_size <<
      "\n nnz in AAT: " << data.nnz_AAT <<
      "\n nnz in L: " << data.nnz_L << 
      "\n solution error: " << data.solution_error << 
      "\n residual error: " << data.residual_error << 
      "\n fill-in: " << data.fill_in_factor << 
      "\n analyse time: " << data.analysis_time<< 
      "\n factorization time: " << data.factorization_time<< 
      "\n solve time: " << data.solve_time << 
      "\n time taken: " << data.time_taken << 
      "\n";
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
