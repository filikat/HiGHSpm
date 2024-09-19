#include "Ipm_aux.h"

#include <fstream>
#include <iostream>

#include "VectorOperations.h"

void scaling2theta(const std::vector<double>& scaling,
                   std::vector<double>& theta) {
  const int dim = scaling.size();
  theta.resize(dim);
  for (int i = 0; i < dim; i++) theta[i] = 1 / scaling[i];
}

Residuals::Residuals(int m, int n)
    : res1(m, 0.0),
      res2(n, 0.0),
      res3(n, 0.0),
      res4(n, 0.0),
      res5(n, 0.0),
      res6(n, 0.0) {}

void Residuals::print(int iter) const {
  std::ofstream out_file;

  char str[50];

  snprintf(str, 50, "res1_%d.txt", iter);
  out_file.open(str);
  for (double d : res1) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res2_%d.txt", iter);
  out_file.open(str);
  for (double d : res2) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res3_%d.txt", iter);
  out_file.open(str);
  for (double d : res3) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res4_%d.txt", iter);
  out_file.open(str);
  for (double d : res4) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res5_%d.txt", iter);
  out_file.open(str);
  for (double d : res5) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "res6_%d.txt", iter);
  out_file.open(str);
  for (double d : res6) {
    out_file << d << '\n';
  }
  out_file.close();
}

bool Residuals::isNaN() const {
  if (isNanVector(res1) || isNanVector(res2) || isNanVector(res3) ||
      isNanVector(res4) || isNanVector(res5) || isNanVector(res6))
    return true;
  return false;
}

bool Residuals::isInf() const {
  if (isInfVector(res1) || isInfVector(res2) || isInfVector(res3) ||
      isInfVector(res4) || isInfVector(res5) || isInfVector(res6))
    return true;
  return false;
}

Iterate::Iterate(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 1.0), xu(n, 1.0), zl(n, 1.0), zu(n, 1.0) {}

bool Iterate::isNaN() const {
  if (isNanVector(x) || isNanVector(xl) || isNanVector(xu) || isNanVector(y) ||
      isNanVector(zl) || isNanVector(zu))
    return true;
  return false;
}

bool Iterate::isInf() const {
  if (isInfVector(x) || isInfVector(xl) || isInfVector(xu) || isInfVector(y) ||
      isInfVector(zl) || isInfVector(zu))
    return true;
  return false;
}

void Iterate::print(int iter) const {
  std::ofstream out_file;

  char str[50];

  snprintf(str, 50, "it_x_%d.txt", iter);
  out_file.open(str);
  for (double d : x) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_y_%d.txt", iter);
  out_file.open(str);
  for (double d : y) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_xl_%d.txt", iter);
  out_file.open(str);
  for (double d : xl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_xu_%d.txt", iter);
  out_file.open(str);
  for (double d : xu) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_zl_%d.txt", iter);
  out_file.open(str);
  for (double d : zl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "it_zu_%d.txt", iter);
  out_file.open(str);
  for (double d : zu) {
    out_file << d << '\n';
  }
  out_file.close();
}

NewtonDir::NewtonDir(int m, int n)
    : x(n, 0.0), y(m, 0.0), xl(n, 0.0), xu(n, 0.0), zl(n, 0.0), zu(n, 0.0) {}

bool NewtonDir::isNaN() const {
  if (isNanVector(x) || isNanVector(xl) || isNanVector(xu) || isNanVector(y) ||
      isNanVector(zl) || isNanVector(zu))
    return true;
  return false;
}

bool NewtonDir::isInf() const {
  if (isInfVector(x) || isInfVector(xl) || isInfVector(xu) || isInfVector(y) ||
      isInfVector(zl) || isInfVector(zu))
    return true;
  return false;
}

void NewtonDir::print(int iter) const {
  std::ofstream out_file;

  char str[50];

  snprintf(str, 50, "D_x_%d.txt", iter);
  out_file.open(str);
  for (double d : x) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_y_%d.txt", iter);
  out_file.open(str);
  for (double d : y) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_xl_%d.txt", iter);
  out_file.open(str);
  for (double d : xl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_xu_%d.txt", iter);
  out_file.open(str);
  for (double d : xu) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_zl_%d.txt", iter);
  out_file.open(str);
  for (double d : zl) {
    out_file << d << '\n';
  }
  out_file.close();

  snprintf(str, 50, "D_zu_%d.txt", iter);
  out_file.open(str);
  for (double d : zu) {
    out_file << d << '\n';
  }
  out_file.close();
}

int computeAThetaAT(const HighsSparseMatrix& matrix,
                    const std::vector<double>& theta, HighsSparseMatrix& AAT,
                    const int max_num_nz) {
  // Create a row-wise copy of the matrix
  HighsSparseMatrix AT = matrix;
  AT.ensureRowwise();

  int AAT_dim = matrix.num_row_;
  AAT.num_col_ = AAT_dim;
  AAT.num_row_ = AAT_dim;
  AAT.start_.resize(AAT_dim + 1, 0);

  std::vector<std::tuple<int, int, double>> non_zero_values;

  // First pass to calculate the number of non-zero elements in each column
  //
  int AAT_num_nz = 0;
  std::vector<double> AAT_col_value(AAT_dim, 0);
  std::vector<int> AAT_col_index(AAT_dim);
  std::vector<bool> AAT_col_in_index(AAT_dim, false);
  for (int iRow = 0; iRow < AAT_dim; iRow++) {
    // Go along the row of A, and then down the columns corresponding
    // to its nonzeros
    int num_col_el = 0;
    for (int iRowEl = AT.start_[iRow]; iRowEl < AT.start_[iRow + 1]; iRowEl++) {
      int iCol = AT.index_[iRowEl];
      const double theta_value = !theta.empty() ? theta[iCol] : 1;
      if (!theta_value) continue;
      const double row_value = theta_value * AT.value_[iRowEl];
      for (int iColEl = matrix.start_[iCol]; iColEl < matrix.start_[iCol + 1];
           iColEl++) {
        int iRow1 = matrix.index_[iColEl];
        if (iRow1 < iRow) continue;
        double term = row_value * matrix.value_[iColEl];
        if (!AAT_col_in_index[iRow1]) {
          // This entry is not yet in the list of possible nonzeros
          AAT_col_in_index[iRow1] = true;
          AAT_col_index[num_col_el++] = iRow1;
          AAT_col_value[iRow1] = term;
        } else {
          // This entry is in the list of possible nonzeros
          AAT_col_value[iRow1] += term;
        }
      }
    }
    for (int iEl = 0; iEl < num_col_el; iEl++) {
      int iCol = AAT_col_index[iEl];
      assert(iCol >= iRow);
      const double value = AAT_col_value[iCol];

      non_zero_values.emplace_back(iRow, iCol, value);
      const int num_new_nz = iRow != iCol ? 2 : 1;
      if (AAT_num_nz + num_new_nz >= max_num_nz)
        return kDecomposerStatusErrorOom;
      AAT.start_[iRow + 1]++;
      if (iRow != iCol) AAT.start_[iCol + 1]++;
      AAT_num_nz += num_new_nz;

      AAT_col_value[iCol] =
          0;  // Not strictly necessary, but simplifies debugging
      AAT_col_in_index[iCol] = false;
    }
  }

  // Prefix sum to get the correct column pointers
  for (int i = 0; i < AAT_dim; ++i) AAT.start_[i + 1] += AAT.start_[i];

  AAT.index_.resize(AAT.start_.back());
  AAT.value_.resize(AAT.start_.back());
  AAT.p_end_ = AAT.start_;
  AAT.p_end_.back() = AAT.index_.size();

  std::vector<int> current_positions = AAT.start_;

  // Second pass to actually fill in the indices and values
  for (const auto& val : non_zero_values) {
    int i = std::get<0>(val);
    int j = std::get<1>(val);
    double dot = std::get<2>(val);

    AAT.index_[current_positions[i]] = j;
    AAT.value_[current_positions[i]] = dot;
    current_positions[i]++;
    AAT.p_end_[i] = current_positions[i];

    if (i != j) {
      AAT.index_[current_positions[j]] = i;
      AAT.value_[current_positions[j]] = dot;
      current_positions[j]++;
      AAT.p_end_[j] = current_positions[j];
    }
  }
  AAT.p_end_.clear();
  return kDecomposerStatusOk;
}

double getWallTime() {
  using namespace std::chrono;
  using wall_clock = std::chrono::high_resolution_clock;
  return duration_cast<duration<double>>(wall_clock::now().time_since_epoch())
      .count();
}

void debug_print(std::string& filestr, const std::vector<int>& data) {
  char filename[100];
  snprintf(filename, 100, "../FactorHiGHS/matlab/%s", filestr.c_str());

  FILE* out_file = fopen(filename, "w+");

  for (int i : data) {
    fprintf(out_file, "%d\n", i);
  }

  fclose(out_file);
}

void debug_print(std::string& filestr, const std::vector<double>& data) {
  char filename[100];
  snprintf(filename, 100, "../FactorHiGHS/matlab/%s", filestr.c_str());

  FILE* out_file = fopen(filename, "w+");

  for (double d : data) {
    fprintf(out_file, "%.10e\n", d);
  }

  fclose(out_file);
}