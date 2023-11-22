#include "SelectMethod.h"

int selectMethod(const HighsSparseMatrix& A, Metis_caller& Metis_data,
                 int& option_nla, int& option_metis) {
  // select NE or AS

  int system_size{};

  if (option_nla == kOptionNlaChoose) option_nla = choose_NE_AS();

  if (option_nla == kOptionNlaAugmented) {
    system_size = A.num_row_ + A.num_col_;
  } else {
    system_size = A.num_row_;
  }

  if (option_metis == kOptionMetisChoose || option_metis == kOptionMetisOn) {
    // try partitioning with 2,4,8 parts to see which one has the smaller Schur
    // complement. If the problem has a specific structure, it may be better to
    // use more parts. In general, 2 will be the best result.

    std::vector<int> parts_to_try = {2, 4, 8};

    int bestSchur{std::numeric_limits<int>::max()};

    for (int parts : parts_to_try) {
      Metis_caller temp_metis(A, option_nla, parts);
      temp_metis.setDebug(false);
      temp_metis.getPartition();
      int metis_status = temp_metis.getPermutation();
      if (metis_status) continue;
      int tempSchur = temp_metis.sizeSchur();
      if (tempSchur < bestSchur) {
        bestSchur = tempSchur;
        Metis_data = std::move(temp_metis);
      };
    }

    if (option_metis == kOptionMetisChoose) {
      if ((double)bestSchur / system_size < kMetisSchurRelativeThreshold &&
          bestSchur < kMetisSchurAbsThreshold) {
        option_metis = kOptionMetisOn;
      } else {
        option_metis = kOptionMetisOff;
      }
    }

    if (bestSchur == std::numeric_limits<int>::max() &&
        option_metis == kOptionMetisOn) {
      // if we couldn't find any acceptable permutation (metis_status was always
      // 1), give an error
      std::cerr << "Error while partitioning the matrix\n";
      return 1;
    }
  }

  return 0;
}

// choice of NE or AS not implemented.
// it may depend on the detection of dense columns
int choose_NE_AS() { return kOptionNlaNewton; }