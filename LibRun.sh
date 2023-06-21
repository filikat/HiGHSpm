highs_root="$HOME/HiGHS"

#highs_build_dir="${highs_root}/test"

highs_build_dir="${highs_root}/build"

#VALGRIND="valgrind"

LD_LIBRARY_PATH=${highs_build_dir}/lib $VALGRIND ./ipm $1
