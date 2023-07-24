highs_root="$HOME/HiGHS"

#highs_build_dir="${highs_root}/test"

highs_build_dir="${highs_root}/build"

#VALGRIND="valgrind"

echo "FRED"

echo "$VALGRIND ./ipm $1 $2"
LD_LIBRARY_PATH=${highs_build_dir}/lib $VALGRIND ./ipm $1 $2

