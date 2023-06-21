CPP_COMPILER="clang++"
highs_root="$HOME/HiGHS"

highs_build="$highs_root/build"

highs_include="$highs_root/src"

highs_library="$highs_build/lib"

#AUX_LIBRARY_FOR_MAC="/Library/Developer/CommandLineTools/SDKs/MacOSX13.1.sdk/usr/lib/"

FLAGS="-std=c++11 -g3 -o ipm"

$CPP_COMPILER \
        $FLAGS \
        -I$highs_build \
        -I$highs_include \
        -L$highs_library -lhighs \
        mainIPM.cpp \
        IPM_caller.cpp \
        IPM_model.cpp \
        SparseMatrix.cpp \
        NormalEquations.cpp \
        ConjugateGradient.cpp \
        VectorOperations.cpp \
        IPM_aux.cpp

	#        -L$AUX_LIBRARY_FOR_MAC -lz \
