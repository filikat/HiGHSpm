
YOUR_CPP_COMPILER="clang++"

YOUR_PATH_TO_INCLUDE="/usr/local/include/highs/"

YOUR_PATH_TO_LIBRARY="/usr/local/lib/"

AUX_LIBRARY_FOR_MAC="/Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/usr/lib/"

FLAGS="-std=c++11 -O3 -o ipm"

$YOUR_CPP_COMPILER \
        $FLAGS \
        -I$YOUR_PATH_TO_INCLUDE \
        -L$YOUR_PATH_TO_LIBRARY -lhighs \
        -L$AUX_LIBRARY_FOR_MAC -lz \
	 mainIPM.cpp \
        IPM_caller.cpp \
        IPM_model.cpp \
        SparseMatrix.cpp \
        NormalEquations.cpp \
        ConjugateGradient.cpp \
        VectorOperations.cpp \
        IPM_aux.cpp
