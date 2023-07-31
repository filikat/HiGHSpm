
YOUR_CPP_COMPILER="clang++"

#AUX_LIBRARY_FOR_MAC="/Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/usr/lib/"  # may be necessary to link libz (-lz)

HIGHS_PATH="$HOME/Documents/HiGHS/"

FLAGS="-std=c++11 -g3 -o ipm"

$YOUR_CPP_COMPILER \
        $FLAGS \
        -I$HIGHS_PATH/build/ \
        -I$HIGHS_PATH/src/ \
        -L$HIGHS_PATH/build/lib/ -lhighs \
        mainIPM.cpp \
        IPM_caller.cpp \
        IPM_model.cpp \
        NormalEquations.cpp \
        ConjugateGradient.cpp \
        VectorOperations.cpp \
        IPM_aux.cpp \
        Direct.cpp
