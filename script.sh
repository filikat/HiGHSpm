
C_COMPILER="clang"
CPP_COMPILER="clang++"

HIGHS_PATH="${HOME}/Documents/HiGHS"
CHOLMOD_PATH="${HOME}/Documents/SuiteSparse"
QDLDL_PATH="${HOME}/Documents/qdldl"
METIS_PATH="${HOME}/Documents/METIS"
LOCAL_PATH="${HOME}/local"

HIGHS_INCLUDE_AND_LIB="-I$HIGHS_PATH/build -I$HIGHS_PATH/src/ -L$HIGHS_PATH/build/lib"
CHOLMOD_INCLUDE_AND_LIB="-I$CHOLMOD_PATH/include -L$CHOLMOD_PATH/lib"
QDLDL_INCLUDE_AND_LIB="-I$QDLDL_PATH/include -L$QDLDL_PATH/build/out"
METIS_INCLUDE_AND_LIB="-I$METIS_PATH/include -L$METIS_PATH/build/libmetis"
LOCAL_INCLUDE_AND_LIB="-I$LOCAL_PATH/include -L$LOCAL_PATH/lib"

LIBS="-lhighs -lcholmod -lqdldl -lmetis -lhsl_ma86 -lfakemetis -lGKlib"

FLAGS="-std=c++11 -o ipm -O3"
OPENMP_FLAGS="-Xclang -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp"

$C_COMPILER -c hsl_ma86_wrapper.c -o hsl_ma86_wrapper.o -I$LOCAL_PATH/include

$CPP_COMPILER \
        $FLAGS \
        $OPENMP_FLAGS \
        $HIGHS_INCLUDE_AND_LIB \
        $CHOLMOD_INCLUDE_AND_LIB \
        $QDLDL_INCLUDE_AND_LIB \
        $METIS_INCLUDE_AND_LIB \
        $LOCAL_INCLUDE_AND_LIB \
        $LIBS \
        mainIPM.cpp \
        IPM_caller.cpp \
        IPM_model.cpp \
        NormalEquations.cpp \
        ConjugateGradient.cpp \
        VectorOperations.cpp \
        IPM_aux.cpp \
        Direct.cpp \
        ExperimentData.cpp \
        Metis_caller.cpp \
        VertexCover.cpp \
        hsl_ma86_wrapper.o \
        -rpath $CHOLMOD_PATH/lib/ \
        -rpath $QDLDL_PATH/build/out/ \
        -rpath $LOCAL_PATH/lib/ \
        -DHAVE_CHOLMOD \
        -DHAVE_MA86


# if shared libraries fail to link:
# DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:${HOME}/local/lib
# export DYLD_LIBRARY_PATH
