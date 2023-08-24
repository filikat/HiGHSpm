#!/bin/bash

YOUR_CPP_COMPILER="g++"

ProtoIPM_HOME="$HOME/ProtoIPM"

# Define directories
#
# Start with location of particular codes
SPRAL_HOME=$ProtoIPM_HOME/spral-2023.03.29
HIGHS_HOME=$HOME/install
MA86_HOME=$ProtoIPM_HOME/hsl_ma86-1.7.3
QDLDL_HOME=$ProtoIPM_HOME/qdldl
CHOLMOD_HOME=$ProtoIPM_HOME/SuiteSparse-7.0.1/CHOLMOD
BOOST_HOME=/usr

# Now the include folders
SPRAL_INCLUDE_DIR=$SPRAL_HOME/include
HIGHS_INCLUDE_DIR=$HIGHS_HOME/include/highs
QDLDL_INCLUDE_DIR=$QDLDL_HOME/include/
CHOLMOD_INCLUDE_DIR=$CHOLMOD_HOME/Include/
MA86_INCLUDE_DIR=$MA86_HOME/include/
BOOST_INCLUDE_DIR=$BOOST_HOME/include

# Now the library folders
LIB_DIR=/usr/local/lib/
LIB64_DIR=/usr/lib/x86_64-linux-gnu/
HIGHS_LIB=$HIGHS_HOME/lib
LIBQDLDL_DIR=$QDLDL_HOME/build/out/
LIBCHOLOMOD_DIR=$CHOLMOD_HOME

MA86WRAPPER_C=hsl_ma86_wrapper.c
#MA86WRAPPER_O=hsl_ma86_wrapper.o

# Define libraries
#LIBS="-lboost_program_options -lhighs -lspral -lblas -llapack -lqdldl -lm -lstdc++ -lgfortran -lz -lcholmod -lhsl_ma86 -lmetis"
LIBS="-lboost_program_options -lhighs -lspral -lblas -llapack -lm -lstdc++ -lgfortran -lz -lmetis"

# Define the compiler flags to use in C and C++.
# 
# When you're doing performance tests, use -O3
COMPILER_FLAGS="-g -fopenmp"
#COMPILER_FLAGS="-O3 -fopenmp -DNDEBUG"

# Possibly define the compiler preprocessor settings
#COMPILER_PREPROCESS="-DHAVE_SPRAL  -DHAVE_QDLDL -DHAVE_CHOLMOD -DHAVE_MA86"
COMPILER_PREPROCESS="-DHAVE_SPRAL"

# Define source file - that might change
#SOURCE_FILE="testSolve.cpp Direct.cpp ExperimentData.cpp VectorOperations.cpp"
SOURCE_FILE="mainIPM.cpp IPM_caller.cpp IPM_model.cpp NormalEquations.cpp ConjugateGradient.cpp Direct.cpp VectorOperations.cpp IPM_aux.cpp ExperimentData.cpp"

OUTPUT_FILE=a.out
rm $OUTPUT_FILE

# Compile
$YOUR_CPP_COMPILER \
    $COMPILER_FLAGS \
    $COMPILER_PREPROCESS \
    $SOURCE_FILE \
    $MA86WRAPPER_O \
    -I$SPRAL_INCLUDE_DIR \
    -I$HIGHS_INCLUDE_DIR \
    -I$QDLDL_INCLUDE_DIR \
    -I$CHOLMOD_INCLUDE_DIR \
    -I$BOOST_INCLUDE_DIR \
    $LIBS \
    -L$LIB_DIR \
    -L$LIB64_DIR \
    -L$LIBCHOLOMOD_DIR \
    -L$HIGHS_LIB \
    -L$LIBQDLDL_DIR -Wl,-rpath=$ProtoIPM_HOME/qdldl/build/out/ \
    -o $OUTPUT_FILE -fno-stack-protector 

