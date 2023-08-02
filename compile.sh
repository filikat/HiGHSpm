#!/bin/bash

YOUR_CPP_COMPILER="g++"

ProtoIPM_HOME="$HOME/ProtoIPM"

# Define directories
#
# Start with location of particular codes
SPRAL_HOME=/usr/local/
HIGHS_HOME=$HOME/install

# Now the include folders
SPRAL_INCLUDE_DIR=$SPRAL_HOME/include
HIGHS_INCLUDE_DIR=$HIGHS_HOME/include/highs

# Now the library folders
LIB_DIR=/usr/local/lib/
LIB64_DIR=/usr/lib/x86_64-linux-gnu/
HIGHS_LIB=$HIGHS_HOME/lib

# Define libraries
LIBS="-lhighs -lspral -lblas -llapack -lm -lstdc++ -lgfortran -lz -lmetis"

# Define the compiler flags to use in C and C++.
# 
# When you're doing performance tests, use -O3
COMPILER_FLAGS="-g -fopenmp "
#COMPILER_FLAGS="-O3 -fopenmp "

# Define source file - that might change
#SOURCE_FILE="testSolve.cpp Direct.cpp ExperimentData.cpp VectorOperations.cpp"
SOURCE_FILE="mainIPM.cpp IPM_caller.cpp IPM_model.cpp NormalEquations.cpp ConjugateGradient.cpp Direct.cpp VectorOperations.cpp IPM_aux.cpp ExperimentData.cpp"

OUTPUT_FILE=a.out
rm $OUTPUT_FILE

# Compilerm
$YOUR_CPP_COMPILER \
    $COMPILER_FLAGS \
    $SOURCE_FILE \
    -I$SPRAL_INCLUDE_DIR \
    -I$HIGHS_INCLUDE_DIR \
    $LIBS \
    -L$LIB_DIR \
    -L$LIB64_DIR \
    -L$HIGHS_LIB \
    -o $OUTPUT_FILE -fno-stack-protector 

