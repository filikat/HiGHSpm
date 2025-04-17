# find stuff (may require to define CMAKE_PREFIX_PATH and/or CMAKE_INCLUDE_PATH)

# BLAS
if(NOT APPLE)
    find_library(BLAS_LIB NAMES blas)
    if(BLAS_LIB)
        message("Found BLAS library at ${BLAS_LIB}")
    else(BLAS_LIB)
        find_library(OPENBLAS_LIB NAMES openblas)
        if(NOT OPENBLAS_LIB)
            message(FATAL_ERROR "No BLAS library found")
        endif(NOT OPENBLAS_LIB)
        message("Found OpenBLAS library at ${OPENBLAS_LIB}")
    endif(BLAS_LIB)
endif(NOT APPLE)

# Metis
find_path(METIS_PATH NAMES "metis.h" REQUIRED)
message("Found Metis header at ${METIS_PATH}")

find_library(METIS_LIB NAMES metis REQUIRED)
message("Found Metis library at ${METIS_LIB}")

# GKlib
find_path(GKLIB_PATH NAMES "GKlib.h" REQUIRED)
message("Found GKlib header at ${GKLIB_PATH}")

find_library(GKLIB_LIB NAMES GKlib REQUIRED)
message("Found GKlib library at ${GKLIB_LIB}")

# HiGHS
find_path(HIGHS_HEADER NAMES "Highs.h" REQUIRED)
message("Found HiGHS header at ${HIGHS_HEADER}")
get_filename_component(HIGHS_PATH ${HIGHS_HEADER} DIRECTORY)

find_library(HIGHS_LIB NAMES highs REQUIRED PATHS "${HIGHS_PATH}/build/lib")
message("Found HiGHS library at ${HIGHS_LIB}")