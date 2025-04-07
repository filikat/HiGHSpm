# set source files
include(sources)
set(sources ${ipm_sources} ${factor_highs_sources} ${util_sources})
set(headers ${ipm_headers} ${factor_highs_headers} ${util_headers})

# create library
add_library(libipm ${sources} ${headers})

# set compiler flags
target_compile_options(libipm PRIVATE -O3 -Wno-deprecated)

include_directories(${PROJECT_SOURCE_DIR}/src)

# find and link stuff:
# (may require to define CMAKE_PREFIX_PATH and/or CMAKE_INCLUDE_PATH)

# BLAS
find_path(BLAS_PATH NAMES "cblas.h" REQUIRED)
message("Found BLAS header at ${BLAS_PATH}")

find_library(BLAS_LIB NAMES blas REQUIRED)
message("Found BLAS library at ${BLAS_LIB}")

target_include_directories(libipm PRIVATE "${BLAS_PATH}")
target_link_libraries(libipm PRIVATE "${BLAS_LIB}")

# Metis
find_path(METIS_PATH NAMES "metis.h" REQUIRED)
message("Found Metis header at ${METIS_PATH}")

find_library(METIS_LIB NAMES metis REQUIRED)
message("Found Metis library at ${METIS_LIB}")

target_include_directories(libipm PRIVATE "${METIS_PATH}")
target_link_libraries(libipm PRIVATE "${METIS_LIB}")

# GKlib
find_path(GKLIB_PATH NAMES "GKlib.h" REQUIRED)
message("Found GKlib header at ${GKLIB_PATH}")

find_library(GKLIB_LIB NAMES GKlib REQUIRED)
message("Found GKlib library at ${GKLIB_LIB}")

target_include_directories(libipm PRIVATE "${GKLIB_PATH}")
target_link_libraries(libipm PRIVATE "${GKLIB_LIB}")

# HiGHS
find_path(HIGHS_HEADER NAMES "Highs.h" REQUIRED)
message("Found HiGHS header at ${HIGHS_HEADER}")

get_filename_component(HIGHS_PATH ${HIGHS_HEADER} DIRECTORY)
find_library(HIGHS_LIB NAMES highs REQUIRED PATHS "${HIGHS_PATH}/build/lib")
message("Found HiGHS library at ${HIGHS_LIB}")

target_include_directories(libipm PUBLIC "${HIGHS_PATH}/src" "${HIGHS_PATH}/build")
target_link_libraries(libipm PUBLIC "${HIGHS_LIB}")