# set source files
include(sources)
set(sources ${ipm_sources} ${factor_highs_sources} ${util_sources})
set(headers ${ipm_headers} ${factor_highs_headers} ${util_headers})

# create library
add_library(highspm ${sources} ${headers})

# set compiler flags
target_compile_options(highspm PRIVATE -O3)

target_include_directories(highspm PUBLIC ${PROJECT_SOURCE_DIR}/src)

# BLAS
if(BLAS_LIB)
    target_link_libraries(highspm PRIVATE "${BLAS_LIB}")
else(BLAS_LIB)
    target_link_libraries(highspm PRIVATE "${OPENBLAS_LIB}")
endif(BLAS_LIB)

# Metis
target_include_directories(highspm PRIVATE "${METIS_PATH}")
target_link_libraries(highspm PRIVATE "${METIS_LIB}")

# GKlib
target_include_directories(highspm PRIVATE "${GKLIB_PATH}")
target_link_libraries(highspm PRIVATE "${GKLIB_LIB}")

# HiGHS
target_include_directories(highspm PRIVATE "${HIGHS_PATH}/src" "${HIGHS_PATH}/build")
target_link_libraries(highspm PRIVATE "${HIGHS_LIB}")