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
if(APPLE)
    target_link_libraries(highspm PRIVATE "-framework Accelerate")
    target_compile_definitions(highspm PRIVATE FRAMEWORK_ACCELERATE)
else(APPLE)
    if(BLAS_LIB)
        target_link_libraries(highspm PRIVATE "${BLAS_LIB}")
    elseif(OPENBLAS_LIB)
        target_link_libraries(highspm PRIVATE "${OPENBLAS_LIB}")
    else()
        message(FATAL_ERROR "No BLAS library available")
    endif(BLAS_LIB)
endif(APPLE)

# Metis
target_include_directories(highspm PRIVATE "${METIS_PATH}")
target_link_libraries(highspm PRIVATE "${METIS_LIB}")

# GKlib
target_include_directories(highspm PRIVATE "${GKLIB_PATH}")
target_link_libraries(highspm PRIVATE "${GKLIB_LIB}")

# HiGHS
target_include_directories(highspm PRIVATE "${HIGHS_PATH}/src" "${HIGHS_PATH}/build")
target_link_libraries(highspm PRIVATE "${HIGHS_LIB}")