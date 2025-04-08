cmake_minimum_required(VERSION 3.15)

project(ProtoIpm VERSION 1.0)

# set cpp version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# include cmake directory in search path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# create main executable
add_executable(ipm src/main.cpp)

# set compiler flags
target_compile_options(ipm PRIVATE -O3 -Wno-deprecated)

add_subdirectory(src)

target_link_libraries(ipm highspm)