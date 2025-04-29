cmake_minimum_required(VERSION 3.15)

project(ProtoIpm VERSION 1.0)

# set cpp version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# include cmake directory in search path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include("find_stuff")

# create main executable
add_executable(ipm src/main.cpp)

# set compiler flags
target_compile_options(ipm PRIVATE -O3 -Wno-deprecated)

add_subdirectory(src)

target_link_libraries(ipm highspm)
target_include_directories(ipm PRIVATE "${HIGHS_PATH}/src" "${HIGHS_PATH}/build")

# flags for address sanitizer
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
#set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")

# executable to run benchmark problems
add_executable(run_test src/mainTest.cpp)
target_compile_options(run_test PRIVATE -O3 -Wno-deprecated)
target_link_libraries(run_test highspm)
target_include_directories(run_test PRIVATE "${HIGHS_PATH}/src" "${HIGHS_PATH}/build")

# tests
enable_testing()

add_test(NAME test_ne COMMAND ipm ../instances/afiro.mps 1 2 0)
set_tests_properties(test_ne PROPERTIES 
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")

add_test(NAME test_as COMMAND ipm ../instances/afiro.mps 0 2 0)
set_tests_properties(test_as PROPERTIES 
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")

add_test(NAME test_ne_crossover COMMAND ipm ../instances/afiro.mps 1 2 1)
set_tests_properties(test_ne_crossover PROPERTIES
                    PASS_REGULAR_EXPRESSION "Status crossover: +optimal")

add_test(NAME test_as_crossover COMMAND ipm ../instances/afiro.mps 0 2 1)
set_tests_properties(test_as_crossover PROPERTIES
                    PASS_REGULAR_EXPRESSION "Status crossover: +optimal")                 