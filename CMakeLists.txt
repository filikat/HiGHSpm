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

# tests
enable_testing()

add_test(NAME test_afiro_ne COMMAND ipm ../instances/afiro.mps 1 2 0)
set_tests_properties(test_afiro_ne PROPERTIES 
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")

add_test(NAME test_afiro_as COMMAND ipm ../instances/afiro.mps 0 2 0)
set_tests_properties(test_afiro_ne PROPERTIES 
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")

add_test(NAME test_afiro_ne_crossover COMMAND ipm ../instances/afiro.mps 1 2 1)
set_tests_properties(test_afiro_ne_crossover PROPERTIES
                    PASS_REGULAR_EXPRESSION "Status crossover: +optimal")

add_test(NAME test_afiro_as_crossover COMMAND ipm ../instances/afiro.mps 0 2 1)
set_tests_properties(test_afiro_as_crossover PROPERTIES
                    PASS_REGULAR_EXPRESSION "Status crossover: +optimal")

add_test(NAME test_cycle COMMAND ipm ../instances/cycle.mps 1 2 0)
set_tests_properties(test_cycle PROPERTIES
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")

add_test(NAME test_pilotwe COMMAND ipm ../instances/pilot.we.mps 1 2 0)
set_tests_properties(test_pilotwe PROPERTIES
                    PASS_REGULAR_EXPRESSION "Primal-dual feasible point found")                    