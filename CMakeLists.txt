cmake_minimum_required(VERSION 3.17)
project(useful_algorithms)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER mpic++)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
find_package(MPI REQUIRED)

include_directories(include)
ENABLE_TESTING()

add_executable(ser_tests
        include/serial/algorithm.hpp
        tests/ser_tests.cpp)

add_executable(par_tests
        include/parallel/algorithm.hpp
        tests/par_tests.cpp)

add_test(NAME serial COMMAND $<TARGET_FILE:ser_tests>)
add_test(NAME parallel COMMAND $<TARGET_FILE:par_tests>)

target_include_directories(par_tests PRIVATE ${MPI_C_INCLUDE_PATH})
target_link_libraries(par_tests PRIVATE ${MPI_C_LIBRARY})