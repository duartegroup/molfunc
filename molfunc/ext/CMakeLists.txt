cmake_minimum_required(VERSION 3.16)
project(molfunc_ext)

set(CMAKE_CXX_STANDARD 11)
include_directories(include/)

add_library(molfunc STATIC
            include/atoms.h src/atoms.cpp
            include/molecules.h src/molecules.cpp
            include/graph.h src/graph.cpp
            include/utils.h src/utils.cpp )

find_package(Catch2 REQUIRED)


add_executable(test
               tests/test_main.cpp
               tests/test_atoms.cpp
               tests/test_molecule.cpp
               tests/test_graphs.cpp)

target_link_libraries(test PRIVATE Catch2::Catch2)
target_link_libraries(test PUBLIC molfunc)