cmake_minimum_required(VERSION 3.16)
project(molfunc_ext)

set(CMAKE_CXX_STANDARD 17)
include_directories(include/)

add_custom_target(
        run ALL
        COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compile_fragments.py ${CMAKE_CURRENT_SOURCE_DIR}
        BYPRODUCTS src/generated_fragments.cpp
        COMMENT "Generating fragment library"
)


add_library(molfunc STATIC
            include/atoms.h src/atoms.cpp
            include/species/molecules.h src/species/molecules.cpp
            include/graph.h src/graph.cpp
            include/utils.h src/utils.cpp
            include/species/fragments.h src/species/fragments.cpp
            include/species/combined.h src/species/combined.cpp
            include/species/species.h src/species/species.cpp
            include/rotation.h src/rotation.cpp
            include/grid.h src/grid.cpp
            include/vector3d.h src/vector3d.cpp
        src/angles.cpp include/angles.h)

target_compile_options(molfunc PRIVATE -O3)

find_package(Catch2 REQUIRED)


add_executable(test
               tests/test_main.cpp
               tests/test_atoms.cpp
               tests/test_molecule.cpp
               tests/test_graphs.cpp
               tests/test_fragments.cpp
               tests/test_utils.cpp
               tests/test_combination.cpp
               tests/test_rotation.cpp
               tests/test_grid.cpp)

target_link_libraries(test PRIVATE Catch2::Catch2)
target_link_libraries(test PUBLIC molfunc)
