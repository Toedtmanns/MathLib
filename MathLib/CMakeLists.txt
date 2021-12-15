cmake_minimum_required(VERSION 3.10)

project(MathLib)

add_library(MathLib SHARED
        ./source/2DVectors.cpp
        ./source/3DVectors.cpp
        ./source/ComplexNumbers.cpp
        ./source/Matrices.cpp
        ./source/Primitives.cpp
        ./source/Utility.cpp)
set_target_properties(MathLib PROPERTIES PUBLIC_HEADER ./include/Maths.h)

option(MATHLIB_EXPORTS "Compile a dynamic library" ON)
configure_file(./include/Maths.h ./include/Maths.h)