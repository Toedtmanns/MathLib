cmake_minimum_required(VERSION 3.12)

project(MathLib)

option(COMPILE_STATIC "Compile MathLib to a static library" ON)

if(COMPILE_STATIC)
    add_library(MathLib STATIC
            ./source/Vectors.cpp
            ./source/ComplexNumbers.cpp
            ./source/Matrices.cpp
            ./source/Primitives.cpp
            ./source/Utility.cpp
            ./source/Geometry.cpp)
    
    add_compile_definitions(MATHLIB_STATIC)
else()
    add_library(MathLib SHARED
            ./source/Vectors.cpp
            ./source/ComplexNumbers.cpp
            ./source/Matrices.cpp
            ./source/Primitives.cpp
            ./source/Utility.cpp
            ./source/Geometry.cpp)
            
    add_compile_definitions(MATHLIB_EXPORTS)
endif()

set_target_properties(MathLib PROPERTIES PUBLIC_HEADER ./include/Maths.h)