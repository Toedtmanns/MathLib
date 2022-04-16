cmake_minimum_required(VERSION 3.12)

project(MathLib)

if(NOT LIBTYPE)
    set(LIBTYPE SHARED)
endif()

if(LIBTYPE STREQUAL "STATIC")
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

target_include_directories(MathLib PUBLIC ./include)