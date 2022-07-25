cmake_minimum_required(VERSION 3.12)

project(MathLib)

set(SOURCES 
    ./source/Vectors.cpp
    ./source/ComplexNumbers.cpp
    ./source/Matrices.cpp
    ./source/Primitives.cpp
    ./source/Utility.cpp
    ./source/Geometry.cpp)

if(NOT LIBTYPE)
    set(LIBTYPE SHARED)
endif()

if(LIBTYPE STREQUAL "STATIC")
    add_library(MathLib STATIC ${SOURCES})
    add_compile_definitions(MATHLIB_STATIC)
else()
    add_library(MathLib SHARED ${SOURCES})     
    add_compile_definitions(MATHLIB_EXPORTS)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
endif()

target_include_directories(MathLib PUBLIC ./include)