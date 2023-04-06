cmake_minimum_required(VERSION 3.23)
project(Scout)

set(CMAKE_CXX_STANDARD 20)

add_library(${PROJECT_NAME}  STATIC
        src/Scout.cpp
        src/Scout.hpp
        src/Parser.cpp
        src/Parser.hpp
        src/Common.hpp
        src/Relation.cpp
        src/Relation.hpp
        src/MatrixOperations.cpp
        src/MatrixOperations.hpp
        src/Scout.hpp)
