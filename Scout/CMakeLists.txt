cmake_minimum_required(VERSION 3.23)
project(Scout)

# ===============================================
# global settings

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


# ===============================================
# options

option(SCOUT_BUILD_SAMPLE "if true, builds the minimal sample" ON)


# ===============================================
# library

add_library(${PROJECT_NAME} STATIC
    src/Scout/Scout.hpp
    src/Scout/Parser.cpp
    src/Scout/Parser.hpp
    src/Scout/Common.hpp
    src/Scout/Relation.cpp
    src/Scout/Relation.hpp
    src/Scout/MatrixOperations.cpp
    src/Scout/MatrixOperations.hpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    src
)

target_compile_options(${PROJECT_NAME} PUBLIC
    # warning + compile settings
    $<$<PLATFORM_ID:Windows>:/MP>
    $<$<PLATFORM_ID:Linux>:-Wall -Wextra>
)


# ===============================================
# sample

if (SCOUT_BUILD_SAMPLE)
    add_subdirectory(sample)
endif()
