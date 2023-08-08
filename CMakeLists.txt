cmake_minimum_required(VERSION 3.5)
project(e3a LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same.")
endif()

include(FetchContent)

find_package(Eigen3)
if(NOT Eigen3_FOUND)
    message("-- Downloading Eigen")
    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG master
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    set(EIGEN_BUILD_DOC OFF)
    set(BUILD_TESTING OFF)
    set(EIGEN_BUILD_PKGCONFIG OFF)
    FetchContent_MakeAvailable(Eigen)
endif(NOT Eigen3_FOUND)

option(LIBIGL_RESTRICTED_TRIANGLE ON)
message("-- Downloading libigl")
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG main
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(libigl)
igl_include(restricted triangle)

message("-- Downloading GMP")
FetchContent_Declare(
    GMP
    GIT_REPOSITORY https://github.com/alisw/GMP
    GIT_TAG master
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(gmp)


add_subdirectory(src)
add_subdirectory(examples)
