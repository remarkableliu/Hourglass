cmake_minimum_required(VERSION 3.28.3)

project(ADAPTIVE_RANGE_FILTER)

option(ADAPT "Build using adaptivity" ON)
option(BUILD_EXAMPLES "Build the examples" ON)
option(BUILD_BENCHMARKS "Build the benchmark targets" ON)

set(CMAKE_CXX_STANDARD 20)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -msse4.2 -Ofast")
endif ()

if (ADAPT)
  ADD_DEFINITIONS(-DADAPT)
endif ()

find_package(OpenSSL REQUIRED)
set(CMAKE_BUILD_TYPE RELEASE)

include_directories(${PROJECT_SOURCE_DIR})

add_library(
  Hourglass
  STATIC
  hourglass.cpp
  hashutil.cc
  )

target_link_libraries(Hourglass OpenSSL::SSL OpenSSL::Crypto)
target_include_directories(Hourglass PRIVATE /usr/lib/openssl)

if (BUILD_EXAMPLES)
    message(STATUS "Building examples")
    add_executable(hourglass example.cpp hashutil.h hashutil.cc)
    target_link_libraries(hourglass OpenSSL::SSL OpenSSL::Crypto)
    target_include_directories(hourglass PRIVATE /usr/lib/openssl)
endif ()

if (BUILD_BENCHMARKS)
    message(STATUS "Building benchmarks")
    add_subdirectory(bench)
endif ()