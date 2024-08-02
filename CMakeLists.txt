cmake_minimum_required(VERSION 3.27.8)

project(ADAPTIVE_RANGE_FILTER)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -msse4.2 -Ofast")
find_package(OpenSSL REQUIRED)
set(CMAKE_BUILD_TYPE RELEASE)
ADD_DEFINITIONS(-DADAPT)

include_directories(${PROJECT_SOURCE_DIR})

add_library(
  Hourglass
  STATIC
  hourglass.cpp
  hashutil.cc
  )

add_executable(hourglass example.cpp hashutil.h hashutil.cc)

target_link_libraries(hourglass OpenSSL::SSL OpenSSL::Crypto)

target_include_directories(hourglass PRIVATE /usr/lib/openssl)

target_link_libraries(Hourglass OpenSSL::SSL OpenSSL::Crypto)

target_include_directories(Hourglass PRIVATE /usr/lib/openssl)

add_subdirectory(bench)

# if (BUILD_EXAMPLES)
#     message(STATUS "Building examples")
#     add_subdirectory(examples)
# endif ()

# if (BUILD_BENCHMARKS)
#     message(STATUS "Building benchmarks")
#     add_subdirectory(bench)
# endif ()