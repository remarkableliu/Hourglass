cmake_minimum_required(VERSION 3.27.8)

project(ADAPTIVE_RANGE_FILTER)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -msse4.2 -Ofast")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -msse")
find_package(OpenSSL REQUIRED)
set(CMAKE_BUILD_TYPE RELEASE)
# add_compile_options(-g)
# add_compile_options(-pg)

include_directories(${PROJECT_SOURCE_DIR})

add_library(
  Hourglass
  STATIC
  hourglass.cpp
  murmur3.c
  hashutil.cc
  )

target_link_libraries(Hourglass OpenSSL::SSL OpenSSL::Crypto)

target_include_directories(Hourglass PRIVATE /usr/lib/openssl)