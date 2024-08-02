cmake_minimum_required(VERSION 3.28.3)
project(rangefilters-bench)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG        v3.0
)
FetchContent_MakeAvailable(argparse)

set(CMAKE_CXX_STANDARD 20)

list(APPEND Targets "hourglass")

function(compile_bench ds)
    if (NOT ds IN_LIST Targets)
        message(WARNING "${ds} is not defined")
        return()
    endif()
    if (ds IN_LIST x86Targets AND NOT CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
        message(WARNING "cannot compile ${ds} on this architecture")
        return()
    endif()
    add_executable(bench_${ds} bench_${ds}.cpp)
    target_link_libraries(bench_${ds} argparse)

    if (ds MATCHES "hourglass")
        target_link_libraries(bench_hourglass Hourglass)
    endif()
endfunction(compile_bench)

foreach(ds ${Targets})
    compile_bench(${ds})
endforeach()

add_executable(workload_gen workload_gen.cpp)
target_link_libraries(workload_gen argparse)