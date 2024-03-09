cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ghost VERSION 0.1)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(BUILD_TEST ON)

set(HEADERS
    include/proto.h
)

set(SOURCES
    src/proto.c
        test/test.c
)

find_package(PkgConfig REQUIRED)
include_directories(include)

add_library(ghost STATIC ${HEADERS} ${SOURCES})

pkg_check_modules(GMP REQUIRED)
target_link_libraries(ghost PRIVATE gmp)
target_include_directories(ghost PRIVATE gmp)

if(BUILD_TEST)
    add_subdirectory(test)
endif()
