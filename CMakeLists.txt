cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ghost VERSION 0.1)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(BUILD_TEST ON)

include(FindPkgConfig)
find_package(PkgConfig REQUIRED)
include_directories(include)

set(HEADERS
        include/ghost-protocol.h
)

set(SOURCES
        src/ghost-protocol.cpp
)

add_library(ghost STATIC ${HEADERS} ${SOURCES})
add_subdirectory(libs/blake3)
add_subdirectory(libs/tinyAES)

pkg_check_modules(GMP REQUIRED)
target_link_libraries(ghost PRIVATE gmp blake3 tiny-aes)
target_include_directories(ghost PRIVATE gmp blake3 tiny-aes)

if(BUILD_TEST)
    add_subdirectory(test)
endif()
