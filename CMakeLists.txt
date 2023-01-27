if (CMAKE_HOST_WIN32)
    # need at least CMake 3.23 for IntelLLVM support of IntelDPCPP package on Windows
    cmake_minimum_required(VERSION 3.23)
else()
    # CMake 3.20.5 is the minimum recommended for IntelLLVM on Linux
    cmake_minimum_required(VERSION 3.20.5)
endif()

set(OBJ lpmcpp)
set(CMAKE_CXX_COMPILER "icpx")

project(${OBJ} LANGUAGES CXX)
find_package(MKL CONFIG REQUIRED)
find_package(OpenMP)

include_directories(include)
include_directories(examples)
file(GLOB SOURCES "src/*.cpp")
add_executable(${OBJ} ${SOURCES})
set(CMAKE_CXX_STANDARD 11)

TARGET_INCLUDE_DIRECTORIES(${OBJ} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
TARGET_COMPILE_OPTIONS(${OBJ} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
TARGET_LINK_libraries(${OBJ} PUBLIC MKL::MKL OpenMP::OpenMP_CXX)
