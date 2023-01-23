cmake_minimum_required(VERSION 3.13)

set(OBJ lpmcpp)
set(CMAKE_CXX_COMPILER "icx")

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
