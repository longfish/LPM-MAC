if (CMAKE_HOST_WIN32)
    # need at least CMake 3.23 for IntelLLVM support of IntelDPCPP package on Windows
    cmake_minimum_required(VERSION 3.23)
else()
    # CMake 3.20.5 is the minimum recommended for IntelLLVM on Linux
    cmake_minimum_required(VERSION 3.20.5)
endif()

set(OBJ lpmcpp)
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)

project(${OBJ} LANGUAGES CXX)
find_package(MKL CONFIG REQUIRED)
find_package(Boost REQUIRED COMPONENTS)
find_package(OpenMP)

include_directories(include)
include_directories(examples)
include_directories( ${BOOST_INCLUDE_DIRS} )
file(GLOB SOURCES "src/*.cpp")
add_executable(${OBJ} ${SOURCES})

TARGET_INCLUDE_DIRECTORIES(${OBJ} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
TARGET_COMPILE_OPTIONS(${OBJ} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
TARGET_LINK_libraries(${OBJ} PUBLIC MKL::MKL OpenMP::OpenMP_CXX)
target_link_libraries(${OBJ} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_REGEX_LIBRARY} )