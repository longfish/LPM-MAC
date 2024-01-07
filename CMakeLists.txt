if (CMAKE_HOST_WIN32)
    # need at least CMake 3.23 for Windows
    cmake_minimum_required(VERSION 3.23)
else()
    cmake_minimum_required(VERSION 3.20.5)
endif()

set(OBJ lpmcpp)
set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_STANDARD 17)

project(${OBJ} LANGUAGES CXX)

# please change this line of code
include_directories(/Users/changyu/Applications/eigen-3.4.0)
find_package(OpenBLAS REQUIRED)
find_package(Boost REQUIRED COMPONENTS)
find_package(OpenMP)

include_directories(include)
include_directories(examples)
include_directories(geometry)
include_directories( ${BOOST_INCLUDE_DIRS} )
file(GLOB SOURCES "src/*.cpp")
add_executable(${OBJ} ${SOURCES})

# Link OpenBLAS and OpenMP libraries
target_link_libraries(${OBJ} PUBLIC OpenBLAS::OpenBLAS OpenMP::OpenMP_CXX)

# link Boost libraries
target_link_libraries(${OBJ} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_REGEX_LIBRARY} )
