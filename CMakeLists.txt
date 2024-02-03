cmake_minimum_required(VERSION 3.20)

set(OBJ lpmcpp)
set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_STANDARD 17)

project(${OBJ} LANGUAGES CXX)

find_package(OpenBLAS REQUIRED)
find_package(Boost REQUIRED COMPONENTS)
find_package(OpenMP)

include_directories(include)
include_directories(examples)
include_directories(geometry)
include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${OpenBLAS_INCLUDE_DIRS} )

file(GLOB SOURCES "src/*.cpp")
add_executable(${OBJ} ${SOURCES})

# Link OpenMP libraries
target_link_libraries(${OBJ} PUBLIC OpenMP::OpenMP_CXX)

# Link OpenBLAS libraries
target_link_libraries(${OBJ} PUBLIC ${OpenBLAS_LIBRARIES})

# link Boost libraries
target_link_libraries(${OBJ} PUBLIC ${Boost_LIBRARIES} )
