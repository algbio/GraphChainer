cmake_minimum_required(VERSION 3.17)
project(GraphAlgChainer)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O9") # Remove this when debugging

include_directories(${LEMON_SOURCE_DIR})
include_directories(${LEMON_BINARY_DIR})
include_directories(includes)


add_subdirectory(external)
add_subdirectory(src)


add_executable(main main.cpp)
target_link_libraries(main lemon)
target_link_libraries(main mpc)