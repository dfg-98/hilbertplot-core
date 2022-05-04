cmake_minimum_required(VERSION 2.8.12)
project(hilbertplot-core)

add_definitions("-std=c++11")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

include_directories(headers)

file(GLOB SRC_FILES    
    "src/*.cpp"
)
file(GLOB HEADER_FILES    
    "headers/*.h"
)

add_library(hilbertplot-core ${SRC_FILES} ${HEADER_FILES})
target_link_libraries(hilbertplot-core ${CONAN_LIBS})