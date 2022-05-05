cmake_minimum_required(VERSION 2.8.12)
project(hilbertplot-core)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 11)

add_definitions("-std=c++11")

include_directories(headers)

file(GLOB SRC_FILES    
    "src/*.cpp"
)
file(GLOB HEADER_FILES    
    "include/*.h"
)

add_library(hilbertplot-core ${SRC_FILES} ${HEADER_FILES})
set_target_properties(hilbertplot-core PROPERTIES PUBLIC_HEADER ${HEADER_FILES}  )
target_link_libraries(hilbertplot-core ${CONAN_LIBS})
target_include_directories(hilbertplot-core PUBLIC include)

install(TARGETS hilbertplot-core)
