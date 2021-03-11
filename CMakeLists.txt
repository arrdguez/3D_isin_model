cmake_minimum_required (VERSION 2.8)

project (project_name CXX)

if (CMAKE_VERSION VERSION_LESS "3.1")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else ()
      set (CMAKE_CXX_STANDARD 11)
endif ()

include_directories($(CMAKE_CURRENT_SOURCE_DIR)/include)

include_directories(${PROJECT_SOURCE_DIR}/include)


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(vector main.cpp ${sources})
