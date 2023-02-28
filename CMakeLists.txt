cmake_minimum_required(VERSION 3.1...3.25)

project(GeomAlgo LANGUAGES CXX)

# Find CGAL and its components
find_package(CGAL QUIET)

if (NOT CGAL_FOUND)
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()


# Boost and its components
find_package( Boost REQUIRED )

if (NOT Boost_FOUND)
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

add_subdirectory(src)
