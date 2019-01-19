cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(Simp)

find_package(PCL 1.5 REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (Simp Simp.cpp)
target_link_libraries (Simp Eigen3::Eigen)
target_link_libraries (Simp ${PCL_LIBRARIES})
