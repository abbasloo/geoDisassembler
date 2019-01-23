cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(Mesh)

find_package(PCL 1.5 REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VTK REQUIRED)

include(${VTK_USE_FILE})
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (Mesh Mesh.cpp)
target_link_libraries (Mesh Eigen3::Eigen)
target_link_libraries (Mesh ${PCL_LIBRARIES})
target_link_libraries (Mesh ${VTK_LIBRARIES})

