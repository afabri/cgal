
cmake_minimum_required(VERSION 3.1...3.23)
project(iglBGL)

include_directories("C:/GitHub/libigl/include")
find_package(CGAL REQUIRED)
find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)

include(CGAL_Eigen3_support)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("main.cpp")
  target_link_libraries(main PUBLIC CGAL::Eigen3_support)
endforeach()
