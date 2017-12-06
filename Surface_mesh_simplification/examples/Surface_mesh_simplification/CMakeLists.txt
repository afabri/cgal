# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Surface_mesh_simplification_Examples )


cmake_minimum_required(VERSION 2.8.11)



# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )

find_package( TBB )
    if( TBB_FOUND )
      include(${TBB_USE_FILE})
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
    endif()

# Boost and its components
find_package( Boost )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()


find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

find_package(METIS)

if(METIS_FOUND)
  include_directories (${METIS_INCLUDE_DIRS})
else()
  message(STATUS "Examples that use the metis library will not be compiled.")
endif()



# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )


create_single_source_cgal_program( "edge_collapse_constrain_sharp_edges.cpp" )
create_single_source_cgal_program( "edge_collapse_constrained_border_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_enriched_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_surface_mesh.cpp" )
create_single_source_cgal_program( "edge_collapse_linear_cell_complex.cpp" )
create_single_source_cgal_program( "edge_collapse_constrained_border_surface_mesh.cpp" )
create_single_source_cgal_program( "edge_collapse_all_short_edges.cpp" )
create_single_source_cgal_program( "edge_collapse_bounded_normal_change.cpp" )

if(OpenMesh_FOUND)
  create_single_source_cgal_program( "edge_collapse_OpenMesh.cpp" )

  target_link_libraries( edge_collapse_OpenMesh PRIVATE ${OPENMESH_LIBRARIES} )
endif()

if(METIS_FOUND)
  create_single_source_cgal_program( "parallel_edge_collapse.cpp" )
  target_link_libraries( parallel_edge_collapse ${METIS_LIBRARIES} )
  
  create_single_source_cgal_program( "parallel_edge_collapse_Polyhedron.cpp" )
  target_link_libraries( parallel_edge_collapse_Polyhedron ${METIS_LIBRARIES} )
endif()
