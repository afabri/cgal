#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_recorder.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;


#define CGAL_SM

#ifdef CGAL_SM
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
#else
typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Surface_mesh;
#endif

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor ;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor ;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef SMS::Edge_profile<Surface_mesh> Profile ;


struct State {
  Point_3 p;
};

struct My_visitor : public SMS::Edge_collapse_visitor_base<Surface_mesh>
{
  State& state;
  Surface_mesh& sm;

  My_visitor(State& state, Surface_mesh& sm)
    : state(state), sm(sm)
  {}
  
  void OnCollapsing(Profile const& profile          
                    ,boost::optional<Point>  placement
                    )
  {
    std::cout << "My_visitor::OnCollapsing()" << std::endl;
    std::cout << profile.v0() << " will collapse into "<< profile.v1() << std::endl; 
  }
  
  // Called AFTER each edge has been collapsed
  void OnCollapsed( Profile const&, vertex_descriptor )
  {}
  
  void OnSplitting(vertex_descriptor v, vertex_descriptor vL, vertex_descriptor vR)
  {
    std::cout << "before split of vertex " << v << std::endl;
    state.p = sm.point(v);
  }

  void OnSplit(halfedge_descriptor hd)
  {
    std::cout << target(hd, sm) << " at " << state.p << " was split into " 
              << hd << " at " <<  sm.point(source(hd,sm)) << " -> " << sm.point(target(hd,sm)) << std::endl;
  }
};


  
int main( int argc, char** argv ) 
{

  Surface_mesh surface_mesh; 

#ifdef CGAL_SM  
  std::ifstream is(argv[1]) ; is >> surface_mesh;
#else
  OpenMesh::IO::read_mesh(surface_mesh, argv[1]);
#endif
  
  if (!CGAL::is_triangle_mesh(surface_mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // In this example, the simplification stops when the number of undirected edges
  // drops below 50% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.5);

  State state;
  My_visitor mvis(state, surface_mesh);

  SMS::Edge_collapse_recorder<Surface_mesh,My_visitor> recorder(surface_mesh);

  int rem = SMS::edge_collapse
    (surface_mesh
     ,stop
     ,CGAL::parameters::get_cost(SMS::Edge_length_cost  <Surface_mesh>())
     .get_placement(SMS::Midpoint_placement<Surface_mesh>())
     .visitor(recorder.visitor(mvis))
     );


#ifdef CGAL_SMS_DUMP_MESH  
  {
    Surface_mesh sm;
    copy_face_graph(surface_mesh,sm);
    std::ofstream out("collapse-0.off");
    out << sm;
  }
#endif

  recorder.undo(mvis);

#ifdef CGAL_SM  
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface_mesh ;
#else  
  surface_mesh.garbage_collection();
  OpenMesh::IO::write_mesh(surface_mesh, "out.off");
#endif

  return EXIT_SUCCESS ;      
}
