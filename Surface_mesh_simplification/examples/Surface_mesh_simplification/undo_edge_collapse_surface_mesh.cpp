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



struct My_visitor : public SMS::Edge_collapse_visitor_base<Surface_mesh>
{
 void OnCollapsing(Profile const&          
                   ,boost::optional<Point>  placement
                   )
  {
    std::cout << "My_visitor::OnCollapsing()" << std::endl;
  }

  void OnSplitting(vertex_descriptor v, vertex_descriptor vL, vertex_descriptor vR)
  {
    std::cout << "before split of vertex " << v << std::endl; 
  }

  void OnSplit(halfedge_descriptor hd)
  {
    std::cout << "new halfedge "<< hd << std::endl;
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

 
  My_visitor mvis;

  SMS::Edge_collapse_recorder<Surface_mesh,My_visitor> recorder(surface_mesh);



  int rem = SMS::edge_collapse
    (surface_mesh
     ,stop
     ,CGAL::parameters::get_cost(SMS::Edge_length_cost  <Surface_mesh>())
     .get_placement(SMS::Midpoint_placement<Surface_mesh>())
     .visitor(recorder.visitor(mvis))
     );

  int i = 0;

#ifdef CGAL_SMS_DUMP_MESH  
  {
    Surface_mesh sm;
    copy_face_graph(surface_mesh,sm);
    std::string fn("collapse-");
    fn += boost::lexical_cast<std::string>(i++) + ".off";
    
    std::ofstream out(fn.c_str());
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
