#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;

//
// Setup an enriched polyhedron type which stores an id() field in the items
//
typedef CGAL::Surface_mesh<Point_3> SM;
typedef CGAL::Face_filtered_graph<SM> Surface_mesh; 

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor ;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor ;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef SMS::Edge_profile<Surface_mesh> Profile ;



int main( int argc, char** argv ) 
{
  SM sm;
  
  std::ifstream is(argv[1]) ; is >> sm ;
  if (!CGAL::is_triangle_mesh(sm)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<face_descriptor> fa;
  fa.push_back(* faces(sm).first);
  Surface_mesh surface_mesh(sm, fa);
  BOOST_FOREACH(edge_descriptor ed, edges(surface_mesh)){
    std::cout << ed << std::endl;
  }
  
  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.1);
  
  // The index maps are not explicitelty passed as in the previous
  // example because the surface mesh items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  int r = SMS::edge_collapse
           (surface_mesh
           ,stop
            ,CGAL::parameters::get_cost     (SMS::Edge_length_cost  <Surface_mesh>())
                              .get_placement(SMS::Midpoint_placement<Surface_mesh>())
           );
  
  
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << sm ;
  
  return EXIT_SUCCESS ;      
}
