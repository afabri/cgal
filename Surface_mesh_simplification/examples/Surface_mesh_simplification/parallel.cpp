
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/lexical_cast.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/partition.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Surface_mesh_simplification/parallel_edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>


typedef CGAL::Simple_cartesian<double> Kernel;

typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh; 

typedef Surface_mesh::Property_map<boost::graph_traits<Surface_mesh>::face_descriptor,std::size_t> CCMap;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;




int main(int argc, char** argv ) 
{
  bool dump = false;

  typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
 
  CGAL::Real_timer t;
  t.start();
  std::ifstream in(argv[1]);


  Surface_mesh sm; 
  in >> sm;
  std::cerr << "Input: #V = "<< num_vertices(sm)  << " #E = "<< num_edges(sm) 
            << " #F = " << num_faces(sm)  <<  " read in " << t.time() << " sec." << std::endl;
  t.reset();
  double ratio = (argc>2)?boost::lexical_cast<double>(argv[2]):0.25;

  typedef Surface_mesh::Property_map<face_descriptor,std::size_t> CCMap;
  Surface_mesh::Property_map<face_descriptor,std::size_t> ccmap 
    = sm.add_property_map<face_descriptor,std::size_t>("f:cc").first;

  std::size_t ncc = 8;
  unsigned int layers = 1;
  bool verbose = false;
  PMP::partition(sm, ccmap, static_cast<int>(ncc));

  std::cerr << "Partition in " << t.time() << " sec."<< std::endl;
  t.reset();
  SMS::LindstromTurk_placement<Surface_mesh> placement;
  SMS::LindstromTurk_cost<Surface_mesh> cost;
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  //SMS::Midpoint_placement<Surface_mesh> placement;
  //SMS::Edge_length_cost<Surface_mesh> cost;
  //SMS::Edge_length_stop_predicate<double> stop(0.01);

  SMS::parallel_edge_collapse(sm, ccmap, placement, stop, cost, ncc, layers, dump, verbose);

  sm.collect_garbage();

  std::cerr << "\nSimplify in " << t.time() << " sec.\n"
            << "Result: #V = " << num_vertices(sm)
            << " #E = " << num_edges(sm) 
            << " #F = " << num_faces(sm) << std::endl;

  t.reset();
  {
    std::ofstream out("out.off");
    out << sm << std::endl;
    out.close();
    std::cerr << "Writing result in " << t.time() << " sec." << std::endl;
  }
  return EXIT_SUCCESS;
}
