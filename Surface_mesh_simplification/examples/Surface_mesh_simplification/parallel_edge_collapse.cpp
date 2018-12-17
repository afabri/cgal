#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>

#include <CGAL/boost/graph/partition.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                                      Kernel;

typedef Kernel::Point_3                                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                                         Triangle_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

// Usage: ./parallel_edge_collapse input_mesh edge_length (default:0.01) number_of_tasks (default: 8)

int main(int argc, char** argv)
{
  typedef boost::graph_traits<Triangle_mesh>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<Triangle_mesh>::face_descriptor face_descriptor;

  std::ifstream in((argc>1) ? argv[1] : "data/elephant.off");
  double edge_length = (argc>2) ? boost::lexical_cast<double>(argv[2]) : 0.01;
  int number_of_parts = (argc>3) ? boost::lexical_cast<int>(argv[3]) : 8;

  CGAL::Real_timer t;
  t.start();

  Triangle_mesh tm;
  in >> tm;


  std::cerr << "Input: #V = "<< num_vertices(tm) << " #E = "<< num_edges(tm)
            << " #F = " << num_faces(tm) << " read in " << t.time() << " sec." << std::endl;
  t.reset();

  // Partition ID map
  Triangle_mesh::Property_map<face_descriptor, std::size_t> fpmap
    = tm.add_property_map<face_descriptor, std::size_t>("f:partition").first;

  // Set some custom options for METIS
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_CONTIG] = 1; // need to have contiguous subdomains
  options[METIS_OPTION_UFACTOR] = 1;

  CGAL::METIS::partition_dual_graph(tm, number_of_parts, CGAL::parameters::METIS_options(&options)
                                                                          .face_partition_id_map(fpmap));
  std::cerr << "Built partition in " << t.time() << " sec."<< std::endl;
  t.reset();

  // Constrained edges map
  Triangle_mesh::Property_map<edge_descriptor, bool> ecmap
    = tm.add_property_map<edge_descriptor, bool>("e:ec", false).first;

  ecmap[*(edges(tm).first)] = true;

  // Simplify the mesh
  SMS::Bounded_normal_change_placement<SMS::LindstromTurk_placement<Triangle_mesh> > placement;
  SMS::Edge_length_cost<Triangle_mesh> cost;
  SMS::Edge_length_stop_predicate<double> stop(edge_length);

  SMS::parallel_edge_collapse(tm, stop, number_of_parts, fpmap,
                              CGAL::parameters::edge_is_constrained_map(ecmap)
                              .get_placement(placement)
                              .get_cost(cost));

  // Specific to the class Surface_mesh, to clean deleted elements
  tm.collect_garbage();

  std::cerr << "\nSimplify in " << t.time() << " sec.\n"
            << "Result: #V = " << num_vertices(tm)
            << " #E = " << num_edges(tm)
            << " #F = " << num_faces(tm) << std::endl;

  t.reset();

  std::ofstream out("simplified_mesh.off");
  out << tm;
  std::cerr << "Writing result in " << t.time() << " sec." << std::endl;

  return EXIT_SUCCESS;
}
