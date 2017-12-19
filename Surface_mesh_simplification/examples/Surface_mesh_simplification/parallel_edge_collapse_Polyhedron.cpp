#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>

#include <CGAL/Polygon_mesh_processing/partition.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>
#include <map>

typedef CGAL::Simple_cartesian<double>                                      Kernel;

typedef Kernel::Point_3                                                     Point_3;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>        Triangle_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

// Usage: ./parallel_edge_collapse input_mesh keep_ratio (default:0.25) number_of_tasks (default: 8)

int main(int argc, char** argv )
{
  typedef boost::graph_traits<Triangle_mesh>::edge_descriptor      edge_descriptor;
  typedef boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

  std::ifstream in((argc>1) ? argv[1] : "data/elephant.off");
  double ratio = (argc>2) ? boost::lexical_cast<double>(argv[2]) : 0.25;
  int ncc = (argc>3) ? boost::lexical_cast<int>(argv[3]) : 8;

  if(!in)
  {
    std::cerr << "Error: could not open input mesh" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Real_timer t;
  t.start();

  Triangle_mesh tm;
  in >> tm;

  CGAL::set_halfedgeds_items_id(pm);

  std::cerr << "Input: #V = "<< num_vertices(tm) << " #E = "<< num_edges(tm)
            << " #F = " << num_faces(tm) << " read in " << t.time() << " sec." << std::endl;
  t.reset();

  std::map<face_descriptor, std::size_t> ccmap;
  boost::associative_property_map<std::map<face_descriptor, std::size_t> > ccpmap(ccmap);

  PMP::partition(tm, ncc /*number of partitions*/, ccpmap);
  std::cerr << "Built partition in " << t.time() << " sec." << std::endl;
  t.reset();

  std::map<edge_descriptor, bool> ecmap;
  boost::associative_property_map<std::map<edge_descriptor, bool> > ecpmap(ecmap);
  BOOST_FOREACH(edge_descriptor ed, edges(tm)) {
    ecmap[ed] = false;
  }

  ecmap[*(edges(tm).first)] = true; // arbitrarily constrain the first edge

  SMS::Bounded_normal_change_placement<SMS::LindstromTurk_placement<Triangle_mesh> > placement;
  SMS::LindstromTurk_cost<Triangle_mesh> cost;
  SMS::Count_ratio_stop_predicate<Triangle_mesh> stop(ratio);

  SMS::parallel_edge_collapse(tm, stop, ccpmap, ncc,
                              CGAL::parameters::edge_is_constrained_map(ecpmap)
                              .get_placement(placement)
                              .get_cost(cost));

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
