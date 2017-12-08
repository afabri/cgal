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
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3>         Surface_mesh;

typedef std::vector<boost::graph_traits<Surface_mesh>::vertices_size_type>  Polygon_3;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

// Usage: ./parallel_edge_collapse input_mesh keep_ratio (default:0.25) number_of_tasks (default: 8)

int main(int argc, char** argv )
{
  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor    vertex_descriptor;
  typedef boost::graph_traits<Surface_mesh>::edge_descriptor      edge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::face_descriptor      face_descriptor;

  CGAL::Real_timer t;
  t.start();

  Surface_mesh sm;
  std::ifstream in((argc>1) ? argv[1] : "data/elephant.ply", std::ios_base::binary);

  std::vector<Point_3> points;
  std::vector<Polygon_3> polygons;
  CGAL::read_PLY(in, points, polygons);
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh (points, polygons, sm);

  int index = 0;

  BOOST_FOREACH(halfedge_descriptor hd, halfedges(sm)) {
    hd->id() = index++;
  }

  index = 0 ;
  BOOST_FOREACH(vertex_descriptor vd, vertices(sm)) {
    vd->id() = index++;
  }

  std::cerr << "Input: #V = "<< num_vertices(sm)  << " #E = "<< num_edges(sm)
            << " #F = " << num_faces(sm)  <<  " read in " << t.time() << " sec." << std::endl;
  t.reset();
  double ratio = (argc>2) ? boost::lexical_cast<double>(argv[2]) : 0.25;

  std::map<face_descriptor,std::size_t> ccmap;
  boost::associative_property_map<std::map<face_descriptor,std::size_t> > ccpmap(ccmap);

  std::map<edge_descriptor,bool> ecmap;
  boost::associative_property_map<std::map<edge_descriptor, bool> > ecpmap(ecmap);
  BOOST_FOREACH(edge_descriptor ed, edges(sm)){
    ecmap[ed] = false;
  }

  ecmap[*(edges(sm).first)] = true;

  int ncc = (argc>3) ? boost::lexical_cast<int>(argv[3]) : 8;

  if(ncc > 1) {
    PMP::partition(sm, ccpmap, ncc);

    std::cerr << "Built partition in " << t.time() << " sec." << std::endl;
    t.reset();
  }

  SMS::Bounded_normal_change_placement<SMS::LindstromTurk_placement<Surface_mesh> > placement;
  SMS::LindstromTurk_cost<Surface_mesh> cost;
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  if(ncc > 1)
  {
    SMS::parallel_edge_collapse(sm, stop, ccpmap, ncc,
                                CGAL::parameters::get_placement(placement)
                                  .edge_is_constrained_map(ecpmap)
                                  .get_cost(cost));
  }
  else
  {
    SMS::edge_collapse(sm, stop,
                       CGAL::parameters::get_placement(placement)
                         .edge_is_constrained_map(ecpmap)
                         .get_cost(cost));
  }

  std::cerr << "\nSimplify in " << t.time() << " sec.\n"
            << "Result: #V = " << num_vertices(sm)
            << " #E = " << num_edges(sm)
            << " #F = " << num_faces(sm) << std::endl;

  t.reset();

  std::ofstream out("outP.ply", std::ios_base::binary);
  CGAL::set_binary_mode(out);
  CGAL::write_PLY(out, sm);
  out.close();

  std::cerr << "Writing result in " << t.time() << " sec." << std::endl;

  return EXIT_SUCCESS;
}
