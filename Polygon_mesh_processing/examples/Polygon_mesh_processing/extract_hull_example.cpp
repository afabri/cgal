#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polygon_mesh_processing/extract_hull.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <iostream>
#include <fstream>
#include <map>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
//typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;


namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/halfcubes-nonmanifold.off";
  std::ifstream input(filename);

  if (!input)
  {
    std::cerr << "Cannot open file " << std::endl;
    return 1;
  }

  std::vector<Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  if (!CGAL::read_OFF(input, points, polygons))
  {
    std::cerr << "Error parsing the OFF file " << std::endl;
    return 1;
  }

  Surface_mesh mesh;

  PMP::orient_polygon_soup(points, polygons);

  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);

  int index = 0;
  std::map<vertex_descriptor,int> vimap;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
    vimap[vd]= index++;
  }
  index = 0;

  std::map<face_descriptor,int> fimap;
  BOOST_FOREACH(face_descriptor fd, faces(mesh)){
    fimap[fd]= index++;
  }
  PMP::extract_hull(mesh,
    params::vertex_index_map(boost::make_assoc_property_map(vimap)).
    face_index_map(boost::make_assoc_property_map(fimap)));

  std::ofstream out("out.off");
  out << mesh;

  std::cerr << "done"<< std::endl;
  return 0;
}
