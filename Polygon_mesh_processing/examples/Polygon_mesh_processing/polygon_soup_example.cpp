#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <vector>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
int main(int argc, char* argv[])
{
#if 0
  CGAL::Polygon_mesh_processing::internal::M m;
  typedef CGAL::Polygon_mesh_processing::internal::M::iterator It;

  std::cout << m << std::endl;
  m.insert(5);
  std::cout << m << std::endl;
  It it = m.begin();
  std::cout << *it << std::endl;
  
  // m.erase(1);
  //std::cout << m << std::endl;
  m.insert(4);
  std::cout << m << std::endl;
 {
    It it = m.begin();
    std::cout << "dereference" << std::endl;
    std::cout << *it << std::endl;
    std::cout << *++it << std::endl;
 }

  m.insert(3);
  std::cout << m << std::endl;
  m.insert(2);
  std::cout << m << std::endl;
  m.insert(1);
  std::cout << m << std::endl;
  {
    It it = m.begin();
    std::cout << *it << std::endl;
    std::cout << *++it << std::endl;
    std::cout << *++it << std::endl;
    std::cout << *++it << std::endl;
  }
  m.erase(4);
  std::cout << m << std::endl;
  m.erase(3);
  std::cout << m << std::endl;
  m.erase(2);
  std::cout << m << std::endl;
  m.erase(1);
  std::cout << m << std::endl;
  m.erase(5);
  std::cout << m << std::endl;
   m.insert(1);
  std::cout << m << std::endl;

#else
  const char* filename = (argc > 1) ? argv[1] : "data/tet-shuffled.off";
  std::ifstream input(filename);

  if (!input)
  {
    std::cerr << "Cannot open file " << std::endl;
    return 1;
  }

  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  if (!CGAL::read_OFF(input, points, polygons))
  {
    std::cerr << "Error parsing the OFF file " << std::endl;
    return 1;
  }
  CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);

  Polyhedron mesh;

  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);

 std::cerr << "After polygon_soup_to_polygon_mesh reading. Continue?" << std::flush;

 //  std::ofstream out("out.off");
 // out << mesh << std::endl;
#endif

  return 0;
}
