#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#define  CGAL_NO_CDT_2_WARNING
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/natural_neighbor_coordinates_3.h>
#include <CGAL/point_generators_3.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

typedef double NT; //Number Type

typedef CGAL::Exact_predicates_exact_constructions_kernel    K;

typedef K::Point_3                                             Point3;
typedef K::Vector_3                                            Vector3;
typedef K::Sphere_3                                            Sphere_3;

typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Dh;

typedef Dh::Facet                                              Facet;
typedef Dh::Vertex_handle                                      Vertex_handle;
typedef Dh::Cell_handle                                        Cell_handle;
typedef Dh::Finite_vertices_iterator                           Finite_vertices_iterator;
typedef Dh::Vertex_iterator                                    Vertex_iterator;
typedef Dh::Facet_circulator                                   Facet_circulator;
typedef Dh::Cell_iterator                                      Cell_iterator;

int main()
{
  Dh triangulation;

  CGAL::Random_points_in_sphere_3<Point3> rpg(1.0, CGAL::Random(0));

  Point3 p, q(0,0,0);
  std::ofstream ps("points.xyz");
  for(int i=0; i < 100; i++){
    Point3 p = *rpg++;

    ps << p << std::endl;
    triangulation.insert(p);
  }

  

  CGAL::natural_neighbor_coordinates_3(triangulation,q);

  return 0;
}
