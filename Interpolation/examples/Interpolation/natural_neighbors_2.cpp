#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K>              Delaunay_triangulation;
typedef Delaunay_triangulation::Vertex_handle Vertex_handle;
typedef std::vector< std::pair< Vertex_handle, K::FT  > > Vertex_handle_coordinate_vector;

int main()
{
  Delaunay_triangulation dt;

  for (int y=0 ; y<3 ; y++)
    for (int x=0 ; x<3 ; x++)
      dt.insert(K::Point_2(x,y));

  //coordinate computation
  K::Point_2 p(1.2, 0.7);
  Vertex_handle_coordinate_vector coords;
  CGAL::Triple<std::back_insert_iterator<Vertex_handle_coordinate_vector>,
               K::FT, bool> result =
      CGAL::natural_neighbors_2(dt, p, std::back_inserter(coords));

  if(!result.third){
    std::cout << "The coordinate computation was not successful."
              << std::endl;
    std::cout << "The point (" <<p << ") lies outside the convex hull."
              << std::endl;
  }

  K::FT  norm = result.second;
  std::cout << "Coordinate computation successful." << std::endl;
  std::cout << "Normalization factor: " <<norm << std::endl;
  std::cout << "done" << std::endl;
  return 0;
}
