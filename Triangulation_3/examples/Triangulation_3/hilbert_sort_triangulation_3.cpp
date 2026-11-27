#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <unordered_map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_3<K>      Triangulation;

typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;


typedef std::pair<Vertex_handle, int>                       Vertex_index_pair;

// Readable property map: pair<Vertex_handle,int>  ->  vertex point
struct Vertex_pair_point_map
{
  typedef Vertex_index_pair                  key_type;
  typedef Point                              value_type;
  typedef const Point&                       reference;
  typedef boost::readable_property_map_tag   category;

  friend reference get(const Vertex_pair_point_map&, const key_type& p)
  {
    return p.first->point();
  }
};

typedef CGAL::Spatial_sort_traits_adapter_3<K, Vertex_pair_point_map> Hilbert_traits;




void reorder(Triangulation& T)
{
  std::unordered_map<Vertex_handle,int> vertex_index;
  std::vector<std::pair<Vertex_handle,int>> vertices;
  vertices.reserve(T.number_of_vertices());
  int index = 0;
  for (auto v = T.finite_vertices_begin(); v != T.finite_vertices_end(); ++v){
    vertex_index[v] = index;
    vertices.push_back(std::make_pair(v, index++));
  }

  CGAL::hilbert_sort(vertices.begin(), vertices.end(),
                     Hilbert_traits(Vertex_pair_point_map()));

}

int main()
{
  // construction from a list of points :
  std::list<Point> L;
  L.push_front(Point(0,0,0));
  L.push_front(Point(1,0,0));
  L.push_front(Point(0,1,0));

  Triangulation T(L.begin(), L.end());

  Triangulation::size_type n = T.number_of_vertices();

  // insertion from a vector :
  std::vector<Point> V(3);
  V[0] = Point(0,0,1);
  V[1] = Point(1,1,1);
  V[2] = Point(2,2,2);

  n = n + T.insert(V.begin(), V.end());

  for(auto v : T.finite_vertex_handles()){
    std::cout << v->point() << std::endl;
  }
  T.reorder_vertices();
  std::cout << "vertices reordered" << std::endl;
  for(auto v : T.finite_vertex_handles()){
    std::cout << v->point() << std::endl;
  }


  assert( T.is_valid() ); // checking validity of T

  return 0;
}
