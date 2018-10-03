#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/natural_neighbor_coordinates_3.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

typedef double NT; //Number Type

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;

typedef K::Point_3                                             Point_3;
typedef K::Sphere_3                                            Sphere_3;

typedef CGAL::Point_set_3<Point_3> Point_set;

typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Dh;

typedef Dh::Facet                                              Facet;
typedef Dh::Vertex_handle                                      Vertex_handle;
typedef Dh::Cell_handle                                        Cell_handle;
typedef Dh::Finite_vertices_iterator                           Finite_vertices_iterator;
typedef Dh::Vertex_iterator                                    Vertex_iterator;
typedef Dh::Facet_circulator                                   Facet_circulator;
typedef Dh::Cell_iterator                                      Cell_iterator;


template <typename V>
struct Value_function
{
  typedef V                                                    argument_type;
  typedef std::pair<double, bool>                              result_type;

  result_type operator()(const argument_type& a) const {
    return result_type(CGAL::sqrt(CGAL::squared_distance(Point_3(0,0,0), a->point())), true);
  }
};


int main()
{
  Point_set input;
  Point_set::Property_map<unsigned char> ired
    = input.add_property_map<unsigned char>("red", 0).first;
  Point_set::Property_map<unsigned char> igreen
    = input.add_property_map<unsigned char>("green", 0).first;
  Point_set::Property_map<unsigned char> iblue
    = input.add_property_map<unsigned char>("blue", 0).first;
  
  Point_set output;
  Point_set::Property_map<unsigned char> ored
    = output.add_property_map<unsigned char>("red", 0).first;
  Point_set::Property_map<unsigned char> ogreen
    = output.add_property_map<unsigned char>("green", 0).first;
  Point_set::Property_map<unsigned char> oblue
    = output.add_property_map<unsigned char>("blue", 0).first;

  Dh triangulation;

  CGAL::Random rng(0);
  CGAL::Random_points_in_sphere_3<Point_3> rpg(1.0,rng);

  int N = 10000;

  input.reserve(N);
  Point_3 p;
  for(int i=0; i < N; i++){
    Point_3 p = *rpg++;
    Point_set::iterator it = input.insert(p);
#ifdef CGAL_CHECK_DISTANCE
    ired[*it] = CGAL::sqrt(CGAL::squared_distance(Point_3(0,0,0), p)) *200;
#endif    
  }
  triangulation.insert(input.points().begin(), input.points().end());

  std::cerr << "Start interpolation" << std::endl;
  Value_function<Vertex_handle> value_function;
  for(int i=0; i < N; i++){
    Point_3 q = Point_3(0,0,0); // *rpg++;

    typedef std::vector<std::pair<Vertex_handle, double> > Coords;
    Coords coords;
    CGAL::Triple<std::back_insert_iterator<Coords>, double, bool> result
      = CGAL::natural_neighbor_coordinates_3(triangulation, q, std::back_inserter(coords));

    return 0;
    
    if(result.third){
      double norm = result.second;
      Point_set::iterator it = output.insert(q);
      double interpolation = CGAL::linear_interpolation(coords.begin(), coords.end(),
                                                        norm, value_function);
      ored[*it] = interpolation* 200;
#ifdef CGAL_CHECK_DISTANCE
      double qdist =  CGAL::sqrt(CGAL::squared_distance(Point_3(0,0,0),q));
      if(qdist < 0.5 * interpolation){
        std::cout << "q = " << q << "  at distance " << qdist << "  and with interpolation " << interpolation << std::endl;
      }
#endif      
    }else{
      // std::cout << "We cannot interpolate " << q << std::endl;
    }
  }

  /*
  {
    std::ofstream ofs("input.ply", std::ios::binary);
    CGAL::set_binary_mode(ofs);
    ofs << input;
  }
  {
    std::ofstream ofs("output.ply", std::ios::binary);
    CGAL::set_binary_mode(ofs);
    ofs << output;
  }
  */
  
  return 0;
}
