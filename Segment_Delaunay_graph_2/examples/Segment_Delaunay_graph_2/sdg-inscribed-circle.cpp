// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

// define the kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Line_2 Line_2;
typedef Kernel::Point_2 Point_2;

#include <CGAL/Polygon_2.h>
#include <CGAL/IO/WKT.h>

typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>

typedef CGAL::Segment_Delaunay_graph_traits_2<Kernel>  Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>             SDG2;


int main()
{
  std::ifstream ifs("data/polygon.wkt");
  assert( ifs );

  Polygon_with_holes pwh;
  CGAL::IO::read_polygon_WKT(ifs, pwh);

  const auto boundary_polygon = pwh.outer_boundary();

  SDG2          sdg;

  sdg.insert_segments(boundary_polygon.edges_begin(), boundary_polygon.edges_end());

  double sd = 0;
  SDG2::Finite_faces_iterator fit = sdg.finite_faces_begin();
  for(; fit !=  sdg.finite_faces_end(); ++fit){
    if(fit->vertex(0)->site().is_segment() && fit->vertex(1)->site().is_segment()  && fit->vertex(2)->site().is_segment()){
      Segment_2 s0 = fit->vertex(0)->site().segment();

      
      Point_2 p = sdg.primal(fit);
      std::cout << p << std::endl;
      if (boundary_polygon.bounded_side(p) == CGAL::ON_BOUNDED_SIDE) {
          double point_segment_sqdist = CGAL::squared_distance(p, s0);
          sd = std::max(point_segment_sqdist, sd);
      }
    }
  }
  std::cout << "radius of largest inscribed circle: " << sqrt(sd) << std::endl;
  return 0;
}
