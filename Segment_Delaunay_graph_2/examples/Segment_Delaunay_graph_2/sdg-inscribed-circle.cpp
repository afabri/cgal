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


int main(int argc, char* argv[])
{
  std::ifstream ifs((argc>1)? argv[1]:"data/polygon.wkt");
  assert( ifs );

  Polygon_with_holes pwh;
  CGAL::IO::read_polygon_WKT(ifs, pwh);

  const auto boundary_polygon = pwh.outer_boundary();

  SDG2          sdg;

  sdg.insert_segments(boundary_polygon.edges_begin(), boundary_polygon.edges_end());

  double sd = 0, sqdist = 0;
  SDG2::Finite_faces_iterator fit = sdg.finite_faces_begin();
  for (; fit != sdg.finite_faces_end(); ++fit) {
      Point_2 pp = sdg.primal(fit);
      std::cout << "Point " << pp << std::endl;
      std::cout << "equidistant to:" << std::endl;
      for (int i = 0; i < 3; ++i) {
          assert(!sdg.is_infinite(fit->vertex(i)));
          if (fit->vertex(i)->site().is_segment()) {
              Segment_2 s = fit->vertex(i)->site().segment();
              double sqdist = CGAL::squared_distance(pp, s);
              std::cout << "  segment " << s << " at distance " << sqrt(sqdist) << std::endl;
          }
          else {
              Point_2 p = fit->vertex(i)->site().point();
              double sqdist = CGAL::squared_distance(pp, p);
              std::cout << "  point " << p << " at distance " << sqrt(sqdist) << std::endl;
          }
      }
      if (boundary_polygon.bounded_side(pp) == CGAL::ON_BOUNDED_SIDE) {
          sd = std::max(sqdist, sd);
      }
      else {
          std::cout << "  but ignored as not on bounded side" << std::endl;
      }
  }
  std::cout << "radius of largest inscribed circle: " << sqrt(sd) << std::endl;
  return 0;
}
