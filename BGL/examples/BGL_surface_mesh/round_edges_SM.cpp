
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor       halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor           edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor           face_descriptor;


int main(int argc, char* argv[])
{
  Mesh sm;

  std::ifstream inm(argv[1]);
  inm >> sm;
  
  std::ifstream ins(argv[2]);
  int i,j;
  std::vector<halfedge_descriptor> hedges;

  while(ins >> i >> j){
    vertex_descriptor vi(i), vj(j);
    std::pair<edge_descriptor,bool> edb = edge(vi,vj,sm);
    assert(edb.second);
    halfedge_descriptor hd = halfedge(edb.first,sm);
    assert(target(hd,sm) == vj);
    hedges.push_back(hd);
  }
  assert(sm.is_valid());
  assert(CGAL::is_triangle_mesh(sm));
  std::cout << num_vertices(sm) << std::endl;
  CGAL::Euler::round_edges(hedges,sm);

  std::cout << num_vertices(sm) << std::endl;
    
  assert(sm.is_valid());
  std::cout << "done" << std::endl;
  return 0;
}
