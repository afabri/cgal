
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

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


template <typename TM>
std::vector<std::pair<typename boost::graph_traits<TM>::vertex_descriptor,
                      typename boost::graph_traits<TM>::face_descriptor> >
round(TM& tm, const std::vector<typename boost::graph_traits<TM>::halfedge_descriptor>& hedges)
{
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;

  assert(tm.is_valid());
  assert(CGAL::is_triangle_mesh(tm));

  std::vector<std::pair<vertex_descriptor,face_descriptor> > res(2*(hedges.size()-1));

  for(std::size_t i = 0; i < hedges.size()-1; ++i)
  {
    res[i] = std::make_pair(target(hedges[i],tm), face(hedges[i],tm));
  }
  CGAL::Euler::round_edges(hedges,tm);
  std::size_t offset = hedges.size() - 1;
  for(std::size_t i = 0; i < offset; ++i){
    halfedge_descriptor hd = hedges[i];
    hd = prev(opposite(hd,tm),tm);
    hd = next(opposite(hd,tm),tm);
    hd = opposite(hd,tm);
    res[i+offset] = std::make_pair(target(hd,tm),face(hd,tm));
  }

  return res;
}


int main(int argc, char* argv[])
{
  Mesh tm;

  std::ifstream inm(argv[1]);
  inm >> tm;
  
  std::ifstream ins(argv[2]);
  int i,j;
  std::vector<halfedge_descriptor> hedges;

  while(ins >> i >> j){
    vertex_descriptor vi(i), vj(j);
    std::pair<edge_descriptor,bool> edb = edge(vi,vj,tm);
    assert(edb.second);
    halfedge_descriptor hd = halfedge(edb.first,tm);
    assert(target(hd,tm) == vj);
    hedges.push_back(hd);
  }

  std::vector<std::pair<boost::graph_traits<Mesh>::vertex_descriptor,
                        boost::graph_traits<Mesh>::face_descriptor> >  vvpairs;
  vvpairs = round<Mesh,Kernel::Vector_3>(tm, hedges);
  assert(tm.is_valid());
  std::cout << "done" << std::endl;
  return 0;
}
