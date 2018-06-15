#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Vector_3                                     Vector;
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

  std::size_t offset = hedges.size();
  if(source(hedges.front(),tm) != target(hedges.back(),tm)) {
    --offset;
  }

  std::vector<std::pair<vertex_descriptor,face_descriptor> > res(2*(offset));
  std::vector<halfedge_descriptor> opp_hedges;
  
  for(std::size_t i = 0; i < offset; ++i)
  {
    res[i] = std::make_pair(target(hedges[i],tm), face(hedges[i],tm));
  }

  if(hedges.size() <= 2) {
    res.clear();
    return res;
  }

  CGAL::Euler::round_edges(hedges,tm);

  for(std::size_t i = 0; i < offset; ++i) {
    halfedge_descriptor hd = hedges[i];
    std::cout << "hd  " << hd << " " << source(hd,tm) << " " << target(hd,tm) << std::endl;
    hd = prev(opposite(hd,tm),tm);
    hd = next(opposite(hd,tm),tm);
    hd = opposite(hd,tm);
    // res[i+offset] = std::make_pair(target(hd,tm),face(hd,tm));
    opp_hedges.push_back(hd);
    std::cout << "opp " << hd << " " << source(hd,tm) << " " << target(hd,tm) << std::endl;
  }

  std::cout << opp_hedges.size() << std::endl;
  std::vector<Point> newpoints;
  BOOST_FOREACH(halfedge_descriptor hd, hedges) {
    Vector v = CGAL::unit_normal(tm.point(source(hd,tm)),
                                 tm.point(target(hd,tm)),
                                 tm.point(target(next(hd,tm),tm)));
    newpoints.push_back(tm.point(target(hd,tm)) + 0.01 * v);
  }

  int i = 0;
  BOOST_FOREACH(halfedge_descriptor hd, hedges) {
    tm.point(target(hd,tm))= newpoints[i++];
  }

  newpoints.clear();
  BOOST_FOREACH(halfedge_descriptor hd, opp_hedges) {
    std::cout << "A" << std::endl;
    std::cout << tm.point(source(hd,tm)) << "|"<< 
                                 tm.point(target(hd,tm)) << "|"<< 
      tm.point(target(next(hd,tm),tm)) << std::endl;
    Vector v = CGAL::unit_normal(tm.point(source(hd,tm)),
                                 tm.point(target(hd,tm)),
                                 tm.point(target(next(hd,tm),tm)));
    newpoints.push_back(tm.point(target(hd,tm)) + 0.01 * v);
  }

  std::cout << opp_hedges.size() << std::endl;

  i = 0;
  BOOST_FOREACH(halfedge_descriptor hd, opp_hedges) {
        std::cout << "set " << hd << " " << source(hd,tm) << " " << target(hd,tm) << std::endl;
    tm.point(target(hd,tm))= newpoints[i++];
  }

  return res;
}

template <typename FeatureEdgeMap>
struct Is_feature {
  Is_feature()
  {}
  
  Is_feature(FeatureEdgeMap m)
    : m_m(m)
  {}
  
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return get(m_m, e);
  }
  FeatureEdgeMap m_m;
};


template <typename Graph, typename Base>
struct Visitor {
  std::vector<typename boost::graph_traits<Base>::vertex_descriptor> poly;
  Graph& g;
  Base& b;

  Visitor(Graph&g, Base& b)
    : g(g),b(b)
  {}
  
    void start_new_polyline()
  {
    std::cout << "start polyline" << std::endl;
    poly.clear();
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor v)
  {
    poly.push_back(v);
  }
  
    void end_polyline()
  {
    std::cout << "end polyline" << std::endl;
    std::vector<typename boost::graph_traits<Base>::halfedge_descriptor> hedges;

    typename std::vector<typename boost::graph_traits<Base>::vertex_descriptor>::iterator it = poly.begin();
    typename boost::graph_traits<Base>::vertex_descriptor vp = *it;
    for(++it; it != poly.end(); ++it) {
      typename boost::graph_traits<Base>::vertex_descriptor vq = *it;
      typename boost::graph_traits<Base>::edge_descriptor ed;
      bool found = false;
      boost::tie(ed,found) = edge(vp,vq,b);
      assert(found);
      hedges.push_back(halfedge(ed,b));
      vp = vq;
    }
    
    std::vector<std::pair<typename boost::graph_traits<Base>::vertex_descriptor,
                          typename boost::graph_traits<Base>::face_descriptor> >  vvpairs;
    vvpairs = round(b, hedges);

  }
};

int main(int argc, char* argv[])
{
  Mesh tm;

  typedef Mesh::Property_map<vertex_descriptor,bool> vFeature;
  typedef Mesh::Property_map<edge_descriptor,bool> eFeature;
  vFeature vfm;
  eFeature efm;
  bool created;
  boost::tie(vfm, created) = tm.add_property_map<vertex_descriptor,bool>("v:feature",false);
  boost::tie(efm, created) = tm.add_property_map<edge_descriptor,bool>("e:feature",false);

  Is_feature<vFeature> isvf(vfm);
  Is_feature<eFeature> isef(efm);

  typedef boost::filtered_graph<Mesh,Is_feature<eFeature>,Is_feature<vFeature> > FG;

  FG fg(tm,isef,isvf);

  std::ifstream inm(argv[1]); // mesh
  inm >> tm;

  std::ifstream ins(argv[2]); // selection
  int i,j;
  std::vector<halfedge_descriptor> hedges;
  std::string ignore;
  std::getline(ins,ignore);
  std::getline(ins,ignore);
  while(ins >> i >> j) {
    vertex_descriptor vi(i), vj(j);
    std::pair<edge_descriptor,bool> edb = edge(vi,vj,tm);
    assert(edb.second);
    halfedge_descriptor hd = halfedge(edb.first,tm);
    assert(target(hd,tm) == vj);
    hedges.push_back(hd);
    put(efm, edge(hd,tm),true);
    put(vfm, vi,true);
    put(vfm, vj,true);
  }

  Visitor<FG,Mesh> vis(fg,tm);
  std::cout << num_vertices(tm) << std::endl;
  CGAL::split_graph_into_polylines(fg,vis);
  std::cout << num_vertices(tm) << std::endl;

  std::ofstream out("out.off");
  out << tm << std::endl;
  std::cout << "done" << std::endl;
  return 0;
}
