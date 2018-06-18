#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <boost/foreach.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

typedef CGAL::Simple_cartesian<double>                       Kernel;

typedef Kernel::FT                                           FT;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Vector_3                                     Vector;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor       halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor           edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor           face_descriptor;

typedef std::pair<halfedge_descriptor, Vector>               Vertex_with_normal;
typedef std::vector<Vertex_with_normal>                      Vertices_with_normals;

template<typename NormalMap>
Vertices_with_normals round(Mesh& tm,
                            NormalMap& vertex_normals,
                            const std::vector<halfedge_descriptor>& hedges,
                            const int layer_n = 4)
{
  std::cout << "round() with hedges size: " << hedges.size() << std::endl;

  Vertices_with_normals vertices_with_normals;

  assert(tm.is_valid());
  assert(CGAL::is_triangle_mesh(tm));

  std::size_t feature_length = hedges.size();
  if(source(hedges.front(),tm) != target(hedges.back(),tm)) {
    --feature_length;
  }

  for(std::size_t i = 0; i < feature_length; ++i)
  {
    const halfedge_descriptor hd = hedges[i];
    vertices_with_normals.push_back(std::make_pair(hd, get(vertex_normals, hd)));
    std::cout << "add : " << target(hd, tm) << std::endl;
  }

  for(int iter = 0; iter < layer_n; ++iter)
  {
    CGAL::Euler::round_edges(hedges, tm);

    // if there is more than one iteration, assign normals to the new points
    for(std::size_t i = 0; i < feature_length; ++i)
    {
      const halfedge_descriptor hd = hedges[i];
      const Vector hd_normal = get(vertex_normals, hd);

      std::size_t next_id = (i == (hedges.size()-1)) ? 0 : i + 1;
      const halfedge_descriptor next_hd = hedges[next_id];
      const halfedge_descriptor opp_hd = opposite(next_hd, tm);
      const Vector opp_hd_normal = get(vertex_normals, opp_hd);

      const FT alpha = FT(iter) / layer_n;
      Vector interpolated_normal = alpha * hd_normal + (1 - alpha) * opp_hd_normal;

      // Won't work if both vertices are opposite...
      // ... but then it shouldn't have been a feature polyline in the first place
      FT l = CGAL::sqrt(interpolated_normal.squared_length());
      assert(l != 0);
      interpolated_normal = interpolated_normal / l;

      halfedge_descriptor new_halfedge = prev(opposite(hd, tm), tm);
      new_halfedge = next(opposite(new_halfedge, tm), tm);
      new_halfedge = next(opposite(new_halfedge, tm), tm);
      new_halfedge = opposite(new_halfedge, tm);

      put(vertex_normals, new_halfedge, interpolated_normal);
      vertices_with_normals.push_back(std::make_pair(new_halfedge, interpolated_normal));
    }
  }

  // move the points
  std::ofstream out_normals("normals.xyz");

  typename Vertices_with_normals::iterator vit = vertices_with_normals.begin(),
                                           end = vertices_with_normals.end();
  for(; vit!=end; ++vit)
  {
    halfedge_descriptor hd = vit->first;
    const Point old_position = tm.point(target(hd,tm));
    tm.point(target(hd,tm)) = old_position + 0.01 * vit->second;
    std::cout << "v: " << target(hd,tm)
              << " old position : " << old_position
              << " new: " << tm.point(target(hd, tm)) << std::endl;

    out_normals << tm.point(target(hd,tm)) << " " << vit->second << std::endl;
  }

  return vertices_with_normals;
}

Vector compute_face_normal(halfedge_descriptor hd,
                           const Mesh& tm)
{
  Point p0 = tm.point(target(hd, tm));
  hd = next(hd, tm);
  Point p1 = tm.point(target(hd, tm));
  hd = next(hd, tm);
  Point p2 = tm.point(target(hd, tm));
  hd = next(hd, tm);

  return CGAL::unit_normal(p1, p2, p0);
}

template<typename Is_feature_edge_map>
boost::optional<Vector> compute_normal(halfedge_descriptor hd,
                                       const Is_feature_edge_map& isef,
                                       const Mesh& tm)
{
  assert(!CGAL::is_border_edge(hd, tm));

  // grab the faces in the sector until we meet another halfedge that is constrained
  halfedge_descriptor done = hd;
  std::list<Vector> incident_faces_normals;

  do
  {
    Vector face_normal = compute_face_normal(hd, tm);
    incident_faces_normals.push_back(face_normal);

    hd = opposite(next(hd, tm), tm);
  }
  while(!isef(edge(hd, tm)));

  // only one feature edge incident to this vertex --> don't care about the normal
  if(hd == done)
    return boost::none;

  assert(incident_faces_normals.size() > 0);

  // compute the average of the normals of these incident faces
  FT x = 0, y = 0, z = 0;

  BOOST_FOREACH(const Vector v, incident_faces_normals)
  {
    x += v.x();
    y += v.y();
    z += v.z();
  }

  std::size_t ifn_size = incident_faces_normals.size();
  return Vector(x / ifn_size, y / ifn_size, z / ifn_size);
}

template<typename Vertex_normal_map, typename Is_feature_edge_map>
void compute_normals(const Mesh& tm,
                     Vertex_normal_map& vertex_normal_map,
                     const Is_feature_edge_map& isef)
{
  std::ofstream out_normals("initial_normals.xyz");

  BOOST_FOREACH(halfedge_descriptor hd, halfedges(tm))
  {
    if(!isef(edge(hd, tm)))
      continue;

    std::cout << hd << "("
              << tm.point(source(hd, tm)) << " ][ "
              << tm.point(target(hd, tm)) << ") is constrained" << std::endl;

    boost::optional<Vector> normal = compute_normal(hd, isef, tm);
    if(normal)
    {
      put(vertex_normal_map, hd, *normal);
      out_normals << tm.point(target(hd, tm)) << " " << *normal << std::endl;
    }
  }
}

template <typename FeatureEdgeMap>
struct Is_feature
{
  Is_feature() {}
  Is_feature(FeatureEdgeMap m) : m_m(m) {}

  template <typename Edge>
  bool operator()(const Edge& e) const { return get(m_m, e); }
  FeatureEdgeMap m_m;
};

template <typename Graph, typename Base, typename NormalMap>
struct Visitor
{
  std::vector<typename boost::graph_traits<Base>::vertex_descriptor> poly;
  Graph& g;
  Base& b;
  NormalMap& normal_map;

  Visitor(Graph&g, Base& b, NormalMap& normal_map)
    : g(g), b(b), normal_map(normal_map)
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

    typename std::vector<vertex_descriptor>::iterator it = poly.begin();
    vertex_descriptor vp = *it;
    for(++it; it != poly.end(); ++it) {
      vertex_descriptor vq = *it;
      edge_descriptor ed;
      bool found = false;
      boost::tie(ed,found) = edge(vp,vq,b);
      assert(found);
      hedges.push_back(halfedge(ed,b));
      vp = vq;
    }

    Vertices_with_normals vertices_with_normals;
    vertices_with_normals = round(b, normal_map, hedges);
  }
};

int main(int argc, char* argv[])
{
  assert(argc > 2);

  Mesh tm;

  typedef Mesh::Property_map<vertex_descriptor, bool> vFeature;
  typedef Mesh::Property_map<edge_descriptor, bool> eFeature;
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
  assert(inm.good());
  inm >> tm;

  std::ifstream ins(argv[2]); // selection
  assert(ins.good());

  int i,j;
  std::vector<halfedge_descriptor> hedges;
  std::string ignore;
  std::getline(ins, ignore);
  std::getline(ins, ignore);

  while(ins >> i >> j) {
    vertex_descriptor vi(i), vj(j);
    std::cout << "constraining " << vi << " " << vj << std::endl;

    std::pair<edge_descriptor,bool> edb = edge(vi,vj,tm);
    assert(edb.second);
    halfedge_descriptor hd = halfedge(edb.first,tm);
    assert(target(hd,tm) == vj);
    hedges.push_back(hd);
    put(efm, edge(hd,tm),true);
    put(vfm, vi,true);
    put(vfm, vj,true);
  }

  typedef CGAL::dynamic_halfedge_property_t<Vector>          Vertex_normal;
  typedef boost::property_map<Mesh, Vertex_normal>::type     Vertex_normal_map;
  Vertex_normal_map vertex_normal_map = get(Vertex_normal(), tm);
  compute_normals(tm, vertex_normal_map, isef);

  Visitor<FG, Mesh, Vertex_normal_map> vis(fg, tm, vertex_normal_map);
  std::cout << num_vertices(tm) << std::endl;
  CGAL::split_graph_into_polylines(fg,vis);
  std::cout << num_vertices(tm) << std::endl;

  std::ofstream out("out.off");
  out << tm << std::endl;

  std::cout << "done" << std::endl;
  return 0;
}
