//#include "tbb/task_group.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;

typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh; 


namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

template <typename G>
class EdgeIterator
  : public boost::iterator_facade<EdgeIterator<G>, typename boost::graph_traits<G>::edge_descriptor, std::forward_iterator_tag> {

public:
 

  friend class boost::iterator_core_access;
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  bool end;
  boost::shared_ptr<std::vector<edge_descriptor> > edges;
  typename std::vector<edge_descriptor>::iterator it;

  EdgeIterator()
    : edges(NULL)
  {
    //    std::cerr << "default construct" << std::endl;
  }

  EdgeIterator(const EdgeIterator& other)
    : end(other.end), edges(other.edges), it(other.it)
  {
    // std::cerr << "copy construct" << std::endl;
  }

  // for end()
  EdgeIterator(int)
    : end(true),edges(NULL)
  {
    // std::cerr << "end construct" << std::endl;
}

  EdgeIterator& 
  operator=(const EdgeIterator other)
  {
    end = other.end;
    edges = other.edges;
    it = other.it;
    return *this;
  }
  
  EdgeIterator(const G& g)
    : end(false)
  {
    //    std::cerr << "G construct" << std::endl;
    edges = boost::shared_ptr<std::vector<edge_descriptor> >(new std::vector<edge_descriptor>());
    typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
    // TODO: replace faces by a functor that writes directly into edges
    //       Or write the iterator so that increment does the work without the vector edges
    std::vector<face_descriptor> faces;
    PMP::connected_component(face(g.cc,g.g),
                             g.g,
                             std::back_inserter(faces),
                             PMP::parameters::edge_is_constrained_map(g.ecmap));

    BOOST_FOREACH(face_descriptor fd, faces){
      BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd,g.g),g.g)){
       halfedge_descriptor hop = opposite(hd,g.g);
       if(is_border(hop,g.g) || get(g.ecmap,edge(hd,g.g)) || (fd < face(hop ,g.g))){
         edges->push_back(edge((hd<hop)?hd:hop, g.g));
        }
      }
    }
    it = edges->begin();
  }
  

  void increment()
  {
    ++it;
  }
  
  // so far only other may be end
  bool equal (const EdgeIterator& other) const
  {
    if(end) std::cerr << "oops"<< std::endl;
    if(other.end){
      return it == edges->end();
    }
    return it == other.it;
  }
  
  edge_descriptor& dereference() const
  {
    return const_cast<edge_descriptor&>(*it);
  }
};


template <typename G, typename ECMap>
class ComponentGraph {
public:
  typedef EdgeIterator<ComponentGraph<G,ECMap> > edge_iterator;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  G& g;
  halfedge_descriptor cc;
  ECMap ecmap;
  int ne;

  ComponentGraph(G& g, ECMap ecmap, halfedge_descriptor cc)
    : g(g), cc(cc), ecmap(ecmap)
  {}
  
  int& num_edges()
  {
    return ne;
  }

  int num_edges() const
  {
    return ne;
  }
}; 


namespace boost {

template <typename G, typename ECMap>
struct graph_traits<ComponentGraph<G,ECMap> > {
  typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<G>::traversal_category traversal_category;
  typedef typename boost::graph_traits<G>::edge_parallel_category edge_parallel_category;
  typedef typename boost::graph_traits<G>::directed_category directed_category;

  typedef typename boost::graph_traits<G>::edges_size_type edges_size_type;
  typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator; // in SMS only for tracing
  typedef typename boost::graph_traits<G>::halfedge_iterator halfedge_iterator; // in order to set bool m_has_border
  typedef typename ComponentGraph<G,ECMap>::edge_iterator edge_iterator;

  static face_descriptor null_face()
  {
    return boost::graph_traits<G>::null_face();
  }
};
} // namespace boost


template <typename G, typename ECMap>
CGAL::Iterator_range<EdgeIterator<ComponentGraph<G,ECMap> > >
edges(const ComponentGraph<G,ECMap>& cg)
{
  return CGAL::make_range(EdgeIterator<ComponentGraph<G,ECMap> >(cg),
                          EdgeIterator<ComponentGraph<G,ECMap> >(0));
}

template <typename G, typename ECMap>
CGAL::Iterator_range<typename boost::graph_traits<G>::halfedge_iterator>
halfedges(const ComponentGraph<G,ECMap>& cg)
{
  return halfedges(cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::face_descriptor
face(typename boost::graph_traits<G>::halfedge_descriptor h,
     const ComponentGraph<G,ECMap>& cg)
{
  return face(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
opposite(typename boost::graph_traits<G>::halfedge_descriptor h,
     const ComponentGraph<G,ECMap>& cg)
{
  return opposite(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::edges_size_type
num_edges(const ComponentGraph<G,ECMap>& cg)
{
  return cg.num_edges();
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
halfedge(typename boost::graph_traits<G>::edge_descriptor e,
     const ComponentGraph<G,ECMap>& cg)
{
  return halfedge(e, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
halfedge(typename boost::graph_traits<G>::vertex_descriptor v,
     const ComponentGraph<G,ECMap>& cg)
{
  return halfedge(v, cg.g);
}

template <typename G, typename ECMap>
std::pair<typename boost::graph_traits<G>::halfedge_descriptor,bool>
halfedge(typename boost::graph_traits<G>::vertex_descriptor v,
         typename boost::graph_traits<G>::vertex_descriptor w,
         const ComponentGraph<G,ECMap>& cg)
{
  return halfedge(v, w, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::edge_descriptor
edge(typename boost::graph_traits<G>::halfedge_descriptor h,
     const ComponentGraph<G,ECMap>& cg)
{
  return edge(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::vertex_descriptor
target(typename boost::graph_traits<G>::halfedge_descriptor h,
     const ComponentGraph<G,ECMap>& cg)
{
  return target(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::vertex_descriptor
source(typename boost::graph_traits<G>::halfedge_descriptor h,
     const ComponentGraph<G,ECMap>& cg)
{
  return source(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::degree_size_type
degree(typename boost::graph_traits<G>::vertex_descriptor v,
     const ComponentGraph<G,ECMap>& cg)
{
  return degree(v, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
prev(typename boost::graph_traits<G>::halfedge_descriptor h,
     const ComponentGraph<G,ECMap>& cg)
{
  return prev(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
next(typename boost::graph_traits<G>::halfedge_descriptor h,
     const ComponentGraph<G,ECMap>& cg)
{
  return next(h, cg.g);
}

template <typename G, typename ECMap>
void
remove_vertex(typename boost::graph_traits<G>::vertex_descriptor v,
     ComponentGraph<G,ECMap>& cg)
{
  remove_vertex(v, cg.g);
}

template <typename G, typename ECMap>
void
remove_edge(typename boost::graph_traits<G>::edge_descriptor e,
     ComponentGraph<G,ECMap>& cg)
{
  return remove_edge(e, cg.g);
}

template <typename G, typename ECMap>
void
remove_face(typename boost::graph_traits<G>::face_descriptor f,
     ComponentGraph<G,ECMap>& cg)
{
  return remove_face(f, cg.g);
}

template <typename G, typename ECMap>
void
set_target(typename boost::graph_traits<G>::halfedge_descriptor h,
           typename boost::graph_traits<G>::vertex_descriptor v,
           const ComponentGraph<G,ECMap>& cg)
{
  return set_target(h,v, cg.g);
}

template <typename G, typename ECMap>
void
set_next(typename boost::graph_traits<G>::halfedge_descriptor h,
           typename boost::graph_traits<G>::halfedge_descriptor n,
           const ComponentGraph<G,ECMap>& cg)
{
  return set_next(h,n, cg.g);
}

template <typename G, typename ECMap>
void
set_halfedge(typename boost::graph_traits<G>::face_descriptor f,
             typename boost::graph_traits<G>::halfedge_descriptor h,
             const ComponentGraph<G,ECMap>& cg)
{
  return set_halfedge(f,h, cg.g);
}
template <typename G, typename ECMap>
void
set_halfedge(typename boost::graph_traits<G>::vertex_descriptor v,
             typename boost::graph_traits<G>::halfedge_descriptor h,
             const ComponentGraph<G,ECMap>& cg)
{
  return set_halfedge(v,h, cg.g);
}

template <typename G, typename ECMap>
void
set_face(typename boost::graph_traits<G>::halfedge_descriptor h,
           typename boost::graph_traits<G>::face_descriptor f,
           const ComponentGraph<G,ECMap>& cg)
{
  return set_face(h,f, cg.g);
}

namespace boost {
template <typename G, typename ECMap>
struct property_map<ComponentGraph<G,ECMap>,CGAL::vertex_point_t>
{
  typedef typename boost::property_map<G,CGAL::vertex_point_t>::type type;
};
}

template <typename G, typename ECMap>
typename boost::property_map<ComponentGraph<G,ECMap>,CGAL::vertex_point_t>::type
get(CGAL::vertex_point_t t, ComponentGraph<G,ECMap>& cg)
{
  return get(t,cg.g);
}



struct Simplify {

  Surface_mesh &sm;

 Simplify(Surface_mesh& sm)
    : sm(sm)
  {}
};

template <typename ZOMap, typename G>
struct Selection_is_constraint_map
{  
  typedef boost::readable_property_map_tag                 category;
  typedef bool                                             value_type;
  typedef bool                                             reference;
  typedef typename boost::graph_traits<G>::edge_descriptor key_type;

  ZOMap zomap;
  const G& g;

  Selection_is_constraint_map(ZOMap zomap, const G& g)
    : zomap(zomap), g(g)
  {}



  friend bool get(const Selection_is_constraint_map& eicm, key_type ed)
  {
    int i = get(eicm.zomap, halfedge(ed, eicm.g));
    return (i ==0) || (i == 1);
  }
  
};


int main(int argc, char** argv ) 
{
  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
 
  std::ifstream is1(argv[1]);
  Surface_mesh sm; 
  is1 >> sm;
  std::cerr << "#V = "<< num_vertices(sm)  << "#E = "<< num_edges(sm) << "  #F = " << num_faces(sm) << std::endl;

  

  typedef Surface_mesh::Property_map<face_descriptor,int> CCMap;
  Surface_mesh::Property_map<face_descriptor,int> ccmap 
    = sm.add_property_map<face_descriptor,int>("f:cc").first;

  typedef Surface_mesh::Property_map<edge_descriptor,bool> ECMap;
  Surface_mesh::Property_map<edge_descriptor,bool> ecmap 
    = sm.add_property_map<edge_descriptor,bool>("e:constrained",false).first;

  typedef Surface_mesh::Property_map<halfedge_descriptor,int> HIMap;
  HIMap himap = sm.add_property_map<halfedge_descriptor,int>("h:index_in_cc",-1).first;

  typedef Selection_is_constraint_map<HIMap, Surface_mesh> SICM;
  SICM sicm(himap,sm);
  
  
  std::ifstream is2(argv[2]);
  std::string line;
  std::size_t id, id2;
  
  if(!std::getline(is2, line)){
    std::cerr << "error in selection: no first line"<< std::endl;
    return 1; 
  }
  if(!std::getline(is2, line)){
    std::cerr << "error in selection: no second line"<< std::endl;
    return 1; 
  }
  if(!std::getline(is2, line)){
    std::cerr << "error in selection: no third line"<< std::endl;
    return 1; 
  }
  std::istringstream edge_line(line);
  while(edge_line >> id >> id2) {
    vertex_descriptor s(id), t(id2);
     halfedge_descriptor hd;
     bool exists;
     boost::tie(hd,exists) = halfedge(s,t,sm);
     if(! exists){
       std::cerr << "error in selection: no edge" << s << " " << t << std::endl;
       return 1; 
     }
     ecmap[edge(hd,sm)] = true;
     himap[hd] = 0;
     himap[opposite(hd,sm)] = 1;
  }
  
  int ncc = PMP::connected_components(sm,
                                      ccmap,
                                      PMP::parameters::edge_is_constrained_map(ecmap));
  std::cout << "#CC = " << ncc << std::endl;


  // For each connected component find one halfedge opposite to a border halfedge
  // as a starting point for each simplification thread
  std::vector<halfedge_descriptor> cc_seed(ncc);
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(sm)){
    int hcc = ccmap[face(hd,sm)];
    int hoppcc = ccmap[face(opposite(hd,sm),sm)];
    if(hcc != hoppcc){
      cc_seed[hcc] = hd;
    }
  }
  typedef ComponentGraph<Surface_mesh,SICM>  Component_graph;
  Component_graph cg(sm, sicm, cc_seed[0]); 

  std::vector<boost::graph_traits<Surface_mesh>::edge_descriptor> cc_edges, expansion;
  int i = 2, count = 0;
  BOOST_FOREACH(boost::graph_traits<Surface_mesh>::edge_descriptor ed, edges(cg)){
    ++count;
    boost::graph_traits<Surface_mesh>::halfedge_descriptor hd, hop;
    hd = halfedge(ed,sm);
    hop = opposite(hd,sm);
    if(ecmap[ed]){
      cc_edges.push_back(ed);
    }else{
      put(himap,hd,i); ++i;
      put(himap,hop,i); ++i;
    }
  }
  std::cerr << cc_edges.size() << "  " << count << std::endl;
  cg.num_edges() = count;

  CGAL::expand_edge_selection(cc_edges,
                              sm,
                              2,
                              ecmap,
                              CGAL::Emptyset_iterator());
  
  std::cerr << "before edge_collapse " << expansion.size()<< std::endl;
  SMS::Count_ratio_stop_predicate<Component_graph> stop(0.25);
  SMS::edge_collapse(cg,
                     stop
                     ,CGAL::parameters::vertex_index_map(get(boost::vertex_index,sm))
                     .halfedge_index_map(himap)
                     .edge_is_constrained_map(ecmap)
                     .vertex_point_map(get(CGAL::vertex_point,sm)));

  sm.collect_garbage();
  std::cerr << "#V = "<< num_vertices(sm)  << "#E = "<< num_edges(sm) << "  #F = " << num_faces(sm) << std::endl;

  std::ofstream out("out.off");
  out << sm << std::endl;
  out.close();


  return EXIT_SUCCESS;
}
