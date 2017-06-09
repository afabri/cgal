// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_COMPONENT_GRAPH_H
#define CGAL_BOOST_COMPONENT_GRAPH_H

namespace CGAL {

  namespace internal {
template <typename G>
class CG_edge_iterator
  : public boost::iterator_facade<CG_edge_iterator<G>, typename boost::graph_traits<G>::edge_descriptor, std::forward_iterator_tag> {

public:
 
  friend class boost::iterator_core_access;
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  bool end;
  boost::shared_ptr<std::vector<edge_descriptor> > edges;
  typename std::vector<edge_descriptor>::iterator it;

  CG_edge_iterator()
    : edges(NULL)
  {
    //    std::cerr << "default construct" << std::endl;
  }

  CG_edge_iterator(const CG_edge_iterator& other)
    : end(other.end), edges(other.edges), it(other.it)
  {
    // std::cerr << "copy construct" << std::endl;
  }

  // for end()
  CG_edge_iterator(int)
    : end(true),edges(NULL)
  {
    // std::cerr << "end construct" << std::endl;
}

  CG_edge_iterator& 
  operator=(const CG_edge_iterator other)
  {
    end = other.end;
    edges = other.edges;
    it = other.it;
    return *this;
  }
  
  CG_edge_iterator(const G& g)
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
  bool equal (const CG_edge_iterator& other) const
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

  } // namespace internal

template <typename G, typename ECMap>
class Component_graph {
public:
  typedef internal::CG_edge_iterator<Component_graph<G,ECMap> > edge_iterator;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  G& g;
  halfedge_descriptor cc;
  ECMap ecmap;
  int ne;

  Component_graph(G& g, ECMap ecmap, halfedge_descriptor cc)
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

} // namespace CGAL


namespace boost {

template <typename G, typename ECMap>
struct graph_traits<CGAL::Component_graph<G,ECMap> > {
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
  typedef typename CGAL::Component_graph<G,ECMap>::edge_iterator edge_iterator;

  static face_descriptor null_face()
  {
    return boost::graph_traits<G>::null_face();
  }
};
} // namespace boost

namespace CGAL {

template <typename G, typename ECMap>
CGAL::Iterator_range<internal::CG_edge_iterator<Component_graph<G,ECMap> > >
edges(const Component_graph<G,ECMap>& cg)
{
  return CGAL::make_range(internal::CG_edge_iterator<Component_graph<G,ECMap> >(cg),
                          internal::CG_edge_iterator<Component_graph<G,ECMap> >(0));
}

template <typename G, typename ECMap>
CGAL::Iterator_range<typename boost::graph_traits<G>::halfedge_iterator>
halfedges(const Component_graph<G,ECMap>& cg)
{
  return halfedges(cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::face_descriptor
face(typename boost::graph_traits<G>::halfedge_descriptor h,
     const Component_graph<G,ECMap>& cg)
{
  return face(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
opposite(typename boost::graph_traits<G>::halfedge_descriptor h,
     const Component_graph<G,ECMap>& cg)
{
  return opposite(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::edges_size_type
num_edges(const Component_graph<G,ECMap>& cg)
{
  return cg.num_edges();
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
halfedge(typename boost::graph_traits<G>::edge_descriptor e,
     const Component_graph<G,ECMap>& cg)
{
  return halfedge(e, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
halfedge(typename boost::graph_traits<G>::vertex_descriptor v,
     const Component_graph<G,ECMap>& cg)
{
  return halfedge(v, cg.g);
}

template <typename G, typename ECMap>
std::pair<typename boost::graph_traits<G>::halfedge_descriptor,bool>
halfedge(typename boost::graph_traits<G>::vertex_descriptor v,
         typename boost::graph_traits<G>::vertex_descriptor w,
         const Component_graph<G,ECMap>& cg)
{
  return halfedge(v, w, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::edge_descriptor
edge(typename boost::graph_traits<G>::halfedge_descriptor h,
     const Component_graph<G,ECMap>& cg)
{
  return edge(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::vertex_descriptor
target(typename boost::graph_traits<G>::halfedge_descriptor h,
     const Component_graph<G,ECMap>& cg)
{
  return target(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::vertex_descriptor
source(typename boost::graph_traits<G>::halfedge_descriptor h,
     const Component_graph<G,ECMap>& cg)
{
  return source(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::degree_size_type
degree(typename boost::graph_traits<G>::vertex_descriptor v,
     const Component_graph<G,ECMap>& cg)
{
  return degree(v, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
prev(typename boost::graph_traits<G>::halfedge_descriptor h,
     const Component_graph<G,ECMap>& cg)
{
  return prev(h, cg.g);
}

template <typename G, typename ECMap>
typename boost::graph_traits<G>::halfedge_descriptor
next(typename boost::graph_traits<G>::halfedge_descriptor h,
     const Component_graph<G,ECMap>& cg)
{
  return next(h, cg.g);
}

template <typename G, typename ECMap>
void
remove_vertex(typename boost::graph_traits<G>::vertex_descriptor v,
     Component_graph<G,ECMap>& cg)
{
  remove_vertex(v, cg.g);
}

template <typename G, typename ECMap>
void
remove_edge(typename boost::graph_traits<G>::edge_descriptor e,
     Component_graph<G,ECMap>& cg)
{
  return remove_edge(e, cg.g);
}

template <typename G, typename ECMap>
void
remove_face(typename boost::graph_traits<G>::face_descriptor f,
     Component_graph<G,ECMap>& cg)
{
  return remove_face(f, cg.g);
}

template <typename G, typename ECMap>
void
set_target(typename boost::graph_traits<G>::halfedge_descriptor h,
           typename boost::graph_traits<G>::vertex_descriptor v,
           const Component_graph<G,ECMap>& cg)
{
  return set_target(h,v, cg.g);
}

template <typename G, typename ECMap>
void
set_next(typename boost::graph_traits<G>::halfedge_descriptor h,
           typename boost::graph_traits<G>::halfedge_descriptor n,
           const Component_graph<G,ECMap>& cg)
{
  return set_next(h,n, cg.g);
}

template <typename G, typename ECMap>
void
set_halfedge(typename boost::graph_traits<G>::face_descriptor f,
             typename boost::graph_traits<G>::halfedge_descriptor h,
             const Component_graph<G,ECMap>& cg)
{
  return set_halfedge(f,h, cg.g);
}
template <typename G, typename ECMap>
void
set_halfedge(typename boost::graph_traits<G>::vertex_descriptor v,
             typename boost::graph_traits<G>::halfedge_descriptor h,
             const Component_graph<G,ECMap>& cg)
{
  return set_halfedge(v,h, cg.g);
}

template <typename G, typename ECMap>
void
set_face(typename boost::graph_traits<G>::halfedge_descriptor h,
           typename boost::graph_traits<G>::face_descriptor f,
           const Component_graph<G,ECMap>& cg)
{
  return set_face(h,f, cg.g);
}

} // namespace CGAL

namespace boost {
template <typename G, typename ECMap>
struct property_map<CGAL::Component_graph<G,ECMap>,CGAL::vertex_point_t>
{
  typedef typename boost::property_map<G,CGAL::vertex_point_t>::type type;
};

} // namespace boost

namespace CGAL {

template <typename G, typename ECMap>
typename boost::property_map<Component_graph<G,ECMap>,CGAL::vertex_point_t>::type
get(CGAL::vertex_point_t t, Component_graph<G,ECMap>& cg)
{
  return get(t,cg.g);
}

} // namespace CGAL

#endif CGAL_BOOST_COMPONENT_GRAPH_H
