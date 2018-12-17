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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_INTERNAL_COMPONENT_GRAPH_H
#define CGAL_BOOST_GRAPH_INTERNAL_COMPONENT_GRAPH_H

#include <tbb/concurrent_vector.h>

namespace CGAL {
namespace internal {

template <typename G_, typename ECMap>
class Component_graph {
public:
  typedef G_ G;

  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<G>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<G>::face_descriptor face_descriptor;
  typedef typename tbb::concurrent_vector<edge_descriptor>::iterator edge_iterator;


  G& g;
  tbb::concurrent_vector<edge_descriptor>& cc;
  ECMap ecmap;
  int ne;

  Component_graph(G& g, ECMap ecmap,  tbb::concurrent_vector<edge_descriptor>&cc)
    : g(g), cc(cc), ecmap(ecmap)
  {}
  
  typename boost::graph_traits<G>::edges_size_type num_edges() const
  {
    return static_cast<typename boost::graph_traits<G>::edges_size_type>(cc.size());
  }
}; 

} } // namespace CGAL::internal


namespace boost {

template <typename G, typename ECMap>
struct graph_traits<CGAL::internal::Component_graph<G,ECMap> > {
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
  typedef typename CGAL::internal::Component_graph<G,ECMap>::edge_iterator edge_iterator;

  static face_descriptor null_face()
  {
    return boost::graph_traits<G>::null_face();
  }
};
} // namespace boost

namespace CGAL {
namespace internal {
template <typename G, typename ECMap>
CGAL::Iterator_range<typename tbb::concurrent_vector<typename boost::graph_traits<G>::edge_descriptor>::iterator>
edges(const Component_graph<G,ECMap>& cg)
{
  return CGAL::make_range(cg.cc.begin(), cg.cc.end());
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
              Component_graph<G,ECMap>& cg,
              bool mark_only = false)
{
  remove_vertex(v, cg.g, mark_only);
}

template <typename G, typename ECMap>
void
remove_edge(typename boost::graph_traits<G>::edge_descriptor e,
     Component_graph<G,ECMap>& cg,
              bool mark_only = false)
{
  return remove_edge(e, cg.g, mark_only);
}

template <typename G, typename ECMap>
void
remove_face(typename boost::graph_traits<G>::face_descriptor f,
            Component_graph<G,ECMap>& cg,
            bool mark_only = false)
{
  return remove_face(f, cg.g, mark_only);
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

} } // namespace CGAL::internal

namespace boost {
template <typename G, typename ECMap>
struct property_map<CGAL::internal::Component_graph<G,ECMap>,CGAL::vertex_point_t>
{
  typedef typename boost::property_map<G,CGAL::vertex_point_t>::type type;
};

} // namespace boost

namespace CGAL {
namespace internal {
template <typename G, typename ECMap>
typename boost::property_map<Component_graph<G,ECMap>,CGAL::vertex_point_t>::type
get(CGAL::vertex_point_t t, Component_graph<G,ECMap>& cg)
{
  return get(t,cg.g);
}

} } // namespace CGAL::internal

#endif // CGAL_BOOST_GRAPH_INTERNAL_COMPONENT_GRAPH_H
