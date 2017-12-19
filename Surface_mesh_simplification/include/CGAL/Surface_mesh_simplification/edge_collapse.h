// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Detail/Edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <CGAL/Surface_mesh_simplification/Detail/parallel_edge_collapse.h>
#endif

#include <CGAL/use.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/named_function_params.h>

#include <CGAL/mutex.h>
#include <CGAL/tags.h>

namespace CGAL {

namespace Surface_mesh_simplification {

template<class TriangleMesh,
         class ShouldStop,
         class ConcurrencyTag,
         class VertexIndexMap,
         class VertexPointMap,
         class EdgeIndexMap,
         class EdgeIsConstrainedMap,
         class GetCost,
         class GetPlacement,
         class Visitor>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  ConcurrencyTag aConcurrency_tag,
                  CGAL_MUTEX* removal_mutex,

                  // optional mesh information policies
                  std::size_t current_num_edges,
                  const VertexIndexMap& aVertex_index_map, // defaults to get(vertex_index,aSurface)
                  const VertexPointMap& aVertex_point_map, // defaults to get(vertex_point,aSurface)
                  const EdgeIndexMap& aEdge_index_map, // defaults to get(edge_index,aSurface)
                  const EdgeIsConstrainedMap& aEdge_is_constrained_map, // defaults to No_constrained_edge_map<TriangleMesh>()

                  // optional strategy policies - defaults to LindstomTurk
                  const GetCost& aGet_cost,
                  const GetPlacement& aGet_placement,
                  Visitor aVisitor);


template<class TriangleMesh,
         class ShouldStop,
         class VertexIndexMap,
         class VertexPointMap,
         class EdgeIndexMap,
         class EdgeIsConstrainedMap,
         class GetCost,
         class GetPlacement,
         class Visitor>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  Sequential_tag,
                  CGAL_MUTEX* removal_mutex,

                  // optional mesh information policies
                  std::size_t current_num_edges,
                  const VertexIndexMap& aVertex_index_map, // defaults to get(vertex_index, aSurface)
                  const VertexPointMap& aVertex_point_map, // defaults to get(vertex_point, aSurface)
                  const EdgeIndexMap& aEdge_index_map, // defaults to get(edge_index, aSurface)
                  const EdgeIsConstrainedMap& aEdge_is_constrained_map, // defaults to No_constrained_edge_map<TriangleMesh>()

                  // optional strategy policies - defaults to LindstomTurk
                  const GetCost& aGet_cost,
                  const GetPlacement& aGet_placement,
                  Visitor aVisitor)
{
  typedef EdgeCollapse<TriangleMesh, ShouldStop, VertexIndexMap, VertexPointMap, EdgeIndexMap,
                       EdgeIsConstrainedMap, GetCost, GetPlacement, Visitor,
                       Edge_profile<TriangleMesh, VertexPointMap> >           Algorithm;

  Algorithm algorithm(aSurface, current_num_edges, aShould_stop, aVertex_index_map,
                      aVertex_point_map, aEdge_index_map, aEdge_is_constrained_map,
                      aGet_cost, aGet_placement, aVisitor);

  return algorithm.run(removal_mutex);
}


template<class TriangleMesh,
         class ShouldStop,
         class VertexIndexMap,
         class VertexPointMap,
         class EdgeIndexMap,
         class EdgeIsConstrainedMap,
         class GetCost,
         class GetPlacement,
         class Visitor>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  Parallel_tag,
                  CGAL_MUTEX* removal_mutex,

                  // optional mesh information policies
                  std::size_t current_num_edges,
                  const VertexIndexMap& aVertex_index_map, // defaults to get(vertex_index,aSurface)
                  const VertexPointMap& aVertex_point_map, // defaults to get(vertex_point,aSurface)
                  const EdgeIndexMap& aEdge_index_map, // defaults to get(edge_index,aSurface)
                  const EdgeIsConstrainedMap& aEdge_is_constrained_map, // defaults to No_constrained_edge_map<TriangleMesh>()

                  // optional strategy policies - defaults to LindstomTurk
                  const GetCost& aGet_cost,
                  const GetPlacement& aGet_placement,
                  Visitor aVisitor)
{
  typedef EdgeCollapse<TriangleMesh, ShouldStop, VertexIndexMap, VertexPointMap, EdgeIndexMap,
                       EdgeIsConstrainedMap, GetCost, GetPlacement, Visitor,
                       CG_Edge_profile<Edge_profile<TriangleMesh, VertexPointMap> > > Algorithm;

  Algorithm algorithm(aSurface, current_num_edges, aShould_stop, aVertex_index_map,
                      aVertex_point_map, aEdge_index_map, aEdge_is_constrained_map,
                      aGet_cost, aGet_placement, aVisitor);

  return algorithm.run(removal_mutex);
}

template <typename G>
struct Dummy_visitor
{
  typedef typename boost::graph_traits<G>::edges_size_type   size_type;

  template<class TriangleMesh>                              void OnStarted(TriangleMesh&) const {}
  template<class TriangleMesh, class Stop, class Size_type> void OnParallelPassFinished(TriangleMesh&, Stop&, Size_type, Size_type) const {}
  template<class TriangleMesh>                              void OnFinished(TriangleMesh&) const {}
  template<class Profile>                                   void OnStopConditionReached(const Profile&) const {}
  template<class Profile, class OFT>                        void OnCollected(const Profile&, const OFT&) const {}
  template<class Profile, class OFT, class Size_type>       void OnSelected(const Profile&, const OFT&, Size_type, Size_type) const {}
  template<class Profile, class OPoint>                     void OnCollapsing(const Profile&, const OPoint&) const {}
  template<class Profile, class VH>                         void OnCollapsed(const Profile&, VH) const {}
  template<class Profile>                                   void OnNonCollapsable(const Profile&) const {}
};

template<class TriangleMesh, class ShouldStop, class NamedParameters>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  Sequential_tag,
                  CGAL_MUTEX* removal_mutex,
                  const NamedParameters& aParams)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  internal_np::graph_visitor_t vis = internal_np::graph_visitor_t();
   return edge_collapse(aSurface, aShould_stop, Sequential_tag(), removal_mutex,
                        choose_param(get_param(aParams, internal_np::current_num_edges), num_edges(aSurface)),
                        choose_const_pmap(get_param(aParams, internal_np::vertex_index), aSurface, boost::vertex_index),
                        choose_pmap(get_param(aParams, internal_np::vertex_point), aSurface, boost::vertex_point),
                        choose_const_pmap(get_param(aParams, internal_np::halfedge_index), aSurface, boost::halfedge_index),
                        choose_param(get_param(aParams, internal_np::edge_is_constrained), No_constrained_edge_map<TriangleMesh>()),
                        choose_param(get_param(aParams, internal_np::get_cost_policy), LindstromTurk_cost<TriangleMesh>()),
                        choose_param(get_param(aParams, internal_np::get_placement_policy), LindstromTurk_placement<TriangleMesh>()),
                        choose_param(get_param(aParams, vis), Dummy_visitor<TriangleMesh>()));
}


template<class TriangleMesh, class ShouldStop, class NamedParameters>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  Parallel_tag,
                  CGAL_MUTEX* removal_mutex,
                  const NamedParameters& aParams)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  internal_np::graph_visitor_t vis = internal_np::graph_visitor_t();
   return edge_collapse(aSurface, aShould_stop, Parallel_tag(), removal_mutex,
                        choose_param(get_param(aParams, internal_np::current_num_edges), num_edges(aSurface)),
                        choose_const_pmap(get_param(aParams, internal_np::vertex_index), aSurface, boost::vertex_index),
                        choose_pmap(get_param(aParams, internal_np::vertex_point), aSurface, boost::vertex_point),
                        choose_const_pmap(get_param(aParams, internal_np::halfedge_index), aSurface, boost::halfedge_index),
                        choose_param(get_param(aParams, internal_np::edge_is_constrained), No_constrained_edge_map<TriangleMesh>()),
                        choose_param(get_param(aParams, internal_np::get_cost_policy), LindstromTurk_cost<TriangleMesh>()),
                        choose_param(get_param(aParams, internal_np::get_placement_policy), LindstromTurk_placement<TriangleMesh>()),
                        choose_param(get_param(aParams, vis), Dummy_visitor<TriangleMesh>()));
}

template<class TriangleMesh, class ShouldStop, class NamedParameters>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  Sequential_tag tag,
                  const NamedParameters& aParams)
{
  return edge_collapse(aSurface, aShould_stop, tag, NULL, aParams);
}

template<class TriangleMesh, class ShouldStop, class NamedParameters>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  const NamedParameters& aParams)
{
  return edge_collapse(aSurface, aShould_stop, Sequential_tag(), NULL, aParams);
}

template<class TriangleMesh, class ShouldStop, class GT, class NamedParameters>
int edge_collapse(TriangleMesh& aSurface,
                  const ShouldStop& aShould_stop,
                  const NamedParameters& aParams)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  internal_np::graph_visitor_t vis = internal_np::graph_visitor_t();
  return edge_collapse(aSurface, aShould_stop,
                       choose_param(get_param(aParams, internal_np::current_num_edges), num_edges(aSurface)),
                       choose_const_pmap(get_param(aParams, internal_np::vertex_index), aSurface, boost::vertex_index),
                       choose_const_pmap(get_param(aParams, internal_np::vertex_point), aSurface, boost::vertex_point),
                       choose_const_pmap(get_param(aParams, internal_np::halfedge_index), aSurface, boost::halfedge_index),
                       choose_param(get_param(aParams, internal_np::edge_is_constrained), No_constrained_edge_map<TriangleMesh>()),
                       choose_param(get_param(aParams, internal_np::get_cost_policy), LindstromTurk_cost<TriangleMesh>()),
                       choose_param(get_param(aParams, internal_np::get_placement_policy), LindstromTurk_placement<TriangleMesh>()),
                       choose_param(get_param(aParams, vis), Dummy_visitor<TriangleMesh>()));
}


template<class TriangleMesh, class ShouldStop>
int edge_collapse(TriangleMesh& aSurface, const ShouldStop& aShould_stop)
{
  return edge_collapse(aSurface, aShould_stop, CGAL::parameters::halfedge_index_map(get(boost::halfedge_index, aSurface)));
}

template<class TriangleMesh, class ShouldStop, class GT>
int edge_collapse(TriangleMesh& aSurface, const ShouldStop& aShould_stop)
{
  return edge_collapse(aSurface, aShould_stop, CGAL::parameters::halfedge_index_map(get(boost::halfedge_index, aSurface)));
}

template<class TriangleMesh, class ShouldStop, class FacePartionMap, class NamedParameters>
int parallel_edge_collapse(TriangleMesh& aSurface,
                           const ShouldStop& aShould_stop,
                           FacePartionMap partition_id_map,
                           int number_of_parts,
                           const NamedParameters& aParams)
{
#ifndef CGAL_LINKED_WITH_TBB
  std::cerr << "CGAL is not linked with TBB, the mesh will be simplified sequentially." << std::endl;
  number_of_parts = 1;
#endif

  if(number_of_parts == 1)
  {
    return edge_collapse(aSurface, aShould_stop, Sequential_tag(), NULL, aParams);
  }
  else
  {
#ifdef CGAL_LINKED_WITH_TBB
    using boost::choose_param;
    using boost::get_param;

    internal_np::graph_visitor_t vis = internal_np::graph_visitor_t();
    return parallel_edge_collapse(aSurface,
                                  aShould_stop,
                                  partition_id_map,
                                  number_of_parts,
                                  choose_param(get_param(aParams, internal_np::edge_is_constrained), No_constrained_edge_map<TriangleMesh>()),
                                  choose_param(get_param(aParams, internal_np::get_placement_policy), LindstromTurk_placement<TriangleMesh>()),
                                  choose_param(get_param(aParams, internal_np::get_cost_policy), LindstromTurk_cost<TriangleMesh>()),
                                  choose_param(get_param(aParams,vis), Dummy_visitor<TriangleMesh>()));
#else
    CGAL_USE(partition_id_map);
#endif
  }
}

} // end namespace Surface_mesh_simplification

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_H
