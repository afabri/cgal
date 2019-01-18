// Copyright (c) 2017  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Andreas Fabri
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_PARALLEL_EDGE_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_PARALLEL_EDGE_COLLAPSE_H 1

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Parallel_stop_predicate_visitor.h>

#include <CGAL/boost/graph/internal/Component_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Real_timer.h>
#include <CGAL/mutex.h>

#include <tbb/task_group.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace CGAL {

namespace Surface_mesh_simplification {

namespace internal {

template <typename ECMap_1, typename ECMap_2>
struct Or_map
{
  typedef boost::readable_property_map_tag                   category;
  typedef bool                                               value_type;
  typedef bool                                               reference;
  typedef typename boost::property_traits<ECMap_1>::key_type key_type;

  ECMap_1 ecm1;
  ECMap_2 ecm2;

  Or_map(const ECMap_1& ecm1, const ECMap_2& ecm2)
    : ecm1(ecm1), ecm2(ecm2)
  {}

  friend bool get(const Or_map& orm, key_type ed)
  {
    return get(orm.ecm1, ed) || get(orm.ecm2, ed);
  }
};


// Adaptor which interprets 0 and 1 as true, and any other index as false
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


// Adapter which interprets char as bool to avoid "performance warning" with VC++
template <typename Char_map>
struct Bool_map
{
  Char_map cm;

  Bool_map(Char_map cm)
    : cm(cm)
  {}

  template <typename K>
  friend bool get(const Bool_map& bm, const K& k)
  {
    return get(bm.cm, k) != 0;
  }
};


template <typename ECMap>
struct Inverted_edge_constraint_map
{
  typedef boost::readable_property_map_tag                 category;
  typedef bool                                             value_type;
  typedef bool                                             reference;
  typedef typename boost::property_traits<ECMap>::key_type key_type;

  ECMap ecmap;

  Inverted_edge_constraint_map(ECMap ecmap)
    : ecmap(ecmap)
  {}

  friend bool get(const Inverted_edge_constraint_map& iecm, key_type ed)
  {
    return ! (get(iecm.ecmap, ed) != 0);
  }
};


template <typename Graph, typename ComponentIdMap, typename EdgeConstraintMap>
struct In_same_component_map
{
  typedef typename boost::property_traits<EdgeConstraintMap>::key_type   key_type;
  typedef typename boost::property_traits<EdgeConstraintMap>::value_type value_type;
  typedef value_type                                                     reference;
  typedef boost::read_write_property_map_tag                             category;

  EdgeConstraintMap ecmap;
  ComponentIdMap ccmap;
  const Graph& g;
  std::size_t id;
  int& count;

  In_same_component_map(const Graph& g, ComponentIdMap ccmap, EdgeConstraintMap ecmap, int id, int& count)
    : ecmap(ecmap), ccmap(ccmap), g(g), id(id), count(count)
  {}

  friend value_type get(const In_same_component_map& iscm, key_type ed)
  {
    return get(iscm.ecmap, ed);
  }

  friend void put(const In_same_component_map& iscm, key_type ed, value_type v)
  {
    if(! get(iscm.ecmap, ed))
    {
      if((! is_border(halfedge(ed, iscm.g), iscm.g) &&
          iscm.ccmap[face(halfedge(ed, iscm.g), iscm.g)] == iscm.id)||
         (! is_border(opposite(halfedge(ed, iscm.g), iscm.g), iscm.g) &&
          iscm.ccmap[face(opposite(halfedge(ed, iscm.g), iscm.g), iscm.g)] == iscm.id))
      {
        put(iscm.ecmap, ed, v);
        ++(iscm.count);
      }
    }
  }
};


// The parallel task
template <typename TriangleMesh, typename Stop, typename CCMap, typename HIMap,
          typename VIpmap, typename ECMap, typename UECMap, typename Placement, typename Cost,
          typename Visitor>
struct Simplify
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor     edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  TriangleMesh& sm;
  Stop stop;
  CCMap ccmap;
  ECMap ecmap;
  UECMap uecmap;
  Placement placement;
  Cost cost;
  HIMap himap;
  VIpmap vim;
  std::vector<int>& buffer_size;
  std::vector<int>& results;
  tbb::concurrent_vector<edge_descriptor>& cc_edges;
  int cc_index;
  unsigned int layers;
  bool dump, verbose, increase;
  CGAL_MUTEX* removal_mutex;
  Visitor& visitor;

  Simplify(TriangleMesh& sm,
           Stop stop,
           CCMap ccmap,
           HIMap himap,
           VIpmap vim,
           ECMap ecmap,
           UECMap uecmap,
           Placement placement,
           Cost cost,
           std::vector<int>& buffer_size,
           std::vector<int>& results,
           tbb::concurrent_vector<edge_descriptor>& cc_edges,
           int cc_index,
           unsigned int layers,
           bool dump,
           bool verbose,
           bool increase,
           CGAL_MUTEX* removal_mutex,
           Visitor& visitor)
    :
      sm(sm),
      stop(stop),
      ccmap(ccmap),
      ecmap(ecmap),
      uecmap(uecmap),
      placement(placement),
      cost(cost),
      himap(himap),
      vim(vim),
      buffer_size(buffer_size),
      results(results),
      cc_edges(cc_edges),
      cc_index(cc_index),
      layers(layers),
      dump(dump),
      verbose(verbose),
      increase(increase),
      removal_mutex(removal_mutex),
      visitor(visitor)
  { }

  void operator()() const
  {
    typedef Selection_is_constraint_map<HIMap, TriangleMesh> SICM;
    typedef CGAL::internal::Component_graph<TriangleMesh,SICM> Component_graph;

    SICM sicm(himap, sm);
    Component_graph cg(sm, sicm, cc_edges);

    std::vector<edge_descriptor> buffer_edges;
    BOOST_FOREACH(edge_descriptor ed, cc_edges)
    {
      halfedge_descriptor hd, hop;
      hd = halfedge(ed, sm);
      hop = opposite(hd, sm);
      if(get(sicm, ed))
      {
        buffer_edges.push_back(ed);
      }
    }

    std::size_t buffer_edges_count = buffer_edges.size();
    std::vector<edge_descriptor> V;
    int number_of_puts = 0;
    In_same_component_map<TriangleMesh, CCMap, ECMap> iscmap(sm, ccmap, ecmap, cc_index, number_of_puts);

    expand_edge_selection(buffer_edges,
                          sm,
                          layers,
                          iscmap,
                          std::back_inserter(V)); // V might contain edges from an adjacent patch

    buffer_edges_count += number_of_puts;
    buffer_size[cc_index] = number_of_puts;

    // A clumsy solution: The thread must constrain also in the neighbor patch as
    // otherwise when we circle around a vertex on the border between two patches
    // may be not constrained, and then have an illegal index
    BOOST_FOREACH(edge_descriptor ed, V)
    {
      put(ecmap, ed, true);
    }

    if(dump)
    {
      std::ostringstream filename;
      filename << "constraints-" << cc_index << ".selection.txt" << std::ends;
      std::ofstream out(filename.str().c_str());
      out << std::endl << std::endl; // edges are on the third line of a CGAL selection file
      BOOST_FOREACH(edge_descriptor ed, V)
      {
        if(get(iscmap,ed))
          out << get(vim, source(ed, sm)) << " " << get(vim, target(ed,sm)) << " ";
      }
      out << std::endl;
    }

    if(verbose)
    {
      std::ostringstream oss;
      oss << "["<< cc_index << "]\t# constrained edges in the partition = "
          << buffer_edges_count << "\n" << std::ends;
      std::cerr << oss.str();
    }

    std::size_t num_unconstrained_edges = cc_edges.size() - buffer_edges_count;

    typedef internal::Or_map<ECMap, UECMap> OrMap;
    OrMap ormap(ecmap, uecmap);

    if(increase)
    {
      Constrained_placement<Placement, OrMap> constrained_placement(ormap, placement);
      results[cc_index]
        = edge_collapse(cg,
                        stop,
                        Parallel_tag(),
                        removal_mutex,
                        parameters::vertex_index_map(vim)
                        .current_num_edges(num_unconstrained_edges - int(double(buffer_edges_count)/0.75))
                        .halfedge_index_map(himap)
                        .get_placement(constrained_placement)
                        .get_cost(cost)
                        .edge_is_constrained_map(ormap)
                        .vertex_point_map(get(CGAL::vertex_point, sm))
                        .visitor(visitor));
    }
    else
    {
      results[cc_index]
        = edge_collapse(cg,
                        stop,
                        Parallel_tag(),
                        removal_mutex,
                        parameters::vertex_index_map(vim)
                        .current_num_edges(num_unconstrained_edges)
                        .halfedge_index_map(himap)
                        .get_placement(placement)
                        .get_cost(cost)
                        .edge_is_constrained_map(ormap)
                        .vertex_point_map(get(CGAL::vertex_point, sm))
                        .visitor(visitor));
    }

    if(verbose)
    {
      std::ostringstream oss;
      oss << "[" << cc_index << "]" << std::endl
          << "\tRemoved " << results[cc_index] << std::endl
          << "\t|| done" << std::endl << std::ends;
      std::cerr << oss.str();
    }
  }
};

} // end namespace internal

template <typename TriangleMesh, typename Stop, typename FacePartitionIDPmap,
          typename VertexIndexMap, typename UECMap, typename Placement, typename Cost,
          typename Visitor>
int parallel_edge_collapse(TriangleMesh& sm,
                           Stop stop,
                           std::size_t number_of_parts,
                           FacePartitionIDPmap partition_id_map,
                           VertexIndexMap vim,
                           UECMap uecmap,
                           Placement placement,
                           Cost cost,
                           Visitor pvis,
                           unsigned int layers = 1,
                           bool increase = true, // whether to increase the buffer by one more layer
                           bool verbose = false,
                           bool dump = false)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edges_size_type size_type;

  Real_timer t;
  Parallel_stop_predicate_visitor<TriangleMesh, Visitor> vis(pvis);


  typedef CGAL::dynamic_halfedge_property_t<int> Halfedge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Halfedge_property_tag>::type HIMap;
  HIMap himap = get(Halfedge_property_tag(), sm);

  typedef CGAL::dynamic_edge_property_t<char> Edge_property_tag;
  typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type ECMap;
  ECMap ecmap = get(Edge_property_tag(), sm);

  std::vector<edge_descriptor> partition_edges;

  std::ofstream out;
  if(dump)
  {
    out.open("partition.selection.txt");
    out << std::endl << std::endl; // edges are on the third line of a CGAL selection file
  }

  // now we have to find the constrained edges and set halfedge indices
  std::vector<int> hindices(number_of_parts, 1);
  BOOST_FOREACH(edge_descriptor ed, edges(sm))
  {
    halfedge_descriptor hd = halfedge(ed, sm);
    halfedge_descriptor hop = opposite(hd, sm);
    if ( is_border(hd, sm) )
      std::swap(hd, hop);

    typename boost::property_traits<FacePartitionIDPmap>::value_type
      pid = get(partition_id_map, face(hd, sm)),
      pid_opp = is_border(hop, sm) ? pid : get(partition_id_map, face(hop, sm));

    if(pid != pid_opp)
    {
      partition_edges.push_back(ed);
      put(ecmap, ed, true);
      put(himap, hd, 0);
      put(himap, hop, 1);

      if(dump)
        out << get(vim, source(ed, sm)) << " " << get(vim, target(ed, sm)) << " ";
    }
    else
    {
      put(himap, hd, ++hindices[pid]);
      put(himap, hop, ++hindices[pid]);
      put(ecmap, ed, false);
    }
  }

  if(dump)
    out << std::endl;

  if(verbose)
  {
    std::cerr << t.time() << " sec.\n";
    t.reset();
  }

  std::vector<tbb::concurrent_vector<edge_descriptor> > cc_edges(number_of_parts);
  BOOST_FOREACH(face_descriptor fd, faces(sm))
  {
    std::size_t partition_id = partition_id_map[fd];
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, sm), sm))
    {
      halfedge_descriptor hop = opposite(hd, sm);
      if(is_border(hop, sm) || get(ecmap, edge(hd, sm)) || (fd < face(hop, sm)))
      {
        cc_edges[partition_id].push_back(edge((hd<hop) ? hd : hop, sm));
      }
    }
  }

  if(verbose)
  {
    for(std::size_t id=0; id<number_of_parts; ++id)
      std::cout << "#edges in component "<< id << ": " << cc_edges[id].size() << std::endl;
  }

  std::vector<int> buffer_size(number_of_parts); // the number of constrained edges inside each component
  std::vector<int> removed(number_of_parts); // the number of collapsed edges in each thread
  if(verbose)
  {
    std::cerr << "#CC = " << number_of_parts << " in " << t.time() << " sec." << std::endl;
    t.reset();
    t.start();
  }

  size_type initial_num_edges = num_edges(sm);

  // Simplify the partition in parallel
  CGAL_MUTEX removal_mutex;
  tbb::task_group tasks;
  for(std::size_t id=0; id<number_of_parts; ++id)
  {
    tasks.run(
      internal::Simplify<TriangleMesh, Stop, FacePartitionIDPmap, HIMap, VertexIndexMap, ECMap, UECMap, Placement, Cost, Visitor>(
        sm, stop, partition_id_map, himap, vim, ecmap, uecmap, placement, cost, buffer_size, removed,
        cc_edges[id], static_cast<int>(id), layers, dump, verbose, increase, &removal_mutex, pvis));
  }
  tasks.wait();

  int result = 0;
  for(std::size_t i=0; i<removed.size(); ++i)
    result += removed[i];

  if(verbose)
  {
    std::cerr << "parallel edge collapse in " << t.time() << " sec." << std::endl;
    t.reset();
  }

  if(dump)
  {
    TriangleMesh sm2;
    copy_face_graph(sm, sm2);
    std::ofstream outi("out-intermediary.off");
    outi << sm2 << std::endl;
  }

  // After the parallel edge_collapse make the same for the buffer
  std::size_t num_unconstrained_edges = 0;
  std::size_t overlapping_unconstrained_edges = 0;
  if(increase)
  {
    // Increase the buffer by one more layer
    BOOST_FOREACH(edge_descriptor ed, edges(sm))
    {
      put(ecmap, ed, 0);
    }

    BOOST_FOREACH(edge_descriptor ed, partition_edges)
    {
      put(ecmap, ed, 1);
    }

    expand_edge_selection(partition_edges, sm, layers+1, ecmap,
                          Counting_output_iterator(&num_unconstrained_edges));
    num_unconstrained_edges += partition_edges.size();

    overlapping_unconstrained_edges = num_unconstrained_edges - partition_edges.size();
    for(std::size_t i=0; i<buffer_size.size(); ++i)
      overlapping_unconstrained_edges -= buffer_size[i];
  }
  else
  {
    num_unconstrained_edges = partition_edges.size();
    for(std::size_t i=0; i<buffer_size.size(); ++i)
    {
      if(verbose)
        std::cout << "# partition edges = "<< partition_edges.size() << " " << buffer_size[i] << std::endl;

      num_unconstrained_edges += buffer_size[i];
    }
  }

  typedef internal::Inverted_edge_constraint_map<ECMap> IECMap;
  typedef internal::Or_map<IECMap,UECMap> OrMap;
  IECMap iecmap(ecmap);
  OrMap ormap(iecmap, uecmap);

  // AF: At this point num_edges(sm) is not the same for a Surface_mesh and a Polyhedron.

  size_type current_num_edges = initial_num_edges - result;
  vis.OnParallelPassFinished(sm, stop, static_cast<size_type>(initial_num_edges), static_cast<size_type>(current_num_edges));

  // Re-number halfedges and indices after the parallel pass
  int index = 0;
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(sm))
  {
    put(himap, hd, index++);
  }

  if(verbose)
    std::cerr << "#removed edges = " << result << std::endl;

  if(increase)
  {
    Constrained_placement <Placement, OrMap> constrained_placement (ormap, placement);
    result += edge_collapse(sm,
                            stop,
                            Sequential_tag(),
                            NULL,
                            parameters::vertex_index_map(vim)
                            .halfedge_index_map(himap)
                            .current_num_edges(current_num_edges)
                            .get_placement(constrained_placement)
                            .get_cost(cost)
                            .edge_is_constrained_map(ormap)
                            .visitor(pvis));
  }
  else
  {
    result += edge_collapse(sm,
                            stop,
                            Sequential_tag(),
                            NULL,
                            parameters::vertex_index_map(vim)
                            .halfedge_index_map(himap)
                            .current_num_edges(current_num_edges)
                            .get_placement(placement)
                            .get_cost(cost)
                            .edge_is_constrained_map(ormap)
                            .visitor(pvis));
  }

  if(verbose)
  {
    std::cerr << "sequential edge collapse on buffer in " << t.time() << " sec." << std::endl;
    std::cerr << "After sequential pass, #removed edges = " << result << std::endl;
  }

  return result;
}

} // end namespace Surface_mesh_simplification
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_PARALLEL_EDGE_COLLAPSE_H
