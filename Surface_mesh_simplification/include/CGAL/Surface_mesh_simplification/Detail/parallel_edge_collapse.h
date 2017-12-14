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

#include <CGAL/boost/graph/Component_graph.h>
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
    return get(bm.cm,k) != 0;
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
    return ! (get(iecm.ecmap, ed)!= 0);
  }
};


template <typename ECMap, typename CCMap, typename G>
struct In_same_component_map
{
  typedef boost::read_write_property_map_tag                 category;
  typedef typename boost::property_traits<ECMap>::value_type value_type;
  typedef value_type                                         reference;
  typedef typename boost::property_traits<ECMap>::key_type   key_type;

  ECMap ecmap;
  CCMap ccmap;
  const G& g;
  std::size_t cc;
  int& count;

  In_same_component_map(ECMap ecmap, CCMap ccmap, const G& g, int cc, int& count)
    : ecmap(ecmap), ccmap(ccmap), g(g), cc(cc), count(count)
  {}

  friend value_type get(const In_same_component_map& iscm, key_type ed)
  {
    return get(iscm.ecmap,ed);
  }

  friend void put(const In_same_component_map& iscm, key_type ed, value_type v)
  {
    if(! get(iscm.ecmap, ed))
    {
      if((! is_border(halfedge(ed,iscm.g),iscm.g) &&
          iscm.ccmap[face(halfedge(ed,iscm.g),iscm.g)] == iscm.cc)||
         (! is_border(opposite(halfedge(ed,iscm.g),iscm.g),iscm.g) &&
          iscm.ccmap[face(opposite(halfedge(ed,iscm.g),iscm.g),iscm.g)] == iscm.cc))
      {
        put(iscm.ecmap, ed, v);
        ++(iscm.count);
      }
    }
  }
};


// The parallel task
template <typename TriangleMesh, typename Placement, typename Cost, typename Stop,
          typename HIMap, typename ECMap, typename UECMap, typename CCMap>
struct Simplify
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor     edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  TriangleMesh& sm;
  HIMap himap;
  ECMap ecmap;
  UECMap uecmap;
  CCMap ccmap;
  Placement placement;
  Cost cost;
  Stop stop;
  std::vector<int>& buffer_size;
  std::vector<int>& results;
  tbb::concurrent_vector<edge_descriptor>& cc_edges;
  int ccindex;
  unsigned int layers;
  bool dump, verbose, increase;
  CGAL_MUTEX* removal_mutex;

  Simplify(TriangleMesh& sm,
           HIMap himap,
           ECMap ecmap,
           UECMap uecmap,
           CCMap ccmap,
           Placement placement,
           Cost cost,
           Stop stop,
           std::vector<int>& buffer_size,
           std::vector<int>& results,
           tbb::concurrent_vector<edge_descriptor>& cc_edges,
           int ccindex,
           unsigned int layers,
           bool dump,
           bool verbose,
           bool increase,
           CGAL_MUTEX* removal_mutex)
    :
      sm(sm),
      himap(himap),
      ecmap(ecmap),
      uecmap(uecmap),
      ccmap(ccmap),
      placement(placement),
      cost(cost),
      stop(stop),
      buffer_size(buffer_size),
      results(results),
      cc_edges(cc_edges),
      ccindex(ccindex),
      layers(layers),
      dump(dump),
      verbose(verbose),
      increase(increase),
      removal_mutex(removal_mutex)
  { }

  void operator()() const
  {
    typedef Selection_is_constraint_map<HIMap, TriangleMesh> SICM;
    typedef Component_graph<TriangleMesh,SICM>               Component_graph;
    SICM sicm(himap,sm);
    Component_graph cg(sm, sicm, cc_edges);

    std::vector<edge_descriptor> buffer_edges;
    int i = 2;
    BOOST_FOREACH(edge_descriptor ed, cc_edges)
    {
      halfedge_descriptor hd, hop;
      hd = halfedge(ed,sm);
      hop = opposite(hd,sm);
      if(get(sicm,ed))
      {
        buffer_edges.push_back(ed);
      }
      else
      {
        put(himap, hd, i++);
        put(himap, hop, i++);
      }
    }

    std::size_t buffer_edges_count = buffer_edges.size();
    std::vector<edge_descriptor> V;
    int number_of_puts = 0;
    In_same_component_map<ECMap,CCMap,TriangleMesh> iscmap(ecmap, ccmap, sm, ccindex, number_of_puts);

    expand_edge_selection(buffer_edges,
                          sm,
                          layers,
                          iscmap,
                          std::back_inserter(V)); // V might contain edges from an adjacent patch

    buffer_edges_count += number_of_puts;
    buffer_size[ccindex] = number_of_puts;

    // A clumsy solution: The thread must constrain also in the neighbor patch as
    // otherwise when we circle around a vertex on the border between two patches
    // may be not constrained, and then have an illegal index
    BOOST_FOREACH(edge_descriptor ed, V)
    {
      put(ecmap, ed, true);
    }

    if(dump)
    {
      typedef typename boost::property_map<TriangleMesh, boost::vertex_index_t>::const_type Indices;
      Indices indices = get(boost::vertex_index, sm);

      std::ostringstream filename;
      filename << "constraints-" << ccindex << ".selection.txt" << std::ends;
      std::ofstream out(filename.str().c_str());
      out << std::endl << std::endl; // edges are on the third line of a CGAL selection file
      BOOST_FOREACH(edge_descriptor ed, V)
      {
        if(get(iscmap,ed))
          out << get(indices, source(ed, sm)) << " " << get(indices, target(ed,sm)) << " ";
      }
      out << std::endl;
    }

    if(verbose)
    {
      std::ostringstream oss;
      oss << "["<< ccindex << "]\t# constrained edges in the partition = "
          << buffer_edges_count << "\n" << std::ends;
      std::cerr << oss.str();
    }

    std::size_t num_unconstrained_edges = cc_edges.size() - buffer_edges_count;

    typedef internal::Or_map<ECMap, UECMap> OrMap;
    OrMap ormap(ecmap, uecmap);

    if(increase)
    {
      Constrained_placement <Placement, OrMap> constrained_placement(ormap, placement);
      results[ccindex]
        = edge_collapse(cg,
                        stop,
                        Parallel_tag(),
                        removal_mutex,
                        parameters::vertex_index_map(get(boost::vertex_index,sm))
                        .current_num_edges(num_unconstrained_edges - int(double(buffer_edges_count)/0.75))
                        .halfedge_index_map(himap)
                        .get_placement(constrained_placement)
                        .get_cost(cost)
                        .edge_is_constrained_map(ormap)
                        .vertex_point_map(get(CGAL::vertex_point,sm)));
    }
    else
    {
      results[ccindex]
        = edge_collapse(cg,
                        stop,
                        Parallel_tag(),
                        removal_mutex,
                        parameters::vertex_index_map(get(boost::vertex_index,sm))
                        .current_num_edges(num_unconstrained_edges)
                        .halfedge_index_map(himap)
                        .get_placement(placement)
                        .get_cost(cost)
                        .edge_is_constrained_map(ormap)
                        .vertex_point_map(get(CGAL::vertex_point,sm)));
    }

    if(verbose)
    {
      std::ostringstream oss;
      oss << "[" << ccindex << "]\tratio = " << stop.ratio() << std::endl
          << "\tRemoved " << results[ccindex] << std::endl
          << "\t|| done" << std::endl << std::ends;
      std::cerr << oss.str();
    }
  }
};

} // end namespace internal

template <typename TriangleMesh, typename Placement, typename CCMap, typename UECMap, typename Stop, typename Cost, typename Visitor>
int parallel_edge_collapse(TriangleMesh& sm,
                           CCMap ccmap, UECMap uecmap,
                           Placement placement, Stop stop, Cost cost,
                           std::size_t ncc,
                           Visitor pvis,
                           unsigned int layers = 1,
                           bool dump = false,
                           bool verbose = false,
                           bool increase = true)
{
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edges_size_type size_type;

  Real_timer t;
  Parallel_stop_predicate_visitor<TriangleMesh, Visitor> vis(pvis);

  typedef CGAL::internal::halfedge_property_t<int> Halfedge_property_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Halfedge_property_tag >::type HIMap;
  HIMap himap = CGAL::internal::add_property(Halfedge_property_tag("h:internal::index_in_cc"), sm);

  typedef CGAL::internal::edge_property_t<char> Edge_property_tag;
  typedef typename CGAL::internal::dynamic_property_map<TriangleMesh, Edge_property_tag >::type ECMap;
  ECMap ecmap = CGAL::internal::add_property(Edge_property_tag("e:internal::constrained"), sm);

  typedef typename boost::property_map<TriangleMesh, boost::vertex_index_t>::const_type Indices;
  Indices indices = get(boost::vertex_index, sm);

  std::vector<edge_descriptor> partition_edges;

  std::ofstream out;
  if(dump)
  {
    out.open("partition.selection.txt");
    out << std::endl << std::endl; // edges are on the third line of a CGAL selection file
  }

  // now we have to find the constrained edges
  BOOST_FOREACH(edge_descriptor ed, edges(sm))
  {
    halfedge_descriptor hd = halfedge(ed, sm);
    halfedge_descriptor hop = opposite(hd, sm);
    put(himap, hd, -1);
    put(himap, hop, -1);
    put(ecmap, ed, false);
    if(is_border(ed, sm))
      continue;

    if(get(ccmap, face(hd, sm)) != get(ccmap, face(hop, sm)))
    {
      partition_edges.push_back(ed);
      put(ecmap, ed, true);
      put(himap, hd, 0);
      put(himap, hop, 1);

      if(dump)
        out << get(indices, source(ed, sm)) << " " << get(indices, target(ed, sm)) << " ";
    }
  }

  if(dump)
    out << std::endl;

  if(verbose)
  {
    std::cerr << t.time() << " sec.\n";
    t.reset();
  }

  std::vector<tbb::concurrent_vector<edge_descriptor> > cc_edges(ncc);
  BOOST_FOREACH(face_descriptor fd, faces(sm))
  {
    std::size_t cc = ccmap[fd]; // this is the partition
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(fd, sm), sm))
    {
      halfedge_descriptor hop = opposite(hd, sm);
      if(is_border(hop, sm) || get(ecmap, edge(hd, sm)) || (fd < face(hop, sm)))
      {
        cc_edges[cc].push_back(edge((hd<hop) ? hd : hop, sm));
      }
    }
  }

  if(verbose)
  {
    for(std::size_t i=0; i<ncc; ++i)
      std::cout << "#edges in component "<< i << ": " << cc_edges[i].size() << std::endl;
  }

  std::vector<int> buffer_size(ncc); // the number of constrained edges inside each component
  std::vector<int> removed(ncc); // the number of collapsed edges in each thread
  if(verbose)
  {
    std::cerr << "#CC = " << ncc << " in " << t.time() << " sec." << std::endl;
    t.reset();
    t.start();
  }

  size_type initial_num_edges = num_edges(sm);
#if 1
  // Simplify the partition in parallel
  CGAL_MUTEX removal_mutex;
  tbb::task_group tasks;
  for(std::size_t i=0; i<ncc; ++i)
  {
    tasks.run(internal::Simplify<TriangleMesh,Placement,Cost,Stop,HIMap,ECMap,UECMap,CCMap>(
            sm, himap, ecmap, uecmap, ccmap, placement, cost, stop, buffer_size, removed,
            cc_edges[i], i, layers, dump, verbose, increase,  &removal_mutex));
  }

  tasks.wait();
#else
  for(std::size_t i=0; i<ncc; ++i)
  {
    tbb::task_group tasks;
    tasks.run(internal::Simplify<TriangleMesh,Placement,Cost,Stop,HIMap,ECMap,UECMap,CCMap>(
                sm, himap, ecmap, uecmap, ccmap, placement, cost, stop, buffer_size, removed,
                cc_edges[i], i, layers, dump, verbose, increase, NULL));
    tasks.wait();
  }
#endif

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

    expand_edge_selection(partition_edges,
                          sm,
                          layers+1,
                          ecmap,
                          Counting_output_iterator(&num_unconstrained_edges));
    num_unconstrained_edges += partition_edges.size();

    overlapping_unconstrained_edges = num_unconstrained_edges - partition_edges.size();
    for(std::size_t i=0; i < buffer_size.size(); ++i)
      overlapping_unconstrained_edges -= buffer_size[i];
  }
  else
  {
    num_unconstrained_edges = partition_edges.size();
    for(std::size_t i=0; i<buffer_size.size(); ++i)
    {
      if(verbose)
      {
        std::cout << "# partition edges = "<< partition_edges.size()
                  <<   " " << buffer_size[i] << std::endl;
      }

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

  {
    int index = 0;
    typedef typename boost::property_map<TriangleMesh,halfedge_index_t>::type HIMap;
    HIMap him = get(halfedge_index, sm);
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(sm))
    {
      put(him,hd,index++);
    }

    index = 0;
    typedef typename boost::property_map<TriangleMesh,boost::vertex_index_t>::type VIMap;
    VIMap vim = get(boost::vertex_index,sm);
    BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd, vertices(sm))
    {
      put(vim,vd,index++);
    }
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
                            parameters::vertex_index_map(get(boost::vertex_index, sm))
                            .current_num_edges(current_num_edges)
                            .get_placement(constrained_placement)
                            .get_cost(cost)
                            .edge_is_constrained_map(ormap));
  }
  else
  {
    result += edge_collapse(sm,
                            stop,
                            Sequential_tag(),
                            NULL,
                            parameters::vertex_index_map(get(boost::vertex_index, sm))
                            .current_num_edges(current_num_edges)
                            .get_placement(placement)
                            .get_cost(cost)
                            .edge_is_constrained_map(ormap));
  }

  CGAL::internal::remove_property(ecmap, sm);
  CGAL::internal::remove_property(himap, sm);

  if(verbose)
  {
    std::cerr << "sequential edge collapse on buffer in " << t.time() << " sec." << std::endl;
    std::cerr << "After sequential pass"<< std::endl;
    std::cerr << "#removed edges = " << result << std::endl;
  }

  return result;
}

} // end namespace Surface_mesh_simplification
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_PARALLEL_EDGE_COLLAPSE_H
