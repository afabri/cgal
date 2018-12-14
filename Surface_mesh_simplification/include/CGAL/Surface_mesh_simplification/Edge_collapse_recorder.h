// Copyright (c) 2018  GeometryFactory (France). All rights reserved.
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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_RECORDER_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_RECORDER_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

namespace CGAL {
  namespace Surface_mesh_simplification {

template <typename Mesh, typename UserVisitor>

class Edge_collapse_recorder {
  typedef typename boost::graph_traits<Mesh>::edges_size_type size_type;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<Mesh, vertex_point_t>::type Vpm;
  typedef typename boost::property_traits<Vpm>::value_type Point_3;

  typedef Edge_profile<Mesh> Profile;

  struct Record {
    vertex_descriptor v0;         // the vertex that gets removed
    vertex_descriptor v1;         // the vertex that survives the collapse_edge()
    vertex_descriptor oppa, oppb; // the vertices opposite to the edge to collapse
    Point_3 p0, p1;               // the original coordinates of the vertices 
  };

public:
  Edge_collapse_recorder(Mesh& sm) 
    : sm(sm) 
  {
    BOOST_FOREACH(vertex_descriptor vd, vertices(sm)){
      v2v[vd]=vd;
    }
  } 

private:
  Mesh& sm;
  std::vector<Record> records;
  boost::unordered_map<vertex_descriptor,vertex_descriptor> v2v;

public:
  
  class Visitor : public Edge_collapse_visitor_base<Mesh>
  {
    UserVisitor visitor;
    Edge_collapse_recorder& recorder;

  public:
    Visitor(UserVisitor visitor, Edge_collapse_recorder& s)
      : visitor(visitor), recorder(s)
    {} 
    
    // Called during the collecting phase for each edge collected.
  void OnCollected( Profile const& profile, boost::optional<double> const& opt)
    {
      visitor.OnCollected(profile,opt);
    }                
    
    // Called during the processing phase for each edge selected.
    // If cost is absent the edge won't be collapsed.
    void OnSelected(Profile const&  profile        
                    ,boost::optional<double> opt
                    ,size_type st1
                    ,size_type st2
                    )
    {
      visitor.OnSelected(profile, opt,st1,st2);
  }                
    
  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
    void OnCollapsing(Profile const& profile         
                      ,boost::optional<Point>  placement
                      )
    {
      visitor.OnCollapsing(profile,placement);
      
      if ( placement ){
        //std::cout << "collapse edge " << profile.v0_v1() << " " << << profile.v0()  << " " << profile.v1() << std::endl;
        //std::cout << profile.v0() << " will get removed "  << profile.p0() << std::endl;
        Record record;
        record.v0 = profile.v0();
        record.v1 = profile.v1();
        record.p0 = profile.p0();
        record.p1 = profile.p1();
        halfedge_descriptor hd = profile.v0_v1();
      
        record.oppa = (! is_border(hd,recorder.sm))
          ? target(next(hd, recorder.sm),recorder.sm)
          : boost::graph_traits<Mesh>::null_vertex();
      
        record.oppb = (! is_border(opposite(hd, recorder.sm),recorder.sm))
          ? target(next(opposite(hd, recorder.sm), recorder.sm),recorder.sm)
          : boost::graph_traits<Mesh>::null_vertex();
      
        recorder.records.push_back(record);
      }
    }                
    
    // Called for each edge which failed the so called link-condition,
    // that is, which cannot be collapsed because doing so would
    // turn the surface mesh into a non-manifold.
    void OnNonCollapsable( Profile const& profile)
    {
      visitor.OnNonCollapsable(profile);
    }                
    
    // Called AFTER each edge has been collapsed
    void OnCollapsed( Profile const& profile, vertex_descriptor vd )
    {
      visitor.OnCollapsed(profile, vd);
    }                
    
  }; // class Visitor 

  
  Visitor visitor(const UserVisitor& user_visitor)
  {
    return Visitor(user_visitor, *this);
  }

  
  std::size_t size() const
  {
    return records.size();
  }

  bool empty() const
  {
    return records.empty();
  }
  
  
  void undo(UserVisitor user_visitor, std::size_t n)
  {
#ifdef CGAL_SMS_DUMP_MESH    
    int i = 1; // 0 can be printed in the .cpp
#endif    
    if(n > records.size()){
      n = records.size();
    }
    while(n-- > 0){
      Record r = records.back();

      user_visitor.OnSplitting(r.v1, r.oppa, r.oppb);
      
      records.pop_back();
      halfedge_descriptor h0, h1;
      if(r.oppa != boost::graph_traits<Mesh>::null_vertex()){
        std::pair<halfedge_descriptor,bool> pa = halfedge(v2v[r.oppa], v2v[r.v1], sm);
        assert(pa.second);
        h0 = pa.first;
      } else {
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(halfedge(v2v[r.v1],sm), sm)){
          if(is_border(hd,sm)){
            h0 = hd;
            break;
          }
        }
      }
      if(r.oppb != boost::graph_traits<Mesh>::null_vertex()){
        std::pair<halfedge_descriptor,bool> pa = halfedge(v2v[r.oppb], v2v[r.v1], sm);
        assert(pa.second);
        h1 = pa.first;
      } else {
        BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(halfedge(v2v[r.v1],sm), sm)){
          if(CGAL::is_border(hd,sm)){
            h1 = hd;
            break;
          }
        }
      }
      
      halfedge_descriptor hnew = Euler::split_vertex(h1, h0, sm);
      assert(target(hnew,sm) == v2v[r.v1]);
      v2v[r.v0] = source(hnew,sm);
      if(r.oppa != boost::graph_traits<Mesh>::null_vertex()){
        Euler::split_face(prev(h0,sm),next(h0, sm),sm);
      }
      if(r.oppb != boost::graph_traits<Mesh>::null_vertex()){
        Euler::split_face(prev(h1,sm),next(h1, sm),sm);
      }
      
      boost::property_map<Mesh,boost::vertex_point_t>::type vpm = get(boost::vertex_point, sm);
      
      put(vpm, v2v[r.v1], r.p1);
      put(vpm, source(hnew, sm), r.p0);

      user_visitor.OnSplit(hnew);
      
#ifdef CGAL_SMS_DUMP_MESH    
      
      Mesh sm2;
      copy_face_graph(sm,sm2);
      std::string fn("collapse-");
      fn += boost::lexical_cast<std::string>(i++) + ".off";
      
      std::ofstream out(fn.c_str());
      out << sm2;
      
#endif
    }
  }

  void undo(const UserVisitor& user_visitor)
  {
    undo(user_visitor, records.size());
  }
  
  
}; // class Edge_collapse_recorder
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif  // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_RECORDER_H
