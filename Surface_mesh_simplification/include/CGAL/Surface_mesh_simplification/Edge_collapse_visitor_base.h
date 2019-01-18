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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_VISITOR_BASE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_VISITOR_BASE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace CGAL {

namespace Surface_mesh_simplification {

template<class TriangleMesh_>
struct Edge_collapse_visitor_base
{
  typedef TriangleMesh_                                                  TriangleMesh;
  typedef Edge_profile<TriangleMesh>                                     Profile;
  typedef boost::graph_traits<TriangleMesh>                              GraphTraits;

  typedef typename GraphTraits::edges_size_type                          size_type;
  typedef typename GraphTraits::vertex_descriptor                        vertex_descriptor;

  typedef typename boost::property_map<TriangleMesh,
                                       CGAL::vertex_point_t>::type       Vertex_point_pmap;
  typedef typename boost::property_traits<Vertex_point_pmap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel                          Kernel;
  typedef typename Kernel::FT                                            FT;

  template <class TM>
  void OnStarted(TM&) {}

  template <typename StopPredicate>
  void OnParallelPassFinished(TriangleMesh&, StopPredicate&, size_type /* initial */ , size_type /* current */) const {}

  template <class TM>
  void OnFinished(TM&) {}
  template <class Profile_>
  void OnStopConditionReached(const Profile_&) {}
  template <class Profile_>
  void OnCollected(const Profile_&, const boost::optional<FT>&) {}
  template <class Profile_>
  void OnSelected(const Profile_&, const boost::optional<FT>&, size_type /* initial */, size_type /* current */) {}
  template <class Profile_>
  void OnCollapsing(const Profile_&, const boost::optional<Point>&) {}
  template <class Profile_>
  void OnCollapsed(const Profile_&, const vertex_descriptor&) {}
  template <class Profile_>
  void OnNonCollapsable(const Profile_&) {}
};

} // end namespace Surface_mesh_simplification

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_EDGE_COLLAPSE_VISITOR_BASE_H
