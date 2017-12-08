// Copyright (c) 2017 GeometryFactory (France). All rights reserved.
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
//
// Author(s)     : Andreas Fabri
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_PARALLEL_STOP_PREDICATE_VISITOR_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_PARALLEL_STOP_PREDICATE_VISITOR_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

namespace CGAL {

namespace Surface_mesh_simplification {

//*******************************************************************************************************************
//                     -=  parallel stopping condition predicate =-
//
// Wraps a basic stop predicate and provides the necessary visitor functions
// that are required in parallel edge collapse.
//
//*******************************************************************************************************************

template<class TriangleMesh, class VisitorBase>
class Parallel_stop_predicate_visitor
  : public VisitorBase
{
public:
  typedef typename VisitorBase::size_type size_type;

  Parallel_stop_predicate_visitor(const VisitorBase& base)
    : VisitorBase(base)
  { }

  template<class Stop_predicate>
  void OnParallelPassFinished(TriangleMesh& tm,
                              Stop_predicate& pred,
                              size_type initial_num_edges,
                              size_type num_current_edges) const
  {
    VisitorBase::OnParallelPassFinished(tm, pred, initial_num_edges, num_current_edges);
  }

  // Specialization for the predicate 'Count_ratio_stop_predicate': adapt the ratio
  // used in the second pass depending on the overall expected ratio and the results
  // of the first pass.
  void OnParallelPassFinished(TriangleMesh& tm,
                              Count_ratio_stop_predicate<TriangleMesh>& pred,
                              size_type initial_num_edges,
                              size_type num_current_edges) const
  {
    // call the wrapped visitor's function
    VisitorBase::OnParallelPassFinished(tm, pred, initial_num_edges, num_current_edges);

    pred.set_ratio((std::min)(1.0, pred.ratio() * ((double)initial_num_edges / (double)num_current_edges)));
  }
};

} // end namespace Surface_mesh_simplification

} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_PARALLEL_STOP_PREDICATE_VISITOR_H
