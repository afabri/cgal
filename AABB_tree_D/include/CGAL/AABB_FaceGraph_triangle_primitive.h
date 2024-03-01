// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
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
//
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_AABB_FACEGRAPH_TRIANGLE_PRIMITIVE_H
#define CGAL_AABB_FACEGRAPH_TRIANGLE_PRIMITIVE_H

#include <CGAL/AABB_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3_property_map.h>
#include <CGAL/Default.h>

namespace CGAL {

/*!
 * \ingroup PkgAABB_tree
 * Primitive type for a facet of a polyhedral surface.
 * It wraps a handle to a facet of a polyhedron to a 3D triangle.
 * The polyhedron from which the primitive is built should not be deleted
 * while the AABB tree holding the primitive is in use.
 *
 * \cgalModels `AABBPrimitiveWithSharedData`
 *
 *\tparam FaceGraph is a \cgal Polyhedron.
 *\tparam VertexPointPMap must be set to `CGAL::Default`
 *        This parameter is useless for the moment and will be useful in an upcoming release of \cgal.
 *\tparam OneFaceGraphPerTree must be set to `CGAL::Default`
 *        This parameter is useless for the moment and will be useful in an upcoming release of \cgal.
 *\tparam CacheDatumTag is either `CGAL::Tag_true` or `CGAL::Tag_false`. In the former case, the datum is stored
 *        in the primitive, while in the latter it is constructed on the fly to reduce the memory footprint.
 *        The default is `CGAL::Tag_false` (datum is not stored).
 *\sa `AABBPrimitive`
 *\sa `AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,ExternalPropertyMaps,CacheDatumTag>`
 *\sa `AABB_HalfedgeGraph_segment_primitive<HalfedgeGraph,OneHalfedgeGraphPerTree,CacheDatumTag>`
 */
template < class FaceGraph,
           class VertexPointPMap = Default,
           class OneFaceGraphPerTree = Default,
           class CacheDatumTag=Tag_false >
class AABB_FaceGraph_triangle_primitive
#ifndef DOXYGEN_RUNNING
: public AABB_primitive<  typename boost::mpl::if_<
                            typename boost::is_const<FaceGraph>::type,
                            typename FaceGraph::Facet_const_handle,
                            typename FaceGraph::Facet_handle
                           >::type,
                         Triangle_from_facet_handle_property_map<FaceGraph>,
                         One_point_from_facet_handle_property_map<FaceGraph>,
                         Tag_true,
                         CacheDatumTag >
#endif
{
  typedef typename boost::mpl::if_<
                          typename boost::is_const<FaceGraph>::type,
                          typename FaceGraph::Facet_const_handle,
                          typename FaceGraph::Facet_handle >::type  Id_;
  typedef Triangle_from_facet_handle_property_map<FaceGraph>  Triangle_property_map;
  typedef One_point_from_facet_handle_property_map<FaceGraph> Point_property_map;

  typedef AABB_primitive< Id_,
                          Triangle_property_map,
                          Point_property_map,
                          Tag_true,
                          CacheDatumTag > Base;

public:
  #ifdef DOXYGEN_RUNNING
  /// \name Types
  /// @{
  /*!
  The point type.
  */
  typedef boost::property_traits< boost::property_map< FaceGraph, vertex_point_t>::type >::value_type Point;
  /*!
  Geometric data type.
  */
  typedef Kernel_traits<Point>::Kernel::Triangle_3 Datum;
  /*!
  Id type.
  */
  typedef boost::graph_traits<FaceGraph>::face_descriptor Id;
  /// @}
  #endif

  // constructors
  /*!
    \tparam Iterator an input iterator with `Id` as value type.
    Constructs a primitive.
  */
  template <class Iterator>
  AABB_FaceGraph_triangle_primitive(Iterator it, FaceGraph& graph)
    : Base( Id_(it),
            Triangle_property_map(&graph),
            Point_property_map(&graph) ){}

  /// For backward-compatibility with AABB_polyhedron_triangle_primitive only.
  /// `Id_` is `Facet_const_handle` if `FaceGraph` is const and `Facet_handle` otherwise.
  AABB_FaceGraph_triangle_primitive(Id_ id)
    : Base( id,
            Triangle_property_map(NULL),
            Point_property_map(NULL) ){}

  static typename Base::Shared_data construct_shared_data( FaceGraph& graph )
  {
    return Base::construct_shared_data(Triangle_property_map(&graph), Point_property_map(&graph));
  }

  /// For backward-compatibility with AABB_polyhedron_triangle_primitive only
  static typename Base::Shared_data construct_shared_data()
  {
    return Base::construct_shared_data(Triangle_property_map(NULL), Point_property_map(NULL));
  }

};

}  // end namespace CGAL

#endif // CGAL_AABB_FACEGRAPH_TRIANGLE_PRIMITIVE_H

