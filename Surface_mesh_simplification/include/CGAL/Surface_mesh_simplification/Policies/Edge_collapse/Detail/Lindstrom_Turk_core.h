// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
// Implementation of the collapsing cost and placement strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H 1

#include <CGAL/license/Surface_mesh_simplification.h>

#include <vector>

#include <CGAL/Cartesian_converter.h>

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_params.h>

namespace CGAL {

namespace Surface_mesh_simplification {

template<class ECM_, class Profile_>
class LindstromTurkCore
{
public:
  typedef ECM_                                                           ECM;
  typedef Profile_                                                       Profile;

  typedef boost::graph_traits<ECM>                                       GraphTraits;

  typedef typename GraphTraits::vertex_descriptor                        vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor                      halfedge_descriptor;

  typedef LindstromTurk_params                                           Params;

  typedef typename Profile::Point                                        Point;

  typedef typename Profile::VertexPointMap                               Vertex_point_pmap;
  typedef typename boost::property_traits<Vertex_point_pmap>::value_type ECM_Point;

  typedef typename Kernel_traits<ECM_Point>::Kernel                      ECM_Kernel;

  typedef typename Kernel_traits<Point>::Kernel                          Kernel;
  typedef typename Kernel::Vector_3                                      Vector;
  typedef typename Kernel::FT                                            FT;

  typedef optional<FT>                                                   Optional_FT;
  typedef optional<Point>                                                Optional_point;
  typedef optional<Vector>                                               Optional_vector;

  typedef MatrixC33<Kernel>                                              Matrix;

  typedef typename Profile::Triangle                                     Triangle;
  typedef typename Profile::vertex_descriptor_vector                     vertex_descriptor_vector;

  typedef typename Profile::Triangle_vector::const_iterator              const_triangle_iterator;
  typedef typename Profile::halfedge_descriptor_vector::const_iterator   const_border_edge_iterator;

public:

  LindstromTurkCore(const Params& aParams, const Profile& aProfile);

  Optional_point compute_placement();
  Optional_FT compute_cost(const Optional_point& p);

private :
  struct Triangle_data
  {
    Triangle_data(const Vector& aNormalV, const FT& aNormalL)
      : NormalV(aNormalV), NormalL(aNormalL)
    { }

    Vector NormalV;
    FT NormalL;
  };

  struct Boundary_data
  {
    Boundary_data(Point s_, Point t_, const Vector& v_, const Vector& n_)
      : s(s_), t(t_), v(v_), n(n_)
    { }

    Point s, t;
    Vector v, n;
  };

  typedef std::vector<Triangle_data>       Triangle_data_vector;
  typedef std::vector<Boundary_data>       Boundary_data_vector;

private :
  void extract_triangle_data();
  void extract_boundary_data();

  void add_boundary_preservation_constraints(const Boundary_data_vector& aBdry);
  void add_volume_preservation_constraints(const Triangle_data_vector& aTriangles);
  void add_boundary_and_volume_optimization_constraints(const Boundary_data_vector& aBdry,
                                                        const Triangle_data_vector& aTriangles);
  void add_shape_optimization_constraints(const vertex_descriptor_vector& aLink);

  FT compute_boundary_cost(const Vector& v, const Boundary_data_vector& aBdry);
  FT compute_volume_cost(const Vector& v, const Triangle_data_vector& aTriangles);
  FT compute_shape_cost(const Point& p, const vertex_descriptor_vector& aLink);

  Point get_point(const vertex_descriptor& v) const
  {
    return convert(get(mProfile.vertex_point_map(), v));
  }

  static Vector point_cross_product(const Point& a, const Point& b)
  {
    return cross_product(a-ORIGIN, b-ORIGIN);
  }

  // This is the (uX)(Xu) product described in the Lindstrom-Turk paper
  static Matrix LT_product(const Vector& u)
  {
    FT a00 = (u.y()*u.y()) + (u.z()*u.z());
    FT a01 = -u.x()*u.y();
    FT a02 = -u.x()*u.z();

    FT a10 = a01;
    FT a11 = (u.x()*u.x()) + (u.z()*u.z());
    FT a12 = - u.y()*u.z();

    FT a20 = a02;
    FT a21 = a12;
    FT a22 = (u.x()*u.x()) + (u.y()*u.y());

    return Matrix(a00,a01,a02,
                  a10,a11,a12,
                  a20,a21,a22);
  }

  static FT big_value() { return static_cast<FT>((std::numeric_limits<double>::max)()); }

  static bool is_finite(const FT& n) { return CGAL_NTS is_finite(n); }
  static bool is_finite(const Point& p) { return is_finite(p.x()) && is_finite(p.y()) && is_finite(p.z()); }
  static bool is_finite(const Vector& v) { return is_finite(v.x()) && is_finite(v.y()) && is_finite(v.z()); }
  static bool is_finite(const Matrix& m) { return is_finite(m.r0()) && is_finite(m.r1()) && is_finite(m.r2()); }

  template<class T>
  static optional<T> filter_infinity(const T& n) { return is_finite(n) ? optional<T>(n) : optional<T>(); }

  ECM& surface() const { return mProfile.surface(); }

private:
  const Params& mParams;
  const Profile& mProfile;

  void add_constraint_if_alpha_compatible(const Vector& Ai, const FT& bi);
  void add_constraint_from_gradient (const Matrix& H, const Vector& c);

private:
  Triangle_data_vector mTriangle_data;
  Boundary_data_vector mBdry_data;

  int mConstraints_n;
  Matrix mConstraints_A;
  Vector mConstraints_b;

  Cartesian_converter<ECM_Kernel, Kernel> convert;

  FT mSquared_cos_alpha;
  FT mSquared_sin_alpha;
};

} // end namespace Surface_mesh_simplification

} // end namespace CGAL

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H
