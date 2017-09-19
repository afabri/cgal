// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/PCA_util.h>

#include <iterator>
#include <vector>
#include <cmath>

namespace CGAL {

namespace internal {
// Fits a line to a 2D triangle set.
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default)

template < typename InputIterator,
           typename Kernel, typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel& k,                   // kernel
			       const CGAL::Dimension_tag<2>& tag,
			       const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename Kernel::Triangle_2  Triangle;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  // compute centroid
  c = centroid(first,beyond,Kernel(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., }};
  assemble_covariance_matrix_2(first,beyond,covariance,c,k,(Triangle*) NULL,tag, diagonalize_traits);
  
  // compute fitting plane
  return fitting_line_2(covariance,c,line,k,diagonalize_traits);
} // end linear_least_squares_fitting_2 for triangle set with 2D tag

template < typename InputIterator,
           typename Kernel,
	   typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
			       const CGAL::Dimension_tag<1>& tag,
			       const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename Kernel::Triangle_2 Triangle;
  typedef typename Kernel::Segment_2  Segment;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Segment> segments;  
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Triangle& t = *it;
    segments.push_back(Segment(t[0],t[1]));
    segments.push_back(Segment(t[1],t[2]));
    segments.push_back(Segment(t[2],t[0]));      
  }    
  
  return linear_least_squares_fitting_2(segments.begin(),segments.end(),line,c,tag,Kernel(),
					diagonalize_traits);
  
} // end linear_least_squares_fitting_2 for triangle set with 1D tag

template < typename InputIterator,
           typename Kernel,
	   typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond, 
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
			       const CGAL::Dimension_tag<0>& tag,
			       const DiagonalizeTraits& diagonalize_traits)
{
  // types

  typedef typename Kernel::Triangle_2 Triangle;
  typedef typename Kernel::Point_2 Point;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);
  
  std::list<Point> points;  
  for(InputIterator it = first;
      it != beyond;
      it++)
  {
    const Triangle& t = *it;
    points.push_back(Point(t[0]));
    points.push_back(Point(t[1]));
    points.push_back(Point(t[2]));      
  }    
  
  return linear_least_squares_fitting_2(points.begin(),points.end(),line,c,tag,Kernel(),
					diagonalize_traits);
  
} // end linear_least_squares_fitting_2 for triangle set with 0D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H
