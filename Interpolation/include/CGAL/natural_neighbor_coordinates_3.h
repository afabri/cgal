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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Raphaelle Chaine

#ifndef CGAL_NATURAL_NEIGHBORS_3_H
#define CGAL_NATURAL_NEIGHBORS_3_H

#include <CGAL/license/Interpolation.h>

#include <CGAL/tags.h>
#include <CGAL/iterator.h>
#include <CGAL/utility.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/Concatenate_iterator.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <iostream> //TO DO : to remove
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace CGAL {

// ====================== Geometric Traits utilities =========================================
// === Declarations

template <class Gt>
typename Gt::FT
signed_area(const typename Gt::Point_3& p, const typename Gt::Point_3& q,
            const typename Gt::Point_3& r, const typename Gt::Point_3& point_of_view,
            const Gt& gt = Gt());

// ====================== Delaunay Triangulation utilities ==========================
// === Declarations

template < class DT>
typename DT::Geom_traits::Point_3
construct_circumcenter(const typename DT::Facet& f,
                       const typename DT::Geom_traits::Point_3& Q,
                       const typename DT::Geom_traits& gt = typename DT::Geom_traits());

// ====================== Natural Neighbors Querries ==========================
// === Definitions

// Given a 3D point Q and a 3D Delaunay triangulation dt,
// the next two functions calculate the natural neighbors and coordinates of Q with regard of dt
//
// OutputIterator has value type
//        std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>
// Result :
// - An OutputIterator providing natural neighbors P_i of Q with unnormalized coordinates a_i associated to them
// - The normalizing coefficient (sum over i of the a_i)
// - A boolean specifying whether the calculation has succeeded or not

template <class Dt, class OutputIterator>
Triple< OutputIterator,  // iterator with value type std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>
        typename Dt::Geom_traits::FT,  // Should provide 0 and 1
        bool >
laplace_natural_neighbor_coordinates_3(const Dt& dt,
                                       const typename Dt::Geom_traits::Point_3& Q,
                                       OutputIterator nn_out,
                                       typename Dt::Geom_traits::FT&  norm_coeff,
                                       const typename Dt::Cell_handle start = CGAL_TYPENAME_DEFAULT_ARG Dt::Cell_handle())
{
  typedef typename Dt::Geom_traits Gt;
  typedef typename Gt::Point_3 Point;
  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Facet Facet;
  typedef typename Dt::Locate_type Locate_type;
  typedef typename Gt::FT Coord_type;

  CGAL_triangulation_precondition (dt.dimension() == 3);

  Locate_type lt;
  int li, lj;
  Cell_handle c = dt.locate( Q, lt, li, lj, start);

  if ( lt == Dt::VERTEX )
  {
    *nn_out++= std::make_pair(c->vertex(li), Coord_type(1));
    return make_triple(nn_out, norm_coeff = Coord_type(1),true);
  }
  else if (dt.is_infinite(c))
  {
    //point outside the convex-hull
    return make_triple(nn_out, Coord_type(1), false);
  }

  std::set<Cell_handle> cells;
  // To replace the forbidden access to the "in conflict" flag :
  // std::find operations on this set
  std::vector<Facet> bound_facets;
  bound_facets.reserve(32);

  // Find the cells in conflict with Q
  dt.find_conflicts(Q, c,
                    std::back_inserter(bound_facets),
                    std::inserter(cells, cells.begin()));

  std::map<Vertex_handle,Coord_type> coordinate;
  typename std::map<Vertex_handle,Coord_type>::iterator coor_it;

  typename std::vector<Facet>::iterator bound_it;
  for (bound_it = bound_facets.begin(); bound_it != bound_facets.end(); ++bound_it)
  {
    //for each facet on the boundary
    Facet f1 = *bound_it;
    Cell_handle cc1 = f1.first;
    if (dt.is_infinite(cc1))
      return make_triple(nn_out, norm_coeff=Coord_type(1), false);//point outside the convex-hull

    CGAL_triangulation_assertion_code(Cell_handle cc2 = cc1->neighbor(f1.second);)
    CGAL_triangulation_assertion(std::find(cells.begin(),cells.end(),cc1) != cells.end());//TODO : Delete
    CGAL_triangulation_assertion(std::find(cells.begin(),cells.end(),cc2) == cells.end());//TODO : Delete

    Point C_1 = construct_circumcenter<Dt>(f1, Q, dt.geom_traits());
    for(int j=1; j<4; j++)
    {
      //for each vertex P of the boundary facet
      Vertex_handle vP = cc1->vertex((f1.second+j)&3);
      Vertex_handle vR = cc1->vertex(dt.next_around_edge(f1.second,(f1.second+j)&3));

      // turn around the oriented edge vR vP
      Cell_handle cc3 = cc1;
      int num_next = dt.next_around_edge((f1.second+j)&3,f1.second);

      Cell_handle next = cc3->neighbor(num_next);
      while (std::find(cells.begin(),cells.end(),next) != cells.end())
      {
        CGAL_triangulation_assertion( next != cc1 );
        cc3 = next;
        num_next = dt.next_around_edge(cc3->index(vR),cc3->index(vP));
        next = cc3->neighbor(num_next);
      }

      Point C_3 = construct_circumcenter<Dt>(Facet(cc3,num_next), Q, dt.geom_traits());
      Point midPQ = midpoint(vP->point(),Q);
      Coord_type coor_add = signed_area<Gt>(C_3,C_1,midPQ, vP->point(), dt.geom_traits());
      ((coor_it = coordinate.find(vP)) == coordinate.end())?
            coordinate[vP] = coor_add : coor_it->second += coor_add; // Replace by a function call
    }
  } //end : for each facet on the boundary

  norm_coeff = 0;
  for (coor_it=coordinate.begin(); coor_it!=coordinate.end(); ++coor_it)
  {
    Coord_type co = coor_it->second /
      (CGAL_NTS sqrt(dt.geom_traits().compute_squared_distance_3_object()(
                                        coor_it->first->point(),Q)));
    *nn_out++ = std::make_pair(coor_it->first,co);
    norm_coeff += co;
  }
  return make_triple(nn_out, norm_coeff, true);
}

template <class Dt, class OutputIterator>
Triple< OutputIterator,  // iterator with value type std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>
        typename Dt::Geom_traits::FT,  // Should provide 0 and 1
        bool >
sibson_natural_neighbor_coordinates_3(const Dt& dt,
                                      const typename Dt::Geom_traits::Point_3& Q,
                                      OutputIterator nn_out,
                                      typename Dt::Geom_traits::FT&  norm_coeff,
                                      const typename Dt::Cell_handle start = CGAL_TYPENAME_DEFAULT_ARG Dt::Cell_handle())
{
  typedef typename Dt::Geom_traits Gt;
  typedef typename Gt::Point_3 Point;
  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Facet Facet;
  typedef typename Dt::Locate_type Locate_type;
  typedef typename Gt::FT Coord_type;

  CGAL_triangulation_precondition (dt.dimension()== 3);

  Locate_type lt;
  int li, lj;
  Cell_handle c = dt.locate( Q, lt, li, lj, start);

  if ( lt == Dt::VERTEX )
  {
    *nn_out++ = std::make_pair(c->vertex(li),Coord_type(1));
    return make_triple(nn_out,norm_coeff=Coord_type(1),true);
  }
  else if (dt.is_infinite(c))
  {
    //point outside the convex-hull
    return make_triple(nn_out, Coord_type(1), false);
  }

  std::set<Cell_handle> cells;
  typename std::set<Cell_handle>::iterator cit;
  // To replace the forbidden access to the "in conflict" flag :
  // std::find operations on this set

  // Find the cells in conflict with Q
  dt.find_conflicts(Q, c,
                    Emptyset_iterator(),
                    std::inserter(cells,cells.begin()));

  std::map<Vertex_handle,Coord_type> coordinate;
  typename std::map<Vertex_handle,Coord_type>::iterator coor_it;

  for (cit = cells.begin(); cit != cells.end(); ++cit)
  {
    // for each cell cc1 in conflict
    Cell_handle cc1 = *cit;
    CGAL_triangulation_assertion(std::find(cells.begin(),cells.end(),cc1)!=cells.end());//TODO : Delete

    if (dt.is_infinite(cc1))
      return make_triple(nn_out,norm_coeff=Coord_type(1), false);//point outside the convex-hull

    typename Dt::Geom_traits::Compute_volume_3 vol =
        dt.geom_traits().compute_volume_3_object();

    Point C1 = dt.dual(cc1);
    for(int i=0; i<4; i++)
    {
      //for each neighboring cell cc2 of cc1
      Cell_handle cc2 = cc1->neighbor(i);
      if(std::find(cells.begin(),cells.end(),cc2) == cells.end())
      {
        // cc2 outside the conflict cavity
        Point C_1 = construct_circumcenter<Dt>(Facet(cc1,i), Q, dt.geom_traits());
        for(int j=1; j<4; j++)
        {
          //for each vertex P of the boundary facet
          Vertex_handle vP = cc1->vertex((i+j)&3);//&3 in place of %4
          Vertex_handle vR = cc1->vertex(dt.next_around_edge(i,(i+j)&3));

          // turn around the oriented edge vR vP
          Cell_handle cc3 = cc1;
          int num_next = dt.next_around_edge((i+j)&3,i);
          Cell_handle next = cc3->neighbor(num_next);

          while (std::find(cells.begin(),cells.end(),next) != cells.end())
          { //next is in conflict
            CGAL_triangulation_assertion( next != cc1 );
            cc3 = next;
            num_next = dt.next_around_edge(cc3->index(vR),cc3->index(vP));
            next = cc3->neighbor(num_next);
          }
          if (dt.is_infinite(cc3))
          {
            //point outside the convex-hull
            return make_triple(nn_out,norm_coeff = Coord_type(1), false);
          }

          Point C3 = dt.dual(cc3);
          Point C_3 = construct_circumcenter<Dt>(Facet(cc3,num_next), Q, dt.geom_traits());
          Point midPQ = midpoint(vP->point(),Q);
          Point midPR = midpoint(vP->point(),vR->point());
          Coord_type coor_add = vol(C_1,C1,midPR,midPQ);
          coor_add -= vol(C_1,C_3,midPR,midPQ);
          coor_add += vol(C3,C_3,midPR,midPQ);
          ((coor_it = coordinate.find(vP)) == coordinate.end())?
                coordinate[vP] = coor_add : coor_it->second += coor_add;// Replace by a function call
        }
      }
      else // cc2 in the conflict cavity
      {
        CGAL_triangulation_assertion(std::find(cells.begin(),cells.end(),cc2)!=cells.end());//TODO : Delete
        if (dt.is_infinite(cc2))
        {
          //point outside the convex-hull
          return make_triple(nn_out,norm_coeff = Coord_type(1), false);
        }

        Point C2 = dt.dual(cc2);
        for(int j=1;j<4;j++)
        {
          //for each vertex P of the internal facet
          Vertex_handle vP=cc1->vertex((i+j)&3);
          Vertex_handle vR=cc1->vertex(dt.next_around_edge(i,(i+j)&3));
          Point midPQ = midpoint(vP->point(),Q);
          Point midPR = midpoint(vP->point(),vR->point());
          Coord_type coor_add = vol(C2,C1,midPR,midPQ);
          ((coor_it=coordinate.find(vP))==coordinate.end())?
                coordinate[vP]=coor_add : coor_it->second+=coor_add;// Replace by a function call
        }
      }
    }
  }

  norm_coeff=0;
  for (coor_it=coordinate.begin(); coor_it!=coordinate.end(); ++coor_it)
  {
    *nn_out++ = std::make_pair(coor_it->first,coor_it->second);
    norm_coeff += coor_it->second;
  }
  return make_triple(nn_out,norm_coeff,true);
}

template <typename Dt, typename InputIterator>
bool is_correct_natural_neighborhood(const Dt& /*dt*/,
                                     const typename Dt::Geom_traits::Point_3&  Q,
                                     InputIterator it_begin, InputIterator it_end,
                                     const typename Dt::Geom_traits::FT&  norm_coeff)
{
  typedef typename Dt::Geom_traits Gt;
  typedef typename Gt::FT Coord_type;
  Coord_type sum_x(0);
  Coord_type sum_y(0);
  Coord_type sum_z(0);
  InputIterator it;
  for(it = it_begin ; it != it_end ; ++it)
  {
    sum_x += it->second*(it->first->point().x());
    sum_y += it->second*(it->first->point().y());
    sum_z += it->second*(it->first->point().z());
  }
  //!!!! to be replaced by a linear combination of points as soon
  // as it is available in the kernel.
  std::cout << sum_x/norm_coeff << " "
            << sum_y/norm_coeff << " "
            << sum_z/norm_coeff << std::endl;
  return ((sum_x == norm_coeff*Q.x()) && (sum_y == norm_coeff*Q.y())
          && (sum_z == norm_coeff*Q.z()));
}


template <class Dt>
void
natural_neighbor_coordinates_3(const Dt& dt,
                               const typename Dt::Geom_traits::Point_3& q)
{
  typedef typename Dt::Geom_traits::Plane_3 Plane_3;
  typedef typename Dt::Geom_traits::Point_3 Point_3;
  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Facet Facet;
  Cell_handle ch = dt.locate(q);

  std::vector<Facet> facets;
  std::set<Vertex_handle> vertices;

  dt.find_conflicts(q, ch, std::back_inserter(facets), Emptyset_iterator());
  BOOST_FOREACH(Facet f, facets){
    for(int i = 0; i < 3; i++){
      int j = Dt::vertex_triple_index(f.second,i);
      vertices.insert(f.first->vertex(j));
    }
  }

  std::vector<Plane_3> planes;
  BOOST_FOREACH(Vertex_handle vh, vertices){
    Point_3 &p = vh->point();
    Plane_3 plane = bisector(p, q);
    planes.push_back(plane);
  }

  Surface_mesh<Point_3> sm;
  std::cout << "construct cell " << std::flush;
  CGAL::halfspace_intersection_with_constructions_3(planes.begin(), planes.end(), sm);
  std::cout << "done" << std::endl;
  CGAL::Polygon_mesh_processing::triangulate_faces(sm);
  std::ofstream mesh("cell.off");
  mesh << sm << std::endl;

  int i = 0;
  BOOST_FOREACH(Vertex_handle vh, vertices){
    Point_3 &p = vh->point();
    std::vector<Vertex_handle> incident;
    std::vector<Plane_3> planes2;
    dt.incident_vertices(vh, std::back_inserter(incident));
    BOOST_FOREACH(Vertex_handle vh, incident){
      if(dt.is_infinite(vh)){
        std::cout << "skip infinite vertex"<< std::endl;
      } else {
        Point_3 &p2 = vh->point();
        Plane_3 plane2 = bisector(p2, p);
        planes2.push_back(plane2);
      }
    }
    planes2.insert(planes2.end(), planes.begin(), planes.end());
    sm.clear();
    std::cout << "Construct cell " << i << std::flush;
    CGAL::halfspace_intersection_with_constructions_3(planes2.begin(), planes2.end(), sm);
    std::cout << "done" << std::endl;
    
    CGAL::Polygon_mesh_processing::triangulate_faces(sm);
    /*
    typedef CGAL::Concatenate_iterator<std::vector<Plane_3>::iterator,std::vector<Plane_3>::iterator> CI;
                                       
    CGAL::halfspace_intersection_3(CI(planes.end(),
                                      planes2.begin(),
                                      planes.begin()),
                                   CI(planes.end(),
                                      planes2.begin(),
                                      planes2.end(), 0),
                                   sm);
    */
    
    std::string fn("cell-");
    fn = fn + boost::lexical_cast<std::string>(i++) + std::string(".off");
    std::ofstream mesh(fn.c_str());
    mesh << sm << std::endl;
  }
}
  
                               
// ====================== Geometric Traits utilities =========================================
// === Definitions

template <class Gt>
typename Gt::FT
signed_area(const typename Gt::Point_3& p, const typename Gt::Point_3& q,
            const typename Gt::Point_3& r,
            const typename Gt::Point_3& point_of_view,
            const Gt& gt /* = Gt() */)
{
  // signed area of the triangle p q r
  return gt.compute_area_3_object()(p,q,r)
      * (gt.orientation_3_object()(p, q, r, point_of_view) == COUNTERCLOCKWISE?+1:-1);
}

// ====================== Delaunay Triangulation utilities ==========================
// === Definitions

template < class DT>
typename DT::Geom_traits::Point_3
construct_circumcenter(const typename DT::Facet& f,
                       const typename DT::Geom_traits::Point_3& Q,
                       const typename DT::Geom_traits& gt /* = typename DT::Geom_traits() */ )
{
  CGAL_triangulation_precondition(//&3 in place of %4
                                  !gt.coplanar_3_object()(
                                        f.first->vertex((f.second+1)&3)->point(),
                                        f.first->vertex((f.second+2)&3)->point(),
                                        f.first->vertex((f.second+3)&3)->point(),
                                        Q));
  // else the facet is not on the enveloppe of the conflict cavity associated to P
  return gt.construct_circumcenter_3_object()(
             f.first->vertex((f.second+1)&3)->point(),
             f.first->vertex((f.second+2)&3)->point(),
             f.first->vertex((f.second+3)&3)->point(),
             Q);
}

} //namespace CGAL

#endif // CGAL_NATURAL_NEIGHBORS_3_H
