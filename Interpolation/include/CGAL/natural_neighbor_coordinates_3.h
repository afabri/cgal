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
#include <CGAL/convex_hull_3.h>

//#define CGAL_NN3_DUMP_OFF
#ifdef CGAL_NN3_DUMP_OFF
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#endif

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

template <typename PolygonMesh>
double
polytope_volume(const PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<PolygonMesh,vertex_point_t>::type VPM;
  VPM vpm = get(vertex_point,pm);
  typedef typename boost::property_traits<VPM>::value_type Point_3;
  Point_3 origin(0, 0, 0);

  typename CGAL::Kernel_traits<Point_3>::Kernel::Compute_volume_3 cv3;

  double volume=0;
  BOOST_FOREACH(face_descriptor f, faces(pm)){
    halfedge_descriptor h = halfedge(f, pm);
    vertex_descriptor v = source(h, pm);
    h = next(h, pm);
    vertex_descriptor tv = target(h, pm);
    while(tv != v){
      volume += to_double(cv3(origin,
                              get(vpm, v),
                              get(vpm, source(h,pm)),
                              get(vpm, target(h,pm))));
      h = next(h, pm);
      tv = target(h, pm);
    }
  }
  return volume;
}

  template <typename P>
  double
  volum(const P& p0, const P& p1, const P& p2, const P& p3)
  {
    static int i2=0;
    std::string fn("tet-");
    fn = fn + boost::lexical_cast<std::string>(i2++) + std::string(".off");
    std::ofstream mesh(fn.c_str());
    mesh << "OFF\n4 4 0\n" << p0  << std::endl<< p1  << std::endl<< p2  << std::endl<< p3  << std::endl;
    mesh << "3 0 1 2\n";
    mesh << "3 1 0 3\n";
    mesh << "3 2 1 3\n";
    mesh << "3 0 2 3\n";
    
    return CGAL::volume(p0,p1,p2,p3);
  }

  
template <class Dt, class OutputIterator>
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_3(const Dt& dt,
                               const typename Dt::Geom_traits::Point_3& q,
                               OutputIterator out,
                               typename Dt::Cell_handle start = typename Dt::Cell_handle())
{
  typedef typename Dt::Geom_traits AK;
  
  typedef typename AK::FT      Coord_type;
  typedef typename AK::Point_3 Point_3;

  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Facet Facet;
  typedef typename Dt::Edge Edge;
  typedef typename Dt::Facet_circulator Facet_circulator;
  typedef typename Dt::Cell_circulator Cell_circulator;

  typename Dt::Locate_type lt;
  int li, lj;
  Cell_handle ch = dt.locate(q, lt, li, lj, start);

  if (lt == Dt::OUTSIDE_AFFINE_HULL || lt == Dt::OUTSIDE_CONVEX_HULL)
  {
    return make_triple(out, Coord_type(1), false);
  }
  
  if (lt == Dt::VERTEX){
    *out++= std::make_pair(ch->vertex(li), Coord_type(1));
    return make_triple(out, Coord_type(1), true);
  }
  
  if((lt == Dt::EDGE) && dt.is_on_convex_hull(Edge(ch,li, lj))){
    std::cout << "to be done: interpolate on edge on convex hull" << std::endl;
    return make_triple(out, Coord_type(1), false);
  }

  if((lt == Dt::FACET) && dt.is_on_convex_hull(Facet(ch,li))){
    std::cout << "to be done: interpolate on face on convex hull" << std::endl;
    return make_triple(out, Coord_type(1), false);     
  }
  
  std::vector<Facet> facets;
  std::set<Cell_handle> cells;
  std::set<Cell_handle> boundary_cells;
  dt.find_conflicts(q, ch, std::back_inserter(facets), std::inserter(cells, cells.end()));

  std::map<Vertex_handle,std::list<Facet> > vertex_facets;
  std::map<Vertex_handle,double> vertex_volume;
  std::set<std::pair<Vertex_handle,Vertex_handle> > boundary_edges;
  
  Point_3 origin(CGAL::ORIGIN);
  
  // For each vertex we collect its umbrella (= incident facets)
  // And we collect all edges on the boundary
  BOOST_FOREACH(Facet f, facets){
    Vertex_handle v0 = f.first->vertex(Dt::vertex_triple_index(f.second,0));
    Vertex_handle v1 = f.first->vertex(Dt::vertex_triple_index(f.second,1));
    Vertex_handle v2 = f.first->vertex(Dt::vertex_triple_index(f.second,2));
    vertex_facets[v0].push_back(f);
    vertex_facets[v1].push_back(f);
    vertex_facets[v2].push_back(f);
    boundary_edges.insert((v0<v1)?std::make_pair(v0,v1):std::make_pair(v1,v0));
    boundary_edges.insert((v1<v2)?std::make_pair(v1,v2):std::make_pair(v2,v1));
    boundary_edges.insert((v0<v2)?std::make_pair(v0,v2):std::make_pair(v2,v0));
  }

  for(std::map<Vertex_handle,std::list<Facet> >::iterator it = vertex_facets.begin();
      it != vertex_facets.end();
      ++it){
    Vertex_handle vh = it->first;
    vertex_volume[vh] = 0;
    std::list<Facet>& list = it->second;
    std::list<Facet> ordered;
    std::map<Vertex_handle, Facet> halfedge_face_map;

    // Order the faces in the umbrella
    BOOST_FOREACH(Facet f, list){
      int vhi=0;
      for(int i=0; i<3;i++){
        if(f.first->vertex(Dt::vertex_triple_index(f.second,i)) == vh){
          vhi = i;
          break;
        }
      }
      
      halfedge_face_map[f.first->vertex(Dt::vertex_triple_index(f.second, Dt::cw(vhi)))] = f; 
    }
    
    Facet f = list.front();
    ordered.push_back(f);
    for(;;){
      int vhi=0;
      for(int i=0; i<3;i++){
        if(f.first->vertex(Dt::vertex_triple_index(f.second,i)) == vh){
          vhi = i;
          break;
        }
      }
      f = halfedge_face_map[f.first->vertex(Dt::vertex_triple_index(f.second, Dt::ccw(vhi)))];
      if(f == ordered.front()){
        break;
      }
      ordered.push_back(f);
    }

    // Compute the volume contribution of the Voronoi face
    {
      double volume = 0;
      std::vector<Point_3> fpoints;
      BOOST_FOREACH (Facet f, ordered){
        fpoints.push_back(circumcenter(q,
                                       f.first->vertex(Dt::vertex_triple_index(f.second,0))->point(),
                                       f.first->vertex(Dt::vertex_triple_index(f.second,1))->point(),
                                       f.first->vertex(Dt::vertex_triple_index(f.second,2))->point()));
      }

      for(int i=1, last= fpoints.size()-1; i < last; ++i){
        double partial_volume =  CGAL::volume(fpoints[0], fpoints[i], fpoints[i+1], vh->point());
        CGAL_assertion(partial_volume > 0);
        volume += partial_volume;
      }

      // todo: check/fix sign
      vertex_volume[vh] -= volume;
    }
      
#ifdef CGAL_NN3_DUMP_OFF
    ordered.push_back(ordered.front());
    std::cout << ordered.size()<< " " ;
    BOOST_FOREACH (Facet f, ordered){
      std::cout << " " << circumcenter(q,
                                       f.first->vertex(Dt::vertex_triple_index(f.second,0))->point(),
                                       f.first->vertex(Dt::vertex_triple_index(f.second,1))->point(),
                                       f.first->vertex(Dt::vertex_triple_index(f.second,2))->point());
    }
    std::cout << std::endl;
#endif
    
    // For each radial edge of the umbrella compute the clipped Voronoi face
    BOOST_FOREACH (Facet f, ordered){
      Cell_handle start = f.first->neighbor(f.second);
      int vhi=0;
      for(int i=0; i<3;i++){
        if(f.first->vertex(Dt::vertex_triple_index(f.second,i)) == vh){
          vhi = i;
          break;
        }
      }
      Vertex_handle vs = f.first->vertex(Dt::vertex_triple_index(f.second, Dt::cw(vhi)));
      int vhs = start->index(vs);
      int vht = start->index(vh);

      std::vector<Point_3> fpoints;
      Point_3 vp = circumcenter(q,
                                f.first->vertex(Dt::vertex_triple_index(f.second,0))->point(),
                                f.first->vertex(Dt::vertex_triple_index(f.second,1))->point(),
                                f.first->vertex(Dt::vertex_triple_index(f.second,2))->point());
      fpoints.push_back(vp);

      Facet_circulator fc = dt.incident_facets(start, vhs, vht), done(fc);
      assert(*fc == dt.mirror_facet(f));
      ++fc;
      do {
        CGAL_assertion(cells.find(fc->first) != cells.end());
        vp = circumcenter(fc->first->vertex(0)->point(),
                          fc->first->vertex(1)->point(),
                          fc->first->vertex(2)->point(),
                          fc->first->vertex(3)->point());
        fpoints.push_back(vp);
        ++fc;
        if(cells.find(fc->first) == cells.end()){
          break;
        }
      }while(fc != done);

      --fc;
      Facet ef = *fc;
      assert(std::find(facets.begin(), facets.end(), ef) != facets.end());
      fpoints.push_back(circumcenter(q,
                                     ef.first->vertex(Dt::vertex_triple_index(ef.second,0))->point(),
                                     ef.first->vertex(Dt::vertex_triple_index(ef.second,1))->point(),
                                     ef.first->vertex(Dt::vertex_triple_index(ef.second,2))->point()));

      // Compute the volume contribution of the Voronoi face to its two Voronoi cells
      {
        double volume = 0;
        for(int i = 1, last=fpoints.size()-1; i < last;i++){
          double partial_volume = CGAL::volume(fpoints[0], fpoints[i], fpoints[i+1], vh->point());
          CGAL_assertion(partial_volume > 0);
          volume += partial_volume;
        }
        
        vertex_volume[vh] += volume;
        
      }
    
#ifdef CGAL_NN3_DUMP_OFF      
      fpoints.push_back(fpoints.front());

      std::cout << fpoints.size();
      for(int i=0; i< fpoints.size(); i++){
        std::cout << " " << fpoints[i];
      }
      std::cout << std::endl;
#endif
    }
  }

  std::set<std::pair<Vertex_handle, Vertex_handle> > treated;
  
  // For each edge of a cell check if it is NOT on the boundary of the conflict zone
  BOOST_FOREACH(Cell_handle c, cells){
    CGAL::array<Vertex_handle,4> vertices;
    for(int i = 0; i < 4; i++){
      vertices[i]=c->vertex(i);
    }
    std::sort(vertices.begin(),vertices.end());
    
    for(int i = 0; i < 3; i++){
      for(int j = i+1; j < 4; j++){
        if(boundary_edges.find(std::make_pair(vertices[i],vertices[j])) == boundary_edges.end()){
          if(! treated.insert(std::make_pair(vertices[i],vertices[j])).second){
            continue; // we treated this edge already in another cell
          }
          int iv = c->index(vertices[i]);
          int jv = c->index(vertices[j]);
          Cell_circulator cc = dt.incident_cells(c, iv, jv), done(cc);
          std::vector<Point_3> fpoints;
          do{
            Point_3 vp = circumcenter(cc->vertex(0)->point(),
                                      cc->vertex(1)->point(),
                                      cc->vertex(2)->point(),
                                      cc->vertex(3)->point());
            fpoints.push_back(vp);
            ++cc;
          }while(cc != done);
          
          // Compute the volume contribution of the Voronoi face to its two Voronoi cells
          {
            double volume_vi = 0;
            double volume_vj = 0;
            for(int ii = 1, last=fpoints.size()-1; ii < last; ii++){
              double partial_volume_vi = CGAL::volume(fpoints[0], fpoints[ii], fpoints[ii+1], vertices[i]->point());
              double partial_volume_vj = CGAL::volume(fpoints[0], fpoints[ii], fpoints[ii+1], vertices[j]->point());
              CGAL_assertion(partial_volume_vi < 0);
              CGAL_assertion(partial_volume_vj > 0);
              volume_vi +=  partial_volume_vi;
              volume_vj +=  partial_volume_vj;
            }
            // todo: check/fix sign
            vertex_volume[vertices[i]] -= volume_vi;
            vertex_volume[vertices[j]] += volume_vj;
          }
          
#ifdef CGAL_NN3_DUMP_OFF
          fpoints.push_back(fpoints.front());
          std::cout << fpoints.size();
          for(int i=0; i< fpoints.size(); i++){
            std::cout << " " << fpoints[i];
          }
          std::cout << std::endl;
#endif
        }
      }
    }
  }

  // Write the coordinates of each vertex of the conflict zone in the output iterator
  // and sum the coordinates up
  double total_volume = 0;
  for(std::map<Vertex_handle,double>::iterator it = vertex_volume.begin();
      it != vertex_volume.end();
      ++it){
    CGAL_assertion(it->second > 0);
    *out++= std::make_pair(it->first, Coord_type(it->second));
    total_volume += it->second;
  }
  
  return make_triple(out, total_volume, true);     
}

  

#if 1
  
template <class Dt, class OutputIterator>
Triple< OutputIterator, typename Dt::Geom_traits::FT, bool >
natural_neighbor_coordinates_3V0(const Dt& dt,
                               const typename Dt::Geom_traits::Point_3& q,
                               OutputIterator out,
                               typename Dt::Cell_handle start = typename Dt::Cell_handle())
{
  typedef typename Dt::Geom_traits AK;
  
  typedef typename AK::FT      Coord_type;
  typedef typename AK::Point_3 Point_3;

  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Facet Facet;
  typedef typename Dt::Edge Edge;

  typename Dt::Locate_type lt;
  int li, lj;
  Cell_handle ch = dt.locate(q, lt, li, lj, start);

  if (lt == Dt::OUTSIDE_AFFINE_HULL || lt == Dt::OUTSIDE_CONVEX_HULL)
  {
    return make_triple(out, Coord_type(1), false);
  }
  
  if (lt == Dt::VERTEX){
    *out++= std::make_pair(ch->vertex(li), Coord_type(1));
    return make_triple(out, Coord_type(1), true);
  }
  
  if((lt == Dt::EDGE) && dt.is_on_convex_hull(Edge(ch,li, lj))){
    std::cout << "to be done: interpolate on edge on convex hull" << std::endl;
    return make_triple(out, Coord_type(1), false);
  }

  if((lt == Dt::FACET) && dt.is_on_convex_hull(Facet(ch,li))){
    std::cout << "to be done: interpolate on face on convex hull" << std::endl;
    return make_triple(out, Coord_type(1), false);     
  }
  
  std::vector<Facet> facets;
  std::set<Vertex_handle> vertices;
  std::vector<Cell_handle> cells;
  std::set<Cell_handle> boundary_cells;
  dt.find_conflicts(q, ch, std::back_inserter(facets), std::back_inserter(cells));
  std::map<Vertex_handle,std::list<Point_3> > voronoi_cell_points;
  std::array<Point_3,3> points;
  std::array<Vertex_handle,3> fvertices;

  double total_volume = 0;
  
  BOOST_FOREACH(Facet f, facets){
    boundary_cells.insert(f.first);
    for(int i = 0; i < 3; i++){
      int j = Dt::vertex_triple_index(f.second,i);
      Vertex_handle v = f.first->vertex(j);
      vertices.insert(v);
      points[i]=v->point();
      fvertices[i] = v;
    }
    Point_3 cc0 = circumcenter(points[0],points[1],points[2],q);
    Point_3 cc1 = circumcenter(points[0],points[1],points[2],f.first->vertex(f.second)->point());
    for(int i =0; i<3; i++){
      voronoi_cell_points[fvertices[i]].push_back(cc0);
      voronoi_cell_points[fvertices[i]].push_back(cc1);
    }
    voronoi_cell_points[f.first->vertex(f.second)].push_back(cc1);
  }

  BOOST_FOREACH(Cell_handle ch, cells){
    if(boundary_cells.find(ch)== boundary_cells.end()){
      Point_3 cc = circumcenter(ch->vertex(0)->point(),
                                ch->vertex(1)->point(),
                                ch->vertex(2)->point(),
                                ch->vertex(3)->point());
      for(int i=0; i<4; i++){
        voronoi_cell_points[ch->vertex(i)].push_back(cc);
      }
    }
  }
  
#ifdef CGAL_NN3_DUMP_OFF  
  std::ofstream pout("points.xyz");
  BOOST_FOREACH(Vertex_handle v, vertices){
    pout << v->point() << std::endl;
  }
  
#endif
  
  int i2 = 0;
  for(std::map<Vertex_handle,std::list<Point_3> >::iterator it = voronoi_cell_points.begin();
      it != voronoi_cell_points.end();
        ++it){
    std::list<Point_3>& list = it->second;

    Surface_mesh<Point_3> sm;
    std::list<Point_3>::iterator it2 = std::unique(list.begin(), list.end());
    convex_hull_3(list.begin(), it2, sm);

    
    double pv = polytope_volume(sm);
    total_volume += pv;
    typename Dt::Geom_traits::FT fpv(pv);
    *out++= std::make_pair(it->first,fpv);
      
#ifdef CGAL_NN3_DUMP_OFF
    CGAL::Polygon_mesh_processing::triangulate_faces(sm);    
    std::string fn("cell-");
    fn = fn + boost::lexical_cast<std::string>(i2++) + std::string(".off");
    std::ofstream mesh(fn.c_str());
    mesh << sm << std::endl;
#endif    
  }

  
  return make_triple(out, typename Dt::Geom_traits::FT(total_volume), true);
}
#endif
 
                               
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
