// Copyright (c) 2016 GeometryFactory (France).
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_EXTRACT_HULL_H
#define CGAL_POLYGON_MESH_PROCESSING_EXTRACT_HULL_H

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>
#include <CGAL/boost/graph/helpers.h>

#include <boost/foreach.hpp>
#include <map>
#include <queue>


namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

// Find any face on the hull
// Algorithm: Shoot a ray from the centroid of any face in direction of its normal.
//            Return the farthest face. Note that the ray might pass through a
//            vertical face so we need potentially shoot again with a perturbed normal
// Don't use an aabb tree as we just want to do 1 shooting, but look at all faces
template <typename Geom_traits, typename TriangleMesh, typename Vpm>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,bool>
face_on_hull(const TriangleMesh& mesh, const Vpm& vpm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_traits<Vpm>::value_type Point_3;
  typedef typename Geom_traits::Vector_3 Vector_3;
  typedef typename Geom_traits::Ray_3 Ray_3;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename Geom_traits::Triangle_3 Triangle_3;

  typedef typename Geom_traits::FT FT;
  Point_3 pmax = get(vpm, target(*(halfedges(mesh).first),mesh));

  // find a point with maximal z-coordinate
  // loop over the halfedges and not the vertices as they may be isolated
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
    Point_3 p = get(vpm,target(hd,mesh));
    if(pmax.z() < p.z()){
      pmax = p;
    }
  }
  
  //  std::cerr << "pmax = "<< pmax << std::endl;
  // find all halfedges where the target point == pmax
  // attention: they may be border halfedges
  std::vector<halfedge_descriptor> he_pmax;
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
    if( get(vpm,target(hd,mesh)) == pmax ){
      he_pmax.push_back(hd);
      }
  }
  
  //std::cerr << "|he_pmax| = " << he_pmax.size() << std::endl;
  
  // among these halfedges find the set of highest halfedges
  halfedge_descriptor emax = he_pmax.front();
  
  std::vector<halfedge_descriptor> he_emax;
  he_emax.push_back(emax);
  // std::cerr << get(vpm,source(emax, mesh)) << "  " << pmax << std::endl;
  for(typename std::vector<halfedge_descriptor>::iterator it = ++(he_pmax.begin());
      it != he_pmax.end();
      ++it){
    //std::cerr << "candidate: "<< get(vpm,source(*it, mesh)) << "  " << pmax << std::endl;
    assert( get(vpm,target(emax, mesh)) == get(vpm,target(*it, mesh)) );
    assert( get(vpm,target(emax, mesh)) == pmax );
  
    if(get(vpm,source(emax,mesh)) == get(vpm,source(*it,mesh))){
      he_emax.push_back(*it);
      //std::cerr << "continue"<< std::endl;
      continue;
    }
    Comparison_result cr = compare_slope(pmax,
                                         get(vpm,source(*it, mesh)),
                                         pmax,
                                         get(vpm,source(emax, mesh)));
    if(cr == EQUAL){
      //std::cerr << "cr == EQUAL" << std::endl;
      // keep only those that at the same time share the segment of emax
    } else if(cr == SMALLER){
      // we found a higher halfedge so reset emax and the collection
      //std::cerr << "cr == SMALLER" << std::endl;
      emax = *it;
      he_emax.clear();
      he_emax.push_back(emax);
    }
  }
  //std::cerr << "|he_emax| = " << he_emax.size() << std::endl;
  
  // he_emax now contains halfedges incident to the highest segment
  // For each of them we look at the incident face
  // and keep the face with the highest 3rd vertex
  emax = he_emax.front();
  for(typename std::vector<halfedge_descriptor>::iterator it = ++(he_emax.begin());
      it != he_emax.end();
      ++it){
    assert( get(vpm,target(*it, mesh)) == pmax );
    if(! is_border(*it,mesh)){
      vertex_descriptor n_emax = target(next(emax,mesh),mesh);
      vertex_descriptor n_it = target(next(*it,mesh),mesh);
      if(compare_slope(pmax,
                       get(vpm, n_it),
                       pmax,
                       get(vpm, n_emax)) == SMALLER){
        emax = *it;
      }
    }
  }

  bool emax_take_opposite = false;
  // We now look at the opposite halfedges of he_max
  // We do not skip the first one this time
  for(typename std::vector<halfedge_descriptor>::iterator it = he_emax.begin();
      it != he_emax.end();
      ++it){
    halfedge_descriptor oit = opposite(*it,mesh);
    halfedge_descriptor oemax = opposite(emax,mesh);
    assert( get(vpm,source(oit, mesh)) == pmax );
    if(! is_border(oit, mesh)){
      vertex_descriptor n_oemax = target(next(oemax,mesh),mesh);
      vertex_descriptor n_oit = target(next(oit,mesh),mesh);
      if(compare_slope(pmax,
                       get(vpm, n_oit),
                       pmax,
                       get(vpm, n_oemax)) == SMALLER){
        emax = *it;
        emax_take_opposite = true;
      }
    }
  }

  if(emax_take_opposite){
    emax = opposite(emax,mesh);
  }
  // As we have a volume the face we have found cannot be vertical
  /*
  std::cerr << "face\n" 
            << get(vpm, target(emax,mesh)) << std::endl
            << get(vpm, target(next(emax,mesh),mesh))  << std::endl
            << get(vpm, target(prev(emax,mesh),mesh))  << std::endl << std::endl;
  */
  Angle a = angle(get(vpm, target(emax,mesh)),
                  get(vpm, target(next(emax,mesh),mesh)),
                  get(vpm, target(prev(emax,mesh),mesh)),
                  Vector_3(0,0,1));

  return std::make_pair(face(emax,mesh),a == ACUTE);

}

} // namespace internal

//named parameters vertex_index, vertex_point_map, face_index, geom_traits
template <typename TriangleMesh, typename NamedParameters>
void extract_hull(TriangleMesh& mesh,
                  const NamedParameters& np)
{
  // extract named parameters
  typedef typename GetGeomTraits<TriangleMesh,
                                 NamedParameters>::type Geom_traits;

  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::type Vpm;

  Vpm vpm = boost::choose_param(get_param(np, vertex_point),
                                get_property_map(vertex_point, mesh));

  typedef typename GetFaceIndexMap<TriangleMesh,
                                   NamedParameters>::type FaceIndexMap;
  typedef typename GetVertexIndexMap<TriangleMesh,
                                     NamedParameters>::type VertexIndexMap;

  FaceIndexMap fim = boost::choose_param(get_param(np, face_index),
                                         get_property_map(face_index, mesh));
  VertexIndexMap vim = boost::choose_param(get_param(np, boost::vertex_index),
                                           get_property_map(boost::vertex_index, mesh));

  // typedefs
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename boost::property_traits<Vpm>::value_type Point_3;

  // for all edges given as a pair of points we store the adjacent border halfedges
  typedef std::multimap<std::pair<Point_3,Point_3>,halfedge_descriptor > Point_pair_2_halfedges;
  Point_pair_2_halfedges point_pair_2_halfedges;

  BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
    if(is_border(hd,mesh)){
      point_pair_2_halfedges.insert(std::make_pair(make_sorted_pair(get(vpm,source(hd,mesh)),
                                                                    get(vpm,target(hd,mesh))),
                                                   hd));
    }
  }

  // assert that there is one or more than two halfedges for each key
  // that is either a border, or a non-manifold situation


  std::vector<std::pair<halfedge_descriptor,halfedge_descriptor> > halfedges_to_stitch;


  // compute for each face in which connected component it is
  /// \todo use bind_maps from corefinement
  std::map<face_descriptor,std::size_t> face_cc_index_map;
  std::size_t nc = connected_components(mesh,
                                        boost::make_assoc_property_map(face_cc_index_map),
                                        parameters::face_index_map(fim));

  #ifdef CGAL_PMP_EXTRACT_HULL_DEBUG
  std::cerr << nc << " connected components\n";
  #endif

  // container to mark which connected components we have already treated in the BFS
  std::vector<bool> cc_is_treated(nc,false);

  std::queue<std::pair<face_descriptor,bool> > Q;
  // find one face on the hull
  face_descriptor fd;
  bool normal_points_outwards;
  boost::tie(fd, normal_points_outwards) = internal::face_on_hull<Geom_traits>(mesh, vpm);

  Q.push(std::make_pair(fd,normal_points_outwards));
  cc_is_treated[face_cc_index_map[fd]] = true;

  // Perform a BFS traversal of the connected components on the hull
  while(! Q.empty()){
    boost::tie(fd, normal_points_outwards) = Q.front();
    Q.pop();
    std::size_t fd_cc_index = face_cc_index_map[fd];

    std::vector<face_descriptor> cc;
    connected_component(fd, mesh, std::back_inserter(cc));
    if(! normal_points_outwards){
      reverse_face_orientations(cc, mesh);
    }
    std::vector<halfedge_descriptor> border;
    border_halfedges(cc, mesh,std::back_inserter(border));


    BOOST_FOREACH(halfedge_descriptor bhd, border){
      Point_3 p = get(vpm, source(bhd,mesh));
      Point_3 q = get(vpm, target(bhd,mesh));
      typename Point_pair_2_halfedges::iterator b,e, b2;
      boost::tie(b,e) = point_pair_2_halfedges.equal_range(make_sorted_pair(p,q));

      b2 = b;
      if(b != e){
        std::vector<halfedge_descriptor> he;
        for(;b!= e; ++b){
          if(bhd != b->second){
            he.push_back(b->second);
          }
        }
        Point_3 p1 = get(vpm,target(next(opposite(bhd,mesh),mesh),mesh));

        halfedge_descriptor nhd = he[0]; // not yet the closest
        for(std::size_t i = 1; i < he.size(); ++i){
          Point_3 p2 = get(vpm,target(next(opposite(nhd,mesh),mesh),mesh));
          Point_3 q2 = get(vpm,target(next(opposite(he[i],mesh),mesh),mesh));
          /// \todo this is not yet a filtered predicate
          if(is_in_interior_of_object<Geom_traits>(q, p, p1, p2,q2)){
            nhd = he[i];
          }
        }

        // What we have so far is outwards oriented
        // Test if the neighbor connected component must get reversed when it comes out of the queue
        normal_points_outwards = (get(vpm,source(bhd,mesh)) == get(vpm,target(nhd,mesh)));

        halfedges_to_stitch.push_back(std::make_pair(bhd, nhd));
        face_descriptor nfd = face(opposite(nhd,mesh),mesh);
        if(! cc_is_treated[face_cc_index_map[nfd]]){
          Q.push(std::make_pair(nfd,normal_points_outwards));
          cc_is_treated[face_cc_index_map[nfd]] = true;
        }
        point_pair_2_halfedges.erase(b2,e);
      }
    }
  }

  stitch_borders(mesh, halfedges_to_stitch);

  std::vector<std::size_t> cc_to_remove;
  for(std::size_t i = 0; i < cc_is_treated.size(); ++i){
    if(! cc_is_treated[i]){
      cc_to_remove.push_back(i);
    }
  }

  remove_connected_components(mesh,
                              cc_to_remove,
                              boost::make_assoc_property_map(face_cc_index_map),
                              parameters::vertex_index_map(vim));

}

template <typename TriangleMesh>
void extract_hull(TriangleMesh& mesh)
{
  extract_hull(mesh, parameters::all_default());
}

} // namespace Polygon_mesh_processing

} // namespace CGAL


#endif  // CGAL_POLYGON_MESH_PROCESSING_EXTRACT_HULL_H
