// Copyright (c) 2019  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_DRACO_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_DRACO_H

#include <CGAL/Simple_cartesian.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <draco/mesh/corner_table.h>
#include <draco/mesh/mesh_misc_functions.h>
#include <boost/iterator/counting_iterator.hpp>
//#include <boost/iterator/iterator_facade.hpp>

namespace CGAL {

  template <typename P>
  struct dMesh {
    typedef P Point_3;
    
  draco::Mesh& mesh;
  draco::CornerTable& ctable;

    dMesh(draco::Mesh& mesh,
          draco::CornerTable& ctable)
      : mesh(mesh), ctable(ctable)
      {}
  };

  template <typename P>
  struct dMesh_vertex_pmap {
    const dMesh<P>& dm;
    const draco::PointAttribute *const att;
    
    dMesh_vertex_pmap(const dMesh<P>& dm)
      : dm(dm), att(dm.mesh.GetNamedAttribute(draco::GeometryAttribute::POSITION))
    {}

    friend typename dMesh<P>::Point_3 get(const dMesh_vertex_pmap<P>& pm, draco::VertexIndex vi)
    {
      draco::AttributeValueIndex i(vi.value()); 
      std::array<float, 3> value;
      pm.att->ConvertValue<float, 3>(i, &value[0]);
      return typename dMesh<P>::Point_3(value[0], value[1], value[2]);
    }
  };

  struct dMesh_vertexi_pmap {
    friend int get(const dMesh_vertexi_pmap& pm, draco::VertexIndex vi)
    {
      return vi.value();
    }
  };
  
  struct dMesh_facei_pmap {
    friend int get(const dMesh_facei_pmap& pm, draco::FaceIndex fi)
    {
      return fi.value();
    }
  };
    
} // namespace CGAL

namespace boost {

  template <typename P>
  struct graph_traits<CGAL::dMesh<P> > {
    typedef typename draco::VertexIndex vertex_descriptor;
    typedef typename draco::FaceIndex face_descriptor;
    
    struct halfedge_descriptor {
      halfedge_descriptor()
      {}
      halfedge_descriptor(const halfedge_descriptor& other)
        : ci(other.ci), border(other.border)
      {}
      
      halfedge_descriptor(draco::CornerIndex ci)
        : ci(ci), border(false)
      {}

      halfedge_descriptor(draco::CornerIndex ci, bool border)
        : ci(ci), border(border)
      {}
      
      bool operator==(const halfedge_descriptor& other) const
      {
        return ci == other.ci && border == other.border;
      }
      
      draco::CornerIndex ci;
      bool border;
    };

    typedef halfedge_descriptor edge_descriptor;  // @todo fix this
    
    typedef boost::counting_iterator<vertex_descriptor, std::random_access_iterator_tag,int> vertex_iterator;
    typedef boost::counting_iterator<face_descriptor, std::random_access_iterator_tag,int> face_iterator;

    // We store an index, and a bool when the iterator is on a border halfedge

    class halfedge_iterator
      : public boost::iterator_facade<halfedge_iterator, halfedge_descriptor, std::bidirectional_iterator_tag>
    {
      halfedge_descriptor hd;
      const CGAL::dMesh<P>* dm;
      
      friend class boost::iterator_core_access;

    public:
      halfedge_iterator()
      {}

      halfedge_iterator(halfedge_descriptor hd, const CGAL::dMesh<P>& dm)
        : hd(hd), dm(&dm)
      {}
    private:
      // when the opposite halfedge is border, we keep the index and set the bool to true
      // when it is not a border we increment the index
      void increment()
      {
        halfedge_descriptor opp = opposite(hd,*dm);
        if(opp.border){
          hd = opp;
        }else{
          hd.ci += 1;
          hd.border = false;
        }
      }
      
      void decrement()
      {}
      
      bool equal(const halfedge_iterator& other) const
      {
        return hd == other.hd;
      }
      
      halfedge_descriptor& dereference() const
      {
        return const_cast<halfedge_descriptor&>(hd);
      }
    };

    
    typedef int vertices_size_type;
    typedef int edges_size_type;
    typedef int halfedges_size_type;
    typedef int faces_size_type;

    
    static vertex_descriptor   null_vertex() { return draco::kInvalidVertexIndex; }
    static face_descriptor     null_face()   { return draco::kInvalidFaceIndex; }
    static halfedge_descriptor null_halfedge()   { return draco::kInvalidCornerIndex; }
  };

  
  template <typename P>
  struct property_map<CGAL::dMesh<P>, CGAL::vertex_point_t >
  {
    typedef CGAL::dMesh_vertex_pmap<P> type;
    typedef type const_type;
  };

  

  template <typename P>
  CGAL::dMesh_vertex_pmap<P>
  get(boost::vertex_point_t, const CGAL::dMesh<P>& dm)
  {
    return CGAL::dMesh_vertex_pmap<P>(dm);
  }

  
  template <typename P>
  CGAL::dMesh_vertexi_pmap
  get(boost::vertex_index_t, const CGAL::dMesh<P>&)
  {
    return CGAL::dMesh_vertexi_pmap();
  }

  
  template <typename P>
  CGAL::dMesh_facei_pmap
  get(boost::face_index_t, const CGAL::dMesh<P>&)
  {
    return CGAL::dMesh_facei_pmap();
  }
} // namespace boost

namespace CGAL {
  
  template <typename P>
  typename boost::graph_traits< dMesh<P> >::vertices_size_type
  num_vertices(const dMesh<P>& dm)
  {
    return dm.ctable.num_vertices();
  }

  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::vertices_size_type
  num_faces(const dMesh<P>& dm)
  {
    return dm.ctable.num_faces();
  }

  
  //  @todo fix for border
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::halfedges_size_type
  num_halfedges(const dMesh<P>& dm)
  {
    return dm.ctable.num_corners();
  }

  
  template <typename P>
  std::pair<typename boost::graph_traits<dMesh<P> >::face_iterator,
            typename boost::graph_traits<dMesh<P> >::face_iterator>
  faces(const dMesh<P>& dm)
{
  return std::make_pair(typename boost::graph_traits<dMesh<P> >::face_iterator(typename boost::graph_traits<dMesh<P> >::face_descriptor(0)),
                        typename boost::graph_traits<dMesh<P> >::face_iterator(typename boost::graph_traits<dMesh<P> >::face_descriptor(dm.ctable.num_faces())));
}


  template <typename P>
  std::pair<typename boost::graph_traits<dMesh<P> >::halfedge_iterator,
            typename boost::graph_traits<dMesh<P> >::halfedge_iterator>
  halfedges(const dMesh<P>& dm)
  { 
  return std::make_pair(typename boost::graph_traits<dMesh<P> >::halfedge_iterator(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor(draco::CornerIndex(0)),dm),
                        typename boost::graph_traits<dMesh<P> >::halfedge_iterator(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor(draco::CornerIndex(dm.ctable.num_corners())),dm));
}

  
  template <typename P>
  std::pair<typename boost::graph_traits<dMesh<P> >::vertex_iterator,
            typename boost::graph_traits<dMesh<P> >::vertex_iterator>
  vertices(const dMesh<P>& dm)
{
  return std::make_pair(typename boost::graph_traits<dMesh<P> >::vertex_iterator(typename boost::graph_traits<dMesh<P> >::vertex_descriptor(0)),
                        typename boost::graph_traits<dMesh<P> >::vertex_iterator(typename boost::graph_traits<dMesh<P> >::vertex_descriptor(dm.ctable.num_vertices())));
}

  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::halfedge_descriptor
  halfedge(typename boost::graph_traits<dMesh<P> >::face_descriptor fd,
           const dMesh<P>& dm)
  {
    return dm.ctable.FirstCorner(fd);
  }

  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::halfedge_descriptor
  halfedge(typename boost::graph_traits<dMesh<P> >::vertex_descriptor vd,
           const dMesh<P>& dm)
  {
    return dm.ctable.LeftMostCorner(vd);
  }

  
 template <typename P>
 typename boost::graph_traits<dMesh<P> >::halfedge_descriptor
 next(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor hd,
      const dMesh<P>& dm)
 {
   typedef typename boost::graph_traits<dMesh<P> >::halfedge_descriptor halfedge_descriptor;
   if(! hd.border){
     return dm.ctable.Next(hd.ci);
   }
    hd = opposite(hd,dm);
    while(face(hd,dm) != boost::graph_traits<dMesh<P> >::null_face()){
      hd = opposite(prev(hd,dm),dm);
    }
    return hd;
 }
  
  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::halfedge_descriptor
  prev(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor hd,
       const dMesh<P>& dm)
  {
    typedef typename boost::graph_traits<dMesh<P> >::halfedge_descriptor halfedge_descriptor;
    if(! hd.border){
      return dm.ctable.Previous(hd.ci);
    }
    hd = opposite(hd,dm);
    while(face(hd,dm) != boost::graph_traits<dMesh<P> >::null_face()){
      hd = opposite(next(hd,dm),dm);
    }
    return hd;
  }

  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::halfedge_descriptor
  opposite(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor hd,
           const dMesh<P>& dm)
  {
    typedef typename boost::graph_traits<dMesh<P> >::halfedge_descriptor halfedge_descriptor;
    
    if(hd.border){
      return halfedge_descriptor(hd.ci,false);
    }
    draco::CornerIndex ci = dm.ctable.Previous(dm.ctable.Opposite(dm.ctable.Next(hd.ci)));
    if(ci == draco::kInvalidCornerIndex){
      return halfedge_descriptor(hd.ci,true);
    }
    return halfedge_descriptor(ci,false);
  }

  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::face_descriptor
  face(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor hd,
       const dMesh<P>& dm)
  {
    if(hd.border){
      return draco::kInvalidFaceIndex;
    }
    return dm.ctable.Face(hd.ci);
  }

  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::vertex_descriptor
  source(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor hd,
         const dMesh<P>& dm)
  {
    if(hd.border){
      return dm.ctable.Vertex(hd.ci);
    }
    return dm.ctable.Vertex(dm.ctable.Previous(hd.ci));
  }

  
  template <typename P>
  typename boost::graph_traits<dMesh<P> >::vertex_descriptor
  target(typename boost::graph_traits<dMesh<P> >::halfedge_descriptor hd,
         const dMesh<P>& dm)
  {
    if(hd.border){
      return dm.ctable.Vertex(dm.ctable.Previous(hd.ci)); 
    }
    return dm.ctable.Vertex(hd.ci);
  }

  
  template <typename P>
  int
  degree(typename boost::graph_traits<dMesh<P> >::vertex_descriptor vd,
         const dMesh<P>& dm)
  {
    return dm.ctable.Valence(vd);
  }
  
}  // namespace CGAL

#endif // CGAL_BOOST_GRAPH_GRAPH_TRAITS_DRACO_H
