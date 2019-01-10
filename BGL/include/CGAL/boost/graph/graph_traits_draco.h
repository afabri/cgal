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

namespace CGAL {

  struct dMesh {
    typedef CGAL::Simple_cartesian<double>::Point_3 Point_3;  // AF: maybe a template parameter?
    
  draco::Mesh& mesh;
  draco::CornerTable& ctable;

    dMesh(draco::Mesh& mesh,
          draco::CornerTable& ctable)
      : mesh(mesh), ctable(ctable)
      {}
  };

  
  struct dMesh_vertex_pmap {
    const dMesh& dm;
    const draco::PointAttribute *const att;
    
    dMesh_vertex_pmap(const dMesh& dm)
      : dm(dm), att(dm.mesh.GetNamedAttribute(draco::GeometryAttribute::POSITION))
    {}

    friend dMesh::Point_3 get(const dMesh_vertex_pmap& pm, draco::VertexIndex vi)
    {
      draco::AttributeValueIndex i(vi.value()); 
      std::array<float, 3> value;
      pm.att->ConvertValue<float, 3>(i, &value[0]);
      return dMesh::Point_3(value[0], value[1], value[2]);
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

  template <>
  struct graph_traits<CGAL::dMesh> {
    typedef draco::VertexIndex vertex_descriptor;
    typedef draco::FaceIndex face_descriptor;
    struct halfedge_descriptor {
      halfedge_descriptor()
      {}
      halfedge_descriptor(const halfedge_descriptor& other)
        : ci(other.ci), border(other.border)
      {}
      
      halfedge_descriptor(draco::CornerIndex ci)
        : ci(ci), border(false)
      {}

      bool operator==(const halfedge_descriptor& other) const
      {
        return ci == other.ci && border == other.border;
      }
      
      draco::CornerIndex ci;
      bool border;
    };

    typedef boost::counting_iterator<vertex_descriptor, std::random_access_iterator_tag,int> vertex_iterator;
    typedef boost::counting_iterator<face_descriptor, std::random_access_iterator_tag,int> face_iterator;

    typedef int vertices_size_type;
    typedef int edges_size_type;
    typedef int halfedges_size_type;
    typedef int faces_size_type;
  };

  
  template <>
  struct property_map<CGAL::dMesh, CGAL::vertex_point_t >
  {
    typedef CGAL::dMesh_vertex_pmap type;
    typedef type const_type;
  };

  
  inline
  CGAL::dMesh_vertex_pmap
  get(boost::vertex_point_t, const CGAL::dMesh& dm)
  {
    return CGAL::dMesh_vertex_pmap(dm);
  }

  
  inline
  CGAL::dMesh_vertexi_pmap
  get(boost::vertex_index_t, const CGAL::dMesh&)
  {
    return CGAL::dMesh_vertexi_pmap();
  }

  
  inline
  CGAL::dMesh_facei_pmap
  get(boost::face_index_t, const CGAL::dMesh&)
  {
    return CGAL::dMesh_facei_pmap();
  }
} // namespace boost

namespace CGAL {
  
  inline
  typename boost::graph_traits<dMesh>::vertices_size_type
  num_vertices(const dMesh& dm)
  {
    return dm.ctable.num_vertices();
  }

  
  inline
  typename boost::graph_traits<dMesh>::vertices_size_type
  num_faces(const dMesh& dm)
  {
    return dm.ctable.num_faces();
  }

  
  //  @todo fix for border
  inline
  typename boost::graph_traits<dMesh>::halfedges_size_type
  num_halfedges(const dMesh& dm)
  {
    return dm.ctable.num_corners();
  }

  
  inline
  std::pair<boost::graph_traits<dMesh>::face_iterator,
            boost::graph_traits<dMesh>::face_iterator>
  faces(const dMesh& dm)
{
  return std::make_pair(boost::graph_traits<dMesh>::face_iterator(boost::graph_traits<dMesh>::face_descriptor(0)),
                        boost::graph_traits<dMesh>::face_iterator(boost::graph_traits<dMesh>::face_descriptor(dm.ctable.num_faces())));
}

  
  inline
  std::pair<boost::graph_traits<dMesh>::vertex_iterator,
            boost::graph_traits<dMesh>::vertex_iterator>
  vertices(const dMesh& dm)
{
  return std::make_pair(boost::graph_traits<dMesh>::vertex_iterator(boost::graph_traits<dMesh>::vertex_descriptor(0)),
                        boost::graph_traits<dMesh>::vertex_iterator(boost::graph_traits<dMesh>::vertex_descriptor(dm.ctable.num_vertices())));
}

  
  inline
  boost::graph_traits<dMesh>::halfedge_descriptor
  halfedge(boost::graph_traits<dMesh>::face_descriptor fd,
           const dMesh& dm)
  {
    return dm.ctable.FirstCorner(fd);
  }

  
  inline
  boost::graph_traits<dMesh>::halfedge_descriptor
  halfedge(boost::graph_traits<dMesh>::vertex_descriptor vd,
           const dMesh& dm)
  {
    return dm.ctable.LeftMostCorner(vd);
  }

  
 inline
  boost::graph_traits<dMesh>::halfedge_descriptor
  next(boost::graph_traits<dMesh>::halfedge_descriptor hd,
       const dMesh& dm)
  {
    return dm.ctable.Next(hd.ci);
  }

  
 inline
  boost::graph_traits<dMesh>::halfedge_descriptor
  prev(boost::graph_traits<dMesh>::halfedge_descriptor hd,
       const dMesh& dm)
  {
    return dm.ctable.Previous(hd.ci);
  }

  
  inline
  boost::graph_traits<dMesh>::halfedge_descriptor
  opposite(boost::graph_traits<dMesh>::halfedge_descriptor hd,
           const dMesh& dm)
  {
    return dm.ctable.Previous(dm.ctable.Opposite(dm.ctable.Next(hd.ci)));
  }

  
 inline
  boost::graph_traits<dMesh>::face_descriptor
  face(boost::graph_traits<dMesh>::halfedge_descriptor hd,
       const dMesh& dm)
  {
    return dm.ctable.Face(hd.ci);
  }

  
  inline
  boost::graph_traits<dMesh>::vertex_descriptor
  source(boost::graph_traits<dMesh>::halfedge_descriptor hd,
         const dMesh& dm)
  {
    return dm.ctable.Vertex(dm.ctable.Previous(hd.ci));
  }

  
  inline
  boost::graph_traits<dMesh>::vertex_descriptor
  target(boost::graph_traits<dMesh>::halfedge_descriptor hd,
         const dMesh& dm)
  {
    return dm.ctable.Vertex(hd.ci);
  }

  
  inline
  int
  degree(boost::graph_traits<dMesh>::vertex_descriptor vd,
         const dMesh& dm)
  {
    return dm.ctable.Valence(vd);
  }
  
}  // namespace CGAL

#endif // CGAL_BOOST_GRAPH_GRAPH_TRAITS_DRACO_H
