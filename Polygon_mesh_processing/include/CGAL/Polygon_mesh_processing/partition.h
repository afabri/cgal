// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_PARTITION_H
#define CGAL_POLYGON_MESH_PROCESSING_PARTITION_H

#include <CGAL/license/Polygon_mesh_processing.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/tuple/tuple.hpp>

#include <metis.h>

#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>

namespace CGAL {

namespace Polygon_mesh_processing {

/// Output a set of partitions
///
/// \param n the number of parts
///
/// \tparam PolygonMesh is a model of the face graph concept.
/// \tparam FacePartitionIDPmap is is a model of `ReadablePropertyMap`
///           with `boost::graph_traits<PolygonMesh>::%face_descriptor`
///           as key type and `boost::face_external_index` as value type.
///
template<typename PolygonMesh, typename FacePartitionIDPmap>
void output_partitions(const PolygonMesh& m,
                       const FacePartitionIDPmap fpmap,
                       const idx_t n)
{
  typedef CGAL::Face_filtered_graph<PolygonMesh>                       Filtered_graph;

  for(int i=0; i<n; ++i)
  {
    std::ostringstream filename;
    filename << "partition_" << i << ".off" << std::ends;
    std::ofstream out(filename.str().c_str());

    Filtered_graph m_part(m, i, fpmap);

    PolygonMesh out_mesh;
    CGAL::copy_face_graph(m_part, out_mesh);

    out << out_mesh;
  }
}

/// Computes a partition of the input mesh into `nparts` roughly equal parts
///
/// \param nparts the number of parts in the final partition
///
/// \tparam PolygonMesh is a model of the `FaceGraph` concept.
/// \tparam FacePartitionIDPmap is is a model of `ReadWritePropertyMap`
///           with `boost::graph_traits<PolygonMesh>::%face_descriptor`
///           as key type and `boost::face_external_index` as value type.
///           The default is `typename boost::property_map<FaceGraph, face_external_index>::%type`.
///
/// \pre `m` is a pure triangular surface mesh: there are no edges
///       without at least one incident face
template<typename PolygonMesh, typename FacePartitionIDPmap>
void partition(const PolygonMesh& m,
               FacePartitionIDPmap partition_id_map,
               int nparts)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_iterator       face_iterator;

  idx_t options[METIS_NOPTIONS];
  options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
  options[METIS_OPTION_CONTIG] = 0;
  options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
//  options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM;
  options[METIS_OPTION_MINCONN] = 0;
  options[METIS_OPTION_NCUTS] = 3;
  options[METIS_OPTION_NITER] = 10;
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
//  options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
  options[METIS_OPTION_SEED] = 12343;
  options[METIS_OPTION_UFACTOR] = 1;

  idx_t nn = num_vertices(m);
  idx_t ne = num_faces(m);
  idx_t d = 3; // number of nodes per element
  idx_t* eptr = new idx_t[ne + 1];
  idx_t* eind = new idx_t[d * ne];

  // fill the adjacency info
  typedef typename boost::property_map<PolygonMesh,
                                       boost::vertex_index_t>::const_type Indices;
  Indices indices = get(boost::vertex_index, m);

  face_iterator fit, fe;
  boost::tie(fit, fe) = faces(m);
  for(int i=0, j=0; fit!=fe; ++fit, ++i)
  {
    eptr[i] = j;

    halfedge_descriptor h = halfedge(*fit, m), done = h;
    do
    {
      vertex_descriptor v = target(h, m);
      assert(j < d * ne);
      eind[j++] = get(indices, v);
      h = next(h, m);
    } while (h != done);

    assert(i < ne);
    eptr[i + 1] = j;
  }

  // a dual edge between elements exists if they share 'nparts' vertices
  idx_t ncommon = 2;

  // either the edgecut or the total communication volume of the dual graph’s partitioning
  idx_t objval;

  // partition info for the nodes
  idx_t* npart = (idx_t*) calloc(nn, sizeof(idx_t));
  assert(npart != NULL);

  // partition info for the elements
  idx_t* epart = (idx_t*) calloc(ne, sizeof(idx_t));
  assert(epart != NULL);

  int ret = METIS_PartMeshDual(&ne, &nn, eptr, eind,
                               NULL /* elements weights*/, NULL /*elements sizes*/,
                               &ncommon, &nparts,
                               NULL /* partitions weights */,
                               options,
                               &objval, epart, npart);

  std::cout << "return: " << ret << " with objval: " << objval << std::endl;
  assert(ret == METIS_OK);

  boost::tie(fit, fe) = faces(m);
  for(int i=0; fit!=fe; ++fit, ++i)
    put(partition_id_map, *fit, epart[i]);

  //  output_partitions(m, partition_id_map, nparts);
}

template<typename PolygonMesh>
void partition(const PolygonMesh& m, int nparts = 3)
{
  typedef typename boost::property_map<PolygonMesh,
                                       boost::face_external_index_t>::type Fpmap;
  Fpmap fpmap = get(boost::face_external_index, m);

  return partition(m, fpmap, nparts);
}

} //end namespace Polygon_mesh_processing

} //end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_PARTITION_H
