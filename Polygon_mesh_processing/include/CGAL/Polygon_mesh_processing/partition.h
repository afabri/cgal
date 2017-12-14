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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_PARTITION_H
#define CGAL_POLYGON_MESH_PROCESSING_PARTITION_H

#include <CGAL/license/Polygon_mesh_processing.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <metis.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/tuple/tuple.hpp>

#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>

namespace CGAL {

namespace Polygon_mesh_processing {

/// Output a set of partitions.
///
/// \param tm a triangle mesh
/// \param nparts the number of parts
/// \param fpmap the property map with the partition indices
///
/// \tparam TriangleMesh must be a model of a `FaceListGraph`, `HalfedgeListGraph`, and \bgllink{VertexListGraph}.
/// \tparam FacePartitionIDPmap is is a model of `ReadablePropertyMap`
///           with `boost::graph_traits<TriangleMesh>::%face_descriptor`
///           as key type and `boost::face_external_index` as value type.
template<typename TriangleMesh, typename FacePartitionIDPmap>
void output_partitions(const TriangleMesh& tm,
                       const idx_t nparts,
                       const FacePartitionIDPmap fpmap)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  typedef CGAL::Face_filtered_graph<TriangleMesh>         Filtered_graph;

  for(int i=0; i<nparts; ++i)
  {
    std::ostringstream filename;
    filename << "partition_" << i << ".off" << std::ends;
    std::ofstream out(filename.str().c_str());

    Filtered_graph m_part(tm, i, fpmap);

    TriangleMesh out_mesh;
    CGAL::copy_face_graph(m_part, out_mesh);

    out << out_mesh;
  }
}

template<typename TriangleMesh, typename FacePartitionIDPmap, typename NamedParameters>
void partition(const TriangleMesh& tm,
               int nparts, FacePartitionIDPmap partition_id_map,
               const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  CGAL_precondition_msg(nparts > 0, ("Partitioning requires a strictly positive number of parts"));

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator       face_iterator;

  using boost::get_param;
  using boost::choose_param;

  idx_t options[METIS_NOPTIONS];
  options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
  options[METIS_OPTION_CONTIG] = 0;
  options[METIS_OPTION_DBGLVL] = 0; // METIS_DBG_INFO;
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

  idx_t nn = static_cast<idx_t>(num_vertices(tm));
  idx_t ne = static_cast<idx_t>(num_faces(tm));
  idx_t d = 3; // number of nodes per element
  idx_t* eptr = new idx_t[ne + 1];
  idx_t* eind = new idx_t[d * ne];

  // fill the adjacency info
  typedef typename boost::property_map<TriangleMesh,
                                       boost::vertex_index_t>::const_type Indices;
  Indices indices = get(boost::vertex_index, tm);

  face_iterator fit, fe;
  boost::tie(fit, fe) = faces(tm);
  for(int i=0, j=0; fit!=fe; ++fit, ++i)
  {
    eptr[i] = j;

    halfedge_descriptor h = halfedge(*fit, tm), done = h;
    do
    {
      vertex_descriptor v = target(h, tm);
      CGAL_assertion(j < d * ne);
      eind[j++] = static_cast<idx_t>(get(indices, v));
      h = next(h, tm);
    } while (h != done);

    CGAL_assertion(i < ne);
    eptr[i + 1] = j;
  }

  // a dual edge between elements exists if they share 'nparts' vertices
  idx_t ncommon = 2;

  // either the edgecut or the total communication volume of the dual graph’s partitioning
  idx_t objval;

  // partition info for the nodes
  idx_t* npart = (idx_t*) calloc(nn, sizeof(idx_t));
  CGAL_assertion(npart != NULL);

  // partition info for the elements
  idx_t* epart = (idx_t*) calloc(ne, sizeof(idx_t));
  CGAL_assertion(epart != NULL);

  int ret = METIS_PartMeshDual(&ne, &nn, eptr, eind,
                               NULL /* elements weights*/, NULL /*elements sizes*/,
                               &ncommon, &nparts,
                               NULL /* partitions weights */,
                               options,
                               &objval, epart, npart);

  std::cout << "return: " << ret << " with objval: " << objval << std::endl;
  CGAL_assertion(ret == METIS_OK);

  boost::tie(fit, fe) = faces(tm);
  for(int i=0; fit!=fe; ++fit, ++i)
    put(partition_id_map, *fit, epart[i]);

  //  output_partitions(m, nparts, partition_id_map);
}

/// Computes a partition of the input mesh into `nparts` parts.
///
/// \param nparts the number of parts in the final partition
///
/// \tparam TriangleMesh is a model of the `FaceListGraph` concept.
/// \tparam NamedParameters a sequence of \ref namedparameters
/// \tparam FacePartitionIDPmap is is a model of `ReadWritePropertyMap`
///           with `boost::graph_traits<TriangleMesh>::%face_descriptor`
///           as key type and `boost::face_external_index` as value type.
///           The default is `typename boost::property_map<FaceGraph, face_external_index>::%type`.
///
/// \param tm a triangle mesh
/// \param np optional \ref namedparameters described below
///
/// \pre `m` is a pure triangular surface mesh: there are no edges
///       without at least one incident face
template<typename TriangleMesh, typename FacePartitionIDPmap>
void partition(const TriangleMesh& tm,
               const int nparts, FacePartitionIDPmap partition_id_map)
{
  return partition(tm, nparts, partition_id_map, CGAL::parameters::all_default());
}

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_PARTITION_H
