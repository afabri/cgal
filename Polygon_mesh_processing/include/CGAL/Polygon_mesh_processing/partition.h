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

// \ingroup PMPMesh_partitioning
//
// Output each part of a partition as a single mesh.
//
// \param tm a triangle mesh
// \param nparts the number of parts
// \param fpmap the property map with the partition indices
// \param filename_base Partitions will be output in `.off` files named
//                      `{filename_base}_[0...nparts].off`
//
// \tparam TriangleMesh must be a model of a `FaceListGraph`, `HalfedgeListGraph`, and \bgllink{VertexListGraph}.
// \tparam FacePartitionIDPmap is a model of `ReadablePropertyMap`
//           with `boost::graph_traits<TriangleMesh>::%face_descriptor`
//           as key type and `boost::graph_traits<Graph>::%faces_size_type` as value type.
template<typename TriangleMesh, typename FacePartitionIDPmap>
void output_partition(const TriangleMesh& tm,
                      const idx_t nparts,
                      const FacePartitionIDPmap fpmap,
                      const std::string filename_base)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));

  typedef CGAL::Face_filtered_graph<TriangleMesh>         Filtered_graph;

  for(int i=0; i<nparts; ++i)
  {
    std::ostringstream filename;
    filename << filename_base << "_" << i << ".off" << std::ends;
    std::ofstream out(filename.str().c_str());

    Filtered_graph m_part(tm, i, fpmap);

    TriangleMesh out_mesh;
    CGAL::copy_face_graph(m_part, out_mesh);

    out << out_mesh;
  }
}

template<typename TriangleMesh, typename FacePartitionIDPmap,
         typename MEDIT_options, typename NamedParameters>
void partition(const TriangleMesh& tm,
               int nparts, FacePartitionIDPmap partition_id_map,
               MEDIT_options options, // pointer to the options array
               const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tm));
  CGAL_precondition_msg(nparts > 1, ("Partitioning requires a number of parts > 1"));

  using boost::choose_param;
  using boost::get_param;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator       face_iterator;

  idx_t nn = static_cast<idx_t>(num_vertices(tm));
  idx_t ne = static_cast<idx_t>(num_faces(tm));
  idx_t d = 3; // number of nodes per element
  idx_t* eptr = new idx_t[ne + 1];
  idx_t* eind = new idx_t[d * ne];

  //Vertex index map
  typedef typename GetVertexIndexMap<TriangleMesh, NamedParameters>::type Indices;
  Indices indices = choose_param(get_param(np, internal_np::vertex_index),
                                 get_const_property_map(boost::vertex_index, tm));

  // fill the adjacency info
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
                               *options,
                               &objval, epart, npart);

  std::cout << "return: " << ret << " with objval: " << objval << std::endl;
  CGAL_assertion(ret == METIS_OK);

  boost::tie(fit, fe) = faces(tm);
  for(int i=0; fit!=fe; ++fit, ++i)
    put(partition_id_map, *fit, epart[i]);
}

template<typename TriangleMesh, typename FacePartitionIDPmap, typename NamedParameters>
void partition(const TriangleMesh& tm,
               int nparts, FacePartitionIDPmap partition_id_map,
               const boost::param_not_found, // no MEDIT options were passed
               const NamedParameters& np)
{
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;
  return partition(tm, nparts, partition_id_map, &options, np);
}

/// \ingroup PMPMesh_partitioning
///
/// Computes a partition of the input mesh into `nparts` parts.
///
/// Property map for `CGAL::vertex_index_t` should be either available
/// as an internal property map to `tm` or provided as \ref pmp_namedparameters "Named Parameters".
///
/// \param tm a triangle mesh
/// \param nparts the number of parts in the final partition
/// \param partition_id_map a property map of type `FacePartitionIDPmap`
/// \param np optional \ref pmp_namedparameters "Named Parameters" described below
///
/// \tparam TriangleMesh is a model of the `FaceListGraph` concept.
/// \tparam FacePartitionIDPmap is a model of `ReadWritePropertyMap`
///           with `boost::graph_traits<TriangleMesh>::%face_descriptor`
///           as key type and `boost::graph_traits<Graph>::%faces_size_type` as value type.
/// \tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// \cgalNamedParamsBegin
///   \cgalParamBegin{vertex_index_map}
///     a property map containing the index of each vertex of `tm` intialized from `0` to `num_vertices(tm)-1`.
///   \cgalParamEnd
///   \cgalParamBegin{METIS_options}
///     is a parameter used in `partition()` to pass options to the METIS mesh
///     partitioner. The many options of METIS are not described here. Instead, users
///     should refer to METIS' <a href="http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf">documentation</a>.
///   \cgalParamEnd
/// \cgalNamedParamsEnd
///
/// \pre `m` is a pure triangular surface mesh: there are no edges
///       without at least one incident face
template<typename TriangleMesh, typename FacePartitionIDPmap, typename NamedParameters>
void partition(const TriangleMesh& tm,
               int nparts, FacePartitionIDPmap partition_id_map,
               const NamedParameters& np)
{
  using boost::get_param;

  return partition(tm, nparts, partition_id_map,
                   get_param(np, internal_np::METIS_options), np);
}

template<typename TriangleMesh, typename FacePartitionIDPmap>
void partition(const TriangleMesh& tm,
               const int nparts, FacePartitionIDPmap partition_id_map)
{
  return partition(tm, nparts, partition_id_map, CGAL::parameters::all_default());
}

} // end namespace Polygon_mesh_processing

} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_PARTITION_H
