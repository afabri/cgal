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
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_BOOST_GRAPH_READ_DRC_H

#include <cinttypes>
#include <fstream>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/named_function_params.h>

#include <draco/compression/decode.h>

namespace CGAL {
/*!
   \ingroup PkgBGLIOFct
    reads the graph `g` in the drc format.

   \cgalNamedParamsBegin
    *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `g`.
    *       If this parameter is omitted, an internal property map for
    *       `CGAL::vertex_point_t` should be available in `FaceGraph`\cgalParamEnd
    * \cgalNamedParamsEnd
    \sa Overloads of this function for specific models of the concept `FaceGraph`.
    \pre The data must represent a 2-manifold
    \attention The graph `g` is not cleared, and the data from the stream are added.
*/
  template <typename FaceGraph, typename NamedParameters>
  bool read_drc(char* fname,
                FaceGraph& g,
                NamedParameters np)
  {
    typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<FaceGraph>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<FaceGraph>::faces_size_type faces_size_type;
    
    typedef typename Polygon_mesh_processing::GetVertexPointMap<FaceGraph, NamedParameters>::type Vpm;
    typedef  typename boost::property_traits<Vpm>::value_type Point_3;
    
    Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_property_map(CGAL::vertex_point, g));
    vertices_size_type nv;
    faces_size_type nf;
  
    std::ifstream input_file(fname, std::ios::binary);

    // Read the file stream into a buffer.
    std::streampos file_size = 0;
    input_file.seekg(0, std::ios::end);
    file_size = input_file.tellg() - file_size;
    input_file.seekg(0, std::ios::beg);
    std::vector<char> data(file_size);
    input_file.read(data.data(), file_size);

    if (data.empty()) {
      std::cerr << "Empty input file.\n";
      return false;
    }

    // Create a draco decoding buffer. Note that no data is copied in this step.
    draco::DecoderBuffer buffer;
    buffer.Init(data.data(), data.size());

    // Decode the input data into a geometry.
    std::unique_ptr<draco::PointCloud> pc;
    draco::Mesh *mesh = nullptr;
    auto type_statusor = draco::Decoder::GetEncodedGeometryType(&buffer);
    if (!type_statusor.ok()) {
      std::cerr << "Decoding failed.\n";
      return false;
    }
    const draco::EncodedGeometryType geom_type = type_statusor.value();
    if (geom_type == draco::TRIANGULAR_MESH) {
      draco::Decoder decoder;
      auto statusor = decoder.DecodeMeshFromBuffer(&buffer);
      if (!statusor.ok()) {
        std::cerr << "Decoding failed.\n";
        return false;
      }
      std::unique_ptr<draco::Mesh> in_mesh = std::move(statusor).value();

      if (in_mesh) {
        mesh = in_mesh.get();
        pc = std::move(in_mesh);
      }
    }else{
      std::cerr << "Not a triangle mesh.\n";
      return false;
    }

    const draco::PointAttribute *const att =  pc->GetNamedAttribute(draco::GeometryAttribute::POSITION);
    if (att == nullptr || att->size() == 0){
      std::cerr << "Missing position attribute.\n";
      return false;
    }
    nv = static_cast<vertices_size_type>(att->size());
    nf = static_cast<faces_size_type>(mesh->num_faces());
    
    std::vector<vertex_descriptor> vertices(nv);
  
    std::array<float, 3> value;
    int ii=0;
    for (draco::AttributeValueIndex i(0); i < static_cast<uint32_t>(att->size()); ++i) {
      att->ConvertValue<float, 3>(i, &value[0]);
      
      vertices[ii] = add_vertex(g);
      put(vpm,vertices[ii],Point_3(value[0], value[1], value[2]));
      ++ii;
    }
    std::vector<vertex_descriptor> face(3);
    for (draco::FaceIndex i(0); i < mesh->num_faces(); ++i) {
      
      for (int j = 0; j < 3; ++j) {
        const draco::PointIndex vert_index = mesh->face(i)[j];
        face[j]= vertices[att->mapped_index(vert_index).value()];  
      }
      Euler::add_face(face,g);
    }
    return true;
  }

   template <typename FaceGraph>
  bool read_drc(char* fname,
                FaceGraph& g)
  {
    return read_drc(fname,g, parameters::all_default());
  }
  
} // namespace CGAL

#endif // CGAL_BOOST_GRAPH_READ_DRC_H
