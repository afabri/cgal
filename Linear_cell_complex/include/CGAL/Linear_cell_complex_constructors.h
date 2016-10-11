// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H
#define CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H 1

#include <CGAL/Combinatorial_map_constructors.h>

#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/Linear_cell_complex_incremental_builder.h>
#include <iostream>
#include <map>
#include <vector>
#include <list>

namespace CGAL {

  template < class LCC >
  void load_off(LCC& alcc, std::istream& in)
  {
    File_header_OFF  m_file_header;
    File_scanner_OFF scanner( in, m_file_header.verbose());
    if ( ! in) return;
    m_file_header = scanner;  // Remember file header after return.

    Linear_cell_complex_incremental_builder_3<LCC> B( alcc);
    B.begin_surface( scanner.size_of_vertices(),
                     scanner.size_of_facets(),
                     scanner.size_of_halfedges());

    typedef typename LCC::Point Point;

    // read in all vertices
    std::size_t  i;
    for ( i = 0; i < scanner.size_of_vertices(); i++) {
      Point p;
      file_scan_vertex( scanner, p);
      B.add_vertex( p);
      scanner.skip_to_next_vertex( i);
    }
    /* TODO rollback
       if ( ! in  || B.error()) {
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
    */

    // read in all facets
    for ( i = 0; i < scanner.size_of_facets(); i++)
    {
      B.begin_facet();
      std::size_t no;
      scanner.scan_facet( no, i);
      /* TODO manage errors
         if( ! in || B.error() || no < 3) {
         if ( scanner.verbose()) {
         std::cerr << " " << std::endl;
         std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
         std::cerr << "operator()(): input error: facet " << i
         << " has less than 3 vertices." << std::endl;
         }
         B.rollback();
         in.clear( std::ios::badbit);
         return;
         } */
      for ( std::size_t j = 0; j < no; j++) {
        std::size_t index;
        scanner.scan_facet_vertex_index( index, i);
        B.add_vertex_to_facet( index);
      }
      B.end_facet();
      scanner.skip_to_next_facet( i);
    }
    /* TODO manage errors
       if ( ! in  || B.error()) {
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
       if ( B.check_unconnected_vertices()) {
       if ( ! B.remove_unconnected_vertices()) {
       if ( scanner.verbose()) {
       std::cerr << " " << std::endl;
       std::cerr << "Polyhedron_scan_OFF<Traits>::" << std::endl;
       std::cerr << "operator()(): input error: cannot "
       "succesfully remove isolated vertices."
       << std::endl;
       }
       B.rollback();
       in.clear( std::ios::badbit);
       return;
       }
       }*/
    B.end_surface();
  }



  /** Export the alcc in off file format. If dimension>2, export all faces but only once.
   */
  template < class LCC >
  void write_off(LCC& alcc, std::ostream& out)
  {
    File_header_OFF header(false);
    header.set_binary(is_binary( out));
    header.set_no_comments(!is_pretty( out));
    File_writer_OFF writer( header);
    writer.header().set_polyhedral_surface(true);
    writer.header().set_halfedges( alcc.number_of_darts());

    // Print header.
    writer.write_header( out,
                         alcc.number_of_vertex_attributes(),
                         alcc.number_of_darts(),
                         alcc.template one_dart_per_cell<2>().size() );

    typedef typename LCC::Vertex_attribute_range::iterator VCI;
    VCI vit, vend = alcc.vertex_attributes().end();
    for ( vit = alcc.vertex_attributes().begin(); vit!=vend; ++vit )
    {
      writer.write_vertex( ::CGAL::to_double( vit->point().x()),
                           ::CGAL::to_double( vit->point().y()),
                           ::CGAL::to_double( vit->point().z()));
    }

    typedef Inverse_index< VCI > Index;
    Index index( alcc.vertex_attributes().begin(),
                 alcc.vertex_attributes().end());
    writer.write_facet_header();

    typename LCC::size_type m = alcc.get_new_mark();

    for ( typename LCC::Dart_range::iterator itall = alcc.darts().begin(),
            itallend = alcc.darts().end(); itall!=itallend; ++itall )
    {
      if ( !alcc.is_marked(itall, m) )
      {
        std::size_t n = 0;
        // First we count the number of vertices of the face.
        for ( typename LCC::template Dart_of_orbit_range<1>::iterator
                itf=alcc.template darts_of_orbit<1>(itall).begin(),
                itfend=alcc.template darts_of_orbit<1>(itall).end();
              itf!=itfend; ++itf, ++n );

        CGAL_assertion( n>=3 );
        writer.write_facet_begin(n);

        // Second we write the indices of vertices.
        for ( typename LCC::template Dart_of_orbit_range<1>::iterator
                itf=alcc.template darts_of_orbit<1>(itall).begin(),
                itfend=alcc.template darts_of_orbit<1>(itall).end();
              itf!=itfend; ++itf )
        {
          // TODO case with index
          writer.write_facet_vertex_index(index[VCI(alcc.vertex_attribute(itf))]);

          for ( typename LCC::template Dart_of_involution_basic_range<1>::iterator
                  itinv=alcc.template darts_of_involution_basic<1>(itf, m).begin(),
                  itinvend=alcc.template darts_of_involution_basic<1>(itf, m).end();
                itinv!=itinvend; ++itinv )
            alcc.mark(itinv, m);
        }
        writer.write_facet_end();
      }
    }
    writer.write_footer();
    alcc.free_mark(m);
  }

} // namespace CGAL

#endif // CGAL_LINEAR_CELL_COMPLEX_CONSTRUCTORS_H //
// EOF //
