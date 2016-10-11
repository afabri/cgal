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
#ifndef CGAL_IMPORT_FROM_PLANE_GRAPH_H
#define CGAL_IMPORT_FROM_PLANE_GRAPH_H 1

#include <CGAL/Combinatorial_map_constructors.h>
#include <iostream>
#include <map>
#include <vector>
#include <list>

namespace CGAL {

  /** @file Linear_cell_complex_constructors.h
   * Some construction operations for a linear cell complex from other
   * CGAL data structures.
   */

  /** Import an embedded plane graph read into a flux into a
   *  linear cell complex.
   * @param alcc the linear cell complex where the graph will be imported.
   * @param ais the istream where read the graph.
   * @return A dart created during the convertion.
   */
  template< class LCC >
  typename LCC::Dart_handle import_from_plane_graph(LCC& alcc,
                                                    std::istream& ais)
  {
    CGAL_static_assertion( LCC::dimension>=2 && LCC::ambient_dimension==2 );
  
    typedef typename LCC::Dart_handle Dart_handle;
    typedef typename LCC::Traits::Direction_2 Direction;
    typedef typename std::list<Dart_handle>::iterator List_iterator;
    typedef typename std::map<Direction, Dart_handle>::iterator LCC_iterator;
  
    // Arrays of vertices
    std::vector< typename LCC::Vertex_attribute_handle > initVertices;
    std::vector< std::list<Dart_handle> > testVertices;

    std::string txt;
    typename LCC::FT x, y;
    Dart_handle d1 = alcc.null_handle, d2 = alcc.null_handle;
    unsigned int v1, v2;
  
    unsigned int nbSommets = 0;
    unsigned int nbAretes = 0;
  
    ais >> nbSommets >> nbAretes;
    while (nbSommets > 0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough vertices."
                  << std::endl;
        return alcc.null_handle;
      }

      ais >> iformat(x) >> iformat(y);
      initVertices.push_back(alcc.create_vertex_attribute
                             (typename LCC::Point(x, y)));
      testVertices.push_back(std::list<Dart_handle>());
      --nbSommets;
    }

    while (nbAretes > 0)
    {
      if (!ais.good())
      {
        std::cout << "Problem: file does not contain enough edges."
                  << std::endl;
        return alcc.null_handle;
      }

      // We read an egde (given by the number of its two vertices).
      ais >> v1 >> v2;
      --nbAretes;

      CGAL_assertion(v1 < initVertices.size());
      CGAL_assertion(v2 < initVertices.size());

      d1 = alcc.create_dart(initVertices[v1]);
      d2 = alcc.create_dart(initVertices[v2]);
      alcc.template link_beta<2>(d1, d2);

      testVertices[v1].push_back(d1);
      testVertices[v2].push_back(d2);
    }

    // LCC associating directions and darts.
    std::map<Direction, Dart_handle> tabDart;
    List_iterator it;
    LCC_iterator  it2;

    Dart_handle first = alcc.null_handle;
    Dart_handle prec = alcc.null_handle;
    typename LCC::Point sommet1, sommet2;
  
    for (unsigned int i = 0; i < initVertices.size(); ++i)
    {
      it = testVertices[i].begin();
      if (it != testVertices[i].end()) // Si la liste n'est pas vide.
      {
        // 1. We insert all the darts and sort them depending on the direction
        tabDart.clear();
      
        sommet1 = alcc.point(*it);
        sommet2 = alcc.point(alcc.beta(*it,2));
      
        tabDart.insert(std::pair<Direction, Dart_handle>
                       (typename LCC::Traits::Construct_direction_2()
                        (typename LCC::Traits::Construct_vector()
                         (sommet1,sommet2)), *it));
      
        ++it;
        while (it != testVertices[i].end())
        {
          sommet2 = alcc.point(alcc.beta(*it,2));
          tabDart.insert(std::pair<Direction, Dart_handle>
                         (typename LCC::Traits::Construct_direction_2()
                          (typename LCC::Traits::Construct_vector()
                           (sommet1,sommet2)), *it));
          ++it;
        }
      
        // 2. We run through the array of darts and 1 links darts.
        it2 = tabDart.begin();
        first = it2->second;
        prec = first;
        ++it2;

        while (it2 != tabDart.end())
        {
          alcc.template link_beta<0>(prec, alcc.beta(it2->second,2));
          prec = it2->second;
          ++it2;
        }
        alcc.template link_beta<0>(prec, alcc.beta(first,2));
      }
    }

    // We return a dart from the imported object.
    return first;
  }

} // namespace CGAL

#endif // CGAL_IMPORT_FROM_PLANE_GRAPH_H //
// EOF //
