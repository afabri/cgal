
#ifndef POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
#define POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H

#include <CGAL/Mesh_3/Detect_features_in_polyhedra.h>

namespace CGAL
{
  template<typename Polyhedron>
  void reset_sharp_edges(Polyhedron* pMesh)
  {
    typename boost::property_map<Polyhedron,halfedge_is_feature_t>::type if_pm;
    BOOST_FOREACH(typename boost::graph_traits<Polyhedron>::edge_descriptor ed, edges(*pMesh))
    {
      put(if_pm,halfedge(ed,*pMesh),false);
      put(if_pm,opposite(halfedge(ed,*pMesh),*pMesh),false);
    }
  }

  template<typename Polyhedron>
  void detect_sharp_edges(Polyhedron* pMesh, const double angle)
  {
    reset_sharp_edges(pMesh);

    // Detect edges in current polyhedron
    typedef typename boost::property_map<Polyhedron,face_patch_id_t>::type PatchID;
    CGAL::Mesh_3::Detect_features_in_polyhedra<Polyhedron,PatchID> 
      detect_features(get(face_patch_id_t(),*pMesh));
    detect_features.detect_sharp_edges(*pMesh, angle);
  }

}//end namespace CGAL

#endif //POLYHEDRON_DEMO_DETECT_SHARP_EDGES_H
