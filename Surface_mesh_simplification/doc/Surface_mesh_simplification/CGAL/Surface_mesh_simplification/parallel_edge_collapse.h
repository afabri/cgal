namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

Simplifies `tm` in-place by collapsing edges, and returns
the number of edges effectively removed.

The function `parallel_edge_collapse()` simplifies in-place a triangulated surface mesh.
A property map that associates a number between `0` and `partition_size-1` to each face
defines a partition of the faces in components. Each component is simplified 
with the sequential algorithm with a layer of edges incident to the boundary
of the components being constrained. The simplification of components is done
in parallel tasks. Once finished 2 layers of edges incident to the boundary
of the components are simplified with the sequential algorithm, while all other
edges are constrained. 

\tparam TriangleMesh must be a model of the concept `EdgeCollapsableSurfaceMesh`

\tparam FacePartionMap must be a model of `ReadablePropertyMap` with the key type
`boost::graph_traits<EdgeCollapsableSurfaceMesh>::%face_descriptor` 
and the value type `boost::graph_traits<EdgeCollapsableSurfaceMesh>::%faces_size_type`.


\param tm  is the mesh to simplify. 

\param should_stop is the stop-condition policy. 
It must be a model of the `StopPredicate` concept. 

\param fpm is the property map that associates to each face the component it is in.

\param partition_size must be the number of different values in `fpm`.

\param np optional sequence of Named Parameters among the ones listed below

\cgalHeading{Named Parameters}

`named_parameters` holds the list of all the additional parameters 
used by the `parallel_edge_collapse()` function (including default parameters). 

The parameters are the same as for `edge_collapse()`.


*/

template<class TriangleMesh, class StopPredicate, class FacePartionMap, class P, class T, class R>
int parallel_edge_collapse ( TriangleMesh& tm
                             , StopPredicate const& should_stop
                             , FacePartionMap fpm
                             , int partition_size
                             , sms_named_params<P,T,R> const& np
) ;

} /* namespace Surface_mesh_simplification */
} /* namespace CGAL */

