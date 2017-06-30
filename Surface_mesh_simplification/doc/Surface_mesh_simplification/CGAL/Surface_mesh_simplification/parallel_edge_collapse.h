namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

Simplifies `surface_mesh` in-place by collapsing edges, and returns
the number of edges effectively removed.

The function `parallel_edge_collapse()` simplifies in-place a triangulated surface mesh.
A property map that associates a number between `0` and `partition_size-1` to each face
defines a partition of the faces in components. Each component is simplified 
with the sequential algorithm with a layer of edges incident to the boundary
of the components being constrained. The simplification of components is done
in parallel tasks. Once finished 4 layers of edges incident to the boundary
of the components are simplified with the sequential algorithm, while all other
edges are constrained. 

\cgalHeading{Non-Named Parameters}

`surface_mesh` is the surface mesh to simplify. 
It must be a `Surface_mesh` or a class derived from it. 

`should_stop` is the stop-condition policy. 
It must be a model of the `StopPredicate` concept. 

`fpm` must be a model of `ReadablePropertyMap` with the key type
`boost::graph_traits<EdgeCollapsableSurfaceMesh>::face_descriptor` 
and the value type `boost::graph_traits<EdgeCollapsableSurfaceMesh>::faces_size_type`.

`partition_size` must be the number of different values in `fpm`.

\cgalHeading{Named Parameters}

`named_parameters` holds the list of all the additional parameters 
used by the `parallel_edge_collapse()` function (including default parameters). 

The parameters are the same as for `edge_collapse()`.


*/

template<class EdgeCollapsableSurfaceMesh,class StopPredicate, class FacePartionMap, class P, class T, class R>
int parallel_edge_collapse ( EdgeCollapsableSurfaceMesh& surface_mesh
                             , StopPredicate const& should_stop
                             , FacePartionMap fpm
                             , int partition_size
                             , sms_named_params<P,T,R> const& named_parameters
) ;

} /* namespace Surface_mesh_simplification */
} /* namespace CGAL */

