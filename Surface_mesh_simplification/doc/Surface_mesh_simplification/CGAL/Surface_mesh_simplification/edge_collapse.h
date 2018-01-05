namespace CGAL {
namespace Surface_mesh_simplification {

/*!
\ingroup PkgSurfaceMeshSimplification

Simplify `tm` in-place by collapsing edges, and returns
the number of edges effectively removed.

@tparam TriangleMesh a model of `EdgeCollapsableSurfaceMesh`
@tparam StopPolicy a model of `StopPredicate`
@tparam NamedParameters a sequence of \ref sms_namedparameters "Named Parameters"

@param tm a triangle mesh
@param should_stop the stop-condition policy
@param np optional sequence of \ref sms_namedparameters "Named Parameters" among the ones listed below

\cgalNamedParamsBegin
  \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of the mesh.
     If this parameter is omitted, an internal property map for
     `CGAL::vertex_point_t` should be available in `PolygonMesh`
  \cgalParamEnd

  \cgalParamBegin{halfedge_index_map} the property map containing an index for each halfedge,
    initialized 0 to `num_halfedges(graph)`
  \cgalParamEnd

  \cgalParamBegin{get_cost}
    The policy which returns the collapse cost for an edge.
  \cgalParamEnd

  \cgalParamBegin{get_placement}
    The policy which returns the placement (position of the replacemet vertex) for an edge.
  \cgalParamEnd

  \cgalParamBegin{edge_is_constrained_map}
    The property map containing the constrained-or-not status of each edge of `pmesh`
   \attention If this parameter is provided, `surface_mesh` must be a model of the
   `EdgeCollapsableSurfaceMeshWithConstraints` concept.
  \cgalParamEnd

  \cgalParamBegin{visitor}
    The visitor that is called by the `edge_collapse` function
    in certain points to allow the user to track the simplification process.
  \cgalParamEnd

  \cgalParamBegin{vertex_index_map}
    is the property map containing the index of each vertex of the input polygon mesh.
    \cgalDebugBegin
    This parameter is only used by debug functions and is usually not needed for users.
    \cgalDebugEnd
  \cgalParamEnd
\cgalNamedParamsEnd

\cgalHeading{Semantics}

The simplification process continues until the `should_stop` predicate returns `true` 
or the surface mesh cannot be simplified any further due to topological constraints. 

`get_cost` and `get_placement` are the policies which control 
the <I>cost-strategy</I>, that is, the order in which edges are collapsed 
and the remaining vertex is re-positioned. 

`visitor` is used to keep track of the simplification process. It has several member functions which 
are called at certain points in the simplification code. 
*/
template<class TriangleMesh, class StopPolicy, class NamedParameters>
int edge_collapse(TriangleMesh& tm,
                  const StopPolicy& should_stop,
                  const NamedParameters& np);

/*!
\ingroup PkgSurfaceMeshSimplification

Simplify the triangulated surface mesh `tm` in-place by iteratively collapsing edges in parallel
and returns the number of edges effectively removed.

The property map `fpm` associates a number between `0` and `partition_size-1` to each face
which defines a partition of the faces in components. Each component is simplified
with the sequential algorithm with a layer of edges incident to the boundary
of the components being constrained. The simplification of components is done
in parallel tasks. Once finished, two layers of edges incident to the boundary
of the components are simplified with the sequential algorithm, while all other
edges are constrained.

\tparam TriangleMesh must be a model of the concept `EdgeCollapsableSurfaceMesh`
\tparam StopPolicy must be a model of the concept `StopPredicate`
\tparam FacePartionMap must be a model of `ReadablePropertyMap` with the key type
`boost::graph_traits<EdgeCollapsableSurfaceMesh>::%face_descriptor`
and the value type `boost::graph_traits<EdgeCollapsableSurfaceMesh>::%faces_size_type`
\tparam NamedParameters a sequence of \ref sms_namedparameters "Named Parameters"

\param tm is the mesh to simplify
\param stop is the stop predicate
\param fpm is the property map that associates to each face the component it is in
\param partition_size must be the number of different values in `fpm`
\param np optional sequence of \ref sms_namedparameters "Named Parameters".
       These named parameters are the same as in `edge_collapse()`. Note that the
       `halfedge_index_map()` parameter is not needed.
*/
template<class TriangleMesh, class Stop, class FacePartionMap, class NamedParameters>
int parallel_edge_collapse(TriangleMesh& tm,
                           const StopPolicy& stop,
                           FacePartionMap fpm,
                           int partition_size,
                           const NamedParameters& np);

} /* namespace Surface_mesh_simplification */
} /* namespace CGAL */

