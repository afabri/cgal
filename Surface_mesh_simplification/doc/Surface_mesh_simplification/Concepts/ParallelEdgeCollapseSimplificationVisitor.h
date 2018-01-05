
/*!
\ingroup PkgSurfaceMeshSimplificationConcepts
\cgalConcept

The concept `ParallelEdgeCollapseSimplificationVisitor` describes the requirements
for the <I>visitor object</I> which is used to track the parallel edge collapse simplification algorithm.

\cgalRefines `EdgeCollapseSimplificationVisitor`

\attention Note that operations such as `EdgeCollapseSimplificationVisitor::On_started()`,
           are called in each parallel thread.

*/
class ParallelEdgeCollapseSimplificationVisitor
{
public:
  /// \name Operations
  /// @{

  /*!
  Called after the parallel pass of the algorithm has finished and before the sequential pass starts.
  */
  void OnParallelPassFinished(TriangleMesh& surface_mesh,
                              Stop_predicate& pred,
                              size_type initial_num_edges,
                              size_type num_current_edges );
  /// @}

}; /* end EdgeCollapseSimplificationVisitor */
