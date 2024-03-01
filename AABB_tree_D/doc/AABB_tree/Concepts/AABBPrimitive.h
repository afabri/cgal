
/*!
\ingroup PkgAABB_treeConcepts
\cgalConcept

The concept `AABBPrimitive` describes the requirements for the primitives stored in the AABB tree data structure. The concept encapsulates a type for the input datum (a geometric object) and an identifier (id) type through which those primitives are referred to. The concept `AABBPrimitive` also refines the concepts DefaultConstructible and Assignable. 

\sa `CGAL::AABB_tree<AABBTraits>` 
\sa `AABBPrimitiveWithSharedData`

### Example ###

The `Primitive` type can be, e.g., a wrapper around a `Handle`. Assume for instance that the input objects are the triangle faces of a mesh stored as a `CGAL::Polyhedron_3`. The `Datum` would be a `Triangle_3` and the `Id` would be a polyhedron `Face_handle`. Method `datum()` can return either a `Triangle_3` constructed on the fly from the face handle or a `Triangle_3` stored internally. This provides a means for the user to trade memory for efficiency. 

\cgalHasModel `CGAL::AABB_primitive<Id,ObjectPropertyMap,PointPropertyMapPolyhedron,Tag_false,CacheDatumTag>`
\cgalHasModel `CGAL::AABB_segment_primitive<GeomTraits,Iterator,CacheDatum>`
\cgalHasModel `CGAL::AABB_triangle_primitive<GeomTraits,Iterator,CacheDatum>`
\cgalHasModel `CGAL::AABB_HalfedgeGraph_segment_primitive<HalfedgeGraph,Tag_false,CacheDatum>`
\cgalHasModel `CGAL::AABB_FaceGraph_triangle_primitive<FaceGraph,Tag_false,CacheDatum>`
*/

class AABBPrimitive {
public:

/// \name Types 
/// @{

/*! 
Point type.
*/ 
typedef unspecified_type Point; 

/*! 
Type of input datum. 
*/ 
typedef unspecified_type Datum; 

/*!
Point reference type returned by the function `point()`. It is convertible to the type `Point`.
 */
typedef unspecified_type Point_reference;

/*!
Datum reference type returned by the function `datum()`. It is convertible to the type `Datum`.
*/
typedef unspecified_type Datum_reference;

/*! 
Type of identifiers through which the input objects are referred to. It must be a model of the concepts DefaultConstructible and Assignable. 
*/ 
typedef unspecified_type Id; 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the datum (geometric object) represented by the primitive. 
*/ 
Datum_reference datum(); 

/*! 
Returns the corresponding identifier. This identifier is only used as a reference for the objects in the output of the `AABB_tree` methods. 
*/ 
Id id(); 

/*! 
Returns a Point located on the geometric object represented by the primitive. This function is used to sort the primitives during the AABB tree construction as well as to construct the search KD-tree internal to the AABB tree used to accelerate distance queries.
*/ 
Point_reference reference_point(); 

/// @}

}; /* end AABBPrimitive */

