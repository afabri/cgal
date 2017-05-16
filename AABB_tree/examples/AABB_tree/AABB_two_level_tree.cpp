#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Primary_AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

namespace CGAL {

template <typename Tree, typename PatchId>
struct Primary_primitive {
  const Tree* secondary_tree;
  PatchId patch_id;

  typedef K::Point_3 Point;
  typedef Bbox_3 Datum;
  typedef const Tree* Id;

  Primary_primitive()
    : secondary_tree(NULL), patch_id()
  {}


  Primary_primitive(const Tree* tree, int)
    : secondary_tree(tree)
  {
    const_cast<Tree*>(secondary_tree)->build();
    secondary_tree->accelerate_distance_queries();
  }

  template <typename Iterator>
  Primary_primitive(const Iterator& it)
    : secondary_tree(it->secondary_tree)
  {}
  
  Id id() const
  { 
    return secondary_tree;
  }

  
  Datum datum() const
  {
    return secondary_tree->bbox();
  }


  Point reference_point() const
  {
    return Point(datum().xmin(),datum().ymin(),datum().zmin());
  }
};

} // namespace CGAL 


Point brute_force(const Point& query, const std::vector<CGAL::Primary_primitive<Tree,std::size_t> >& primary_primitives)
{
  Point p = primary_primitives[0].secondary_tree->closest_point(query);
  for(int i = 1; i < primary_primitives.size(); i++){
    Point p2 = primary_primitives[i].secondary_tree->closest_point(query);
    if(CGAL::squared_distance(query,p2) < CGAL::squared_distance(query,p)){
      p = p2;
    }
  }
  return p;
}

typedef CGAL::Primary_primitive<Tree,std::size_t> Primary_primitive;
typedef CGAL::Primary_AABB_traits<K, Primary_primitive> Primary_AABB_traits;
typedef CGAL::AABB_tree<Primary_AABB_traits> Primary_tree;

int main(int argc, char* argv[])
{

  std::vector<Primary_primitive> primary_primitives;
  std::list<Polyhedron> polyhedra;

  for(int i =1; i < argc; i++){
    std::ifstream in(argv[i]);
    Polyhedron polyhedron;
    in >> polyhedron;
    polyhedra.push_back(polyhedron);
  }
  BOOST_FOREACH(Polyhedron& polyhedron, polyhedra){
    primary_primitives.push_back(Primary_primitive(new Tree(faces(polyhedron).first, 
                                                            faces(polyhedron).second,
                                                            polyhedron),0));
  }
    Primary_tree ptree(primary_primitives.begin(), primary_primitives.end());
    // query point
    Point query(0.0, 0.0, 0.0);

    Point res = ptree.closest_point(query);
#if COMPARE_RESULT
    Point bres = brute_force(query, primary_primitives);
    if(res != bres){
      std::cerr << "Error: " << res << " " << bres << std::endl;
    }
#endif
    K::Ray_3 ray_query(K::Point_3(1.0, 0.0, 0.0), K::Point_3(0.0, 1.0, 0.0));
    std::cout << ptree.number_of_intersected_primitives(ray_query) << std::endl;

    return EXIT_SUCCESS;
}
