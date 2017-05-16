#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Primary_AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Mesh_3/experimental/AABB_filtered_projection_traits.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Timer.h>

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

  template <typename PatchId>
  struct Get_patch_id_from_primary_primitive {
    typedef void key_type;
    typedef PatchId value_type;
    typedef const value_type& reference;
    typedef boost::readable_property_map_tag                category;
    template <typename K>
    friend value_type get(const Get_patch_id_from_primary_primitive& pm, const K& k)
    {
      return k.patch_id;
    }
  };

template <typename Tree, typename PatchId>
struct Primary_primitive {
  typedef Primary_primitive<Tree,PatchId> Self;
  const Tree* secondary_tree;
  PatchId patch_id;

  typedef K::Point_3 Point;
  typedef Bbox_3 Datum;
  typedef Self Id;

  Primary_primitive()
    : secondary_tree(NULL), patch_id()
  {}


  Primary_primitive(const Tree* tree, PatchId pid)
    : secondary_tree(tree), patch_id(pid)
  {
    const_cast<Tree*>(secondary_tree)->build();
    secondary_tree->accelerate_distance_queries();
  }

  template <typename Iterator>
  Primary_primitive(const Iterator& it)
    : secondary_tree(it->secondary_tree), patch_id(it->patch_id)
  {}
  
  Id id() const
  { 
    return *this;
  }

  
  Datum datum() const
  {
    return secondary_tree->bbox();
  }


  Point reference_point() const
  {
    return secondary_tree->any_reference_point_and_id().first;
  }
};

} // namespace CGAL 


Point brute_force(const Point& query,
                  const std::vector<CGAL::Primary_primitive<Tree,std::size_t> >& primary_primitives,
                  const std::set<std::size_t> ids)
{
  Point p = primary_primitives[0].secondary_tree->closest_point(query);
  for(int i = 1; i < primary_primitives.size(); i++){
    if(ids.find(i) == ids.end()){
      Point p2 = primary_primitives[i].secondary_tree->closest_point(query);
      if(CGAL::compare_distance(query,p2,p) == CGAL::SMALLER){
        p = p2;
      }
    }
  }
  return p;
}

typedef CGAL::Primary_primitive<Tree,std::size_t> Primary_primitive;
typedef CGAL::Primary_AABB_traits<K, Primary_primitive> Primary_AABB_traits;
typedef CGAL::AABB_tree<Primary_AABB_traits> Primary_tree;

int main(int argc, char* argv[])
{

  std::cerr.precision(17);
  std::vector<Point> queries, result;
  {
    std::ifstream in(argv[1]);
    CGAL::read_xyz_points(in, std::back_inserter(queries));
  }
  result.resize(queries.size());

  std::set<std::size_t> ids;
  {
    std::ifstream in(argv[2]);
    std::size_t id;
    while(in >> id){
      ids.insert(id);
    }
  }

  std::vector<Primary_primitive> primary_primitives;
  primary_primitives.reserve(argc-3);
  std::list<Polyhedron> polyhedra;

  for(int i=3; i < argc; i++){
    std::ifstream in(argv[i]);
    Polyhedron polyhedron;
    in >> polyhedron;
    polyhedra.push_back(polyhedron);

    primary_primitives.push_back(Primary_primitive(new Tree(faces(polyhedra.back()).first, 
                                                            faces(polyhedra.back()).second,
                                                            polyhedra.back()),i-3));
  }
    Primary_tree ptree(primary_primitives.begin(), primary_primitives.end());

    CGAL::Timer t;
    t.start();
    for(int i = 0; i < queries.size(); i++){
      Point& query = queries[i];
      
      CGAL::Mesh_3::Filtered_projection_traits<
        Primary_AABB_traits,
        CGAL::Get_patch_id_from_primary_primitive<std::size_t>
        > projection_traits(ids.begin(), ids.end(),
                            ptree.traits());
      ptree.traversal(query, projection_traits);
      if(projection_traits.found()) {
        result[i] = projection_traits.closest_point(); 
      }else{
        std::cerr << "ERROR: no closest point" <<std::endl;
        return EXIT_FAILURE;
      }
    }
    std::cerr << t.time() << " sec."<< std::endl;


#if 1
    t.reset();
    t.start();
    for(int i = 0; i < queries.size(); i++){
      Point& query = queries[i];
      Point bres = brute_force(query, primary_primitives, ids);
      if(result[i] != bres){
        std::cerr << "Error for query point " << query << std::endl 
                  << "  TL: " << result[i]<<  " " << CGAL::squared_distance(query,result[i]) << std::endl 
                  << "  BF: " << bres <<  " " << CGAL::squared_distance(query,bres) << std::endl;
      }
    }
    std::cerr << t.time() << " sec."<< std::endl;
#endif
    
    K::Ray_3 ray_query(K::Point_3(1.0, 0.0, 0.0), K::Point_3(0.0, 1.0, 0.0));
    ptree.first_intersection(ray_query);
    return EXIT_SUCCESS;
}
