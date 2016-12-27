#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>
#include <CGAL/Random.h>
#include <CGAL/boost/graph/helpers.h>
#include <iostream>
#include <fstream>
#include <map>
#include <queue>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Ray_3 Ray_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

std::pair<Point_3,Point_3> canonical(const Point_3& p, const Point_3& q)
{
  if(p < q){
    return std::make_pair(p,q);
  }
  return std::make_pair(q,p);
}


// Find any face on the hull
// Algorithm: Shoot a ray from the centroid of any face in direction of its normal.
//            Return the farthest face. Note that the ray might pass through a
//            vertical face so we need potentially shoot again with a perturbed normal
// Don't use an aabb tree as we just want to do 1 shooting, but look at all faces
std::pair<face_descriptor,bool>
 face_on_hull(Surface_mesh mesh)
{
  CGAL::Random rng;
  boost::property_map<Surface_mesh,CGAL::vertex_point_t>::type vpm = get(CGAL::vertex_point,mesh);

  face_descriptor fd = *faces(mesh).first;
  halfedge_descriptor hd = halfedge(fd,mesh);
  Vector_3 n = PMP::compute_face_normal(fd,mesh);
  Point_3 centroid = CGAL::centroid(get(vpm, source(hd,mesh)),
                                    get(vpm, target(hd,mesh)),
                                    get(vpm, target(next(hd,mesh),mesh)));


  face_descriptor farthest_fd;
  bool normal_points_outwards = true;

  while(true){
    double sd = 0;
    farthest_fd = fd;
    Ray_3 ray(centroid,n);
    bool perturb = false;

    for(face_descriptor fd2 : faces(mesh)){
      if(fd == fd2){
        continue;
      }
      halfedge_descriptor hd2 = halfedge(fd2,mesh);
      Triangle_3 triangle(get(vpm, source(hd2,mesh)),
                          get(vpm, target(hd2,mesh)),
                          get(vpm, target(next(hd2,mesh),mesh)));

      CGAL::cpp11::result_of<K::Intersect_3(Triangle_3, Ray_3)>::type res = CGAL::intersection(ray,triangle);
      if(! res){
        continue;
      }
      if (const Point_3 *p = boost::get<Point_3>(&*res)){ 
        double sd2 = CGAL::squared_distance(*p,centroid);
        if(sd2 > sd){
          sd = sd2;
          farthest_fd = fd2;
          continue;
        }
      } else if (const Segment_3 *s = boost::get<Segment_3>(&*res)){
        perturb = true;
        break;
      }
    }
    if(perturb){
      n = Vector_3(n.x()+rng.get_double(0.99,1.01),
                   n.y()+rng.get_double(0.99,1.01),
                   n.z()+rng.get_double(0.99,1.01));
    }else{
      break;
    }
  }
  if(farthest_fd != fd){
    halfedge_descriptor hd = halfedge(farthest_fd,mesh);
    normal_points_outwards = CGAL::orientation(get(vpm, source(hd,mesh)),
                                               get(vpm, target(hd,mesh)),
                                               get(vpm, target(next(hd,mesh),mesh)),
                                               centroid) == CGAL::NEGATIVE;
  }

  return std::make_pair(farthest_fd, normal_points_outwards);
}

      
int main(int argc, char* argv[])
{ 
  const char* filename = (argc > 1) ? argv[1] : "data/tet-shuffled.off";
  std::ifstream input(filename);

  if (!input)
  {
    std::cerr << "Cannot open file " << std::endl;
    return 1;
  }

  std::vector<K::Point_3> points;
  std::vector< std::vector<std::size_t> > polygons;
  if (!CGAL::read_OFF(input, points, polygons))
  {
    std::cerr << "Error parsing the OFF file " << std::endl;
    return 1;
  }

  Surface_mesh mesh;

  PMP::orient_polygon_soup(points, polygons);

  PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
  std::size_t rem = CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);
  std::cerr << "Removed " << rem << " isolated vertices" << std::endl;

  // for all edges given as a pair of points we store the adjacent border halfedges
  typedef std::multimap<std::pair<Point_3,Point_3>,halfedge_descriptor > Point_pair_2_halfedges;
  Point_pair_2_halfedges point_pair_2_halfedges;

  boost::property_map<Surface_mesh,CGAL::vertex_point_t>::type vpm = get(CGAL::vertex_point,mesh);
  for(halfedge_descriptor hd : halfedges(mesh)){
    if(is_border(hd,mesh)){
      point_pair_2_halfedges.insert(std::make_pair(canonical(get(vpm,source(hd,mesh)),
                                                             get(vpm,target(hd,mesh))),
                                                   hd));
    }
  }
  
  // assert that there is one or more than two halfedges for each key
  // that is either a border, or a non-manifold situation


  std::vector<std::pair<halfedge_descriptor,halfedge_descriptor> > halfedges_to_stitch;

  // compute for each face in which connected component it is
  std::map<face_descriptor,int> face_cc_index_map;
  int nc = PMP::connected_components(mesh, boost::make_assoc_property_map(face_cc_index_map));

  std::cerr << nc << " connected components"<< std::endl;

  // container to mark which connected components we have already treated in the BFS
  std::vector<bool> cc_is_treated(nc,false);
 
  std::queue<std::pair<face_descriptor,bool> > Q;
  // find one face on the hull
  face_descriptor fd;
  bool normal_points_outwards;
  boost::tie(fd, normal_points_outwards) = face_on_hull(mesh);

  Q.push(std::make_pair(fd,normal_points_outwards));
  cc_is_treated[face_cc_index_map[fd]] = true;
  
  // Perform a BFS traversal of the connected components on the hull
  while(! Q.empty()){
    boost::tie(fd, normal_points_outwards) = Q.front();
    Q.pop();
    int fd_cc_index = face_cc_index_map[fd];

    std::vector<face_descriptor> cc;
    PMP::connected_component(fd, mesh, std::back_inserter(cc));
    if(! normal_points_outwards){
      PMP::reverse_face_orientations(cc, mesh);
    }
    std::vector<halfedge_descriptor> border;
    PMP::border_halfedges(cc, mesh,std::back_inserter(border));
  

    for(halfedge_descriptor bhd : border){
      Point_3 p = get(vpm, source(bhd,mesh));
      Point_3 q = get(vpm, target(bhd,mesh));
      Point_pair_2_halfedges::iterator b,e, b2;
      boost::tie(b,e) = point_pair_2_halfedges.equal_range(canonical(p,q));

      b2 = b;
      if(b != e){
        std::vector<halfedge_descriptor> he;
        for(;b!= e; ++b){
          if(bhd != b->second){
            he.push_back(b->second);
          }
        }
        Point_3 p1 = get(vpm,target(next(opposite(bhd,mesh),mesh),mesh));

        halfedge_descriptor nhd = he[0]; // not yet the closest
        for(int i = 1; i < he.size(); ++i){
          Point_3 p2 = get(vpm,target(next(opposite(nhd,mesh),mesh),mesh));
          Point_3 q2 = get(vpm,target(next(opposite(he[i],mesh),mesh),mesh));
          if(CGAL::is_in_interior_of_object<K>(q, p, p1, p2,q2)){
            nhd = he[i];
          }
        }
                
        // What we have so far is outwards oriented
        // Test if the neighbor connected component must get reversed when it comes out of the queue
        normal_points_outwards = (get(vpm,source(bhd,mesh)) == get(vpm,target(nhd,mesh)));

        halfedges_to_stitch.push_back(std::make_pair(bhd, nhd));
        face_descriptor nfd = face(opposite(nhd,mesh),mesh);
        if(! cc_is_treated[face_cc_index_map[nfd]]){
          Q.push(std::make_pair(nfd,normal_points_outwards));
          cc_is_treated[face_cc_index_map[nfd]] = true;
        }
        point_pair_2_halfedges.erase(b2,e);
      }
    }
  }

  PMP::stitch_borders(mesh, halfedges_to_stitch);
  
  std::vector<int> cc_to_remove;
  for(int i = 0; i < cc_is_treated.size(); ++i){
    if(! cc_is_treated[i]){
      cc_to_remove.push_back(i);
    }
  }

  CGAL::Polygon_mesh_processing::remove_connected_components(mesh,
                                                             cc_to_remove,
                                                             boost::make_assoc_property_map(face_cc_index_map));

  std::ofstream out("out.off");
  out << mesh;

  std::cerr << "done"<< std::endl;
  return 0;
}
