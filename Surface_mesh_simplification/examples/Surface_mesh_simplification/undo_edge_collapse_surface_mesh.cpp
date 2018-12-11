#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

//
// Setup an enriched polyhedron type which stores an id() field in the items
//

//#define CGAL_SM

#ifdef CGAL_SM
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
#else
typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Surface_mesh;
#endif

typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor ;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor ;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef SMS::Edge_profile<Surface_mesh> Profile ;



// The following is a Visitor that keeps track of the simplification process.


struct Record {
  vertex_descriptor v0;         // the vertex that gets removed
  vertex_descriptor v1;         // the vertex that survives the collapse_edge()
  vertex_descriptor oppa, oppb; // the vertices opposite to the edge to collapse
  Point_3 p0, p1;               // the original coordinates of the vertices 
};

  
struct Stats
{
  Stats(const Surface_mesh& sm) 
    : sm(sm) 
  {
    BOOST_FOREACH(vertex_descriptor vd, vertices(sm)){
      v2v[vd]=vd;
    }
  } 

  const Surface_mesh& sm;
  std::vector<Record> records;
  std::map<vertex_descriptor,vertex_descriptor> v2v;
} ;

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh>
{
  My_visitor( Stats* s) : stats(s){} 

  // Called during the collecting phase for each edge collected.
  void OnCollected( Profile const&, boost::optional<double> const& )
  {}                
  
  // Called during the processing phase for each edge selected.
  // If cost is absent the edge won't be collapsed.
  void OnSelected(Profile const&          
                 ,boost::optional<double>
                 ,std::size_t
                 ,std::size_t
                 )
  {}                
  
  // Called during the processing phase for each edge being collapsed.
  // If placement is absent the edge is left uncollapsed.
  void OnCollapsing(Profile const& profile         
                   ,boost::optional<Point>  placement
                   )
  {
    if ( placement ){
      //std::cout << "collapse edge " << profile.v0_v1() << " " << << profile.v0()  << " " << profile.v1() << std::endl;
      //std::cout << profile.v0() << " will get removed "  << profile.p0() << std::endl;
      Record record;
      record.v0 = profile.v0();
      record.v1 = profile.v1();
      record.p0 = profile.p0();
      record.p1 = profile.p1();
      halfedge_descriptor hd = profile.v0_v1();
      
      record.oppa = (! CGAL::is_border(hd,stats->sm))
        ? target(next(hd, stats->sm),stats->sm)
        : boost::graph_traits<Surface_mesh>::null_vertex();

      record.oppb = (! CGAL::is_border(opposite(hd, stats->sm),stats->sm))
        ? target(next(opposite(hd, stats->sm), stats->sm),stats->sm)
        : boost::graph_traits<Surface_mesh>::null_vertex();
      
      stats->records.push_back(record);
    }
  }                
  
  // Called for each edge which failed the so called link-condition,
  // that is, which cannot be collapsed because doing so would
  // turn the surface mesh into a non-manifold.
  void OnNonCollapsable( Profile const& )
  {}                
  
  // Called AFTER each edge has been collapsed
  void OnCollapsed( Profile const& prof, vertex_descriptor vd )
  {}                
  
  Stats* stats ;
} ;


int main( int argc, char** argv ) 
{

  Surface_mesh surface_mesh; 

#ifdef CGAL_SM  
  std::ifstream is(argv[1]) ; is >> surface_mesh;
#else
  OpenMesh::IO::read_mesh(surface_mesh, argv[1]);
#endif
  
  if (!CGAL::is_triangle_mesh(surface_mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // In this example, the simplification stops when the number of undirected edges
  // drops below 50% of the initial count
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.5);

 
  Stats stats(surface_mesh) ;
  
  My_visitor vis(&stats) ;

  // The index maps are not explicitelty passed as in the previous
  // example because the surface mesh items have a proper id() field.
  // On the other hand, we pass here explicit cost and placement
  // function which differ from the default policies, ommited in
  // the previous example.
  int rem = SMS::edge_collapse
    (surface_mesh
     ,stop
     ,CGAL::parameters::get_cost(SMS::Edge_length_cost  <Surface_mesh>())
     .get_placement(SMS::Midpoint_placement<Surface_mesh>())
     .visitor(vis)
     );

  int i = 0;

#ifdef DUMP_MESH  
  {
    Surface_mesh sm;
    copy_face_graph(surface_mesh,sm);
    std::string fn("collapse-");
    fn += boost::lexical_cast<std::string>(i++) + ".off";
    
    std::ofstream out(fn.c_str());
    out << sm;
  }
#endif
  
  while(! stats.records.empty()){
    //    std::cout << "undo collapse" << std::endl;
    Record r = stats.records.back();

    stats.records.pop_back();
    halfedge_descriptor h0, h1;
    if(r.oppa != boost::graph_traits<Surface_mesh>::null_vertex()){
      std::pair<halfedge_descriptor,bool> pa = halfedge(stats.v2v[r.oppa], stats.v2v[r.v1], surface_mesh);
      assert(pa.second);
      h0 = pa.first;
    } else {
      BOOST_FOREACH(halfedge_descriptor hd, CGAL::halfedges_around_target(halfedge(stats.v2v[r.v1],surface_mesh), surface_mesh)){
        if(CGAL::is_border(hd,surface_mesh)){
          h0 = hd;
          break;
        }
      }
    }
    if(r.oppb != boost::graph_traits<Surface_mesh>::null_vertex()){
      std::pair<halfedge_descriptor,bool> pa = halfedge(stats.v2v[r.oppb], stats.v2v[r.v1], surface_mesh);
      assert(pa.second);
      h1 = pa.first;
    } else {
      BOOST_FOREACH(halfedge_descriptor hd, CGAL::halfedges_around_target(halfedge(stats.v2v[r.v1],surface_mesh), surface_mesh)){
        if(CGAL::is_border(hd,surface_mesh)){
          h1 = hd;
          break;
        }
      }
    }

    halfedge_descriptor hnew = CGAL::Euler::split_vertex(h1, h0, surface_mesh);
    assert(target(hnew,surface_mesh) == stats.v2v[r.v1]);
    stats.v2v[r.v0] = source(hnew,surface_mesh);
    if(r.oppa != boost::graph_traits<Surface_mesh>::null_vertex()){
      CGAL::Euler::split_face(prev(h0,surface_mesh),next(h0, surface_mesh),surface_mesh);
    }
    if(r.oppb != boost::graph_traits<Surface_mesh>::null_vertex()){
      CGAL::Euler::split_face(prev(h1,surface_mesh),next(h1, surface_mesh),surface_mesh);
    }

    
    boost::property_map<Surface_mesh,boost::vertex_point_t>::type vpm = get(boost::vertex_point, surface_mesh);
    
    
    put(vpm, stats.v2v[r.v1], r.p1);
    put(vpm, source(hnew,surface_mesh), r.p0);

#ifdef DUMP_MESH    
    {
      Surface_mesh sm;
      copy_face_graph(surface_mesh,sm);
      std::string fn("collapse-");
      fn += boost::lexical_cast<std::string>(i++) + ".off";
      
      std::ofstream out(fn.c_str());
      out << sm;
    }
#endif
  }

#ifdef CGAL_SM  
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface_mesh ;
#else  
  surface_mesh.garbage_collection();
  OpenMesh::IO::write_mesh(surface_mesh, "out.off");
#endif

  return EXIT_SUCCESS ;      
}
