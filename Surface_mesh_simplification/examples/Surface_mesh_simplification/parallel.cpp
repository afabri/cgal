
#define INCREASE

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/partition.h>

#include <CGAL/Surface_mesh_simplification/parallel_edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

#include <CGAL/Real_timer.h>

#include <boost/lexical_cast.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;

typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh; 
typedef Surface_mesh::Property_map<boost::graph_traits<Surface_mesh>::halfedge_descriptor,int> HIMap;


typedef Surface_mesh::Property_map<boost::graph_traits<Surface_mesh>::edge_descriptor,char> ECMap;
typedef Surface_mesh::Property_map<boost::graph_traits<Surface_mesh>::face_descriptor,std::size_t> CCMap;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;



template <typename TriangleMesh, typename Placement, typename CCMap, typename Stop, typename Cost>
int parallel_edge_collapse(TriangleMesh& sm, CCMap ccmap, Placement placement, Stop stop, Cost cost, std::size_t ncc, bool dump = false)
{
  typedef boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
 
  CGAL::Real_timer t;


  TriangleMesh::Property_map<edge_descriptor,char> ecmap 
    = sm.add_property_map<edge_descriptor,char>("e:constrained",false).first;

  
  HIMap himap = sm.add_property_map<halfedge_descriptor,int>("h:index_in_cc",-1).first;
  
  std::vector<edge_descriptor> partition_edges;
 
  std::ofstream out;
  if(dump){
    out.open("partition.selection.txt");
    out << std::endl << std::endl;
  }
  // now we have to find the constrained edges
  BOOST_FOREACH(edge_descriptor ed, edges(sm)){
    if(is_border(ed,sm)){
      continue;
    }
    halfedge_descriptor hd = halfedge(ed,sm);
    halfedge_descriptor hop = opposite(hd,sm);
    if (ccmap[face(hd,sm)] !=ccmap[face(hop,sm)]) {
      partition_edges.push_back(ed);
      ecmap[ed] = true;
      himap[hd] = 0;
      himap[hop] = 1;
      if(dump){
        out << int(source(ed,sm)) << " " << int(target(ed,sm)) <<  " ";
      }
    }
  }
  if(dump){
    out << std::endl;
  }

  
  std::cerr << t.time() << " sec.\n";
  t.reset();

  std::vector<tbb::concurrent_vector<edge_descriptor> > cc_edges(ncc);
  tbb::parallel_for(tbb::blocked_range<size_t>( 0, num_faces(sm)),
                    Collect_edges<TriangleMesh>(cc_edges, ncc, ecmap, ccmap, sm));
 
  for(int i=0; i < ncc; i++){
    std::cout << "#edges in component "<< i << ": " << cc_edges[i].size() << std::endl;
  }
 
  std::vector<int> buffer_size(ncc); // the number of constrained edges inside each component
  std::cerr << "#CC = " << ncc << "  in " << t.time() << " sec." << std::endl;
  t.reset();

#if 1
  // Simplify the partition in parallel
  tbb::task_group tasks;
  for(int i = 0; i < ncc; i++){
    tasks.run(Simplify<TriangleMesh,Placement,Cost>(sm, himap, ecmap, ccmap, placement, cost,buffer_size, cc_edges[i], i, dump));
  }

  tasks.wait();
#else
  for(int i = 0; i < ncc; i++){
    tbb::task_group tasks;
    tasks.run(Simplify<TriangleMesh,Placement,Cost>(sm, himap, ecmap, ccmap, placement, cost, buffer_size, cc_edges[i], i,dump));
    tasks.wait();
  }
#endif

  std::cerr << "parallel edge collapse in " << t.time() << " sec." << std::endl;
  t.reset();

  if(dump){
    TriangleMesh sm2;
    std::cerr << "before copy_face_graph()"<< std::endl;
    CGAL::copy_face_graph(sm,sm2);
    std::cerr << "after copy_face_graph()"<< std::endl;
    std::ofstream outi("out-intermediary.off");
    outi << sm2 << std::endl;
  }

#if 1
  // After the parallel edge_collapse make the same for the buffer
#ifdef INCREASE
  // Increase the buffer by one more layer
  std::size_t num_constrained_edges = 0;
  BOOST_FOREACH(edge_descriptor ed, edges(sm)){
    ecmap[ed]=0;
  }

  BOOST_FOREACH(edge_descriptor ed, partition_edges){
    ecmap[ed]=1;
  }
  CGAL::expand_edge_selection(partition_edges,
                              sm,
                              3,
                              ecmap,
                              CGAL::Counting_output_iterator(&num_constrained_edges));
  num_constrained_edges += partition_edges.size();
#else
  std::size_t num_constrained_edges = partition_edges.size();
  for(int i =0; i < buffer_size.size(); i++){
    std::cout << "# partition edges = "<<   partition_edges.size() <<   " " << buffer_size[i] << std::endl;
    num_constrained_edges += buffer_size[i];
  }
#endif

  typedef Inverted_edge_constraint_map<ECMap> IECMap;
  IECMap iecmap(ecmap);
  
  num_constrained_edges = num_edges(sm) - num_constrained_edges; 

  
#ifdef INCREASE
  SMS::Constrained_placement <Placement, IECMap> constrained_placement (iecmap);

#endif
  double uc_ratio = 0.3;
  double cc_ratio = double(num_constrained_edges)/double(num_edges(sm));
  double ratio = uc_ratio + (1.0 - uc_ratio)*cc_ratio;
  assert(ratio < 1.0);
   
  std::cout << num_constrained_edges << " constrained; cc_ratio: " << cc_ratio << std::endl;
  std::cout << "ratio : " << ratio << std::endl;
  
  stop.second_pass(0,0);
  
#if 1
  SMS::edge_collapse(sm,
                     stop,
                     CGAL::Sequential_tag(),
                     CGAL::parameters::vertex_index_map(get(boost::vertex_index,sm))
                     .get_placement(constrained_placement)
                     .get_cost(cost)
                     .edge_is_constrained_map(iecmap)
);
#endif
  std::cerr << "sequential edge collapse on buffer in " << t.time() << " sec." << std::endl;

#endif
  return 0;
}


int main(int argc, char** argv ) 
{
  bool dump = false;

  typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
 
  CGAL::Real_timer t;
  t.start();
  std::ifstream is1(argv[1]);
  Surface_mesh sm; 
  is1 >> sm;
  std::cerr << "Input: #V = "<< num_vertices(sm)  << " #E = "<< num_edges(sm) << "  #F = " << num_faces(sm) << std::endl;

  typedef Surface_mesh::Property_map<face_descriptor,std::size_t> CCMap;
  Surface_mesh::Property_map<face_descriptor,std::size_t> ccmap 
    = sm.add_property_map<face_descriptor,std::size_t>("f:cc").first;

  std::size_t ncc = 8;
  PMP::partition(sm, ccmap, static_cast<int>(ncc));

  SMS::LindstromTurk_placement<Surface_mesh> placement;
  SMS::LindstromTurk_cost<Surface_mesh> cost;
  double ratio = 0.25;
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  //SMS::Midpoint_placement<Surface_mesh> placement;
  //SMS::Edge_length_cost<Surface_mesh> cost;
  //SMS::Edge_length_stop_predicate<double> stop(0.01);

  SMS::parallel_edge_collapse(sm, ccmap, placement, stop, cost, ncc, dump);


  std::cerr << "collect garbage" << std::endl;

  sm.collect_garbage();
  std::cerr << "\nResult: #V = "<< num_vertices(sm)  << " #E = "<< num_edges(sm) 
            << "  #F = " << num_faces(sm) << "  in " << t.time() << " sec." << std::endl;

  {
    std::ofstream out("out.off");
    out << sm << std::endl;
    out.close();
    std::cerr << "Writing result in " << t.time() << " sec." << std::endl;
  }
  return EXIT_SUCCESS;
}
