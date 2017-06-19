
#define INCREASE

bool dump = false;

#include "tbb/task_group.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/partition.h>
#include <CGAL/boost/graph/Component_graph.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>

#include <CGAL/Real_timer.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
typedef CGAL::Simple_cartesian<double> Kernel;

typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh; 
typedef Surface_mesh::Property_map<boost::graph_traits<Surface_mesh>::halfedge_descriptor,int> HIMap;


typedef Surface_mesh::Property_map<boost::graph_traits<Surface_mesh>::edge_descriptor,char> ECMap;
typedef Surface_mesh::Property_map<boost::graph_traits<Surface_mesh>::face_descriptor,std::size_t> CCMap;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;


// Adaptor which interprets 0 and 1 as true, and any other index as false
template <typename ZOMap, typename G>
struct Selection_is_constraint_map
{  
  typedef boost::readable_property_map_tag                 category;
  typedef bool                                             value_type;
  typedef bool                                             reference;
  typedef typename boost::graph_traits<G>::edge_descriptor key_type;

  ZOMap zomap;
  const G& g;

  Selection_is_constraint_map(ZOMap zomap, const G& g)
    : zomap(zomap), g(g)
  {}



  friend bool get(const Selection_is_constraint_map& eicm, key_type ed)
  {
    int i = get(eicm.zomap, halfedge(ed, eicm.g));
    return (i ==0) || (i == 1);
  }
  
};


// Adapter which interprets char as bool to avoid "performance warning" with VC++
template <typename Char_map>
struct Bool_map
{
  Char_map cm;

  Bool_map(Char_map cm)
    : cm(cm)
  {}

  template <typename K>
  friend bool get(const Bool_map& bm, const K& k)
  {
    return get(bm.cm,k) != 0;
  }
};


template <typename ECMap>
struct Inverted_edge_constraint_map
{  
  typedef boost::readable_property_map_tag                 category;
  typedef bool                                             value_type;
  typedef bool                                             reference;
  typedef typename boost::property_traits<ECMap>::key_type key_type;

  ECMap ecmap;

  Inverted_edge_constraint_map(ECMap ecmap)
    : ecmap(ecmap)
  {}


  friend bool get(const Inverted_edge_constraint_map& iecm, key_type ed)
  {
    return ! (get(iecm.ecmap, ed)!= 0);
  }
  
};


template <typename ECMap, typename CCMap, typename G>
struct In_same_component_map
{  
  typedef boost::read_write_property_map_tag                 category;
  typedef typename boost::property_traits<ECMap>::value_type value_type;
  typedef value_type                                         reference;
  typedef typename boost::property_traits<ECMap>::key_type   key_type;

  ECMap ecmap;
  CCMap ccmap;
  const G& g;
  int cc;
  mutable int& count;

  In_same_component_map(ECMap ecmap, CCMap ccmap, const G& g, int cc, int& count)
    : ecmap(ecmap), ccmap(ccmap), g(g), cc(cc), count(count)
  {}


  friend value_type get(const In_same_component_map& iscm, key_type ed)
  {
    return get(iscm.ecmap,ed);
  }
  
  
  friend void put(const In_same_component_map& iscm, key_type ed, value_type v)
  {
    if( (! is_border(halfedge(ed,iscm.g),iscm.g) && iscm.ccmap[face(halfedge(ed,iscm.g),iscm.g)]== iscm.cc)||
        (! is_border(opposite(halfedge(ed,iscm.g),iscm.g),iscm.g) && iscm.ccmap[face(opposite(halfedge(ed,iscm.g),iscm.g),iscm.g)]== iscm.cc) ){
      if(! get(iscm.ecmap, ed)){
        put(iscm.ecmap, ed, v);
        ++(iscm.count);
      }
    }
  }
  
};


// The parallel task
struct Simplify {

  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
  Surface_mesh &sm;
  HIMap himap;
  ECMap ecmap;
  CCMap ccmap;
  std::vector<int>& buffer_size;
  halfedge_descriptor hd;
  int ccindex;

  Simplify(Surface_mesh& sm,
           HIMap himap,
           ECMap ecmap,
           CCMap ccmap,
           std::vector<int>& buffer_size,
           halfedge_descriptor hd, int ccindex)
    : sm(sm), himap(himap), ecmap(ecmap), ccmap(ccmap), buffer_size(buffer_size), hd(hd), ccindex(ccindex)
  {}


  void operator()() const
  {
    typedef Selection_is_constraint_map<HIMap, Surface_mesh> SICM;
    typedef CGAL::Component_graph<Surface_mesh,SICM>  Component_graph;
    SICM sicm(himap,sm);
    Component_graph cg(sm, sicm, hd); 

    std::vector<edge_descriptor> cc_edges;
    int i = 2, count = 0;
    BOOST_FOREACH(edge_descriptor ed, edges(cg)){
      ++count;
      halfedge_descriptor hd, hop;
      hd = halfedge(ed,sm);
      hop = opposite(hd,sm);
      if(get(sicm,ed)){
        cc_edges.push_back(ed);
      }else{
        put(himap,hd,i); ++i;
        put(himap,hop,i); ++i;
      }
    }
    cg.num_edges() = count;

    std::size_t cc_edges_count = cc_edges.size();

    std::vector<boost::graph_traits<Surface_mesh>::edge_descriptor> V;
    int number_of_puts = 0;
    In_same_component_map<ECMap,CCMap,Surface_mesh> iscmap(ecmap,ccmap,sm,ccindex, number_of_puts);

    CGAL::expand_edge_selection(cc_edges,
                                sm,
                                2,
                                iscmap,
                                std::back_inserter(V));
    cc_edges_count += number_of_puts;
    buffer_size[ccindex] = number_of_puts;

    if(dump){
      std::ofstream out(std::string("constraints-")+boost::lexical_cast<std::string>(ccindex)+".selection.txt");
      out << std::endl << std::endl;
      BOOST_FOREACH(edge_descriptor ed, V){
        if(get(iscmap,ed)){
          out << int(source(ed,sm)) << " " << int(target(ed,sm)) <<  " ";
        }
      }
      out << std::endl;
    }

    std::cerr << "# constrained edges in the partition = " << cc_edges_count << std::endl;
    double uc_ratio = 0.25;
    double cc_ratio = double(cc_edges_count)/double(count);
    double ratio = uc_ratio + (1.0 - uc_ratio)*cc_ratio;
    std::cerr << cc_ratio << " " << ratio << std::endl;
    assert(ratio < 1.0);

#ifdef INCREASE  
    SMS::Constrained_placement < SMS::LindstromTurk_placement<Component_graph>, ECMap> placement (ecmap);
#else
    SMS::LindstromTurk_placement<Component_graph> placement;
#endif
    SMS::Count_ratio_stop_predicate<Component_graph> stop(ratio);
    Bool_map<ECMap> becmap(ecmap);
    SMS::edge_collapse(cg,
                       stop
                       ,CGAL::parameters::vertex_index_map(get(boost::vertex_index,sm))
                       .halfedge_index_map(himap)
                       .get_placement(placement)
                       .edge_is_constrained_map(becmap)
                       .vertex_point_map(get(CGAL::vertex_point,sm)));

    
    
    std::cerr << "|| done" << std::endl;
  }
};


int main(int argc, char** argv ) 
{
  typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
  typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor;
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

  Surface_mesh::Property_map<edge_descriptor,char> ecmap 
    = sm.add_property_map<edge_descriptor,char>("e:constrained",false).first;

  
  HIMap himap = sm.add_property_map<halfedge_descriptor,int>("h:index_in_cc",-1).first;
  
  std::vector<edge_descriptor> partition_edges;
  std::size_t ncc = 0;
  if(argc>2){
    std::ifstream is2(argv[2]);
    std::string line;
    int id, id2;
    
    if(!std::getline(is2, line)){
      std::cerr << "error in selection: no first line"<< std::endl;
      return 1; 
    }
    if(!std::getline(is2, line)){
      std::cerr << "error in selection: no second line"<< std::endl;
      return 1; 
    }
    if(!std::getline(is2, line)){
      std::cerr << "error in selection: no third line"<< std::endl;
      return 1; 
    }
    std::istringstream edge_line(line);
    while(edge_line >> id >> id2) {
      vertex_descriptor s(id), t(id2);
      halfedge_descriptor hd;
      bool exists;
      boost::tie(hd,exists) = halfedge(s,t,sm);
      if(! exists){
        std::cerr << "error in selection: no edge" << s << " " << t << std::endl;
        return 1; 
      }
      partition_edges.push_back(edge(hd,sm));
      ecmap[edge(hd,sm)] = true;
      himap[hd] = 0;
      himap[opposite(hd,sm)] = 1;
    }

    ncc = PMP::connected_components(sm,
                                    ccmap,
                                    PMP::parameters::edge_is_constrained_map(ecmap));
 
  } else {
    ncc = 8;
    PMP::partition(sm, ccmap, static_cast<int>(ncc));

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
  }
  
  std::cerr << t.time() << " sec.\n";
  t.reset();


  // For each connected component find one halfedge on the "border"
  // as a starting point for each simplification thread
  std::vector<halfedge_descriptor> cc_seed(ncc);
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(sm)){
    if(is_border(hd,sm)){
      continue;
    }
    if(is_border(opposite(hd,sm),sm)){
      cc_seed[ccmap[face(hd,sm)]] = hd;
    }else{
      std::size_t hcc = ccmap[face(hd,sm)];
      std::size_t hoppcc = ccmap[face(opposite(hd,sm),sm)];
      if(hcc != hoppcc){
        cc_seed[hcc] = hd;
      }
    }
  } 
  
  std::vector<int> buffer_size(ncc); // the number of constrained edges inside each component
  std::cerr << "#CC = " << ncc << "  in " << t.time() << " sec." << std::endl;
  t.reset();

#if 1
  // Simplify the partition in parallel
  tbb::task_group tasks;
  for(int i = 0; i < ncc; i++){
    tasks.run(Simplify(sm, himap, ecmap, ccmap, buffer_size, cc_seed[i], i));
  }

  tasks.wait();
#else
  for(int i = 0; i < ncc; i++){
    tbb::task_group tasks;
    tasks.run(Simplify(sm, himap, ecmap, ccmap, buffer_size, cc_seed[i], i));
    tasks.wait();
  }
#endif

  std::cerr << "parallel edge collapse in " << t.time() << " sec." << std::endl;
  t.reset();

  if(dump){
    Surface_mesh sm2;
    CGAL::copy_face_graph(sm,sm2);
    
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
  SMS::Constrained_placement <SMS::LindstromTurk_placement<Surface_mesh>, IECMap> placement (iecmap);
#else
  SMS::LindstromTurk_placement<Surface_mesh> placement;
#endif
  double uc_ratio = 0.3;
  double cc_ratio = double(num_constrained_edges)/double(num_edges(sm));
  double ratio = uc_ratio + (1.0 - uc_ratio)*cc_ratio;
  assert(ratio < 1.0);
   
  std::cout << num_constrained_edges << " constrained; cc_ratio: " << cc_ratio << std::endl;
  std::cout << "ratio : " << ratio << std::endl;

  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);
  SMS::edge_collapse(sm,
                     stop
                     ,CGAL::parameters::vertex_index_map(get(boost::vertex_index,sm))
                     .get_placement(placement)
                     .edge_is_constrained_map(iecmap));

  std::cerr << "sequential edge collapse on buffer in " << t.time() << " sec." << std::endl;
  t.reset();
  
#endif
  std::cerr << "collect garbage" << std::endl;

  sm.collect_garbage();
  std::cerr << "\nResult: #V = "<< num_vertices(sm)  << " #E = "<< num_edges(sm) 
            << "  #F = " << num_faces(sm) << "  in " << t.time() << " sec." << std::endl;

  std::ofstream out("out.off");
  out << sm << std::endl;
  out.close();
  std::cerr << "Writing result in " << t.time() << " sec." << std::endl;

  return EXIT_SUCCESS;
}
