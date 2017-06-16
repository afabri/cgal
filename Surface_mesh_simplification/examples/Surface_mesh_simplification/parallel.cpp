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
    return ! get(iecm.ecmap, ed);
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

  In_same_component_map(ECMap ecmap, CCMap ccmap, const G& g, int cc)
    : ecmap(ecmap), ccmap(ccmap), g(g), cc(cc)
  {}


  friend value_type get(const In_same_component_map& iscm, key_type ed)
  {
    return get(iscm.ecmap,ed);
  }
  
  
  friend void put(const In_same_component_map& iscm, key_type ed, value_type v)
  {
    if( (! is_border(halfedge(ed,iscm.g),iscm.g) && iscm.ccmap[face(halfedge(ed,iscm.g),iscm.g)]== iscm.cc)||
        (! is_border(opposite(halfedge(ed,iscm.g),iscm.g),iscm.g) && iscm.ccmap[face(opposite(halfedge(ed,iscm.g),iscm.g),iscm.g)]== iscm.cc) ){
      put(iscm.ecmap, ed, v);
    }
  }
  
};


// The parallel task
struct Simplify {

  Surface_mesh &sm;
  HIMap himap;
  ECMap ecmap;
  CCMap ccmap;
  boost::graph_traits<Surface_mesh>::halfedge_descriptor hd;
  int ccindex;

  Simplify(Surface_mesh& sm,
           HIMap himap,
           ECMap ecmap,
           CCMap ccmap,
           boost::graph_traits<Surface_mesh>::halfedge_descriptor hd, int ccindex)
    : sm(sm), himap(himap), ecmap(ecmap), ccmap(ccmap), hd(hd), ccindex(ccindex)
  {}


  void operator()() const
  {
    typedef Selection_is_constraint_map<HIMap, Surface_mesh> SICM;
    typedef CGAL::Component_graph<Surface_mesh,SICM>  Component_graph;
    SICM sicm(himap,sm);
    Component_graph cg(sm, sicm, hd); 

    std::vector<boost::graph_traits<Surface_mesh>::edge_descriptor> cc_edges;
    int i = 2, count = 0;
    BOOST_FOREACH(boost::graph_traits<Surface_mesh>::edge_descriptor ed, edges(cg)){
      ++count;
      boost::graph_traits<Surface_mesh>::halfedge_descriptor hd, hop;
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
    std::cerr << "cc_edges_count = " << cc_edges_count << std::endl;

    std::vector<boost::graph_traits<Surface_mesh>::edge_descriptor> V;
    In_same_component_map<ECMap,CCMap,Surface_mesh> iscmap(ecmap,ccmap,sm,ccindex);

    CGAL::expand_edge_selection(cc_edges,
                                sm,
                                2,
                                iscmap,
                                std::back_inserter(V));
    cc_edges_count += V.size();
    std::cerr << "cc_edges_count = " << cc_edges_count << std::endl;
    double uc_ratio = 0.25;
    double cc_ratio = double(cc_edges_count)/double(count);
    double ratio = uc_ratio + (1.0 - uc_ratio)*cc_ratio;
    std::cerr << cc_ratio << " " << ratio << std::endl;
    assert(ratio < 1.0);
  
#if 1
   SMS::Constrained_placement
     <SMS::LindstromTurk_placement<Component_graph>,
      ECMap>
          placement (ecmap);
    SMS::Count_ratio_stop_predicate<Component_graph> stop(ratio);
    SMS::edge_collapse(cg,
                       stop
                       ,CGAL::parameters::vertex_index_map(get(boost::vertex_index,sm))
                       .halfedge_index_map(himap)
                       .get_placement(placement)
                       .edge_is_constrained_map(ecmap)
                       .vertex_point_map(get(CGAL::vertex_point,sm)));
#endif
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
  
  int ncc = 0;
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
      ecmap[edge(hd,sm)] = true;
      himap[hd] = 0;
      himap[opposite(hd,sm)] = 1;
    }

    ncc = PMP::connected_components(sm,
                                    ccmap,
                                    PMP::parameters::edge_is_constrained_map(ecmap));
 
  } else {
    ncc = 8;
    PMP::partition(sm, ccmap, ncc);

    // now we have to find the constrained edges
    BOOST_FOREACH(halfedge_descriptor hd, halfedges(sm)){
      if(is_border(edge(hd,sm),sm)){
        continue;
      }
      if (ccmap[face(hd,sm)] !=ccmap[face(opposite(hd,sm),sm)]) {
        ecmap[edge(hd,sm)] = true;
        himap[hd] = 0;
        himap[opposite(hd,sm)] = 1;
      }
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
      int hcc = ccmap[face(hd,sm)];
      int hoppcc = ccmap[face(opposite(hd,sm),sm)];
      if(hcc != hoppcc){
        cc_seed[hcc] = hd;
      }
    }
  } 
  std::cerr << "#CC = " << ncc << "  in " << t.time() << " sec." << std::endl;
  t.reset();
#if 1
  tbb::task_group tasks;
  for(int i = 0; i < ncc; i++){
    tasks.run(Simplify(sm, himap, ecmap, ccmap, cc_seed[i],i));
  }

  tasks.wait();
#else
  for(int i = 0; i < ncc; i++){
    tbb::task_group tasks;
    tasks.run(Simplify(sm, himap, ecmap, ccmap, cc_seed[i],i));
    tasks.wait();
  }
#endif

  std::cerr << "parallel edge collapse in " << t.time() << " sec." << std::endl;
  t.reset();
  
 

#if 1

 typedef Inverted_edge_constraint_map<ECMap> IECMap;
 IECMap iecmap(ecmap);

 
  double icc_count = 0;
  BOOST_FOREACH(edge_descriptor ed, edges(sm)){
    if(get(iecmap,ed)){
      icc_count += 1;
    }
  }

  SMS::Constrained_placement
    <SMS::LindstromTurk_placement<Surface_mesh>,
     IECMap>
    placement (iecmap);
  
  double uc_ratio = 0.25;
  double cc_ratio = double(icc_count)/double(sm.number_of_edges());
  double ratio = uc_ratio + (1.0 - uc_ratio)*cc_ratio;
  assert(ratio < 1.0);
   
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.25);
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
