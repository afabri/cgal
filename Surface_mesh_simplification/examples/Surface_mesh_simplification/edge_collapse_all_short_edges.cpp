#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>                          Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  if(argc<3)
  {
    std::cerr << "Usage: " << argv[0] << " input.off minimal_edge_length [out.off]\n";
    return EXIT_FAILURE;
  }

  Surface_mesh sm;
  std::ifstream is(argv[1]);
  is >> sm;
  double threshold = atof(argv[2]);

  if (!CGAL::is_triangle_mesh(sm))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  int r = SMS::edge_collapse(sm,
                             CGAL::Surface_mesh_simplification::Edge_length_stop_predicate<double>(threshold),
                             CGAL::Sequential_tag(),
                             CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index, sm))
                                              .halfedge_index_map(get(CGAL::halfedge_external_index, sm))
                                              .get_cost(SMS::Edge_length_cost<Surface_mesh>())
                                              .get_placement(SMS::Midpoint_placement<Surface_mesh>()));

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
            << (sm.size_of_halfedges()/2) << " final edges.\n";

  std::ofstream os(argc > 3 ? argv[3] : "out.off");
  os << sm;

  return EXIT_SUCCESS;
}
