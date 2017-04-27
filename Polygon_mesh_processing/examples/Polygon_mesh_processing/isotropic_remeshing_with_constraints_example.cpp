#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/cube_quad.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh)) {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  double target_edge_length = 1;

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  Mesh::Property_map<edge_descriptor,bool> constrained;
  constrained = mesh.add_property_map<edge_descriptor,bool>("e:constrained", false).first;

  vertex_descriptor vd = *vertices(mesh).first;
  BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(halfedge(vd,mesh),mesh)){
    constrained[edge(hd,mesh)] = true;
  }

  PMP::triangulate_faces(mesh);

  PMP::isotropic_remeshing(
      faces(mesh),
      target_edge_length,
      mesh,
      PMP::parameters::number_of_iterations(3)
      .number_of_relaxation_steps(3)
      .protect_constraints(true)
      .edge_is_constrained_map(constrained)
      );

  std::ofstream out("remeshed.off");
  out << mesh << std::endl;
  std::cout << "Remeshing done."   << " (" << num_faces(mesh) << " faces)..." << std::endl;

  return 0;
}
