#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/partition.h>

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  typedef CGAL::Simple_cartesian<double>                           K;

#if 1
  // with a polyhedron mesh
  typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>  PM;
  PM pm;
  std::ifstream in((argc>1)?argv[1]:"data/blobby.off");
  in >> pm;
  CGAL::set_halfedgeds_items_id(pm);
  PMP::partition(pm);
#else
  // with a surface mesh, a property map, and a custom number of iterations
  typedef CGAL::Surface_mesh<K::Point_3> SM;
  SM sm;
  CGAL::read_off(sm, (argc>1)?argv[1]:"data/blobby.off");

  typedef SM::Property_map<SM::Face_index, std::size_t> Facet_int_map;
  Facet_int_map partition_property_map =
    sm.add_property_map<SM::Face_index, std::size_t>("f:pid").first;

  PMP::partition(sm, partition_property_map, 8);
#endif
  return 0;
}

