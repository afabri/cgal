#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/partition.h>

namespace PMP = CGAL::Polygon_mesh_processing;

int main()
{
  typedef CGAL::Simple_cartesian<double>                           K;
  typedef CGAL::Surface_mesh<K::Point_3>                           SM;
  typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>  PM;

  // with a polyhedron mesh
  PM pm;
  std::ifstream in("data/blobby.off");
  in >> pm;
  CGAL::set_halfedgeds_items_id(pm);
  PMP::partition(pm);

  // with a surface mesh, a property map, and a custom number of iterations
  SM sm;
  CGAL::read_off(sm, "data/blobby.off");

  typedef SM::Property_map<SM::Face_index, std::size_t> Facet_int_map;
  Facet_int_map partition_property_map =
    sm.add_property_map<SM::Face_index, std::size_t>("f:pid").first;

  PMP::partition(sm, partition_property_map, 5);
}

