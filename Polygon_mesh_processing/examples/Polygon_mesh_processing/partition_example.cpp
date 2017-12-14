#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Polygon_mesh_processing/partition.h>

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  typedef CGAL::Simple_cartesian<double>                             K;

  // with a polyhedron mesh
  {
    typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>  PM;
    PM pm;
    std::ifstream in((argc>1) ? argv[1] : "data/blobby.off");
    in >> pm;
    CGAL::set_halfedgeds_items_id(pm);

    typedef boost::property_map<PM, boost::face_external_index_t>::type Face_id_map;
    Face_id_map partition_id_map = get(boost::face_external_index, pm);

    PMP::partition(pm, 8 /*number of partitions*/, partition_id_map);
  }

  // with a surface mesh
  {
    typedef CGAL::Surface_mesh<K::Point_3>                           SM;
    SM sm;
    std::ifstream in((argc>1) ? argv[1] : "data/blobby.off");
    CGAL::read_off(in, sm);

    typedef SM::Property_map<SM::Face_index, std::size_t>            Face_id_map;
    Face_id_map partition_id_map = sm.add_property_map<SM::Face_index, std::size_t>("f:pid").first;

    PMP::partition(sm, 8 /*number of partitions*/, partition_id_map);
  }

  return 0;
}

