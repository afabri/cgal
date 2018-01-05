#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/partition.h>

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  typedef CGAL::Simple_cartesian<double>                             K;

  std::ifstream in((argc>1) ? argv[1] : "data/mech-holes-shark.off");
  if(!in)
  {
    std::cerr << "Error: could not read input file" << std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Surface_mesh<K::Point_3>                           SM;
  SM sm;
  CGAL::read_off(in, sm);

  // The face <--> partition_id property map
  typedef SM::Property_map<SM::Face_index, std::size_t>            Face_id_map;
  Face_id_map partition_id_map = sm.add_property_map<SM::Face_index, std::size_t>("f:pid").first;

  // Partition the mesh and output its parts
  PMP::partition(sm, 8 /*number of partitions*/, partition_id_map);

  // Extract the part n°4 of the partition into a new, independent mesh
  typedef CGAL::Face_filtered_graph<SM>                       Filtered_graph;
  Filtered_graph filtered_sm(sm, 4 /*id of th part*/, partition_id_map);
  SM part_sm;
  CGAL::copy_face_graph(filtered_sm, part_sm);

  std::ofstream out("sm_part_4.off");
  CGAL::write_off(out, part_sm);

  return EXIT_SUCCESS;
}

