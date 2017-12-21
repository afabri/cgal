#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/partition.h>

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  typedef CGAL::Simple_cartesian<double>                             K;

  std::ifstream in((argc>1) ? argv[1] : "data/mech-holes-shark.off");
  if(!in)
  {
    std::cerr << "Error: could not read file " << in << std::endl;
    return EXIT_FAILURE;
  }

  // With a polyhedron mesh
  {
    typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3>  PM;
    PM pm;
    in >> pm; // read the mesh
    CGAL::set_halfedgeds_items_id(pm);

    // The face <--> partition_id property map
    typedef boost::property_map<PM, boost::face_external_index_t>::type Face_id_map;
    Face_id_map partition_id_map = get(boost::face_external_index, pm);

    // Partition the mesh and output its parts
    PMP::partition(pm, 8 /*number of partitions*/, partition_id_map);

    // Output all partitions to files named "polyhedron_partition_[0...7].off"
    PMP::output_partition(pm, 8, partition_id_map, "polyhedron_partition");
  }

  // With a surface mesh and custom METIS options
  {
    typedef CGAL::Surface_mesh<K::Point_3>                           SM;
    SM sm;
    in.clear();
    in.seekg(0, std::ios::beg); // reset the reader to the beginning of the file
    CGAL::read_off(in, sm);

    // The face <--> partition_id property map
    typedef SM::Property_map<SM::Face_index, std::size_t>            Face_id_map;
    Face_id_map partition_id_map = sm.add_property_map<SM::Face_index, std::size_t>("f:pid").first;

    // Set some custom options for METIS
    idx_t options[METIS_NOPTIONS];

    // Set all options to default ahead of manually editing some values
    METIS_SetDefaultOptions(options);

    // See METIS documentation for details on the many options
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    options[METIS_OPTION_NCUTS] = 3;
    options[METIS_OPTION_NITER] = 10;
    options[METIS_OPTION_SEED] = 12345;
    options[METIS_OPTION_MINCONN] = 0;
    options[METIS_OPTION_CONTIG] = 0;
    options[METIS_OPTION_UFACTOR] = 1;

    // Partition the mesh and output its parts
    PMP::partition(sm, 8 /*number of partitions*/, partition_id_map, CGAL::parameters::METIS_options(&options));

    // Extract the part n°4 of the partition into a new, independent mesh
    typedef CGAL::Face_filtered_graph<SM>                       Filtered_graph;
    Filtered_graph filtered_sm(sm, 4 /*id of th part*/, partition_id_map);
    SM part_sm;
    CGAL::copy_face_graph(filtered_sm, part_sm);
  }

  return EXIT_SUCCESS;
}

