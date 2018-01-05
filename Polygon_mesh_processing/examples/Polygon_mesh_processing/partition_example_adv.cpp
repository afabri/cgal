#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
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

  typedef CGAL::Polyhedron_3<K>  PM;
  PM pm;
  in >> pm; // read the mesh

  // The face <--> partition_id property map  
  typedef CGAL::dynamic_face_property_t<int> Face_property_tag;
  typedef typename boost::property_map<PM, Face_property_tag>::type Face_id_map;
  Face_id_map partition_id_map = get(Face_property_tag(), pm);

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
  PMP::partition(pm, 8 /*number of partitions*/,
                 partition_id_map,
                 CGAL::parameters::METIS_options(&options).
                 vertex_index_map(get(boost::vertex_external_index, pm)) );

  return EXIT_SUCCESS;
}

