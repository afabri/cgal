#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/read_drc.h>
#include <draco/mesh/corner_table.h>
#include <draco/mesh/mesh_misc_functions.h>

#include <CGAL/boost/graph/graph_traits_draco.h>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

int main(int argc, char **argv)
{
  std::ifstream input_file(argv[1], std::ios::binary);

  // Read the file stream into a buffer.
  std::streampos file_size = 0;
  input_file.seekg(0, std::ios::end);
  file_size = input_file.tellg() - file_size;
  input_file.seekg(0, std::ios::beg);
  std::vector<char> data(file_size);
  input_file.read(data.data(), file_size);

  if (data.empty()) {
    std::cerr << "Empty input file.\n";
    return false;
  }

  // Create a draco decoding buffer. Note that no data is copied in this step.
  draco::DecoderBuffer buffer;
  buffer.Init(data.data(), data.size());

  // Decode the input data into a geometry.
  std::unique_ptr<draco::PointCloud> pc;
  draco::Mesh *mesh = nullptr;
  auto type_statusor = draco::Decoder::GetEncodedGeometryType(&buffer);
  if (!type_statusor.ok()) {
    std::cerr << "Decoding failed.\n";
    return false;
  }
  const draco::EncodedGeometryType geom_type = type_statusor.value();
  if (geom_type == draco::TRIANGULAR_MESH) {
    draco::Decoder decoder;
    auto statusor = decoder.DecodeMeshFromBuffer(&buffer);
    if (!statusor.ok()) {
      std::cerr << "Decoding failed.\n";
      return false;
    }
    std::unique_ptr<draco::Mesh> in_mesh = std::move(statusor).value();

    if (in_mesh) {
      mesh = in_mesh.get();
      pc = std::move(in_mesh);
    }
  }else{
    std::cerr << "Not a triangle mesh.\n";
    return false;
  }


  std::unique_ptr<draco::CornerTable> in_ct = draco::CreateCornerTableFromPositionAttribute(mesh);
  draco::CornerTable& ct = *(in_ct.get());

  CGAL::dMesh dm(*mesh,ct);
  typedef boost::graph_traits<CGAL::dMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<CGAL::dMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<CGAL::dMesh>::vertex_descriptor vertex_descriptor;
  std::cout << "# faces = " << num_vertices(dm) << std::endl;
  BOOST_FOREACH(face_descriptor fd, faces(dm)){
    std::cout << fd << std::endl;
    halfedge_descriptor hd = halfedge(fd,dm);
    vertex_descriptor vd = target(hd,dm);
    halfedge_descriptor hd2 = halfedge(vd,dm);
    vertex_descriptor vd2 = target(hd2,dm);
    assert(vd == vd2);

    hd2 = opposite(hd,dm);
    assert(hd == opposite(hd2,dm));

    boost::property_map<CGAL::dMesh,boost::vertex_point_t>::type vpm = get(boost::vertex_point,dm);
    BOOST_FOREACH(halfedge_descriptor hd, CGAL::halfedges_around_face(halfedge(fd,dm),dm)){
      std::cout << "hd: " << hd.ci << "  " << get(vpm, target(hd,dm))<< std::endl;
    }
  }

  
  /*
    Mesh sm;
    CGAL::read_drc(argv[1],sm);
    std::cout << sm << std::endl;
  */
  return 0;
}
