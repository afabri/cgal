
int CUTOFF_FACTOR;

#define BID_PARALLEL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Real_timer.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  CGAL::Real_timer timer;
  
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";

  
  CUTOFF_FACTOR = (argc > 3)? atoi(argv[3]) : 2;
  
  std::ifstream input(filename1);

  Mesh mesh1, mesh2;
  if (!input || !(input >> mesh1))
  {
    std::cerr << "First mesh is not a valid off file." << std::endl;
    return 1;
  }
  input.close();
  input.open(filename2);
  if (!input || !(input >> mesh2))
  {
    std::cerr << "Second mesh is not a valid off file." << std::endl;
    return 1;
  }

  std::cout << "Number of vertices before corefinement "
            << num_vertices(mesh1) << " and "
            << num_vertices(mesh2) << "\n";

  timer.start();
  for(int i =0; i < 1; ++i){
    PMP::corefine(mesh1,mesh2);
  }
  timer.stop();
  std::cout << timer.time() << " sec."<< std::endl;
  
  std::cout << "Number of vertices after corefinement "
            << num_vertices(mesh1) << " and "
            << num_vertices(mesh2) << "\n";
  /*
  std::ofstream output("mesh1_refined.off");
  output.precision(17);
  output << mesh1;
  output.close();
  output.open("mesh2_refined.off");
  output << mesh2;
  */
  return 0;
}
