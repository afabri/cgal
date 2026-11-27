#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/IO/File_medit.h>
#include <CGAL/Real_timer.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/sphere.mesh";
  const double target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1;

  std::cout << "using target_edge_length = " << target_edge_length << "\n";

  Remeshing_triangulation tr;

  std::ifstream is(filename, std::ios_base::in);
  CGAL::IO::read_MEDIT(is, tr);

  CGAL::Real_timer timer;
  timer.start();
  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length, CGAL::parameters::number_of_iterations(5));
  timer.stop();
  std::cout << "Done in " << timer.time() << "\n";

  std::ofstream os("after_remeshing.mesh");
  CGAL::IO::write_MEDIT(os, tr);

  return EXIT_SUCCESS;
}
