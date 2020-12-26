#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;

int main(int argc, char* argv[])
{
  Polyhedron P;

  P.make_triangle(Point(-1,1,1), Point(0,0,1), Point(-1,-1,1));

  return EXIT_SUCCESS;
}