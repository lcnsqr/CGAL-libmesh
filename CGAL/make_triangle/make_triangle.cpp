#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Mesh_polyhedron_3<Kernel>::type                Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  Polyhedron polyhedron;

  polyhedron.make_triangle(Point(-1,-1,-1),Point(0,0,-1),Point(-1,1,-1));
  polyhedron.make_triangle(Point(-1,1,-1),Point(0,0,-1),Point(1,1,-1));
  polyhedron.make_triangle(Point(1,1,-1),Point(0,0,-1),Point(1,-1,-1));
  polyhedron.make_triangle(Point(1,-1,-1),Point(0,0,-1),Point(-1,-1,-1));
  polyhedron.make_triangle(Point(-1,-1,-1),Point(0,-1,0),Point(1,-1,-1));
  polyhedron.make_triangle(Point(1,-1,-1),Point(0,-1,0),Point(1,-1,1));
  polyhedron.make_triangle(Point(1,-1,1),Point(0,-1,0),Point(-1,-1,1));
  polyhedron.make_triangle(Point(-1,-1,1),Point(0,-1,0),Point(-1,-1,-1));
  polyhedron.make_triangle(Point(1,-1,-1),Point(1,0,0),Point(1,1,-1));
  polyhedron.make_triangle(Point(1,1,-1),Point(1,0,0),Point(1,1,1));
  polyhedron.make_triangle(Point(1,1,1),Point(1,0,0),Point(1,-1,1));
  polyhedron.make_triangle(Point(1,-1,1),Point(1,0,0),Point(1,-1,-1));
  polyhedron.make_triangle(Point(1,1,-1),Point(0,1,0),Point(-1,1,-1));
  polyhedron.make_triangle(Point(-1,1,-1),Point(0,1,0),Point(-1,1,1));
  polyhedron.make_triangle(Point(-1,1,1),Point(0,1,0),Point(1,1,1));
  polyhedron.make_triangle(Point(1,1,1),Point(0,1,0),Point(1,1,-1));
  polyhedron.make_triangle(Point(-1,1,-1),Point(-1,0,0),Point(-1,-1,-1));
  polyhedron.make_triangle(Point(-1,-1,-1),Point(-1,0,0),Point(-1,-1,1));
  polyhedron.make_triangle(Point(-1,-1,1),Point(-1,0,0),Point(-1,1,1));
  polyhedron.make_triangle(Point(-1,1,1),Point(-1,0,0),Point(-1,1,-1));
  polyhedron.make_triangle(Point(-1,-1,1),Point(0,0,1),Point(1,-1,1));
  polyhedron.make_triangle(Point(1,-1,1),Point(0,0,1),Point(1,1,1));
  polyhedron.make_triangle(Point(1,1,1),Point(0,0,1),Point(-1,1,1));
  polyhedron.make_triangle(Point(-1,1,1),Point(0,0,1),Point(-1,-1,1));

  if (!CGAL::is_triangle_mesh(polyhedron)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain
  Mesh_domain domain(polyhedron);

  // Get sharp features
  domain.detect_features();

  // Mesh criteria
  Mesh_criteria criteria(edge_size = 0.1,
                         facet_angle = 25, facet_size = 0.2, facet_distance = 0.02,
                         cell_radius_edge_ratio = 3, cell_size = 0.2);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  // Output
  std::ofstream file("out.mesh");
  c3t3.output_to_medit(file);

  std::cout << "Sucesso." << std::endl;

  return EXIT_SUCCESS;
}
