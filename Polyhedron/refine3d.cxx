#include "refine3d.hxx"

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries.
  LibMeshInit init (argc, argv);

  // Dimension
  int dim = 3;

  // Create a mesh with user-defined dimension.
  // Number of elements
  int ps = 1;

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());

  // Use the MeshTools::Generation mesh generator to create a uniform
  // grid on the square [-1,1]^D.  We instruct the mesh generator
  // to build a mesh of TET4 elements in 3D.

  Real halfwidth = 1.;
  Real halfheight = 1.;

  MeshTools::Generation::build_cube (mesh,
                                     ps,
                                     (dim>1) ? ps : 0,
                                     (dim>2) ? ps : 0,
                                     -1., 1.,
                                     -halfwidth, halfwidth,
                                     -halfheight, halfheight,
                                     TET4);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Todos os nós que estão no contorno
  std::unordered_set<dof_id_type> boundary_nodes;
  boundary_nodes = libMesh::MeshTools::find_boundary_nodes(mesh);

  // Poliedro CGAL
  Polyhedron polyhedron;

  // Passar por todos os elementos ativos da malha original
  MeshBase::const_element_iterator el = mesh.active_elements_begin(), end_el = mesh.active_elements_end();
  for ( ; el != end_el ; el++ )
	{
		const Elem * elem = *el;

    // Ignorar se não for um elemento de contorno
    if ( ! elem->on_boundary() ) continue;

    //libMesh::out << "Face externa do elemento de contorno " << elem->id() << ":" << std::endl;

    // Os três vértices de uma face de contorno
    // Kernel::Point_3 é o tipo ponto 3D da CGAL
    std::array<Kernel::Point_3, 3> vertices;

    // Passar por todos os nós do elemento
    uint v = 0;
    for (uint nnid=0; nnid < elem->n_nodes(); nnid++)
		{
      // Obter o ID do nó
      const Node * nptr = elem->node_ptr(nnid);
      uint nid = nptr->id();

      // Verificar se nó está na face externa do elemento
      if (boundary_nodes.find(nid) != boundary_nodes.end()){

        // Coordenadas do nó de contorno
        libMesh::Point p = *nptr;
        //libMesh::out << "Nó de contorno (" << nid << "): " << p(0) << ", " << p(1) << ", " << p(2) << std::endl;

        // Atualizar vértice correspondente
        vertices[v] = Kernel::Point_3(p(0), p(1), p(2));
        v++;

      }

		}

    // Adicionar triângulo de contorno ao poliedro CGAL
		polyhedron.make_triangle(vertices[0], vertices[1], vertices[2]);
	}

  libMesh::out << "Salvando cube.e" << std::endl;
  mesh.write("cube.e");

  if (!CGAL::is_triangle_mesh(polyhedron)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
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

  std::cout << "number_of_facets: " << c3t3.number_of_facets() << std::endl;
  std::cout << "number_of_edges: " << c3t3.number_of_edges() << std::endl;
  std::cout << "number_of_corners: " << c3t3.number_of_corners() << std::endl;

  // Output
  libMesh::out << "Salvando cube_refined.mesh" << std::endl;
  std::ofstream file("cube_refined.mesh");
  c3t3.output_to_medit(file);


  // All done.
  return 0;
}
