// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>
#include <unordered_set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/boundary_info.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries, like in example 2.
  LibMeshInit init (argc, argv);

  // Brief message to the user regarding the program name
  // and command line arguments.
  libMesh::out << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

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

	/*
  // Iterate through nodes
  libMesh::out << "Todos os nós da malha:" << std::endl;
  for (auto & node : mesh.local_node_ptr_range())
  {
    libMesh::out << "ID: " << node->id() << ": ";
    libMesh::out << node->slice(0) << " " << node->slice(1) << " " << node->slice(2) << std::endl;
  }
	*/

  // Boundary nodes
  libMesh::out << "Nós na fronteira:" << std::endl;
  std::unordered_set<dof_id_type> block_boundary_nodes;
  block_boundary_nodes = libMesh::MeshTools::find_boundary_nodes(mesh);
  for (auto & id : block_boundary_nodes)
  {
    libMesh::out << "ID: " << id << std::endl;

		Point p = mesh.node_ref(id);
		libMesh::out << p(0) << ", " << p(1) << ", " << p(2) << std::endl;
  }

	/*
  const BoundaryInfo & boundary_info = mesh.get_boundary_info();
  boundary_info.print_summary();

  for (auto & id : boundary_info.get_boundary_ids ())
  {
    libMesh::out << "ID: " << id << ", nome: " << boundary_info.get_sideset_name(id) << std::endl;
  }

  // Element iterator
  MeshBase::const_element_iterator el = mesh.active_elements_begin(), end_el = mesh.active_elements_end();
  for ( ; el != end_el ; el++ )
	{
		const Elem * elem = *el;
    for (uint nnid=0; nnid < elem->n_nodes(); nnid++)
		{
      const Node * nptr = elem->node_ptr(nnid);
      uint nid = nptr->id();

			// Fills a user-provided std::vector with the boundary ids associated with Node nptr
      std::vector<boundary_id_type> vec;
      boundary_info.boundary_ids( nptr, vec );
			if ( ! vec.size() ) continue;

			Point p = *nptr;
			libMesh::out << p(0) << ", " << p(1) << ", " << p(2) << std::endl;
			// Inserção na malha CGAL
			//CPoint cpt = tf->transf( CPoint( p(0), p(1) ) );
			//cdt.insert( cpt );

		}
	}
	*/

  mesh.write("teste.mesh");

  // All done.
  return 0;
}
