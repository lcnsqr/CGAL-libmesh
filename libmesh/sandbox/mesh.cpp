// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"

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
  // to build a mesh of Hex27
  // elements in 3D.  Building these higher-order elements allows
  // us to use higher-order approximation, as in example 3.

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

  // Iterate through nodes
  for (auto & node : mesh.local_node_ptr_range())
  {
    node->print_info();
  }

  mesh.write("teste.mesh");

  // All done.
  return 0;
}
