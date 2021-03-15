#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_solver_package.h"


// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the EquationSystems object and the
// name of the system we are assembling as input.  From the
// EquationSystems object we have access to the Mesh and
// other objects we might need.
void assemble_poisson(libMesh::EquationSystems & es, const std::string & system_name);

// Exact solution function prototype.
libMesh::Real exact_solution (const libMesh::Real x, const libMesh::Real y, const libMesh::Real z = 0.);

// Define a wrapper for exact_solution that will be needed below
void exact_solution_wrapper (libMesh::DenseVector<libMesh::Number> & output, const libMesh::Point & p, const libMesh::Real);

void poisson(libMesh::Mesh & mesh);
