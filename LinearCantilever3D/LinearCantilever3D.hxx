// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1> Systems Example 6 - 3D Linear Elastic Cantilever </h1>
// \author David Knezevic
// \date 2012
//
// This is a 3D version of systems_of_equations_ex4. The weak form PDE for
// equilibrium elasticity is:
//
//     \int_\Omega Sigma_ij v_i,j = \int_\Omega f_i v_i + \int_\Gamma g_i v_i ds,
//
// for all admissible test functions v, where:
//  * Sigma is the stress tensor, which for linear elasticity is
//    given by Sigma_ij = C_ijkl u_k,l.
//  * f is a body load.
//  * g is a surface traction on the surface \Gamma.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/utility.h"

// boundary IDs
#define BOUNDARY_ID_MIN_Z 0
#define BOUNDARY_ID_MIN_Y 1
#define BOUNDARY_ID_MAX_X 2
#define BOUNDARY_ID_MAX_Y 3
#define BOUNDARY_ID_MIN_X 4
#define BOUNDARY_ID_MAX_Z 5
#define NODE_BOUNDARY_ID 10
#define EDGE_BOUNDARY_ID 20
#define PUSH_BOUNDARY_ID 30

#ifdef LIBMESH_HAVE_PETSC
// This class allows us to set the solver and preconditioner
// to be appropriate for linear elasticity.
class PetscSolverConfiguration : public libMesh::SolverConfiguration
{
public:

  PetscSolverConfiguration(libMesh::PetscLinearSolver<libMesh::Number> & petsc_linear_solver) :
    _petsc_linear_solver(petsc_linear_solver)
  {
  }

  virtual void configure_solver();

  // The linear solver object that we are configuring
  libMesh::PetscLinearSolver<libMesh::Number> & _petsc_linear_solver;

};
#endif

class LinearElasticity : public libMesh::System::Assembly
{
private:
  libMesh::EquationSystems & es;

public:

  LinearElasticity (libMesh::EquationSystems & es_in) :
    es(es_in)
  {}

  /**
   * Kronecker delta function.
   */
  libMesh::Real kronecker_delta(unsigned int i,
                       unsigned int j);

  /**
   * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
   */
  libMesh::Real elasticity_tensor(unsigned int i,
                         unsigned int j,
                         unsigned int k,
                         unsigned int l);

  /**
   * Assemble the system matrix and right-hand side vector.
   */
  void assemble();

  // Post-process the solution to compute stresses
  void compute_stresses();
};

