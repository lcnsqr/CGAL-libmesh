#include "poisson.hxx"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void poisson(Mesh & mesh)
{

  // FE order
  std::string order = "FIRST";

  // FE Family
  std::string family = "LAGRANGE";

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Poisson"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Poisson");


  // Add the variable "u" to "Poisson".  "u"
  // will be approximated using second-order approximation.
  unsigned int u_var = system.add_variable("u",
                                           Utility::string_to_enum<Order>   (order),
                                           Utility::string_to_enum<FEFamily>(family));

  // Give the system a pointer to the matrix assembly
  // function.
  system.attach_assemble_function (assemble_poisson);

  // Construct a Dirichlet boundary condition object

  // Indicate which boundary IDs we impose the BC on
  // We either build a line, a square or a cube, and
  // here we indicate the boundaries IDs in each case
  BoundaryInfo & boundary_info = mesh.get_boundary_info();
  std::set<boundary_id_type> boundary_ids;
  boundary_ids = boundary_info.get_boundary_ids();

  // Create a vector storing the variable numbers which the BC applies to
  std::vector<unsigned int> variables(1);
  variables[0] = u_var;

  // Create an AnalyticFunction object that we use to project the BC
  // This function just calls the function exact_solution via exact_solution_wrapper
  AnalyticFunction<> exact_solution_object(exact_solution_wrapper);

#ifdef LIBMESH_ENABLE_DIRICHLET
  // In general, when reusing a system-indexed exact solution, we want
  // to use the default system-ordering constructor for
  // DirichletBoundary, so we demonstrate that here.  In this case,
  // though, we have only one variable, so system- and local-
  // orderings are the same.
  DirichletBoundary dirichlet_bc
    (boundary_ids, variables, exact_solution_object);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
#endif

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();
  mesh.print_info();

  // Solve the system "Poisson", just like example 2.
  system.solve();

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO (mesh).write_equation_systems("solution.e", equation_systems);
#endif

}

// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions.
void assemble_poisson(EquationSystems & es,
                      const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Poisson");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the LinearImplicitSystem we are solving
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Poisson");

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

  // A 5th order Gauss quadrature rule for numerical integration.
  QGauss qrule (dim, FIFTH);

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);

  // Declare a special finite element object for
  // boundary integration.
  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim-1, FIFTH);

  // Tell the finite element object to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  // We begin with the element Jacobian * quadrature weight at each
  // integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The physical XY locations of the quadrature points on the element.
  // These might be useful for evaluating spatially varying material
  // properties at the quadrature points.
  const std::vector<Point> & q_point = fe->get_xyz();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe". More detail is in example 3.
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // The global system matrix
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.  See example 3 for a discussion of the
  // element iterators.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Cache the number of degrees of freedom on this element, for
      // use as a loop bound later.  We use cast_int to explicitly
      // convert from size() (which may be 64-bit) to unsigned int
      // (which may be 32-bit but which is definitely enough to count
      // *local* degrees of freedom.
      const unsigned int n_dofs =
        cast_int<unsigned int>(dof_indices.size());

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // With one variable, we should have the same number of degrees
      // of freedom as shape functions.
      libmesh_assert_equal_to (n_dofs, phi.size());

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);

      Fe.resize (n_dofs);

      // Now we will build the element matrix.  This involves
      // a double loop to integrate the test functions (i) against
      // the trial functions (j).
      //
      // We have split the numeric integration into two loops
      // so that we can log the matrix and right-hand-side
      // computation separately.

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        for (unsigned int i=0; i != n_dofs; i++)
          for (unsigned int j=0; j != n_dofs; j++)
            Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

      // Now we build the element right-hand-side contribution.
      // This involves a single loop in which we integrate the
      // "forcing function" in the PDE against the test functions.

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // fxy is the forcing function for the Poisson equation.
          // In this case we set fxy to be a finite difference
          // Laplacian approximation to the (known) exact solution.
          //
          // We will use the second-order accurate FD Laplacian
          // approximation, which in 3D on a structured grid is
          //
          // u_xx + u_yy = (u(i-1,j) + u(i+1,j) +
          //                u(i,j-1) + u(i,j+1) +
          //                -4*u(i,j))/h^2
          //
          // Since the value of the forcing function depends only
          // on the location of the quadrature point (q_point[qp])
          // we will compute it here, outside of the i-loop
          const Real x = q_point[qp](0);
#if LIBMESH_DIM > 1
          const Real y = q_point[qp](1);
#else
          const Real y = 0.;
#endif
#if LIBMESH_DIM > 2
          const Real z = q_point[qp](2);
#else
          const Real z = 0.;
#endif
          const Real eps = 1.e-3;

          const Real uxx = (exact_solution(x-eps, y, z) +
                            exact_solution(x+eps, y, z) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          const Real uyy = (exact_solution(x, y-eps, z) +
                            exact_solution(x, y+eps, z) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          const Real uzz = (exact_solution(x, y, z-eps) +
                            exact_solution(x, y, z+eps) +
                            -2.*exact_solution(x, y, z))/eps/eps;

          Real fxy;
          if (dim==1)
            {
              // In 1D, compute the rhs by differentiating the
              // exact solution twice.
              const Real pi = libMesh::pi;
              fxy = (0.25*pi*pi)*sin(.5*pi*x);
            }
          else
            {
              fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
            }

          // Add the RHS contribution
          for (unsigned int i=0; i != n_dofs; i++)
            Fe(i) += JxW[qp]*fxy*phi[i][qp];
        }

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      // Also, note that here we call heterogenously_constrain_element_matrix_and_vector
      // to impose a inhomogeneous Dirichlet boundary conditions.
      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.

      matrix.add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

}

// Define a wrapper for exact_solution that will be needed below
void exact_solution_wrapper (DenseVector<Number> & output,
                             const Point & p,
                             const Real)
{
  output(0) = exact_solution(p(0),
                             (LIBMESH_DIM>1)?p(1):0,
                             (LIBMESH_DIM>2)?p(2):0);
}

/**
 * This is the exact solution that
 * we are trying to obtain.  We will solve
 *
 * - (u_xx + u_yy) = f
 *
 * and take a finite difference approximation using this
 * function to get f.  This is the well-known "method of
 * manufactured solutions".
 */
Real exact_solution (const Real x,
                     const Real y,
                     const Real z)
{
  static const Real pi = acos(-1.);

  return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
}
