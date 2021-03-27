//#include "Triangulation.hxx"
#include "LinearCantilever3D.hxx"

// libMesh namespace
using namespace libMesh;

// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // Create a 3D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  MeshTools::Generation::build_cube (mesh,
                                     24,
                                     24,
                                     3,
                                     -.5, .5,
                                     -.5, .5,
                                     -.05, .05,
                                     TET4);

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Let's add some node and edge boundary conditions
  // Each processor should know about each boundary condition it can
  // see, so we loop over all elements, not just local elements.
  for (const auto & elem : mesh.element_ptr_range())
    {

      for (auto s : elem->side_index_range())
        {
          if (mesh.get_boundary_info().has_boundary_id(elem, s, BOUNDARY_ID_MAX_Z))
            {
              // Identificar se há pontos de impacto na face
              int contact_points = 0;
              // Ponteiro para a face do elemento
              std::unique_ptr<Elem> side = elem->side_ptr(s);
              // Percorrer cada nó na face
              for (auto n : side->node_index_range())
              {
                // Ponto (coordenadas) referente ao nó
                Point p = side->node_ref(n);
                if ( sqrt(pow(p(0), 2.)+pow(p(1), 2.)) < .1 )
                {
                  contact_points++;
                }
              }
              if ( contact_points == side->n_nodes() ){
                // Incluir face como parte da área de impacto
                mesh.get_boundary_info().add_side(elem, s, PUSH_BOUNDARY_ID);
              }
            }
        }

    }

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

#ifdef LIBMESH_HAVE_PETSC
  // Attach a SolverConfiguration object to system.linear_solver
  PetscLinearSolver<Number> * petsc_linear_solver =
    cast_ptr<PetscLinearSolver<Number>*>(system.get_linear_solver());
  libmesh_assert(petsc_linear_solver);
  PetscSolverConfiguration petsc_solver_config(*petsc_linear_solver);
  petsc_linear_solver->set_solver_configuration(petsc_solver_config);
#endif

  LinearElasticity le(equation_systems);
  system.attach_assemble_object(le);

#ifdef LIBMESH_ENABLE_DIRICHLET
  // Add three displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
  unsigned int v_var = system.add_variable("v", FIRST, LAGRANGE);
  unsigned int w_var = system.add_variable("w", FIRST, LAGRANGE);

  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(BOUNDARY_ID_MIN_X);
  boundary_ids.insert(BOUNDARY_ID_MAX_X);
  //boundary_ids.insert(NODE_BOUNDARY_ID);
  //boundary_ids.insert(EDGE_BOUNDARY_ID);

  // Create a vector storing the variable numbers which the BC applies to
  std::vector<unsigned int> variables;
  variables.push_back(u_var);
  variables.push_back(v_var);
  variables.push_back(w_var);

  // Create a ZeroFunction to initialize dirichlet_bc
  ZeroFunction<> zf;

  // Most DirichletBoundary users will want to supply a "locally
  // indexed" functor
  DirichletBoundary dirichlet_bc(boundary_ids, variables, zf,
                                 LOCAL_VARIABLE_ORDER);

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
#endif // LIBMESH_ENABLE_DIRICHLET

  // Also, initialize an ExplicitSystem to store stresses
  ExplicitSystem & stress_system =
    equation_systems.add_system<ExplicitSystem> ("StressSystem");

  stress_system.add_variable("sigma_00", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_01", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_02", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_11", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_12", FIRST, L2_LAGRANGE);
  stress_system.add_variable("sigma_22", FIRST, L2_LAGRANGE);
  stress_system.add_variable("vonMises", FIRST, L2_LAGRANGE);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Post-process the solution to compute the stresses
  le.compute_stresses();

  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API

  // Use single precision in this case (reduces the size of the exodus file)
  ExodusII_IO exo_io(mesh, /*single_precision=*/true);
  exo_io.write_discontinuous_exodusII("displacement_and_stress.exo", equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}
