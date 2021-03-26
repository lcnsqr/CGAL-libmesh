#include "LinearCantilever3D.hxx"

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifdef LIBMESH_HAVE_PETSC
void PetscSolverConfiguration::configure_solver()
{
  PetscErrorCode ierr = 0;
  ierr = KSPSetType (_petsc_linear_solver.ksp(), const_cast<KSPType>(KSPCG));
  CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);

  ierr = PCSetType (_petsc_linear_solver.pc(), const_cast<PCType>(PCBJACOBI));
  CHKERRABORT(_petsc_linear_solver.comm().get(), ierr);
}
#endif

Real LinearElasticity::kronecker_delta(unsigned int i,
                     unsigned int j)
{
  return i == j ? 1. : 0.;
}

/**
 * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
 */
Real LinearElasticity::elasticity_tensor(unsigned int i,
                       unsigned int j,
                       unsigned int k,
                       unsigned int l)
{
  // Hard code material parameters for the sake of simplicity
  const Real poisson_ratio = 0.3;
  const Real young_modulus = 1.;

  // Define the Lame constants
  const Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
  const Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));

  return lambda_1 * kronecker_delta(i, j) * kronecker_delta(k, l) +
    lambda_2 * (kronecker_delta(i, k) * kronecker_delta(j, l) + kronecker_delta(i, l) * kronecker_delta(j, k));
}

/**
 * Assemble the system matrix and right-hand side vector.
 */
void LinearElasticity::assemble()
{
  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

  const unsigned int u_var = system.variable_number ("u");

  const DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[3][3] =
    {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
    };

  DenseVector<Number> Fe;

  DenseSubVector<Number> Fe_var[3] =
    {DenseSubVector<Number>(Fe),
     DenseSubVector<Number>(Fe),
     DenseSubVector<Number>(Fe)};

  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_var(3);

  SparseMatrix<Number> & matrix = system.get_system_matrix();

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      dof_map.dof_indices (elem, dof_indices);
      for (unsigned int var=0; var<3; var++)
        dof_map.dof_indices (elem, dof_indices_var[var], var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      for (unsigned int var_i=0; var_i<3; var_i++)
        for (unsigned int var_j=0; var_j<3; var_j++)
          Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);

      Fe.resize (n_dofs);
      for (unsigned int var=0; var<3; var++)
        Fe_var[var].reposition (var*n_var_dofs, n_var_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // assemble \int_Omega C_ijkl u_k,l v_i,j \dx
          for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
            for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
              for (unsigned int i=0; i<3; i++)
                for (unsigned int j=0; j<3; j++)
                  for (unsigned int k=0; k<3; k++)
                    for (unsigned int l=0; l<3; l++)
                      Ke_var[i][k](dof_i,dof_j) +=
                        JxW[qp] * elasticity_tensor(i,j,k,l) * dphi[dof_j][qp](l) * dphi[dof_i][qp](j);

          // assemble \int_Omega f_i v_i \dx
          VectorValue<Number> f_vec(0., -1., 0.);
          //for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
          //  for (unsigned int i=0; i<3; i++)
          //    Fe_var[i](dof_i) += JxW[qp] * (f_vec(i) * phi[dof_i][qp]);
        }

      // assemble \int_\Gamma g_i v_i \ds
      VectorValue<Number> g_vec(0., 0., -1.);
      {
        for (auto side : elem->side_index_range())
          if (elem->neighbor_ptr(side) == nullptr)
            {
              const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
              const std::vector<Real> & JxW_face = fe_face->get_JxW();

              fe_face->reinit(elem, side);

              // Apply a traction
              for (unsigned int qp=0; qp<qface.n_points(); qp++)
                if (mesh.get_boundary_info().has_boundary_id(elem, side, PUSH_BOUNDARY_ID))
                  for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
                    for (unsigned int i=0; i<3; i++)
                      Fe_var[i](dof_i) += JxW_face[qp] * (g_vec(i) * phi_face[dof_i][qp]);
            }
      }

      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      matrix.add_matrix         (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
}

// Post-process the solution to compute stresses
void LinearElasticity::compute_stresses()
{
  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

  unsigned int displacement_vars[3];
  displacement_vars[0] = system.variable_number ("u");
  displacement_vars[1] = system.variable_number ("v");
  displacement_vars[2] = system.variable_number ("w");
  const unsigned int u_var = system.variable_number ("u");

  const DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem & stress_system = es.get_system<ExplicitSystem>("StressSystem");
  const DofMap & stress_dof_map = stress_system.get_dof_map();
  unsigned int sigma_vars[6];
  sigma_vars[0] = stress_system.variable_number ("sigma_00");
  sigma_vars[1] = stress_system.variable_number ("sigma_01");
  sigma_vars[2] = stress_system.variable_number ("sigma_02");
  sigma_vars[3] = stress_system.variable_number ("sigma_11");
  sigma_vars[4] = stress_system.variable_number ("sigma_12");
  sigma_vars[5] = stress_system.variable_number ("sigma_22");
  unsigned int vonMises_var = stress_system.variable_number ("vonMises");

  // Storage for the stress dof indices on each element
  std::vector<std::vector<dof_id_type>> dof_indices_var(system.n_vars());
  std::vector<dof_id_type> stress_dof_indices_var;
  std::vector<dof_id_type> vonmises_dof_indices_var;

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      for (unsigned int var=0; var<3; var++)
        dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);

      const unsigned int n_var_dofs = dof_indices_var[0].size();

      fe->reinit (elem);

      std::vector<TensorValue<Number>> stress_tensor_qp(qrule.n_points());
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Row is variable u1, u2, or u3, column is x, y, or z
          TensorValue<Number> grad_u;
          for (unsigned int var_i=0; var_i<3; var_i++)
            for (unsigned int var_j=0; var_j<3; var_j++)
              for (unsigned int j=0; j<n_var_dofs; j++)
                grad_u(var_i,var_j) += dphi[j][qp](var_j) * system.current_solution(dof_indices_var[var_i][j]);

          for (unsigned int var_i=0; var_i<3; var_i++)
            for (unsigned int var_j=0; var_j<3; var_j++)
              for (unsigned int k=0; k<3; k++)
                for (unsigned int l=0; l<3; l++)
                  stress_tensor_qp[qp](var_i,var_j) += elasticity_tensor(var_i,var_j,k,l) * grad_u(k,l);
        }

      stress_dof_map.dof_indices (elem, vonmises_dof_indices_var, vonMises_var);
      std::vector<TensorValue<Number>> elem_sigma_vec(vonmises_dof_indices_var.size());

      // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
      // Note that this gives a discontinuous stress plot on element boundaries, which is
      // appropriate. We then also get the von Mises stress from the projected stress tensor.
      unsigned int stress_var_index = 0;
      for (unsigned int var_i=0; var_i<3; var_i++)
        for (unsigned int var_j=var_i; var_j<3; var_j++)
          {
            stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[stress_var_index]);

            const unsigned int n_proj_dofs = stress_dof_indices_var.size();

            DenseMatrix<Real> Me(n_proj_dofs, n_proj_dofs);
            for (unsigned int qp=0; qp<qrule.n_points(); qp++)
              {
                for(unsigned int i=0; i<n_proj_dofs; i++)
                  for(unsigned int j=0; j<n_proj_dofs; j++)
                    {
                      Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
                    }
              }

            DenseVector<Number> Fe(n_proj_dofs);
            for (unsigned int qp=0; qp<qrule.n_points(); qp++)
              for(unsigned int i=0; i<n_proj_dofs; i++)
                {
                  Fe(i) += JxW[qp] * stress_tensor_qp[qp](var_i,var_j) * phi[i][qp];
                }

            DenseVector<Number> projected_data;
            Me.cholesky_solve(Fe, projected_data);

            for(unsigned int index=0; index<n_proj_dofs; index++)
              {
                dof_id_type dof_index = stress_dof_indices_var[index];
                if ((stress_system.solution->first_local_index() <= dof_index) &&
                    (dof_index < stress_system.solution->last_local_index()))
                  stress_system.solution->set(dof_index, projected_data(index));

                elem_sigma_vec[index](var_i,var_j) = projected_data(index);
              }

            stress_var_index++;
          }

      for (std::size_t index=0; index<elem_sigma_vec.size(); index++)
        {
          elem_sigma_vec[index](1,0) = elem_sigma_vec[index](0,1);
          elem_sigma_vec[index](2,0) = elem_sigma_vec[index](0,2);
          elem_sigma_vec[index](2,1) = elem_sigma_vec[index](1,2);

          // Get the von Mises stress from the projected stress tensor
          Number vonMises_value = std::sqrt(0.5*(Utility::pow<2>(elem_sigma_vec[index](0,0) - elem_sigma_vec[index](1,1)) +
                                                 Utility::pow<2>(elem_sigma_vec[index](1,1) - elem_sigma_vec[index](2,2)) +
                                                 Utility::pow<2>(elem_sigma_vec[index](2,2) - elem_sigma_vec[index](0,0)) +
                                                 6.*(Utility::pow<2>(elem_sigma_vec[index](0,1)) +
                                                     Utility::pow<2>(elem_sigma_vec[index](1,2)) +
                                                     Utility::pow<2>(elem_sigma_vec[index](2,0)))));

          dof_id_type dof_index = vonmises_dof_indices_var[index];

          if ((stress_system.solution->first_local_index() <= dof_index) &&
              (dof_index < stress_system.solution->last_local_index()))
            stress_system.solution->set(dof_index, vonMises_value);
        }
    }

  // Should call close and update when we set vector entries directly
  stress_system.solution->close();
  stress_system.update();
}
