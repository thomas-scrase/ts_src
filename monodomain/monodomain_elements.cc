#include "monodomain_elements.h"

// Essentially the same as the implementation in GeneralisedAdvectionDiffusionEquations but with
// zero wind, a capacitance scaling term on the time-derivative,
// and the (analogous numbner to the) Peclet number scales the source term
template<unsigned DIM>
void MonodomainEquations::fill_in_generic_residual_contribution_monodomain(Vector<double>& residuals,
                                                                           DenseMatrix<double>& jacobian,
                                                                           DenseMatrix<double>& mass_matrix,
                                                                           unsigned flag)
{
 // Find out how many nodes there are
 const unsigned n_node = nnode();

 // Get the nodal index at which the unknown is stored
 const unsigned u_nodal_index = u_index_monodomain();

 // Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

 // Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 // Set the Vector to hold local coordinates
 Vector<double> s(DIM);

 // Get Peclet number
 const double peclet = omega();

 // Get the Peclet*Strouhal number
 const double peclet_st = omega() * gamma();

 // Integers used to store the local equation number and local unknown
 // indices for the residuals and jacobians
 int local_eqn = 0, local_unknown = 0;

 // Loop over the integration points
 for (unsigned ipt = 0; ipt < n_intpt; ipt++)
 {
  // Assign values of s
  for (unsigned i = 0; i < DIM; i++) s[i] = integral_pt()->knot(ipt, i);

  // Get the integral weight
  double w = integral_pt()->weight(ipt);

  // Call the derivatives of the shape and test functions
  double J = dshape_and_dtest_eulerian_at_knot_monodomain(ipt, psi, dpsidx, test, dtestdx);

  // Premultiply the weights and the Jacobian
  double W = w * J;

  // Calculate local values of the solution and its derivatives
  // Allocate
  double interpolated_u = 0.0;
  double dudt = 0.0;
  Vector<double> interpolated_x(DIM, 0.0);
  Vector<double> interpolated_dudx(DIM, 0.0);
  Vector<double> mesh_velocity(DIM, 0.0);


  // Calculate function value and derivatives:
  //-----------------------------------------
  // Loop over nodes
  for (unsigned l = 0; l < n_node; l++)
  {
   // Get the value at the node
   double u_value = raw_nodal_value(l, u_nodal_index);
   interpolated_u += u_value * psi(l);
   dudt += du_dt_monodomain(l) * psi(l);
   // Loop over directions
   for (unsigned j = 0; j < DIM; j++)
   {
    interpolated_x[j] += raw_nodal_position(l, j) * psi(l);
    interpolated_dudx[j] += u_value * dpsidx(l, j);
   }
  }

  // Mesh velocity?
  if (!ALE_is_disabled)
  {
   for (unsigned l = 0; l < n_node; l++)
   {
    for (unsigned j = 0; j < DIM; j++)
    {
     mesh_velocity[j] += raw_dnodal_position_dt(l, j) * psi(l);
    }
   }
  }


  // Get source function
  //-------------------
  double source;
  get_source_monodomain(ipt, interpolated_x, source);

  // Get diffusivity tensor
  DenseMatrix<double> D(DIM, DIM);
  get_diff_monodomain(ipt, s, interpolated_x, D);

  double capacitance;
  get_capacitance_monodomain(ipt, interpolated_x, capacitance);

  // Assemble residuals and Jacobian
  //--------------------------------

  // Loop over the test functions
  for (unsigned l = 0; l < n_node; l++)
  {
   // Set the local equation number
   local_eqn = nodal_local_eqn(l, u_nodal_index);

   /*IF it's not a boundary condition*/
   if (local_eqn >= 0)
   {
    // Add body force/source term and time derivative
    residuals[local_eqn] -= (peclet_st * capacitance * dudt + peclet * source) * test(l) * W;

    // The Generalised Advection Diffusion bit itself
    for (unsigned k = 0; k < DIM; k++)
    {
     // Terms that multiply the test function
     // divergence-free wind
     double tmp = 0.0;
     // If the mesh is moving need to subtract the mesh velocity
     if (!ALE_is_disabled)
     {
      tmp -= peclet_st * capacitance * mesh_velocity[k];
     }
     tmp *= interpolated_dudx[k];

     // Terms that multiply the derivative of the test function
     // Advective term
     double tmp2 = 0.0;
     // Now the diuffusive term
     for (unsigned j = 0; j < DIM; j++)
     {
      tmp2 += D(k, j) * interpolated_dudx[j];
     }
     // Now construct the contribution to the residuals
     residuals[local_eqn] -= (tmp * test(l) + tmp2 * dtestdx(l, k)) * W;
    }

    // Calculate the jacobian
    //-----------------------
    if (flag)
    {
     // Loop over the velocity shape functions again
     for (unsigned l2 = 0; l2 < n_node; l2++)
     {
      // Set the number of the unknown
      local_unknown = nodal_local_eqn(l2, u_nodal_index);

      // If at a non-zero degree of freedom add in the entry
      if (local_unknown >= 0)
      {
       // Mass matrix term
       jacobian(local_eqn, local_unknown) -=
         peclet_st * capacitance * test(l) * psi(l2) *
         node_pt(l2)->time_stepper_pt()->weight(1, 0) * W;

       // Add the mass matrix term
       if (flag == 2)
       {
        mass_matrix(local_eqn, local_unknown) +=
          peclet_st * capacitance * test(l) * psi(l2) * W;
       }

       // Add contribution to Elemental Matrix
       for (unsigned k = 0; k < DIM; k++)
       {
        // Temporary term used in assembly
        double tmp = 0.0;
        if (!ALE_is_disabled)
        {
         tmp -= peclet_st * capacitance * mesh_velocity[k];
        }
        tmp *= dpsidx(l2, k);

        double tmp2 = 0.0;
        // Now the diffusive term
        for (unsigned j = 0; j < DIM; j++)
        {
         tmp2 += D(k, j) * dpsidx(l2, j);
        }

        // Now assemble Jacobian term
        jacobian(local_eqn, local_unknown) -=
          (tmp * test(l) + tmp2 * dtestdx(l, k)) * W;
       }
      }
     }
    }
   }
  }
 } // End of loop over integration points
}