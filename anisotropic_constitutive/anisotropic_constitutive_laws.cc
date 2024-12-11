#include "anisotropic_constitutive_laws.h"
#include "../generic/elements.h"


namespace oomph
{
//=========================================================================
/// Calculate the derivatives of the contravariant
/// 2nd Piola Kirchhoff stress tensor with respect to the deformed metric
/// tensor. Arguments are the
/// covariant undeformed and deformed metric tensor and the
/// matrix in which to return the derivatives of the stress tensor
/// The default implementation uses finite differences, but can be
/// overloaded for constitutive laws in which an analytic formulation
/// is possible.
//==========================================================================
void AnisotropicConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
 const DenseMatrix<double>& g,
 const DenseMatrix<double>& G,
 const DenseMatrix<double>& sigma,
 const Vector<Vector<double>>& a,
 RankFourTensor<double>& d_sigma_dG,
 const bool& symmetrize_tensor)
{
 // Initial error checking
#ifdef PARANOID
 // Test that the matrices are of the same dimension
 if (!are_matrices_of_equal_dimensions(g, G))
 {
   throw OomphLibError("Matrices passed are not of equal dimension",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
 }

 // Test that there are enough PVA
 if(!verify_number_of_pva(a, N_Principal_Vectors_Of_Anisotropy))
 {
  throw OomphLibError("Number of principal vectors of anisotropy is incorrect",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 
 // Find the dimension of the matrix (assuming that it's square)
 const unsigned dim = G.ncol();

 // Find the dimension
 // FD step
 const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

 // Advanced metric tensor
 DenseMatrix<double> G_pls(dim, dim);
 DenseMatrix<double> sigma_pls(dim, dim);

 // Copy across the original value
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   G_pls(i, j) = G(i, j);
  }
 }

 // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
 // NOTE: We exploit the symmetry of the stress and metric tensors
 //       by incrementing G(i,j) and G(j,i) simultaenously and
 //       only fill in the "upper" triangles without copying things
 //       across the lower triangle. This is taken into account
 //       in the solid mechanics codes.
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = i; j < dim; j++)
  {
   G_pls(i, j) += eps_fd;
   G_pls(j, i) = G_pls(i, j);

   // Get advanced stress
   this->calculate_second_piola_kirchhoff_stress(g, G_pls, a, sigma_pls);

   for (unsigned ii = 0; ii < dim; ii++)
   {
    for (unsigned jj = ii; jj < dim; jj++)
    {
     d_sigma_dG(ii, jj, i, j) =
      (sigma_pls(ii, jj) - sigma(ii, jj)) / eps_fd;
    }
   }

   // Reset
   G_pls(i, j) = G(i, j);
   G_pls(j, i) = G(j, i);
  }
 }

 // If we are symmetrising the tensor, do so
 if (symmetrize_tensor)
 {
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < i; j++)
   {
    for (unsigned ii = 0; ii < dim; ii++)
    {
     for (unsigned jj = 0; jj < ii; jj++)
     {
      d_sigma_dG(ii, jj, i, j) = d_sigma_dG(jj, ii, j, i);
     }
    }
   }
  }
 }
}

//=========================================================================
/// Calculate the derivatives of the contravariant
/// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
/// with respect to the deformed metric tensor.
/// Also return the derivatives of the determinant of the
/// deformed metric tensor with respect to the deformed metric tensor.
/// This form is appropriate
/// for truly-incompressible materials.
/// The default implementation uses finite differences for the
/// derivatives that depend on the constitutive law, but not
/// for the derivatives of the determinant, which are generic.
//========================================================================
void AnisotropicConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
 const DenseMatrix<double>& g,
 const DenseMatrix<double>& G,
 const DenseMatrix<double>& sigma,
 const double& detG,
 const double& interpolated_solid_p,
 const Vector<Vector<double>>& a,
 RankFourTensor<double>& d_sigma_dG,
 DenseMatrix<double>& d_detG_dG,
 const bool& symmetrize_tensor)
{
 // Initial error checking
#ifdef PARANOID
 // Test that the matrices are of the same dimension
 if (!are_matrices_of_equal_dimensions(g, G))
 {
  throw OomphLibError("Matrices passed are not of equal dimension",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }

 // Test that there are enough PVA
 if(!verify_number_of_pva(a, N_Principal_Vectors_Of_Anisotropy))
 {
  throw OomphLibError("Number of principal vectors of anisotropy is incorrect",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Find the dimension of the matrix (assuming that it's square)
 const unsigned dim = G.ncol();

 // FD step
 const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

 // Advanced metric tensor etc.
 DenseMatrix<double> G_pls(dim, dim);
 DenseMatrix<double> sigma_dev_pls(dim, dim);
 DenseMatrix<double> Gup_pls(dim, dim);
 double detG_pls;

 // Copy across
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   G_pls(i, j) = G(i, j);
  }
 }

 // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
 // NOTE: We exploit the symmetry of the stress and metric tensors
 //       by incrementing G(i,j) and G(j,i) simultaenously and
 //       only fill in the "upper" triangles without copying things
 //       across the lower triangle. This is taken into account
 //       in the remaining code further below.
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = i; j < dim; j++)
  {
   G_pls(i, j) += eps_fd;
   G_pls(j, i) = G_pls(i, j);

   // Get advanced stress
   this->calculate_second_piola_kirchhoff_stress(
    g, G_pls, a, sigma_dev_pls, Gup_pls, detG_pls);


   // Derivative of determinant of deformed metric tensor
   d_detG_dG(i, j) = (detG_pls - detG) / eps_fd;

   // Derivatives of deviatoric stress and "upper" deformed metric
   // tensor
   for (unsigned ii = 0; ii < dim; ii++)
   {
    for (unsigned jj = ii; jj < dim; jj++)
    {
     d_sigma_dG(ii, jj, i, j) =
      (sigma_dev_pls(ii, jj) - interpolated_solid_p * Gup_pls(ii, jj) -
       sigma(ii, jj)) /
      eps_fd;
    }
   }

   // Reset
   G_pls(i, j) = G(i, j);
   G_pls(j, i) = G(j, i);
  }
 }

 // If we are symmetrising the tensor, do so
 if (symmetrize_tensor)
 {
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < i; j++)
   {
    d_detG_dG(i, j) = d_detG_dG(j, i);

    for (unsigned ii = 0; ii < dim; ii++)
    {
     for (unsigned jj = 0; jj < ii; jj++)
     {
      d_sigma_dG(ii, jj, i, j) = d_sigma_dG(jj, ii, j, i);
     }
    }
   }
  }
 }
}

//========================================================================
/// Calculate the derivatives of the contravariant
/// 2nd Piola Kirchoff stress tensor with respect to the deformed metric
/// tensor. Also return the derivatives of the generalised dilatation,
/// \f$ d, \f$ with respect to the deformed metric tensor.
/// This form is appropriate
/// for near-incompressible materials.
/// The default implementation uses finite differences.
//=======================================================================
void AnisotropicConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
 const DenseMatrix<double>& g,
 const DenseMatrix<double>& G,
 const DenseMatrix<double>& sigma,
 const double& gen_dil,
 const double& inv_kappa,
 const double& interpolated_solid_p,
 const Vector<Vector<double>>& a,
 RankFourTensor<double>& d_sigma_dG,
 DenseMatrix<double>& d_gen_dil_dG,
 const bool& symmetrize_tensor)
{
 // Initial error checking
#ifdef PARANOID
 // Test that the matrices are of the same dimension
 if (!are_matrices_of_equal_dimensions(g, G))
 {
  throw OomphLibError("Matrices passed are not of equal dimension",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }

 // Test that there are enough PVA
 if(!verify_number_of_pva(a, N_Principal_Vectors_Of_Anisotropy))
 {
  throw OomphLibError("Number of principal vectors of anisotropy is incorrect",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Find the dimension of the matrix (assuming that it's square)
 const unsigned dim = G.ncol();

 // FD step
 const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

 // Advanced metric tensor etc
 DenseMatrix<double> G_pls(dim, dim);
 DenseMatrix<double> sigma_dev_pls(dim, dim);
 DenseMatrix<double> Gup_pls(dim, dim);
 double gen_dil_pls;
 double inv_kappa_pls;

 // Copy across
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   G_pls(i, j) = G(i, j);
  }
 }

 // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
 // NOTE: We exploit the symmetry of the stress and metric tensors
 //       by incrementing G(i,j) and G(j,i) simultaenously and
 //       only fill in the "upper" triangles without copying things
 //       across the lower triangle. This is taken into account
 //       in the remaining code further below.
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = i; j < dim; j++)
  {
   G_pls(i, j) += eps_fd;
   G_pls(j, i) = G_pls(i, j);

   // Get advanced stress
   this->calculate_second_piola_kirchhoff_stress(
    g, G_pls, a, sigma_dev_pls, Gup_pls, gen_dil_pls, inv_kappa_pls);

   // Derivative of generalised dilatation
   d_gen_dil_dG(i, j) = (gen_dil_pls - gen_dil) / eps_fd;

   // Derivatives of deviatoric stress and "upper" deformed metric
   // tensor
   for (unsigned ii = 0; ii < dim; ii++)
   {
    for (unsigned jj = ii; jj < dim; jj++)
    {
     d_sigma_dG(ii, jj, i, j) =
      (sigma_dev_pls(ii, jj) - interpolated_solid_p * Gup_pls(ii, jj) -
       sigma(ii, jj)) /
      eps_fd;
    }
   }

   // Reset
   G_pls(i, j) = G(i, j);
   G_pls(j, i) = G(j, i);
  }
 }

 // If we are symmetrising the tensor, do so
 if (symmetrize_tensor)
 {
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < i; j++)
   {
    d_gen_dil_dG(i, j) = d_gen_dil_dG(j, i);

    for (unsigned ii = 0; ii < dim; ii++)
    {
     for (unsigned jj = 0; jj < ii; jj++)
     {
      d_sigma_dG(ii, jj, i, j) = d_sigma_dG(jj, ii, j, i);
     }
    }
   }
  }
 }
}

//========================================================================
/// Calculate the contravariant 2nd Piola Kirchhoff
/// stress tensor. Arguments are the
/// covariant undeformed and deformed metric tensor and the
/// matrix in which to return the stress tensor.
/// Uses correct 3D invariants for 2D (plane strain) problems.
//=======================================================================
void AnisotropicStrainEnergyFunctionConstitutiveLaw::
 calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                         const DenseMatrix<double>& G,
                                         const Vector<Vector<double>>& a,
                                         DenseMatrix<double>& sigma)
{
// Error checking
#ifdef PARANOID
 error_checking_in_input(g, G, sigma);

 // Test that there are enough PVA
 if(!verify_number_of_pva(a, N_Principal_Vectors_Of_Anisotropy))
 {
  throw OomphLibError("Number of principal vectors of anisotropy is incorrect",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Find the dimension of the problem
 const unsigned dim = g.nrow();

#ifdef PARANOID
 if (dim == 1)
 {
  std::string function_name =
   "AnisotropicStrainEnergyFunctionConstitutiveLaw::";
  function_name += "calculate_second_piola_kirchhoff_stress()";

  throw OomphLibError("Check constitutive equations carefully when dim=1",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Calculate the contravariant undeformed and deformed metric tensors
 // and get the determinants of the metric tensors
 DenseMatrix<double> gup(dim), Gup(dim);
 double detg = calculate_contravariant(g, gup);
 double detG = calculate_contravariant(G, Gup);
 
 // Find the number of additional strain invariants introduced by the PVA
 const unsigned n_ai = Strain_Energy_Function_pt->get_n_additional_strain_invariants();

 // Allocate storage for the strain invariants, we need the 3 isotropic ones plus however many the strain energy function
 // model creates.
 Vector<double> I(3 + n_ai, 0.0);

 // The third strain invaraint is the volumetric change
 I[2] = detG / detg;
 // The first and second are a bit more complex --- see G&Z
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   I[0] += gup(i, j) * G(i, j);
   I[1] += g(i, j) * Gup(i, j);
  }
 }

 // If 2D we assume plane strain: In this case the 3D tensors have
 // a 1 on the diagonal and zeroes in the off-diagonals of their
 // third rows and columns. Only effect: Increase the first two
 // invariants by one; rest of the computation can just be performed
 // over the 2d set of coordinates.
 if (dim == 2)
 {
  I[0] += 1.0;
  I[1] += 1.0;
 }

 // Second strain invariant is multiplied by the third.
 I[1] *= I[2];

 // Create the additional strain invariants which arise due to the PVA
 Vector<double> I_add(n_ai, 0.0);
 // Create the derivatives of the additional strain invariants with respect to G
 Vector<DenseMatrix<double>> dI_add_dG(n_ai, DenseMatrix<double>(dim));
 Strain_Energy_Function_pt->I(g, gup, G, Gup, detg, detG, a, I_add, dI_add_dG);
 for(unsigned i=0;i<n_ai; i++)
 {
  I[i+3] = I_add[i];
 }

 // Calculate the derivatives of the strain energy function wrt the
 // strain invariants
 Vector<double> dWdI(3 + n_ai, 0.0);
 Strain_Energy_Function_pt->derivatives(I, dWdI);

 // Only bother to compute the tensor B^{ij} (Green & Zerna notation)
 // if the derivative wrt the second strain invariant is non-zero
 DenseMatrix<double> Bup(dim, dim, 0.0);
 if (std::fabs(dWdI[1]) > 0.0)
 {
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < dim; j++)
   {
    Bup(i, j) = I[0] * gup(i, j);
    for (unsigned r = 0; r < dim; r++)
    {
     for (unsigned s = 0; s < dim; s++)
     {
      Bup(i, j) -= gup(i, r) * gup(j, s) * G(r, s);
     }
    }
   }
  }
 }

 // Now set the values of the functions phi, psi and p (Green & Zerna
 // notation) Note that the Green & Zerna stress \tau^{ij} is
 // s^{ij}/sqrt(I[2]), where s^{ij} is the desired second Piola-Kirchhoff
 // stress tensor so we multiply their constants by sqrt(I[2])
 double phi = 2.0 * dWdI[0];
 double psi = 2.0 * dWdI[1];
 double p = 2.0 * dWdI[2] * I[2];

 // Put it all together to get the stress
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   sigma(i, j) = phi * gup(i, j) + psi * Bup(i, j) + p * Gup(i, j);
  }
 }

 // Add the contributions from the additional strain invariants from the PVA
 for(unsigned k = 0; k < n_ai; k++)
 {
  // if (!(std::fabs(dWdI[3 + k]) > 0.0)) continue;
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < dim; j++)
   {
    // We multiply by 2.0 because d/dsigma(i,j) = 2.0d/dG(i,j)
    sigma(i, j) += dWdI[3 + k] * 2.0 * dI_add_dG[k](i,j);
   }
  }
 }
}

//===========================================================================
/// Calculate the deviatoric part
/// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
/// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
/// Also return the contravariant deformed metric
/// tensor and the determinant of the deformed metric tensor.
/// Uses correct 3D invariants for 2D (plane strain) problems.
/// This is the version for the pure incompressible formulation.
//============================================================================
void AnisotropicStrainEnergyFunctionConstitutiveLaw::
 calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                         const DenseMatrix<double>& G,
                                         const Vector<Vector<double>>& a,
                                         DenseMatrix<double>& sigma_dev,
                                         DenseMatrix<double>& Gup,
                                         double& detG)
{
// Error checking
#ifdef PARANOID
 error_checking_in_input(g, G, sigma);

 // Test that there are enough PVA
 if(!verify_number_of_pva(a, N_Principal_Vectors_Of_Anisotropy))
 {
  throw OomphLibError("Number of principal vectors of anisotropy is incorrect",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Find the dimension of the problem
 const unsigned dim = g.nrow();

#ifdef PARANOID
 if (dim == 1)
 {
  std::string function_name =
   "AnisotropicStrainEnergyFunctionConstitutiveLaw::";
  function_name += "calculate_second_piola_kirchhoff_stress()";

  throw OomphLibError("Check constitutive equations carefully when dim=1",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Calculate the contravariant undeformed and deformed metric tensors
 // and get the determinants of the metric tensors
 DenseMatrix<double> gup(dim);
 // Don't need this determinant
 (void)calculate_contravariant(g, gup);
 detG = calculate_contravariant(G, Gup);
 
 // Find the number of additional strain invariants introduced by the PVA
 const unsigned n_ai = Strain_Energy_Function_pt->get_n_additional_strain_invariants();

 // Allocate storage for the strain invariants, we need the 3 isotropic ones plus however many the strain energy function
 // model creates.
 Vector<double> I(3 + n_ai, 0.0);

 // The third strain invaraint is the volumetric change
 I[2] = 1.0;
 // The first and second are a bit more complex --- see G&Z
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   I[0] += gup(i, j) * G(i, j);
   I[1] += g(i, j) * Gup(i, j);
  }
 }

 // If 2D we assume plane strain: In this case the 3D tensors have
 // a 1 on the diagonal and zeroes in the off-diagonals of their
 // third rows and columns. Only effect: Increase the first two
 // invariants by one; rest of the computation can just be performed
 // over the 2d set of coordinates.
 if (dim == 2)
 {
  I[0] += 1.0;
  I[1] += 1.0;
 }

 // Create the additional strain invariants which arise due to the PVA
 Vector<double> I_add(n_ai, 0.0);
 // Create the derivatives of the additional strain invariants with respect to G
 Vector<DenseMatrix<double>> dI_add_dG(n_ai, DenseMatrix<double>(dim));
 Strain_Energy_Function_pt->I(g, gup, G, Gup, detG, detG, a, I_add, dI_add_dG);
 for(unsigned i=0;i<n_ai; i++)
 {
  I[i+3] = I_add[i];
 }

 // Calculate the derivatives of the strain energy function wrt the
 // strain invariants
 Vector<double> dWdI(3 + n_ai, 0.0);
 Strain_Energy_Function_pt->derivatives(I, dWdI);

 // Only bother to compute the tensor B^{ij} (Green & Zerna notation)
 // if the derivative wrt the second strain invariant is non-zero
 DenseMatrix<double> Bup(dim, dim, 0.0);
 if (std::fabs(dWdI[1]) > 0.0)
 {
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < dim; j++)
   {
    Bup(i, j) = I[0] * gup(i, j);
    for (unsigned r = 0; r < dim; r++)
    {
     for (unsigned s = 0; s < dim; s++)
     {
      Bup(i, j) -= gup(i, r) * gup(j, s) * G(r, s);
     }
    }
   }
  }
 }

 // Now set the values of the functions phi, psi and p (Green & Zerna
 // notation) Note that the Green & Zerna stress \tau^{ij} is
 // s^{ij}/sqrt(I[2]), where s^{ij} is the desired second Piola-Kirchhoff
 // stress tensor so we multiply their constants by sqrt(I[2])
 double phi = 2.0 * dWdI[0];
 double psi = 2.0 * dWdI[1];

 // Put it all together to get the stress
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   sigma_dev(i, j) = phi * gup(i, j) + psi * Bup(i, j);
  }
 }

 // Add the contributions from the additional strain invariants from the PVA
 for(unsigned k = 0; k < n_ai; k++)
 {
  // if (!(std::fabs(dWdI[3 + k]) > 0.0)) continue;
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < dim; j++)
   {
    // We multiply by 2.0 because d/dsigma(i,j) = 2.0d/dG(i,j)
    sigma_dev(i, j) += dWdI[3 + k] * 2.0 * dI_add_dG[k](i,j);
   }
  }
 }

 // Remove the hydrostatic stress - trace/dim
 double K = 0.0;
 for(unsigned i=0; i<dim; i++)
 {
  K += sigma_dev(i,i);
 }
 K /= dim;
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   sigma_dev(i, j) -= K * Gup(i, j);
  }
 }
}

//===========================================================================
/// Calculate the deviatoric part of the contravariant
/// 2nd Piola Kirchoff stress tensor. Also return the contravariant
/// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
/// the inverse of the bulk modulus \f$ \kappa\f$.
/// Uses correct 3D invariants for 2D (plane strain) problems.
/// This is the version for the near-incompressible formulation.
//===========================================================================
void AnisotropicStrainEnergyFunctionConstitutiveLaw::
  calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                          const DenseMatrix<double>& G,
                                          const Vector<Vector<double>>& a,
                                          DenseMatrix<double>& sigma_dev,
                                          DenseMatrix<double>& Gup,
                                          double& gen_dil,
                                          double& inv_kappa)
{
// Error checking
#ifdef PARANOID
 error_checking_in_input(g, G, sigma);

 // Test that there are enough PVA
 if(!verify_number_of_pva(a, N_Principal_Vectors_Of_Anisotropy))
 {
  throw OomphLibError("Number of principal vectors of anisotropy is incorrect",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Find the dimension of the problem
 const unsigned dim = g.nrow();

#ifdef PARANOID
 if (dim == 1)
 {
  std::string function_name =
   "AnisotropicStrainEnergyFunctionConstitutiveLaw::";
  function_name += "calculate_second_piola_kirchhoff_stress()";

  throw OomphLibError("Check constitutive equations carefully when dim=1",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }
#endif

 // Calculate the contravariant undeformed and deformed metric tensors
 // and get the determinants of the metric tensors
 DenseMatrix<double> gup(dim);
 double detg = calculate_contravariant(g, gup);
 double detG = calculate_contravariant(G, Gup);
 
 // Find the number of additional strain invariants introduced by the PVA
 const unsigned n_ai = Strain_Energy_Function_pt->get_n_additional_strain_invariants();

 // Allocate storage for the strain invariants, we need the 3 isotropic ones plus however many the strain energy function
 // model creates.
 Vector<double> I(3 + n_ai, 0.0);
 // The third strain invaraint is the volumetric change
 I[2] = detG / detg;
 // The first and second are a bit more complex --- see G&Z
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   I[0] += gup(i, j) * G(i, j);
   I[1] += g(i, j) * Gup(i, j);
  }
 }

 // If 2D we assume plane strain: In this case the 3D tensors have
 // a 1 on the diagonal and zeroes in the off-diagonals of their
 // third rows and columns. Only effect: Increase the first two
 // invariants by one; rest of the computation can just be performed
 // over the 2d set of coordinates.
 if (dim == 2)
 {
  I[0] += 1.0;
  I[1] += 1.0;
 }

 I[1] *= I[2];

 // Create the additional strain invariants which arise due to the PVA
 Vector<double> I_add(n_ai, 0.0);
 // Create the derivatives of the additional strain invariants with respect to G
 Vector<DenseMatrix<double>> dI_add_dG(n_ai, DenseMatrix<double>(dim));
 Strain_Energy_Function_pt->I(g, gup, G, Gup, detG, detG, a, I_add, dI_add_dG);
 for(unsigned i=0;i<n_ai; i++)
 {
  I[i+3] = I_add[i];
 }

 // Calculate the derivatives of the strain energy function wrt the
 // strain invariants
 Vector<double> dWdI(3 + n_ai, 0.0);
 Strain_Energy_Function_pt->derivatives(I, dWdI);

 // Only bother to compute the tensor B^{ij} (Green & Zerna notation)
 // if the derivative wrt the second strain invariant is non-zero
 DenseMatrix<double> Bup(dim, dim, 0.0);
 if (std::fabs(dWdI[1]) > 0.0)
 {
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < dim; j++)
   {
    Bup(i, j) = I[0] * gup(i, j);
    for (unsigned r = 0; r < dim; r++)
    {
     for (unsigned s = 0; s < dim; s++)
     {
      Bup(i, j) -= gup(i, r) * gup(j, s) * G(r, s);
     }
    }
   }
  }
 }

 // Now set the values of the functions phi and psi (Green & Zerna notation)
 // but multiplied by sqrt(I[2]) to recover the second Piola-Kirchhoff stress
 double phi = 2.0 * dWdI[0];
 double psi = 2.0 * dWdI[1];

 // Put it all together to get the stress
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   sigma_dev(i, j) = phi * gup(i, j) + psi * Bup(i, j);
  }
 }

 // Add the contributions from the additional strain invariants from the PVA
 for(unsigned k = 0; k < n_ai; k++)
 {
  // if (!(std::fabs(dWdI[3 + k]) > 0.0)) continue;
  for (unsigned i = 0; i < dim; i++)
  {
   for (unsigned j = 0; j < dim; j++)
   {
    // We multiply by 2.0 because d/dsigma(i,j) = 2.0d/dG(i,j)
    sigma_dev(i, j) += dWdI[3 + k] * 2.0 * dI_add_dG[k](i,j);
   }
  }
 }

 // Calculate the hydrostatic stress - trace/dim
 double K = 0.0;
 for(unsigned i=0; i<dim; i++)
 {
  K += sigma_dev(i,i);
 }
 K /= dim;

 // Choose inverse kappa to be one...
 inv_kappa = 1.0;

 //...then the generalised dilation is the same as p  in Green & Zerna's
 // notation, but multiplied by sqrt(I[2]) with the addition of the
 // terms that are subtracted to make the other part of the stress
 // deviatoric
 gen_dil = 2.0 * dWdI[2] * I[2] + K;

 // Calculate the deviatoric part of the stress by subtracting
 // the computed trace/dim
 for (unsigned i = 0; i < dim; i++)
 {
  for (unsigned j = 0; j < dim; j++)
  {
   sigma_dev(i, j) -= K * Gup(i, j);
  }
 }
}




} // End namespace