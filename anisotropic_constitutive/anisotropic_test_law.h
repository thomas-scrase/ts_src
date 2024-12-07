#ifndef ANISOTROPIC_TEST_LAW_HEADER
#define ANISOTROPIC_TEST_LAW_HEADER

#include "anisotropic_constitutive.h"

namespace oomph
{
// A nonsense strain energy function to illustrate how to use anisotropic strain energy functions
class AnisotropicTestStrainEnergy : AnisotropicStrainEnergyFunction
{
public:
 AnisotropicTestStrainEnergy(double* c1_pt, double* c1_pt, double* c3_pt) : N_Principal_Vectors_Of_Anisotropy(1),
                                                                            N_Additional_Strain_Invariants(1), // One additional strain invariant arises from the fibre orientation
                                                                            C1_pt(c1_pt),
                                                                            C2_pt(c2_pt),
                                                                            C3_pt(c3_pt)
 {
  
 }

 // Tell the constitutive law how to compute the additional strain invariants and the derivatives of the strain invariants with respect to the deformed metric tensor
 void I(const DenseMatrix<double>& g,
        const DenseMatrix<double>& g_up,
        const DenseMatrix<double>& G,
        const DenseMatrix<double>& G_up,
        const double& detg,
        const double& getG,
        const Vector<Vector<double>>& a,
        Vector<double>& I,
        Vector<DenseMatrix<double>>& dIdG)
 {
  const unsigned dim = g.ncol();
  // The only additional strain invariant this model uses is a_iG_{ij}a_j
  I[0] = 0.0;
  for(unsigned i=0; i<dim; i++)
  {
   for(unsigned j=0; j<dim; j++)
   {
    I[0] += a[0][i] * G(i,j) * a[0][j];

    dIdG[0](i,j) = a[0][i] * a[0][j];
   }
  }
 }

 // Compute the strain energy
 double W(const Vector<double>& I)
 {
  return (*C1_pt)*(I[0] - 3.0) + (*C2_pt) * (I[1] - 3.0)
         + (*C3_pt) * I[3];
 }

 /// Return the derivatives of the strain energy function with
 /// respect to the strain invariants
 void derivatives(Vector<double>& I, Vector<double>& dWdI)
 {
   dWdI[0] = (*C1_pt);
   dWdI[1] = (*C2_pt);
   dWdI[2] = 0.0;
   dWdI[3] = (*C3_pt);
 }

 bool requires_incompressibility_constraint()
 {
  return true;
 }

private:
 double* C1_pt;
 double* C2_pt;
 double* C3_pt;
};


// A nonsense strain energy function to illustrate how to use anisotropic strain energy functions with active stress
class AnisotropicActiveStressTestStrainEnergy : AnisotropicStrainEnergyFunction
{
public:
 AnisotropicActiveStressTestStrainEnergy(double* c1_pt, double* c1_pt, double* c3_pt) : N_Principal_Vectors_Of_Anisotropy(2),
                                                                                        N_Additional_Strain_Invariants(2), // A second additional strain invariant is included as a dummy invariant to store the active stress along the fibres
                                                                                        C1_pt(c1_pt),
                                                                                        C2_pt(c2_pt),
                                                                                        C3_pt(c3_pt)
 {
 }

 // Tell the constitutive law how to compute the additional strain invariants
 void I(const DenseMatrix<double>& g,
        const DenseMatrix<double>& g_up,
        const DenseMatrix<double>& G,
        const DenseMatrix<double>& G_up,
        const double& detg,
        const double& getG,
        const Vector<Vector<double>>& a,
        Vector<double>& I,
        Vector<DenseMatrix<double>>& dIdG)
 {
  const unsigned dim = g.ncol();
  // The only additional strain invariant this model uses is a_iG_{ij}a_j
  I[0] = 0.0;
  // We use a dummy strain invariant to encode the active stress along the fibres - we get the stress along the fibres from the first element in the dummy additional principal vector of anisotropy
  I[1] = a[1][0];
  for(unsigned i=0; i<dim; i++)
  {
   for(unsigned j=0; j<dim; j++)
   {
    I[0] += a[0][i] * G(i,j) * a[0][j];
    dIdG[0](i,j) = a[0][i] * a[0][j];

    // The second "strain invariant" is actually just used to store the active stress along the fibre direction
    // so it has zero derivatives with respect to the deformed metric tensor
    dIdG[1](i,j) = 0.0;
   }
  }
 }

 // Compute the strain energy - adds an additional term which represents the active stress in the direction of the fibres
 double W(const Vector<double>& I)
 {
  return (*C1_pt)*(I[0] - 3.0) + (*C2_pt) * (I[1] - 3.0)
         + (*C3_pt + I[4]) * I[3];
 }

 /// Return the derivatives of the strain energy function with
 /// respect to the strain invariants
 /// - adds an additional term which represents the active stress in the direction of the fibres
 void derivatives(Vector<double>& I, Vector<double>& dWdI)
 {
   dWdI[0] = (*C1_pt);
   dWdI[1] = (*C2_pt);
   dWdI[2] = 0.0;
   dWdI[3] = (*C3_pt) + I[4];
   dWdI[4] = 0.0; // This is a dummy strain invariant so it should not be included in the fill-in of the 2nd Piola-Kirchhoff stress tensor
 }

 bool requires_incompressibility_constraint()
 {
  return true;
 }

private:
 double* C1_pt;
 double* C2_pt;
 double* C3_pt;
};

};

#endif