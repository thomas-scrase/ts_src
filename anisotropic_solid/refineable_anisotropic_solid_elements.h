#ifndef REFINEABLE_ANISOTROPIC_SOLID_ELEMENTS_H
#define REFINEABLE_ANISOTROPIC_SOLID_ELEMENTS_H

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "anisotropic_solid_elements.h"
#include "../solid/refineable_solid_elements.h"


// Implement

namespace oomph
{
template<unsigned DIM>
class RefineableAnisotropicPVDEquations : public virtual AnisotropicPVDEquations<DIM>,
                                          public virtual RefineablePVDEquations<DIM>
{
public:
 RefineableAnisotropicPVDEquations() : 
  PVDEquations<DIM>(),
  RefineableElement(),
  RefineableSolidElement(),
  ElementWithZ2ErrorEstimator(),
  RefineablePVDEquations<DIM>(),
  AnisotropicPVDEquationsBase<DIM>(),
  AnisotropicPVDEquations<DIM>()
 {
 }

 void further_build()
 {
  // Call the refineable PVD further build
  RefineablePVDEquations<DIM>::further_build();

  // Complete the further build by getting the principal vectors function pointer from father element
  AnisotropicPVDEquations<DIM>* cast_father_element_pt = dynamic_cast<AnisotropicPVDEquations<DIM>*>(this->father_element_pt());
  this->Principal_vectors_of_anisotropy_fct_pt = cast_father_element_pt->principal_vectors_of_anisotropy_fct_pt();
  this->Anisotropic_constitutive_law_pt = cast_father_element_pt->anisotropic_constitutive_law_pt();
 }

protected:
 // Override the fill in to contain the necessary steps of including the anisotropic components
 virtual void fill_in_generic_contribution_to_residuals_pvd(Vector<double>& residuals,
                                                            DenseMatrix<double>& jacobian,
                                                            const unsigned& flag) override;
};


// Implement RefineableAnisotropicQPVDElement
template<unsigned DIM, unsigned NNODE_1D>
class RefineableAnisotropicQPVDElement : public virtual RefineableQPVDElement<DIM, NNODE_1D>,
                                         public virtual RefineableAnisotropicPVDEquations<DIM>,
                                         public virtual RefineableSolidQElement<DIM>
{
public:
 /// Constructor:
 RefineableAnisotropicQPVDElement() :
  SolidQElement<DIM, NNODE_1D>(),
  PVDEquations<DIM>(),
  QPVDElement<DIM, NNODE_1D>(),
  RefineablePVDEquations<DIM>(),
  RefineableSolidQElement<DIM>(),
  RefineableElement(),
  RefineableSolidElement(),
  ElementWithZ2ErrorEstimator(),
  AnisotropicPVDEquationsBase<DIM>(),
  AnisotropicPVDEquations<DIM>(),
  RefineableAnisotropicPVDEquations<DIM>()
 {
 }

protected:
 // Override the fill in to contain the necessary steps of including the anisotropic components
 virtual void fill_in_generic_contribution_to_residuals_pvd(Vector<double>& residuals,
                                                            DenseMatrix<double>& jacobian,
                                                            const unsigned& flag) override
 {
  RefineableAnisotropicPVDEquations<DIM>::fill_in_generic_contribution_to_residuals_pvd(residuals, jacobian, flag);
 }
};

//==============================================================
/// FaceGeometry of the 2D RefineableAnisotropicQPVDElement elements
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<RefineableAnisotropicQPVDElement<2, NNODE_1D>>
  : public virtual SolidQElement<1, NNODE_1D>
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
};

//==============================================================
/// FaceGeometry of the FaceGeometry of the 2D RefineableAnisotropicQPVDElement
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<FaceGeometry<RefineableAnisotropicQPVDElement<2, NNODE_1D>>>
  : public virtual PointElement
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : PointElement() {}
};


//==============================================================
/// FaceGeometry of the 3D RefineableAnisotropicQPVDElement elements
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<RefineableAnisotropicQPVDElement<3, NNODE_1D>>
  : public virtual SolidQElement<2, NNODE_1D>
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<2, NNODE_1D>() {}
};

//==============================================================
/// FaceGeometry of the FaceGeometry of the 3D RefineableAnisotropicQPVDElement
//==============================================================
template<unsigned NNODE_1D>
class FaceGeometry<FaceGeometry<RefineableAnisotropicQPVDElement<3, NNODE_1D>>>
  : public virtual SolidQElement<1, NNODE_1D>
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
};


// Implement RefineableAnisotropicPVDEquationsWithPressure and face elements
template<unsigned DIM>
class RefineableAnisotropicPVDEquationsWithPressure : public virtual RefineablePVDEquationsWithPressure<DIM>,
                                                      public virtual AnisotropicPVDEquationsWithPressure<DIM>
{
 RefineableAnisotropicPVDEquationsWithPressure() :
  PVDEquationsBase<DIM>(),
  PVDEquationsWithPressure<DIM>(),
  RefineableElement(),
  RefineableSolidElement(),
  ElementWithZ2ErrorEstimator(),
  AnisotropicPVDEquationsBase<DIM>(),
  AnisotropicPVDEquationsWithPressure<DIM>(),
  RefineablePVDEquations<DIM>()
 {
 }

 void further_build()
 {
  // Call the refineable PVD further build
  RefineablePVDEquations<DIM>::further_build();

  // Complete the further build by getting the principal vectors function pointer from father element
  AnisotropicPVDEquations<DIM>* cast_father_element_pt = dynamic_cast<AnisotropicPVDEquations<DIM>*>(this->father_element_pt());
  this->Principal_vectors_of_anisotropy_fct_pt = cast_father_element_pt->principal_vectors_of_anisotropy_fct_pt();
  this->Anisotropic_constitutive_law_pt = cast_father_element_pt->anisotropic_constitutive_law_pt();
 }

protected:
 void fill_in_generic_residual_contribution_pvd_with_pressure(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix,
  const unsigned& flag) override;
};

// Implement RefineableAnisotropicQPVDElementWithPressure and face elements
template<unsigned DIM>
class RefineableAnisotropicQPVDElementWithPressure : public virtual RefineableQPVDElementWithPressure<DIM>,
                                                     public virtual RefineableAnisotropicPVDEquationsWithPressure<DIM>
{
public:
 RefineableAnisotropicQPVDElementWithPressure() : 
  QPVDElementWithPressure<DIM>(),
  RefineableElement(),
  RefineableSolidElement(),
  RefineablePVDEquationsWithPressure<DIM>(),
  RefineableSolidQElement<DIM>(),
  PVDEquationsWithPressure<DIM>(),
  AnisotropicPVDEquationsBase<DIM>(),
  AnisotropicPVDEquationsWithPressure<DIM>(),
  RefineableAnisotropicPVDEquationsWithPressure<DIM>()
 {
 }


protected:
 void fill_in_generic_residual_contribution_pvd_with_pressure(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix,
  const unsigned& flag) override
 {
  RefineableAnisotropicPVDEquationsWithPressure<DIM>::fill_in_generic_residual_contribution_pvd_with_pressure(residuals, jacobian, mass_matrix, flag);
 }
};


//======================================================================
/// FaceGeometry of the 2D RefineableAnisotropicQPVDElementWithPressure
//=======================================================================
template<>
class FaceGeometry<RefineableAnisotropicQPVDElementWithPressure<2>>
 : public virtual SolidQElement<1, 3>
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<1, 3>() {}
};


//===========================================================================
/// FaceGeometry of the FaceGeometry of the 2D
/// RefineableAnisotropicQPVDElementWithPressure
//============================================================================
template<>
class FaceGeometry<FaceGeometry<RefineableAnisotropicQPVDElementWithPressure<2>>>
 : public virtual PointElement
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : PointElement() {}
};

//=========================================================================
/// FaceGeometry of the 3D RefineableAnisotropicQPVDElementWithPressure
//========================================================================
template<>
class FaceGeometry<RefineableAnisotropicQPVDElementWithPressure<3>>
 : public virtual SolidQElement<2, 3>
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<2, 3>() {}
};


//========================================================================
/// FaceGeometry of the FaceGeometry of the 3D
/// RefineableAnisotropicQPVDElementWithPressure
//==========================================================================
template<>
class FaceGeometry<FaceGeometry<RefineableAnisotropicQPVDElementWithPressure<3>>>
 : public virtual SolidQElement<1, 3>
{
public:
 // Make sure that we call the constructor of the SolidQElement
 // Only the Intel compiler seems to need this!
 FaceGeometry() : SolidQElement<1, 3>() {}
};



// Implement RefineableAnisotropicQPVDElementWithContinuousPressure and face elements

} // End namespace

#endif 