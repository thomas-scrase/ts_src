#include "refineable_anisotropic_solid_elements.h"

namespace oomph
{
// Implement the refineable anisotropic fill in procedure
template<unsigned DIM>
void RefineableAnisotropicPVDEquations<DIM>::fill_in_generic_contribution_to_residuals_pvd(Vector<double>& residuals,
                                                                                           DenseMatrix<double>& jacobian,
                                                                                           const unsigned& flag)
{

}

template class RefineableAnisotropicPVDEquations<2>;
template class RefineableAnisotropicPVDEquations<3>;
} // End namespace