
#ifndef MONODOMAIN_ELEMENTS_HEADER
#define MONODOMAIN_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
  #include <oomph-lib-config.h>
#endif

#include "generic.h"

namespace oomph
{

//=============================================================
/// A class for all elements that solve the Monodomain
/// equations.
/// \f[ \omega\left(\gamma C_m\frac{\partial V}{\partial t} + I\right) = \frac{\partial}{\partial x_i}\left(D_{ij}\frac{\partial V}{\partial x_k}\right) \f]
/// \omega = \frac{[I][X]^2}{[D][U]} analogous to the peclet number in advection diffusion problems
/// \gamma = \frac{[C][U]}{[T][I]} analogous to the strouhal number in advection diffusion problems
/// where:
///  C_m is the non-dimensionalised membrane capacitance per unit volume
///  V is the non-dimensionalised transmembrane potential
///  t is the non-dimensionalised time
///  I is the non-dimensionalised transmembrane current per unit volume
///  x_i is the non-dimensionalised i-th spatial coordinate
///  D_{ij} is the non-dimensionalised conductivity tensor

///  [C] is the characteristic capacitance per unit volume of the problem
///  [X] is the characteristic spatial length-scale of the problem
///  [D] is the characteristic conductance of the problem
///  [T] is the characteristic time-scale of the problem
///  [U] is the characteristic transmembrane potential of the problem
///  [I] is the characteristic transmembrane current per unit volume of the problem

// Based heavily on the GeneralisedAdvectionDiffusionEquations etc classes

/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template<unsigned DIM>
class MonodomainEquations : FiniteElement
{
public:
 /// Function pointer to source function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*MonodomainSourceFctPt)(const Vector<double>& x,
                                       double& f);

 /// Funciton pointer to a diffusivity function
 typedef void (*MonodomainDiffFctPt)(const Vector<double>& x,
                                     DenseMatrix<double>& D);

 /// Function pointer to capacitance function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*MonodomainCapacitanceFctPt)(const Vector<double>& x,
                                            double& C);

 /// Constructor: Initialise the Source_fct_pt, Diff_fct_pt, and Capacitance_fct_pt to null
 /// and set the (pointer to) Omega number and Gamma number to default;
 MonodomainEquations() : Source_fct_pt(0),
                         Diff_fct_pt(0),
                         Capacitance_fct_pt(0),
                         ALE_is_disabled(false)
 {
  Omega_pt = &Default_omega_number;
  Gamma_pt = &Default_gamma_number;
 }

 /// Broken copy constructor
 MonodomainEquations(const MonodomainEquations& dummy) = delete;

 /// Broken assignment operator
 void operator=(const MonodomainEquations&) = delete;

 /// Return the index at which the unknown value
 /// is stored. The default value, 0, is appropriate for single-physics
 /// problems, when there is only one variable, the value that satisfies
 /// the monodomain equation.
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
 virtual inline unsigned u_index_monodomain() const
 {
  return 0;
 }

 /// du/dt at local node n.
 /// Uses suitably interpolated value for hanging nodes.
 /// Virtual to allow for overriding in operator splitting methods.
 /// In operator splitting we want to be able to grab other values
 /// as a substitute for the values at the previous time-step.
 virtual double du_dt_monodomain(const unsigned& n) const
 {
  // Get the data's timestepper
  TimeStepper* time_stepper_pt = this->node_pt(n)->time_stepper_pt();

  // Initialise dudt
  double dudt = 0.0;
  // Loop over the timesteps, if there is a non Steady timestepper
  if (!time_stepper_pt->is_steady())
  {
   // Find the index at which the variable is stored
   const unsigned u_nodal_index = u_index_monodomain();

   // Number of timsteps (past & present)
   const unsigned n_time = time_stepper_pt->ntstorage();

   for (unsigned t = 0; t < n_time; t++)
   {
    dudt += time_stepper_pt->weight(1, t) * nodal_value(t, n, u_nodal_index);
   }
  }
  return dudt;
 }

 /// Disable ALE, i.e. assert the mesh is not moving -- you do this
 /// at your own risk!
 void disable_ALE()
 {
  ALE_is_disabled = true;
 }

 /// (Re-)enable ALE, i.e. take possible mesh motion into account
 /// when evaluating the time-derivative. Note: By default, ALE is
 /// enabled, at the expense of possibly creating unnecessary work
 /// in problems where the mesh is, in fact, stationary.
 void enable_ALE()
 {
  ALE_is_disabled = false;
 }

 /// Output with default number of plot points
 void output(std::ostream& outfile)
 {
  unsigned nplot = 5;
  output(outfile, nplot);
 }

 /// Output FE representation of soln: x,y,u or x,y,z,u at
 /// nplot^DIM plot points
 void output(std::ostream& outfile, const unsigned& nplot);


 /// C_style output with default number of plot points
 void output(FILE* file_pt)
 {
  unsigned n_plot = 5;
  output(file_pt, n_plot);
 }

 /// C-style output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned& n_plot);


 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at nplot^DIM plot points
 void output_fct(std::ostream& outfile,
                 const unsigned& nplot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at
 /// nplot^DIM plot points (dummy time-dependent version to
 /// keep intel compiler happy)
 virtual void output_fct(
   std::ostream& outfile,
   const unsigned& nplot,
   const double& time,
   FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
 {
  throw OomphLibError("There is no time-dependent output_fct() for "
                      "Monodomain elements",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
 }

 /// Get error against and norm of exact solution
 void compute_error(std::ostream& outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error,
                    double& norm);


 /// Dummy, time dependent error checker
 void compute_error(std::ostream& outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time,
                    double& error,
                    double& norm)
 {
  throw OomphLibError(
    "No time-dependent compute_error() for Advection Diffusion elements",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
 }

 /// Integrate the concentration over the element
 double integrate_u();

 /// Access function: Pointer to source function
 MonodomainSourceFctPt& source_fct_pt()
 {
  return Source_fct_pt;
 }

 /// Access function: Pointer to source function. Const version
 MonodomainSourceFctPt source_fct_pt() const
 {
  return Source_fct_pt;
 }

 /// Access function: Pointer to diffusion  function
 MonodomainDiffFctPt& diff_fct_pt()
 {
  return Diff_fct_pt;
 }

 /// Access function: Pointer to diffusion function. Const version
 MonodomainDiffFctPt diff_fct_pt() const
 {
  return Diff_fct_pt;
 }

 /// Access function: Pointer to diffusion  function
 MonodomainCapacitanceFctPt& capacitance_fct_pt()
 {
  return Capacitance_fct_pt;
 }

 /// Access function: Pointer to diffusion function. Const version
 MonodomainCapacitanceFctPt capacitance_fct_pt() const
 {
  return Capacitance_fct_pt;
 }

 const double& omega() const
 {
  return *Omega_pt;
 }

 double*& omega_pt()
 {
  return Omega_pt;
 }

 const double& gamma() const
 {
  return *Gamma_pt;
 }

 double*& gamma_pt()
 {
  return Gamma_pt;
 }

 /// Get source term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the source function might be determined by
 /// another system of equations
 inline virtual void get_source_monodomain(const unsigned& ipt,
                                           const Vector<double>& x,
                                           double& source) const
 {
  // If no source function has been set, return zero
  if (Source_fct_pt == 0)
  {
   source = 0.0;
  }
  else
  {
   // Get source strength
   (*Source_fct_pt)(x, source);
  }
 }

 /// Get diffusivity tensor at (Eulerian) position
 /// x and/or local coordinate s.
 /// This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the wind function might be determined by
 /// another system of equations
 inline virtual void get_diff_monodomain(const unsigned& ipt,
                                         const Vector<double>& s,
                                         const Vector<double>& x,
                                         DenseMatrix<double>& D) const
 {
  // If no wind function has been set, return identity
  if (Diff_fct_pt == 0)
  {
   for (unsigned i = 0; i < DIM; i++)
   {
    for (unsigned j = 0; j < DIM; j++)
    {
     if (i == j)
     {
      D(i, j) = 1.0;
     }
     else
     {
      D(i, j) = 0.0;
     }
    }
   }
  }
  else
  {
   // Get diffusivity tensor
   (*Diff_fct_pt)(x, D);
  }
 }

 /// Get capacitance term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the capacitance function might be determined by
 /// another system of equations
 inline virtual void get_capacitance_monodomain(const unsigned& ipt,
                                                const Vector<double>& x,
                                                double& capacitance) const
 {
  // If no source function has been set, return zero
  if (Capacitance_fct_pt == 0)
  {
   capacitance = 0.0;
  }
  else
  {
   // Get source strength
   (*Capacitance_fct_pt)(x, capacitance);
  }
 }

 /// Get flux: \f$\mbox{flux}[i] = \mbox{d}u / \mbox{d}x_i \f$
 void get_flux(const Vector<double>& s, Vector<double>& flux) const
 {
  // Find out how many nodes there are in the element
  unsigned n_node = nnode();

  // Get the nodal index at which the unknown is stored
  unsigned u_nodal_index = u_index_monodomain();

  // Set up memory for the shape and test functions
  Shape psi(n_node);
  DShape dpsidx(n_node, DIM);

  // Call the derivatives of the shape and test functions
  dshape_eulerian(s, psi, dpsidx);

  // Initialise to zero
  for (unsigned j = 0; j < DIM; j++)
  {
   flux[j] = 0.0;
  }

  // Loop over nodes
  for (unsigned l = 0; l < n_node; l++)
  {
   // Loop over derivative directions
   for (unsigned j = 0; j < DIM; j++)
   {
    flux[j] += nodal_value(l, u_nodal_index) * dpsidx(l, j);
   }
  }
 }

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double>& residuals)
 {
  // Call the generic residuals function with flag set to 0 and using
  // a dummy matrix
  fill_in_generic_residual_contribution_monodomain(residuals,
                                                   GeneralisedElement::Dummy_matrix,
                                                   GeneralisedElement::Dummy_matrix,
                                                   0);
 }

 /// Add the element's contribution to its residual vector and
 /// the element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                       DenseMatrix<double>& jacobian)
 {
  // Call the generic routine with the flag set to 1
  fill_in_generic_residual_contribution_monodomain(residuals,
                                                      jacobian,
                                                      GeneralisedElement::Dummy_matrix,
                                                      1);
 }

 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double>& residuals,
                                                       DenseMatrix<double>& jacobian,
                                                       DenseMatrix<double>& mass_matrix)
 {
  // Call the generic routine with the flag set to 2
  fill_in_generic_residual_contribution_monodomain(residuals,
                                                      jacobian,
                                                      mass_matrix,
                                                      2);
 }
 

 /// Return FE representation of function value u(s) at local coordinate s
 inline double interpolated_u_monodomain(const Vector<double>& s) const
 {
  // Find number of nodes
  unsigned n_node = nnode();

  // Get the nodal index at which the unknown is stored
  unsigned u_nodal_index = u_index_monodomain();

  // Local shape function
  Shape psi(n_node);

  // Find values of shape function
  shape(s, psi);

  // Initialise value of u
  double interpolated_u = 0.0;

  // Loop over the local nodes and sum
  for (unsigned l = 0; l < n_node; l++)
  {
   interpolated_u += nodal_value(l, u_nodal_index) * psi[l];
  }

  return (interpolated_u);
 }

 /// Self-test: Return 0 for OK
 unsigned self_test();


protected:
 /// Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_monodomain(const Vector<double>& s,
                                                     Shape& psi,
                                                     DShape& dpsidx,
                                                     Shape& test,
                                                     DShape& dtestdx) const = 0;

 /// Shape/test functions and derivs w.r.t. to global coords at
 /// integration point ipt; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_at_knot_monodomain(const unsigned& ipt,
                                                             Shape& psi,
                                                             DShape& dpsidx,
                                                             Shape& test,
                                                             DShape& dtestdx) const = 0;

 /// Add the element's contribution to its residual vector only
 /// (if flag=and/or element  Jacobian matrix
 virtual void fill_in_generic_residual_contribution_monodomain(Vector<double>& residuals,
                                                               DenseMatrix<double>& jacobian,
                                                               DenseMatrix<double>& mass_matrix,
                                                               unsigned flag);

 /// Pointer to global Omega (Peclet) number
 double* Omega_pt;

 /// Pointer to global Gamma (Strouhal) number
 double* Gamma_pt;

 /// Pointer to source function:
 MonodomainSourceFctPt Source_fct_pt;

 /// Pointer to diffusivity funciton
 MonodomainDiffFctPt Diff_fct_pt;

 /// Pointer to capacitance function
 MonodomainCapacitanceFctPt Capacitance_fct_pt

 /// Boolean flag to indicate if ALE formulation is disabled when
 /// time-derivatives are computed. Only set to false if you're sure
 /// that the mesh is stationary.
 bool ALE_is_disabled;

private:
 /// Static default value for the Omega (Peclet) number
 static double Default_omega_number;

 /// Static default value for the Gamma (Strouhal) number
 static double Default_gamma_number;
};



// Q monodomain element class
template<unsigned DIM, unsigned NNODE>
class QMonodomainElement : public virtual MonodomainEquations<DIM>,
                           public virtual QElement<DIM, NNODE>
{

};


// T monodomain element class
template<unsigned DIM, unsigned NNODE>
class TMonodomainElement : public virtual MonodomainEquations<DIM>,
                           public virtual TElement<DIM, NNODE>
{

};


} // End namespace
#endif