#include "definitions.h"

CustomWeakFormPoisson::CustomWeakFormPoisson(const std::string& mat_motor, double eps_motor, 
                                             const std::string& mat_air, double eps_air, bool is_matfree) : WeakForm<double>(1)
{
  this->is_matfree = is_matfree;

  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, new Hermes1DFunction<double>(eps_motor), mat_motor));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, new Hermes1DFunction<double>(eps_air), mat_air));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, new Hermes1DFunction<double>(eps_motor), mat_motor));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, new Hermes1DFunction<double>(eps_air), mat_air));
}

