// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Wed 16 Oct 2019 22:56:04

#include "CMSSMSemiAnalytic_semi_analytic_susy_scale_constraint.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"

#include <cmath>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
#define EXTRAPARAMETER(p) model->get_##p()
#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETAPARAMETER1(l,p) beta_functions_##l##L.get_##p()
#define BETA(p) beta_##p
#define BETA1(l,p) beta_##l##L_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME CMSSMSemiAnalytic<Semi_analytic>

CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::CMSSMSemiAnalytic_susy_scale_constraint(
   CMSSMSemiAnalytic<Semi_analytic>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

double CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::get_scale() const
{
   if (!scale_getter) {
      return scale;
   }
   return scale_getter();

}

void CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto MuInput = INPUTPARAMETER(MuInput);

   MODEL->set_Mu(Re(MuInput));



}

   const CMSSMSemiAnalytic_input_parameters& CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::get_input_parameters() const
{
   check_model_ptr();

   return model->get_input();
}

void CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::set_model(Model* model_)
{
   model = cast_model<CMSSMSemiAnalytic<Semi_analytic>*>(model_);
}

void CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
}

void CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::initialize()
{
   check_model_ptr();

   const auto Azero = INPUTPARAMETER(Azero);
   const auto m12 = INPUTPARAMETER(m12);

   initial_scale_guess = Sqrt(-3*Azero*m12 + Sqr(Azero) + 14*Sqr(m12));

   scale = initial_scale_guess;
}

void CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::update_scale()
{
   check_model_ptr();



}

void CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>::check_model_ptr() const
{
   if (!model)
      throw SetupError("CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
