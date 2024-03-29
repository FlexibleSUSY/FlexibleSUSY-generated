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


#include "SplitMSSM_two_scale_susy_scale_constraint.hpp"
#include "SplitMSSM_two_scale_model.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
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
#define MODELCLASSNAME SplitMSSM<Two_scale>

SplitMSSM_susy_scale_constraint<Two_scale>::SplitMSSM_susy_scale_constraint(
   SplitMSSM<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void SplitMSSM_susy_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   
   MODEL->solve_ewsb();

   // calculate SM-like Higgs pole mass
   // for usage in MW calculation at low-energy scale
   {
      auto tmp = *MODEL;
      tmp.do_force_output(true); // enforce calculation of pole masses
      tmp.solve_ewsb();
      tmp.calculate_Mhh_pole();
      MODEL->get_physical().Mhh = tmp.get_physical().Mhh;
   }

}

double SplitMSSM_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double SplitMSSM_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const SplitMSSM_input_parameters& SplitMSSM_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   check_model_ptr();

   return model->get_input();
}

SplitMSSM<Two_scale>* SplitMSSM_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void SplitMSSM_susy_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<SplitMSSM<Two_scale>*>(model_);
}

void SplitMSSM_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& SplitMSSM_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void SplitMSSM_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
}

void SplitMSSM_susy_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto MEWSB = INPUTPARAMETER(MEWSB);

   initial_scale_guess = MEWSB;

   scale = initial_scale_guess;
}

void SplitMSSM_susy_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto MEWSB = INPUTPARAMETER(MEWSB);

   scale = MEWSB;


}

void SplitMSSM_susy_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("SplitMSSM_susy_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
