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

// File generated at Tue 8 Sep 2015 13:06:08

#include "NUTSMSSM_two_scale_susy_scale_constraint.hpp"
#include "NUTSMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"

#include <cassert>
#include <cmath>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole Electroweak_constants::MZ
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME NUTSMSSM<Two_scale>

NUTSMSSM_susy_scale_constraint<Two_scale>::NUTSMSSM_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

NUTSMSSM_susy_scale_constraint<Two_scale>::NUTSMSSM_susy_scale_constraint(
   NUTSMSSM<Two_scale>* model_)
   : Constraint<Two_scale>()
   , model(model_)
{
   initialize();
}

NUTSMSSM_susy_scale_constraint<Two_scale>::~NUTSMSSM_susy_scale_constraint()
{
}

void NUTSMSSM_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: NUTSMSSM_susy_scale_constraint::apply():"
          " model pointer must not be zero");

   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto LambdaSInput = INPUTPARAMETER(LambdaSInput);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   MODEL->set_vS(Re(LambdaSInput/Lambdax));
   MODEL->solve_ewsb();

}

double NUTSMSSM_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double NUTSMSSM_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const NUTSMSSM_input_parameters& NUTSMSSM_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   assert(model && "Error: NUTSMSSM_susy_scale_constraint::"
          "get_input_parameters(): model pointer is zero.");

   return model->get_input();
}

NUTSMSSM<Two_scale>* NUTSMSSM_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void NUTSMSSM_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<NUTSMSSM<Two_scale>*>(model_);
}

void NUTSMSSM_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void NUTSMSSM_susy_scale_constraint<Two_scale>::initialize()
{
   assert(model && "NUTSMSSM_susy_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   const auto m0 = INPUTPARAMETER(m0);
   const auto m12 = INPUTPARAMETER(m12);

   initial_scale_guess = Sqrt(Sqr(m0) + 4*Sqr(m12));

   scale = initial_scale_guess;
}

void NUTSMSSM_susy_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "NUTSMSSM_susy_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const auto ZU = MODELPARAMETER(ZU);
   const auto MSu = MODELPARAMETER(MSu);

   scale = Sqrt(Power(MSu(0),Sqr(Abs(ZU(0,2))) + Sqr(Abs(ZU(0,5))))*Power(MSu(1
      ),Sqr(Abs(ZU(1,2))) + Sqr(Abs(ZU(1,5))))*Power(MSu(2),Sqr(Abs(ZU(2,2))) +
      Sqr(Abs(ZU(2,5))))*Power(MSu(3),Sqr(Abs(ZU(3,2))) + Sqr(Abs(ZU(3,5))))*Power
      (MSu(4),Sqr(Abs(ZU(4,2))) + Sqr(Abs(ZU(4,5))))*Power(MSu(5),Sqr(Abs(ZU(5,2))
      ) + Sqr(Abs(ZU(5,5)))));


}

} // namespace flexiblesusy
