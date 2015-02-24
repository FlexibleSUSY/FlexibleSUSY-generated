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

// File generated at Tue 24 Feb 2015 17:52:39

#include "MSSM_two_scale_susy_scale_constraint.hpp"
#include "MSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"

#include <cassert>
#include <cmath>

namespace flexiblesusy {

#define INPUTPARAMETER(p) inputPars.p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define SM(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME MSSM<Two_scale>

MSSM_susy_scale_constraint<Two_scale>::MSSM_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , inputPars()
{
}

MSSM_susy_scale_constraint<Two_scale>::MSSM_susy_scale_constraint(
   MSSM<Two_scale>* model_,
   const MSSM_input_parameters& inputPars_)
   : Constraint<Two_scale>()
   , model(model_)
   , inputPars(inputPars_)
{
   initialize();
}

MSSM_susy_scale_constraint<Two_scale>::~MSSM_susy_scale_constraint()
{
}

void MSSM_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: MSSM_susy_scale_constraint::apply():"
          " model pointer must not be zero");

   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints


   // the parameters, which are fixed by the EWSB eqs., will now be
   // defined at this scale (at the EWSB loop level defined in the
   // model)
   model->solve_ewsb();
}

double MSSM_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MSSM_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const MSSM_input_parameters& MSSM_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   return inputPars;
}

MSSM<Two_scale>* MSSM_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void MSSM_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<MSSM<Two_scale>*>(model_);
}

void MSSM_susy_scale_constraint<Two_scale>::set_input_parameters(const MSSM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void MSSM_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void MSSM_susy_scale_constraint<Two_scale>::initialize()
{
   assert(model && "MSSM_susy_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = SM(MZ);

   scale = initial_scale_guess;
}

void MSSM_susy_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "MSSM_susy_scale_constraint<Two_scale>::"
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
