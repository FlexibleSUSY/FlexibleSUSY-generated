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

// File generated at Tue 24 Feb 2015 17:31:35

#include "MRSSM_two_scale_susy_scale_constraint.hpp"
#include "MRSSM_two_scale_model.hpp"
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
#define MODELCLASSNAME MRSSM<Two_scale>

MRSSM_susy_scale_constraint<Two_scale>::MRSSM_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , inputPars()
{
}

MRSSM_susy_scale_constraint<Two_scale>::MRSSM_susy_scale_constraint(
   MRSSM<Two_scale>* model_,
   const MRSSM_input_parameters& inputPars_)
   : Constraint<Two_scale>()
   , model(model_)
   , inputPars(inputPars_)
{
   initialize();
}

MRSSM_susy_scale_constraint<Two_scale>::~MRSSM_susy_scale_constraint()
{
}

void MRSSM_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: MRSSM_susy_scale_constraint::apply():"
          " model pointer must not be zero");

   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto LamTDInput = INPUTPARAMETER(LamTDInput);
   const auto LamTUInput = INPUTPARAMETER(LamTUInput);
   const auto LamSDInput = INPUTPARAMETER(LamSDInput);
   const auto LamSUInput = INPUTPARAMETER(LamSUInput);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto MuDInput = INPUTPARAMETER(MuDInput);
   const auto MuUInput = INPUTPARAMETER(MuUInput);
   const auto vTInput = INPUTPARAMETER(vTInput);
   const auto vSInput = INPUTPARAMETER(vSInput);
   const auto BMuInput = INPUTPARAMETER(BMuInput);
   const auto BMuDInput = INPUTPARAMETER(BMuDInput);
   const auto BMuUInput = INPUTPARAMETER(BMuUInput);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto moc2Input = INPUTPARAMETER(moc2Input);
   const auto mRd2Input = INPUTPARAMETER(mRd2Input);
   const auto mRu2Input = INPUTPARAMETER(mRu2Input);
   const auto MDBSInput = INPUTPARAMETER(MDBSInput);
   const auto MDWBTInput = INPUTPARAMETER(MDWBTInput);
   const auto MDGocInput = INPUTPARAMETER(MDGocInput);

   MODEL->set_LamTD(LamTDInput);
   MODEL->set_LamTU(LamTUInput);
   MODEL->set_LamSD(LamSDInput);
   MODEL->set_LamSU(LamSUInput);
   MODEL->set_Mu(MuInput);
   MODEL->set_MuD(MuDInput);
   MODEL->set_MuU(MuUInput);
   MODEL->set_vT(vTInput);
   MODEL->set_vS(vSInput);
   MODEL->set_BMu(BMuInput);
   MODEL->set_BMuD(BMuDInput);
   MODEL->set_BMuU(BMuUInput);
   MODEL->set_mq2(mq2Input);
   MODEL->set_ml2(ml2Input);
   MODEL->set_md2(md2Input);
   MODEL->set_mu2(mu2Input);
   MODEL->set_me2(me2Input);
   MODEL->set_moc2(moc2Input);
   MODEL->set_mRd2(mRd2Input);
   MODEL->set_mRu2(mRu2Input);
   MODEL->set_MDBS(MDBSInput);
   MODEL->set_MDWBT(MDWBTInput);
   MODEL->set_MDGoc(MDGocInput);


   // the parameters, which are fixed by the EWSB eqs., will now be
   // defined at this scale (at the EWSB loop level defined in the
   // model)
   model->solve_ewsb();
}

double MRSSM_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MRSSM_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const MRSSM_input_parameters& MRSSM_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   return inputPars;
}

MRSSM<Two_scale>* MRSSM_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void MRSSM_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<MRSSM<Two_scale>*>(model_);
}

void MRSSM_susy_scale_constraint<Two_scale>::set_input_parameters(const MRSSM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void MRSSM_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void MRSSM_susy_scale_constraint<Two_scale>::initialize()
{
   assert(model && "MRSSM_susy_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = 1000;

   scale = initial_scale_guess;
}

void MRSSM_susy_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "MRSSM_susy_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   scale = 1000;


}

} // namespace flexiblesusy
