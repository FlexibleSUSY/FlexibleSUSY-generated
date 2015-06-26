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

// File generated at Fri 26 Jun 2015 18:59:20

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

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME MRSSM<Two_scale>

MRSSM_susy_scale_constraint<Two_scale>::MRSSM_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

MRSSM_susy_scale_constraint<Two_scale>::MRSSM_susy_scale_constraint(
   MRSSM<Two_scale>* model_)
   : Constraint<Two_scale>()
   , model(model_)
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

   MODEL->set_LamTD(Re(LamTDInput));
   MODEL->set_LamTU(Re(LamTUInput));
   MODEL->set_LamSD(Re(LamSDInput));
   MODEL->set_LamSU(Re(LamSUInput));
   MODEL->set_Mu(Re(MuInput));
   MODEL->set_MuD(Re(MuDInput));
   MODEL->set_MuU(Re(MuUInput));
   MODEL->set_vT(Re(vTInput));
   MODEL->set_vS(Re(vSInput));
   MODEL->set_BMu(Re(BMuInput));
   MODEL->set_BMuD(Re(BMuDInput));
   MODEL->set_BMuU(Re(BMuUInput));
   MODEL->set_mq2((mq2Input).real());
   MODEL->set_ml2((ml2Input).real());
   MODEL->set_md2((md2Input).real());
   MODEL->set_mu2((mu2Input).real());
   MODEL->set_me2((me2Input).real());
   MODEL->set_moc2(Re(moc2Input));
   MODEL->set_mRd2(Re(mRd2Input));
   MODEL->set_mRu2(Re(mRu2Input));
   MODEL->set_MDBS(Re(MDBSInput));
   MODEL->set_MDWBT(Re(MDWBTInput));
   MODEL->set_MDGoc(Re(MDGocInput));


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
   assert(model && "Error: MRSSM_susy_scale_constraint::"
          "get_input_parameters(): model pointer is zero.");

   return model->get_input();
}

MRSSM<Two_scale>* MRSSM_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void MRSSM_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<MRSSM<Two_scale>*>(model_);
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
