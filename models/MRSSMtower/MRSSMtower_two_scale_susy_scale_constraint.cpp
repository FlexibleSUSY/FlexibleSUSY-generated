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

// File generated at Thu 15 Dec 2016 12:39:52

#include "MRSSMtower_two_scale_susy_scale_constraint.hpp"
#include "MRSSMtower_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"

#include <cassert>
#include <cmath>

namespace flexiblesusy {

#define DERIVEDPARAMETER(p) model->p()
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
#define MODELCLASSNAME MRSSMtower<Two_scale>

MRSSMtower_susy_scale_constraint<Two_scale>::MRSSMtower_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , qedqcd()
{
}

MRSSMtower_susy_scale_constraint<Two_scale>::MRSSMtower_susy_scale_constraint(
   MRSSMtower<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : Constraint<Two_scale>()
   , model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

MRSSMtower_susy_scale_constraint<Two_scale>::~MRSSMtower_susy_scale_constraint()
{
}

void MRSSMtower_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: MRSSMtower_susy_scale_constraint::apply():"
          " model pointer must not be zero");



   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto mS2Input = INPUTPARAMETER(mS2Input);
   const auto mT2Input = INPUTPARAMETER(mT2Input);
   const auto moc2Input = INPUTPARAMETER(moc2Input);
   const auto mRd2Input = INPUTPARAMETER(mRd2Input);
   const auto mRu2Input = INPUTPARAMETER(mRu2Input);
   const auto BMuInput = INPUTPARAMETER(BMuInput);
   const auto LamSDInput = INPUTPARAMETER(LamSDInput);
   const auto LamSUInput = INPUTPARAMETER(LamSUInput);
   const auto LamTDInput = INPUTPARAMETER(LamTDInput);
   const auto LamTUInput = INPUTPARAMETER(LamTUInput);
   const auto MDBSInput = INPUTPARAMETER(MDBSInput);
   const auto MDGocInput = INPUTPARAMETER(MDGocInput);
   const auto MDWBTInput = INPUTPARAMETER(MDWBTInput);
   const auto MuDInput = INPUTPARAMETER(MuDInput);
   const auto MuUInput = INPUTPARAMETER(MuUInput);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);

   MODEL->set_vu(Re(TanBeta*Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_vd(Re(Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_mq2((mq2Input).real());
   MODEL->set_mu2((mu2Input).real());
   MODEL->set_md2((md2Input).real());
   MODEL->set_ml2((ml2Input).real());
   MODEL->set_me2((me2Input).real());
   MODEL->set_mS2(Re(mS2Input));
   MODEL->set_mT2(Re(mT2Input));
   MODEL->set_moc2(Re(moc2Input));
   MODEL->set_mRd2(Re(mRd2Input));
   MODEL->set_mRu2(Re(mRu2Input));
   MODEL->set_BMu(Re(BMuInput));
   MODEL->set_LamSD(Re(LamSDInput));
   MODEL->set_LamSU(Re(LamSUInput));
   MODEL->set_LamTD(Re(LamTDInput));
   MODEL->set_LamTU(Re(LamTUInput));
   MODEL->set_MDBS(Re(MDBSInput));
   MODEL->set_MDGoc(Re(MDGocInput));
   MODEL->set_MDWBT(Re(MDWBTInput));
   MODEL->set_MuD(Re(MuDInput));
   MODEL->set_MuU(Re(MuUInput));
   MODEL->set_Mu(Re(0));
   MODEL->set_BMuD(Re(0));
   MODEL->set_BMuU(Re(0));
   MODEL->solve_ewsb();


}

double MRSSMtower_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MRSSMtower_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const MRSSMtower_input_parameters& MRSSMtower_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   assert(model && "Error: MRSSMtower_susy_scale_constraint::"
          "get_input_parameters(): model pointer is zero.");

   return model->get_input();
}

MRSSMtower<Two_scale>* MRSSMtower_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void MRSSMtower_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<MRSSMtower<Two_scale>*>(model_);
}

void MRSSMtower_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& MRSSMtower_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void MRSSMtower_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
   qedqcd = softsusy::QedQcd();
}

void MRSSMtower_susy_scale_constraint<Two_scale>::initialize()
{
   assert(model && "MRSSMtower_susy_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   const auto MS = INPUTPARAMETER(MS);

   initial_scale_guess = MS;

   scale = initial_scale_guess;
}

void MRSSMtower_susy_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "MRSSMtower_susy_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const auto MS = INPUTPARAMETER(MS);

   scale = MS;


}

} // namespace flexiblesusy
