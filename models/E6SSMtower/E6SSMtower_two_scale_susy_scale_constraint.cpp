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

// File generated at Thu 15 Dec 2016 12:44:07

#include "E6SSMtower_two_scale_susy_scale_constraint.hpp"
#include "E6SSMtower_two_scale_model.hpp"
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
#define MODELCLASSNAME E6SSMtower<Two_scale>

E6SSMtower_susy_scale_constraint<Two_scale>::E6SSMtower_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , qedqcd()
{
}

E6SSMtower_susy_scale_constraint<Two_scale>::E6SSMtower_susy_scale_constraint(
   E6SSMtower<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : Constraint<Two_scale>()
   , model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

E6SSMtower_susy_scale_constraint<Two_scale>::~E6SSMtower_susy_scale_constraint()
{
}

void E6SSMtower_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: E6SSMtower_susy_scale_constraint::apply():"
          " model pointer must not be zero");



   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto gNInput = INPUTPARAMETER(gNInput);
   const auto M1Input = INPUTPARAMETER(M1Input);
   const auto M2Input = INPUTPARAMETER(M2Input);
   const auto M3Input = INPUTPARAMETER(M3Input);
   const auto M4Input = INPUTPARAMETER(M4Input);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto mDx2Input = INPUTPARAMETER(mDx2Input);
   const auto mDxbar2Input = INPUTPARAMETER(mDxbar2Input);
   const auto mH1I2Input = INPUTPARAMETER(mH1I2Input);
   const auto mH2I2Input = INPUTPARAMETER(mH2I2Input);
   const auto msI2Input = INPUTPARAMETER(msI2Input);
   const auto mHp2Input = INPUTPARAMETER(mHp2Input);
   const auto mHpbar2Input = INPUTPARAMETER(mHpbar2Input);
   const auto BMuPrInput = INPUTPARAMETER(BMuPrInput);
   const auto MuPrInput = INPUTPARAMETER(MuPrInput);
   const auto KappaInput = INPUTPARAMETER(KappaInput);
   const auto Lambda12Input = INPUTPARAMETER(Lambda12Input);
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);
   const auto AKappaInput = INPUTPARAMETER(AKappaInput);
   const auto ALambda12Input = INPUTPARAMETER(ALambda12Input);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto mAInput = INPUTPARAMETER(mAInput);
   const auto AuInput = INPUTPARAMETER(AuInput);
   const auto AdInput = INPUTPARAMETER(AdInput);
   const auto AeInput = INPUTPARAMETER(AeInput);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto vs = MODELPARAMETER(vs);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_vu(Re(TanBeta*Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_vd(Re(Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_gN(Re(gNInput));
   MODEL->set_MassB(Re(M1Input));
   MODEL->set_MassWB(Re(M2Input));
   MODEL->set_MassG(Re(M3Input));
   MODEL->set_MassBp(Re(M4Input));
   MODEL->set_mq2((mq2Input).real());
   MODEL->set_mu2((mu2Input).real());
   MODEL->set_md2((md2Input).real());
   MODEL->set_ml2((ml2Input).real());
   MODEL->set_me2((me2Input).real());
   MODEL->set_mDx2((mDx2Input).real());
   MODEL->set_mDxbar2((mDxbar2Input).real());
   MODEL->set_mH1I2((mH1I2Input).real());
   MODEL->set_mH2I2((mH2I2Input).real());
   MODEL->set_msI2((msI2Input).real());
   MODEL->set_mHp2(Re(mHp2Input));
   MODEL->set_mHpbar2(Re(mHpbar2Input));
   MODEL->set_BMuPr(Re(BMuPrInput));
   MODEL->set_MuPr(Re(MuPrInput));
   MODEL->set_Kappa((KappaInput).real());
   MODEL->set_Lambda12((Lambda12Input).real());
   MODEL->set_Lambdax(Re(LambdaInput));
   MODEL->set_TKappa((AKappaInput*KappaInput).real());
   MODEL->set_TLambda12((ALambda12Input*Lambda12Input).real());
   MODEL->set_vs(Re((1.4142135623730951*MuInput)/LambdaInput));
   MODEL->set_TLambdax(Re((1.4142135623730951*Sqr(mAInput))/((1/TanBeta +
      TanBeta)*vs)));
   MODEL->set_TYu((AuInput*Yu).real());
   MODEL->set_TYd((AdInput*Yd).real());
   MODEL->set_TYe((AeInput*Ye).real());
   MODEL->solve_ewsb();


}

double E6SSMtower_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double E6SSMtower_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const E6SSMtower_input_parameters& E6SSMtower_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   assert(model && "Error: E6SSMtower_susy_scale_constraint::"
          "get_input_parameters(): model pointer is zero.");

   return model->get_input();
}

E6SSMtower<Two_scale>* E6SSMtower_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void E6SSMtower_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<E6SSMtower<Two_scale>*>(model_);
}

void E6SSMtower_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& E6SSMtower_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void E6SSMtower_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
   qedqcd = softsusy::QedQcd();
}

void E6SSMtower_susy_scale_constraint<Two_scale>::initialize()
{
   assert(model && "E6SSMtower_susy_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   initial_scale_guess = MSUSY;

   scale = initial_scale_guess;
}

void E6SSMtower_susy_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "E6SSMtower_susy_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   scale = MSUSY;


}

} // namespace flexiblesusy
