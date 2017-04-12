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

// File generated at Wed 12 Apr 2017 12:39:07

#include "lowNMSSM_two_scale_susy_scale_constraint.hpp"
#include "lowNMSSM_two_scale_model.hpp"
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
#define MODELCLASSNAME lowNMSSM<Two_scale>

lowNMSSM_susy_scale_constraint<Two_scale>::lowNMSSM_susy_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , qedqcd()
{
}

lowNMSSM_susy_scale_constraint<Two_scale>::lowNMSSM_susy_scale_constraint(
   lowNMSSM<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : Constraint<Two_scale>()
   , model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

lowNMSSM_susy_scale_constraint<Two_scale>::~lowNMSSM_susy_scale_constraint()
{
}

void lowNMSSM_susy_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: lowNMSSM_susy_scale_constraint::apply():"
          " model pointer must not be zero");



   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto KappaInput = INPUTPARAMETER(KappaInput);
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);
   const auto AKappaInput = INPUTPARAMETER(AKappaInput);
   const auto ALambdaInput = INPUTPARAMETER(ALambdaInput);
   const auto MuEffInput = INPUTPARAMETER(MuEffInput);
   const auto M1Input = INPUTPARAMETER(M1Input);
   const auto M2Input = INPUTPARAMETER(M2Input);
   const auto M3Input = INPUTPARAMETER(M3Input);
   const auto me1Input = INPUTPARAMETER(me1Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto me3Input = INPUTPARAMETER(me3Input);
   const auto ml1Input = INPUTPARAMETER(ml1Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto ml3Input = INPUTPARAMETER(ml3Input);
   const auto md1Input = INPUTPARAMETER(md1Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto md3Input = INPUTPARAMETER(md3Input);
   const auto mu1Input = INPUTPARAMETER(mu1Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto mu3Input = INPUTPARAMETER(mu3Input);
   const auto mq1Input = INPUTPARAMETER(mq1Input);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto mq3Input = INPUTPARAMETER(mq3Input);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto AbInput = INPUTPARAMETER(AbInput);
   const auto ATauInput = INPUTPARAMETER(ATauInput);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_vd(Re(Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_vu(Re(TanBeta*Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_Kappa(Re(KappaInput));
   MODEL->set_Lambdax(Re(LambdaInput));
   MODEL->set_TKappa(Re(AKappaInput*KappaInput));
   MODEL->set_TLambdax(Re(ALambdaInput*LambdaInput));
   MODEL->set_vS(Re((1.4142135623730951*MuEffInput)/LambdaInput));
   MODEL->set_MassB(Re(M1Input));
   MODEL->set_MassWB(Re(M2Input));
   MODEL->set_MassG(Re(M3Input));
   MODEL->set_me2(0,0,Re(Sqr(me1Input)));
   MODEL->set_me2(1,1,Re(Sqr(me2Input)));
   MODEL->set_me2(2,2,Re(Sqr(me3Input)));
   MODEL->set_ml2(0,0,Re(Sqr(ml1Input)));
   MODEL->set_ml2(1,1,Re(Sqr(ml2Input)));
   MODEL->set_ml2(2,2,Re(Sqr(ml3Input)));
   MODEL->set_md2(0,0,Re(Sqr(md1Input)));
   MODEL->set_md2(1,1,Re(Sqr(md2Input)));
   MODEL->set_md2(2,2,Re(Sqr(md3Input)));
   MODEL->set_mu2(0,0,Re(Sqr(mu1Input)));
   MODEL->set_mu2(1,1,Re(Sqr(mu2Input)));
   MODEL->set_mu2(2,2,Re(Sqr(mu3Input)));
   MODEL->set_mq2(0,0,Re(Sqr(mq1Input)));
   MODEL->set_mq2(1,1,Re(Sqr(mq2Input)));
   MODEL->set_mq2(2,2,Re(Sqr(mq3Input)));
   MODEL->set_me2(0,0,Re(Sqr(me1Input)));
   MODEL->set_me2(1,1,Re(Sqr(me2Input)));
   MODEL->set_TYu(2,2,Re(AtInput*Yu(2,2)));
   MODEL->set_TYd(2,2,Re(AbInput*Yd(2,2)));
   MODEL->set_TYe(2,2,Re(ATauInput*Ye(2,2)));
   MODEL->solve_ewsb();


}

double lowNMSSM_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double lowNMSSM_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const lowNMSSM_input_parameters& lowNMSSM_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   assert(model && "Error: lowNMSSM_susy_scale_constraint::"
          "get_input_parameters(): model pointer is zero.");

   return model->get_input();
}

lowNMSSM<Two_scale>* lowNMSSM_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void lowNMSSM_susy_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<lowNMSSM<Two_scale>*>(model_);
}

void lowNMSSM_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& lowNMSSM_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void lowNMSSM_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
   qedqcd = softsusy::QedQcd();
}

void lowNMSSM_susy_scale_constraint<Two_scale>::initialize()
{
   assert(model && "lowNMSSM_susy_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   const auto Qin = INPUTPARAMETER(Qin);

   initial_scale_guess = Qin;

   scale = initial_scale_guess;
}

void lowNMSSM_susy_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "lowNMSSM_susy_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const auto Qin = INPUTPARAMETER(Qin);

   scale = Qin;


}

} // namespace flexiblesusy
