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

// File generated at Sun 26 Aug 2018 14:36:00

#include "lowNMSSMTanBetaAtMZ_two_scale_susy_scale_constraint.hpp"
#include "lowNMSSMTanBetaAtMZ_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
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
#define MODELCLASSNAME lowNMSSMTanBetaAtMZ<Two_scale>

lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::lowNMSSMTanBetaAtMZ_susy_scale_constraint(
   lowNMSSMTanBetaAtMZ<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
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
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

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

double lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const lowNMSSMTanBetaAtMZ_input_parameters& lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   check_model_ptr();

   return model->get_input();
}

lowNMSSMTanBetaAtMZ<Two_scale>* lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<lowNMSSMTanBetaAtMZ<Two_scale>*>(model_);
}

void lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
}

void lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   initial_scale_guess = Qin;

   scale = initial_scale_guess;
}

void lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto Qin = INPUTPARAMETER(Qin);

   scale = Qin;


}

void lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("lowNMSSMTanBetaAtMZ_susy_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
