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


#include "NMSSMEFTHiggs_two_scale_susy_scale_constraint.hpp"
#include "NMSSMEFTHiggs_two_scale_model.hpp"
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
#define MODELCLASSNAME NMSSMEFTHiggs<Two_scale>

NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::NMSSMEFTHiggs_susy_scale_constraint(
   NMSSMEFTHiggs<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

void NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   model->calculate_DRbar_masses();
   update_scale();

   // apply user-defined susy scale constraints
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto M1Input = INPUTPARAMETER(M1Input);
   const auto M2Input = INPUTPARAMETER(M2Input);
   const auto M3Input = INPUTPARAMETER(M3Input);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto AuInput = INPUTPARAMETER(AuInput);
   const auto AdInput = INPUTPARAMETER(AdInput);
   const auto AeInput = INPUTPARAMETER(AeInput);
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto KappaInput = INPUTPARAMETER(KappaInput);
   const auto AKappaInput = INPUTPARAMETER(AKappaInput);
   const auto ALambdaInput = INPUTPARAMETER(ALambdaInput);
   const auto vd = MODELPARAMETER(vd);
   const auto vu = MODELPARAMETER(vu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_vu(Re(TanBeta*Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_vd(Re(Sqrt((Sqr(vd) + Sqr(vu))/(1 + Sqr(TanBeta)))));
   MODEL->set_MassB(Re(M1Input));
   MODEL->set_MassWB(Re(M2Input));
   MODEL->set_MassG(Re(M3Input));
   MODEL->set_mq2((mq2Input).real());
   MODEL->set_mu2((mu2Input).real());
   MODEL->set_md2((md2Input).real());
   MODEL->set_ml2((ml2Input).real());
   MODEL->set_me2((me2Input).real());
   MODEL->set_TYu(((AuInput).cwiseProduct(Yu)).real());
   MODEL->set_TYd(((AdInput).cwiseProduct(Yd)).real());
   MODEL->set_TYe(((AeInput).cwiseProduct(Ye)).real());
   MODEL->set_vS(Re((1.4142135623730951*MuInput)/LambdaInput));
   MODEL->set_Kappa(Re(KappaInput));
   MODEL->set_Lambdax(Re(LambdaInput));
   MODEL->set_TKappa(Re(AKappaInput*KappaInput));
   MODEL->set_TLambdax(Re(ALambdaInput*LambdaInput));
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

double NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const NMSSMEFTHiggs_input_parameters& NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::get_input_parameters() const
{
   check_model_ptr();

   return model->get_input();
}

NMSSMEFTHiggs<Two_scale>* NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<NMSSMEFTHiggs<Two_scale>*>(model_);
}

void NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
   qedqcd = softsusy::QedQcd();
}

void NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   initial_scale_guess = MSUSY;

   scale = initial_scale_guess;
}

void NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   scale = MSUSY;


}

void NMSSMEFTHiggs_susy_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("NMSSMEFTHiggs_susy_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
