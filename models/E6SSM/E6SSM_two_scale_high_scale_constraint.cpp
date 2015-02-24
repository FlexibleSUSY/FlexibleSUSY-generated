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

// File generated at Tue 24 Feb 2015 17:35:55

#include "E6SSM_two_scale_high_scale_constraint.hpp"
#include "E6SSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "numerics.hpp"

#include <cassert>
#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define INPUTPARAMETER(p) inputPars.p
#define MODELPARAMETER(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define SM(p) Electroweak_constants::p
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define MODEL model
#define MODELCLASSNAME E6SSM<Two_scale>

E6SSM_high_scale_constraint<Two_scale>::E6SSM_high_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , inputPars()
{
}

E6SSM_high_scale_constraint<Two_scale>::E6SSM_high_scale_constraint(
   E6SSM<Two_scale>* model_,
   const E6SSM_input_parameters& inputPars_)
   : Constraint<Two_scale>()
   , model(model_)
   , inputPars(inputPars_)
{
   initialize();
}

E6SSM_high_scale_constraint<Two_scale>::~E6SSM_high_scale_constraint()
{
}

void E6SSM_high_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: E6SSM_high_scale_constraint::apply():"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("E6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(1.0);
   }
   if (std::fabs(model->get_g2()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("E6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(1.0);
   }
   if (std::fabs(model->get_g3()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("E6SSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(1.0);
   }

   update_scale();

   const auto Azero = INPUTPARAMETER(Azero);
   const auto m0 = INPUTPARAMETER(m0);
   const auto KappaInput = INPUTPARAMETER(KappaInput);
   const auto Lambda12Input = INPUTPARAMETER(Lambda12Input);
   const auto LambdaInput = INPUTPARAMETER(LambdaInput);
   const auto m12 = INPUTPARAMETER(m12);
   const auto g1 = MODELPARAMETER(g1);
   const auto Ye = MODELPARAMETER(Ye);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);

   MODEL->set_gN(g1);
   MODEL->set_TYe(Azero*Ye);
   MODEL->set_TYd(Azero*Yd);
   MODEL->set_TYu(Azero*Yu);
   MODEL->set_mq2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_ml2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_md2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mu2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_me2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mDx2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mDxbar2(Sqr(m0)*UNITMATRIX(3));
   MODEL->set_mH1I2(Sqr(m0)*UNITMATRIX(2));
   MODEL->set_mH2I2(Sqr(m0)*UNITMATRIX(2));
   MODEL->set_msI2(Sqr(m0)*UNITMATRIX(2));
   MODEL->set_mHp2(Sqr(m0));
   MODEL->set_mHpbar2(Sqr(m0));
   MODEL->set_Kappa(KappaInput*UNITMATRIX(3));
   MODEL->set_Lambda12(Lambda12Input*UNITMATRIX(2));
   MODEL->set_Lambdax(LambdaInput);
   MODEL->set_TKappa(Azero*KappaInput*UNITMATRIX(3));
   MODEL->set_TLambda12(Azero*Lambda12Input*UNITMATRIX(2));
   MODEL->set_TLambdax(Azero*LambdaInput);
   MODEL->set_MassB(m12);
   MODEL->set_MassWB(m12);
   MODEL->set_MassG(m12);
   MODEL->set_MassBp(m12);

   {
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto g3 = MODELPARAMETER(g3);
      const auto gN = MODELPARAMETER(gN);
      const auto Yd = MODELPARAMETER(Yd);
      const auto Ye = MODELPARAMETER(Ye);
      const auto Kappa = MODELPARAMETER(Kappa);
      const auto Lambda12 = MODELPARAMETER(Lambda12);
      const auto Lambdax = MODELPARAMETER(Lambdax);
      const auto Yu = MODELPARAMETER(Yu);

      if (MaxAbsValue(g1) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("g1", MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("g1");
      if (MaxAbsValue(g2) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("g2", MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("g2");
      if (MaxAbsValue(g3) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("g3", MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("g3");
      if (MaxAbsValue(gN) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("gN", MaxAbsValue(gN), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("gN");
      if (MaxAbsValue(Yd) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Yd", MaxAbsValue(Yd), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Yd");
      if (MaxAbsValue(Ye) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Ye", MaxAbsValue(Ye), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Ye");
      if (MaxAbsValue(Kappa) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Kappa", MaxAbsValue(Kappa), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Kappa");
      if (MaxAbsValue(Lambda12) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Lambda12", MaxAbsValue(Lambda12), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Lambda12");
      if (MaxAbsValue(Lambdax) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Lambdax", MaxAbsValue(Lambdax), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Lambdax");
      if (MaxAbsValue(Yu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Yu", MaxAbsValue(Yu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Yu");

   }
}

double E6SSM_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double E6SSM_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void E6SSM_high_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<E6SSM<Two_scale>*>(model_);
}

void E6SSM_high_scale_constraint<Two_scale>::set_input_parameters(const E6SSM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void E6SSM_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void E6SSM_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void E6SSM_high_scale_constraint<Two_scale>::initialize()
{
   assert(model && "E6SSM_high_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = 1.e16;

   scale = initial_scale_guess;
}

void E6SSM_high_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "E6SSM_high_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const double currentScale = model->get_scale();
   const E6SSM_soft_parameters beta_functions(model->calc_beta());

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto beta_g1 = BETAPARAMETER(g1);
   const auto beta_g2 = BETAPARAMETER(g2);

   scale = currentScale*exp((-g1 + g2)/(BETA(g1) - BETA(g2)));


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("E6SSM_high_scale_constraint<Two_scale>: Overflow error"
            " during calculation of high scale: " << strerror(errno) << '\n'
            << "   current scale = " << currentScale << '\n'
            << "   new scale = " << scale << '\n'
            << "   resetting scale to " << get_initial_scale_guess());
#endif
      scale = get_initial_scale_guess();
      errno = 0;
   }

   if (scale < 10000.) {
   #ifdef ENABLE_VERBOSE
      WARNING("scale < 10000.");
   #endif
      scale = 10000.;
   }
   if (scale > 5.e17) {
   #ifdef ENABLE_VERBOSE
      WARNING("scale > 5.e17");
   #endif
      scale = 5.e17;
   }

}

} // namespace flexiblesusy
