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

// File generated at Tue 24 Feb 2015 17:31:34

#include "MRSSM_two_scale_high_scale_constraint.hpp"
#include "MRSSM_two_scale_model.hpp"
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
#define MODELCLASSNAME MRSSM<Two_scale>

MRSSM_high_scale_constraint<Two_scale>::MRSSM_high_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , inputPars()
{
}

MRSSM_high_scale_constraint<Two_scale>::MRSSM_high_scale_constraint(
   MRSSM<Two_scale>* model_,
   const MRSSM_input_parameters& inputPars_)
   : Constraint<Two_scale>()
   , model(model_)
   , inputPars(inputPars_)
{
   initialize();
}

MRSSM_high_scale_constraint<Two_scale>::~MRSSM_high_scale_constraint()
{
}

void MRSSM_high_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: MRSSM_high_scale_constraint::apply():"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("MRSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(1.0);
   }
   if (std::fabs(model->get_g2()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("MRSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(1.0);
   }
   if (std::fabs(model->get_g3()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("MRSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(1.0);
   }

   update_scale();


   {
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto g3 = MODELPARAMETER(g3);
      const auto Yd = MODELPARAMETER(Yd);
      const auto Ye = MODELPARAMETER(Ye);
      const auto LamTD = MODELPARAMETER(LamTD);
      const auto LamTU = MODELPARAMETER(LamTU);
      const auto LamSD = MODELPARAMETER(LamSD);
      const auto LamSU = MODELPARAMETER(LamSU);
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
      if (MaxAbsValue(Yd) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Yd", MaxAbsValue(Yd), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Yd");
      if (MaxAbsValue(Ye) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Ye", MaxAbsValue(Ye), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Ye");
      if (MaxAbsValue(LamTD) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("LamTD", MaxAbsValue(LamTD), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("LamTD");
      if (MaxAbsValue(LamTU) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("LamTU", MaxAbsValue(LamTU), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("LamTU");
      if (MaxAbsValue(LamSD) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("LamSD", MaxAbsValue(LamSD), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("LamSD");
      if (MaxAbsValue(LamSU) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("LamSU", MaxAbsValue(LamSU), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("LamSU");
      if (MaxAbsValue(Yu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Yu", MaxAbsValue(Yu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Yu");

   }
}

double MRSSM_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MRSSM_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void MRSSM_high_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<MRSSM<Two_scale>*>(model_);
}

void MRSSM_high_scale_constraint<Two_scale>::set_input_parameters(const MRSSM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void MRSSM_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void MRSSM_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void MRSSM_high_scale_constraint<Two_scale>::initialize()
{
   assert(model && "MRSSM_high_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   initial_scale_guess = 2.e16;

   scale = initial_scale_guess;
}

void MRSSM_high_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "MRSSM_high_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const double currentScale = model->get_scale();
   const MRSSM_soft_parameters beta_functions(model->calc_beta());

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto beta_g1 = BETAPARAMETER(g1);
   const auto beta_g2 = BETAPARAMETER(g2);

   scale = currentScale*exp((-g1 + g2)/(BETA(g1) - BETA(g2)));


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("MRSSM_high_scale_constraint<Two_scale>: Overflow error"
            " during calculation of high scale: " << strerror(errno) << '\n'
            << "   current scale = " << currentScale << '\n'
            << "   new scale = " << scale << '\n'
            << "   resetting scale to " << get_initial_scale_guess());
#endif
      scale = get_initial_scale_guess();
      errno = 0;
   }


}

} // namespace flexiblesusy
