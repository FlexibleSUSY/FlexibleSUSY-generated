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

// File generated at Tue 24 Feb 2015 17:52:38

#include "MSSM_two_scale_high_scale_constraint.hpp"
#include "MSSM_two_scale_model.hpp"
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
#define MODELCLASSNAME MSSM<Two_scale>

MSSM_high_scale_constraint<Two_scale>::MSSM_high_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
   , inputPars()
{
}

MSSM_high_scale_constraint<Two_scale>::MSSM_high_scale_constraint(
   MSSM<Two_scale>* model_,
   const MSSM_input_parameters& inputPars_)
   : Constraint<Two_scale>()
   , model(model_)
   , inputPars(inputPars_)
{
   initialize();
}

MSSM_high_scale_constraint<Two_scale>::~MSSM_high_scale_constraint()
{
}

void MSSM_high_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: MSSM_high_scale_constraint::apply():"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("MSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(1.0);
   }
   if (std::fabs(model->get_g2()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("MSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(1.0);
   }
   if (std::fabs(model->get_g3()) > 3.0) {
#ifdef ENABLE_VERBOSE
      ERROR("MSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(1.0);
   }

   update_scale();

   const auto Aeij = INPUTPARAMETER(Aeij);
   const auto Adij = INPUTPARAMETER(Adij);
   const auto Auij = INPUTPARAMETER(Auij);
   const auto mHd2IN = INPUTPARAMETER(mHd2IN);
   const auto mHu2IN = INPUTPARAMETER(mHu2IN);
   const auto mq2Input = INPUTPARAMETER(mq2Input);
   const auto ml2Input = INPUTPARAMETER(ml2Input);
   const auto md2Input = INPUTPARAMETER(md2Input);
   const auto mu2Input = INPUTPARAMETER(mu2Input);
   const auto me2Input = INPUTPARAMETER(me2Input);
   const auto MassBInput = INPUTPARAMETER(MassBInput);
   const auto MassWBInput = INPUTPARAMETER(MassWBInput);
   const auto MassGInput = INPUTPARAMETER(MassGInput);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_TYe(Aeij*Ye);
   MODEL->set_TYd(Adij*Ye);
   MODEL->set_TYu(Auij*Ye);
   MODEL->set_mHd2(mHd2IN);
   MODEL->set_mHu2(mHu2IN);
   MODEL->set_mq2(mq2Input);
   MODEL->set_ml2(ml2Input);
   MODEL->set_md2(md2Input);
   MODEL->set_mu2(mu2Input);
   MODEL->set_me2(me2Input);
   MODEL->set_MassB(MassBInput);
   MODEL->set_MassWB(MassWBInput);
   MODEL->set_MassG(MassGInput);

   {
      const auto g1 = MODELPARAMETER(g1);
      const auto g2 = MODELPARAMETER(g2);
      const auto g3 = MODELPARAMETER(g3);
      const auto Yd = MODELPARAMETER(Yd);
      const auto Ye = MODELPARAMETER(Ye);
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
      if (MaxAbsValue(Yu) > 3.5449077018110318)
         model->get_problems().flag_non_perturbative_parameter_warning("Yu", MaxAbsValue(Yu), model->get_scale(), 3.5449077018110318);
      else
         model->get_problems().unflag_non_perturbative_parameter_warning("Yu");

   }
}

double MSSM_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double MSSM_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void MSSM_high_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<MSSM<Two_scale>*>(model_);
}

void MSSM_high_scale_constraint<Two_scale>::set_input_parameters(const MSSM_input_parameters& inputPars_)
{
   inputPars = inputPars_;
}

void MSSM_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void MSSM_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void MSSM_high_scale_constraint<Two_scale>::initialize()
{
   assert(model && "MSSM_high_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   const auto Qin = INPUTPARAMETER(Qin);

   initial_scale_guess = Qin;

   scale = initial_scale_guess;
}

void MSSM_high_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "MSSM_high_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const double currentScale = model->get_scale();
   const MSSM_soft_parameters beta_functions(model->calc_beta());

   const auto Qin = INPUTPARAMETER(Qin);

   scale = Qin;


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("MSSM_high_scale_constraint<Two_scale>: Overflow error"
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
