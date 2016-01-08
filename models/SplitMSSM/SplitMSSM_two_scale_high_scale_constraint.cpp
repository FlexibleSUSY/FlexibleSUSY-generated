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

// File generated at Fri 8 Jan 2016 11:56:25

#include "SplitMSSM_two_scale_high_scale_constraint.hpp"
#include "SplitMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#include <cassert>
#include <cmath>
#include <cerrno>
#include <cstring>

namespace flexiblesusy {

#define INPUTPARAMETER(p) model->get_input().p
#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole Electroweak_constants::MZ
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME SplitMSSM<Two_scale>

SplitMSSM_high_scale_constraint<Two_scale>::SplitMSSM_high_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , initial_scale_guess(0.)
   , model(0)
{
}

SplitMSSM_high_scale_constraint<Two_scale>::SplitMSSM_high_scale_constraint(
   SplitMSSM<Two_scale>* model_)
   : Constraint<Two_scale>()
   , model(model_)
{
   initialize();
}

SplitMSSM_high_scale_constraint<Two_scale>::~SplitMSSM_high_scale_constraint()
{
}

void SplitMSSM_high_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: SplitMSSM_high_scale_constraint::apply():"
          " model pointer must not be zero");

   if (std::fabs(model->get_g1()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("SplitMSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g1 = " << model->get_g1());
#endif
      model->set_g1(3.54491);
   }
   if (std::fabs(model->get_g2()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("SplitMSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g2 = " << model->get_g2());
#endif
      model->set_g2(3.54491);
   }
   if (std::fabs(model->get_g3()) > 3.54491) {
#ifdef ENABLE_VERBOSE
      ERROR("SplitMSSM_high_scale_constraint: Non-perturbative gauge "
            "coupling g3 = " << model->get_g3());
#endif
      model->set_g3(3.54491);
   }

   update_scale();

   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto mAInput = INPUTPARAMETER(mAInput);
   const auto msq2 = INPUTPARAMETER(msq2);
   const auto msu2 = INPUTPARAMETER(msu2);
   const auto msd2 = INPUTPARAMETER(msd2);
   const auto mse2 = INPUTPARAMETER(mse2);
   const auto msl2 = INPUTPARAMETER(msl2);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);

   MODEL->set_Lambdax(Re(0.25*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)
      )) + If(IsCloseRel(msq2(2,2),msu2(2,2),0.01), 0.000053468657576480914*Sqr(g3
      )*(-18 - Power(AtInput - Mu/TanBeta,4)/Abs(msq2(2,2)*msu2(2,2)) + (12*Sqr(
      AtInput - Mu/TanBeta))/Power(Abs(msq2(2,2)*msu2(2,2)),0.5) + 6*Log(msq2(2,2)
      /Sqr(SCALE))*(4 + Power(AtInput - Mu/TanBeta,4)/Abs(msq2(2,2)*msu2(2,2)) - (
      12*Sqr(AtInput - Mu/TanBeta))/Power(Abs(msq2(2,2)*msu2(2,2)),0.5)) - 36*Sqr(
      Log(msq2(2,2)/Sqr(SCALE))))*Power(Yu(2,2),4), -0.0003208119454588855*Sqr(g3)
      *(3 + 4*Log(Power(Abs(msq2(2,2)/msu2(2,2)),0.5)) - 4*(1 + 3*Log(Power(Abs(
      msq2(2,2)/msu2(2,2)),0.5)))*Log(msq2(2,2)/Sqr(SCALE)) + 8*Sqr(Log(Power(Abs(
      msq2(2,2)/msu2(2,2)),0.5))) + (Sqr(AtInput - Mu/TanBeta)*((12*Power(Abs(msq2
      (2,2)/msu2(2,2)),0.5)*Log(Power(Abs(msq2(2,2)/msu2(2,2)),0.5))*(-1 + 2*Log(
      msq2(2,2)/Sqr(SCALE))))/(-1 + Abs(msq2(2,2)/msu2(2,2))) - (16*(-2 + Abs(msq2
      (2,2)/msu2(2,2)))*Power(Abs(msq2(2,2)/msu2(2,2)),0.5)*Sqr(Log(Power(Abs(msq2
      (2,2)/msu2(2,2)),0.5))))/Sqr(-1 + Abs(msq2(2,2)/msu2(2,2)))))/Power(Abs(msq2
      (2,2)*msu2(2,2)),0.5) + (Power(AtInput - Mu/TanBeta,4)*((6*Abs(msq2(2,2)
      /msu2(2,2))*(5 + Abs(msq2(2,2)/msu2(2,2)))*Log(Power(Abs(msq2(2,2)/msu2(2,2)
      ),0.5)))/Power(-1 + Abs(msq2(2,2)/msu2(2,2)),3) - (10*Abs(msq2(2,2)/msu2(2,2
      )))/Sqr(-1 + Abs(msq2(2,2)/msu2(2,2))) + (12*Abs(msq2(2,2)/msu2(2,2))*(1 - (
      (1 + Abs(msq2(2,2)/msu2(2,2)))*Log(Power(Abs(msq2(2,2)/msu2(2,2)),0.5)))/(-1
      + Abs(msq2(2,2)/msu2(2,2))))*Log(msq2(2,2)/Sqr(SCALE)))/Sqr(-1 + Abs(msq2(2
      ,2)/msu2(2,2))) + (4*Abs(msq2(2,2)/msu2(2,2))*(-5 - 4*Abs(msq2(2,2)/msu2(2,2
      )) + Sqr(Abs(msq2(2,2)/msu2(2,2))))*Sqr(Log(Power(Abs(msq2(2,2)/msu2(2,2)),
      0.5))))/Power(-1 + Abs(msq2(2,2)/msu2(2,2)),4)))/Abs(msq2(2,2)*msu2(2,2)) +
      6*Sqr(Log(msq2(2,2)/Sqr(SCALE))))*Power(Yu(2,2),4))*UnitStep(-2 + THRESHOLD)
      + UnitStep(-1 + THRESHOLD)*(0.006332573977646111*(-0.09*Power(g1,4) - 0.3*
      Sqr(g1)*Sqr(g2) - Power(g2,4)*(0.75 - 0.16666666666666666*Sqr(Cos(2*ArcTan(
      TanBeta))))) + 0.006332573977646111*(0.00020833333333333335*Log(Sqr(mAInput)
      /Sqr(SCALE))*(261*Power(g1,4) + 1325*Power(g2,4) + 630*Sqr(g1)*Sqr(g2) - 4*
      Cos(4*ArcTan(TanBeta))*(9*Power(g1,4) + 175*Power(g2,4) + 90*Sqr(g1)*Sqr(g2)
      ) - 9*Cos(8*ArcTan(TanBeta))*Sqr(3*Power(g1,2) + 5*Power(g2,2))) +
      0.0033333333333333335*(6*Power(g1,4)*(Log(msd2(0,0)/Sqr(SCALE)) + Log(msd2(1
      ,1)/Sqr(SCALE)) + Log(msd2(2,2)/Sqr(SCALE))) + 18*Power(g1,4)*(Log(mse2(0,0)
      /Sqr(SCALE)) + Log(mse2(1,1)/Sqr(SCALE)) + Log(mse2(2,2)/Sqr(SCALE))) + (9*
      Power(g1,4) + 25*Power(g2,4))*(Log(msl2(0,0)/Sqr(SCALE)) + Log(msl2(1,1)/Sqr
      (SCALE)) + Log(msl2(2,2)/Sqr(SCALE))) + 3*(Power(g1,4) + 25*Power(g2,4))*(
      Log(msq2(0,0)/Sqr(SCALE)) + Log(msq2(1,1)/Sqr(SCALE)) + Log(msq2(2,2)/Sqr(
      SCALE))) + 24*Power(g1,4)*(Log(msu2(0,0)/Sqr(SCALE)) + Log(msu2(1,1)/Sqr(
      SCALE)) + Log(msu2(2,2)/Sqr(SCALE))))*Sqr(Cos(2*ArcTan(TanBeta))) - 0.1875*
      Sqr(0.6*Power(g1,2) + Power(g2,2))*Sqr(Sin(4*ArcTan(TanBeta))) + 3*Log(msu2(
      2,2)/Sqr(SCALE))*Sqr(Yu(2,2))*(0.4*Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(Yu(2
      ,2))) + 3*Log(msq2(2,2)/Sqr(SCALE))*Sqr(Yu(2,2))*(0.5*Cos(2*ArcTan(TanBeta))
      *(-0.2*Sqr(g1) + Sqr(g2)) + Sqr(Yu(2,2))) + (6*Sqr(AtInput - Mu/TanBeta)*
      Power(Yu(2,2),4)*(TCF(1)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))) - (
      0.08333333333333333*Sqr(AtInput - Mu/TanBeta)*TCF(2)(Sqrt(Abs(msq2(2,2)/msu2
      (2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2))) +
      (0.75*Cos(2*ArcTan(TanBeta))*Sqr(AtInput - Mu/TanBeta)*Sqr(Yu(2,2))*(0.6*Sqr
      (g1)*TCF(3)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))) + Sqr(g2)*TCF(4)(Sqrt(Abs(msq2(2
      ,2)/msu2(2,2))))))/Sqrt(Abs(msq2(2,2)*msu2(2,2))) - (0.25*(0.6*Sqr(g1) + Sqr
      (g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(AtInput - Mu/TanBeta)*Sqr(Yu(2,2))*TCF
      (5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))))));
   MODEL->set_gYu(Re((0.7745966692414834*g1*TanBeta)/Sqrt(1 + Sqr(TanBeta)) + (
      0.004905190710809969*g1*TanBeta*(0.1*(2*Log(mse2(0,0)/Sqr(SCALE)) + 2*Log(
      mse2(1,1)/Sqr(SCALE)) + 2*Log(mse2(2,2)/Sqr(SCALE)) + Log(msl2(0,0)/Sqr(
      SCALE)) + Log(msl2(1,1)/Sqr(SCALE)) + Log(msl2(2,2)/Sqr(SCALE)))*Sqr(g1) +
      0.03333333333333333*(2*Log(msd2(0,0)/Sqr(SCALE)) + 2*Log(msd2(1,1)/Sqr(SCALE
      )) + 2*Log(msd2(2,2)/Sqr(SCALE)) + Log(msq2(0,0)/Sqr(SCALE)) + Log(msq2(1,1)
      /Sqr(SCALE)) + Log(msq2(2,2)/Sqr(SCALE)) + 8*Log(msu2(0,0)/Sqr(SCALE)) + 8*
      Log(msu2(1,1)/Sqr(SCALE)) + 8*Log(msu2(2,2)/Sqr(SCALE)))*Sqr(g1) + 0.0375*
      Sqr(g1)*(-44 + 7/(1 + Sqr(TanBeta))) + 0.1875*Sqr(g2)*(-2 + 7/(1 + Sqr(
      TanBeta))) + 0.025*Log(Sqr(mAInput)/Sqr(SCALE))*(4*Sqr(g1) - (9*(Sqr(g1) + 5
      *Sqr(g2)))/(1 + Sqr(TanBeta))) + (2.25*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (0.25*(7*Log(msq2(2,2)/Sqr(SCALE)) - 13*Log(msu2(2,2)/Sqr(SCALE))
      )*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))*UnitStep(-1 + THRESHOLD))
      /Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_gYd(Re((0.7745966692414834*g1)/Sqrt(1 + Sqr(TanBeta)) + (
      0.004905190710809969*g1*(0.1*(2*Log(mse2(0,0)/Sqr(SCALE)) + 2*Log(mse2(1,1)
      /Sqr(SCALE)) + 2*Log(mse2(2,2)/Sqr(SCALE)) + Log(msl2(0,0)/Sqr(SCALE)) + Log
      (msl2(1,1)/Sqr(SCALE)) + Log(msl2(2,2)/Sqr(SCALE)))*Sqr(g1) +
      0.03333333333333333*(2*Log(msd2(0,0)/Sqr(SCALE)) + 2*Log(msd2(1,1)/Sqr(SCALE
      )) + 2*Log(msd2(2,2)/Sqr(SCALE)) + Log(msq2(0,0)/Sqr(SCALE)) + Log(msq2(1,1)
      /Sqr(SCALE)) + Log(msq2(2,2)/Sqr(SCALE)) + 8*Log(msu2(0,0)/Sqr(SCALE)) + 8*
      Log(msu2(1,1)/Sqr(SCALE)) + 8*Log(msu2(2,2)/Sqr(SCALE)))*Sqr(g1) + 0.0375*
      Sqr(g1)*(-44 + (7*Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.1875*Sqr(g2)*(-2 + (
      7*Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.025*Log(Sqr(mAInput)/Sqr(SCALE))*(4*
      Sqr(g1) - (9*(Sqr(g1) + 5*Sqr(g2))*Sqr(TanBeta))/(1 + Sqr(TanBeta))))*
      UnitStep(-1 + THRESHOLD))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_g2u(Re((g2*TanBeta)/Sqrt(1 + Sqr(TanBeta)) + (
      0.006332573977646111*g2*TanBeta*(0.16666666666666666*(Log(msl2(0,0)/Sqr(
      SCALE)) + Log(msl2(1,1)/Sqr(SCALE)) + Log(msl2(2,2)/Sqr(SCALE)))*Sqr(g2) +
      0.5*(Log(msq2(0,0)/Sqr(SCALE)) + Log(msq2(1,1)/Sqr(SCALE)) + Log(msq2(2,2)
      /Sqr(SCALE)))*Sqr(g2) - Sqr(g2)*(0.6666666666666666 + 0.6875/(1 + Sqr(
      TanBeta))) + 0.0375*Sqr(g1)*(-2 + 7/(1 + Sqr(TanBeta))) +
      0.008333333333333333*Log(Sqr(mAInput)/Sqr(SCALE))*(20*Sqr(g2) + (3*(-9*Sqr(
      g1) + 35*Sqr(g2)))/(1 + Sqr(TanBeta))) + (2.25*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      )))/Sqr(TanBeta) - (0.75*(3*Log(msq2(2,2)/Sqr(SCALE)) - Log(msu2(2,2)/Sqr(
      SCALE)))*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))*UnitStep(-1 +
      THRESHOLD))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_g2d(Re(g2/Sqrt(1 + Sqr(TanBeta)) + (0.006332573977646111*g2*(
      0.16666666666666666*Log(msl2(0,0)/Sqr(SCALE))*Sqr(g2) + 0.16666666666666666*
      Log(msl2(1,1)/Sqr(SCALE))*Sqr(g2) + 0.16666666666666666*Log(msl2(2,2)/Sqr(
      SCALE))*Sqr(g2) + 0.5*Log(msq2(0,0)/Sqr(SCALE))*Sqr(g2) + 0.5*Log(msq2(1,1)
      /Sqr(SCALE))*Sqr(g2) + 0.5*Log(msq2(2,2)/Sqr(SCALE))*Sqr(g2) - Sqr(g2)*(
      0.6666666666666666 + (0.6875*Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.0375*Sqr(
      g1)*(-2 + (7*Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.008333333333333333*Log(
      Sqr(mAInput)/Sqr(SCALE))*(20*Sqr(g2) + (3*(-9*Sqr(g1) + 35*Sqr(g2))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta))))*UnitStep(-1 + THRESHOLD))/Sqrt(1 + Sqr(
      TanBeta))));


   check_non_perturbative();
}

bool SplitMSSM_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto gYd = MODELPARAMETER(gYd);
   const auto g2d = MODELPARAMETER(g2d);
   const auto gYu = MODELPARAMETER(gYu);
   const auto g2u = MODELPARAMETER(g2u);

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("g1", MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("g1");
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("g2", MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("g2");
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("g3", MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("g3");
   }
   if (MaxAbsValue(Lambdax) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("Lambdax", MaxAbsValue(Lambdax), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("Lambdax");
   }
   if (MaxAbsValue(Yu) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("Yu", MaxAbsValue(Yu), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("Yu");
   }
   if (MaxAbsValue(Yd) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("Yd", MaxAbsValue(Yd), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("Yd");
   }
   if (MaxAbsValue(Ye) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("Ye", MaxAbsValue(Ye), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("Ye");
   }
   if (MaxAbsValue(gYd) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("gYd", MaxAbsValue(gYd), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("gYd");
   }
   if (MaxAbsValue(g2d) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("g2d", MaxAbsValue(g2d), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("g2d");
   }
   if (MaxAbsValue(gYu) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("gYu", MaxAbsValue(gYu), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("gYu");
   }
   if (MaxAbsValue(g2u) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter("g2u", MaxAbsValue(g2u), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter("g2u");
   }


   return problem;
}

double SplitMSSM_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double SplitMSSM_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const SplitMSSM_input_parameters& SplitMSSM_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

SplitMSSM<Two_scale>* SplitMSSM_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void SplitMSSM_high_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<SplitMSSM<Two_scale>*>(model_);
}

void SplitMSSM_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void SplitMSSM_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = NULL;
}

void SplitMSSM_high_scale_constraint<Two_scale>::initialize()
{
   assert(model && "SplitMSSM_high_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   initial_scale_guess = MSUSY;

   scale = initial_scale_guess;
}

void SplitMSSM_high_scale_constraint<Two_scale>::update_scale()
{
   assert(model && "SplitMSSM_high_scale_constraint<Two_scale>::"
          "update_scale(): model pointer is zero.");

   const double currentScale = model->get_scale();
   const SplitMSSM_soft_parameters beta_functions(model->calc_beta());

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   scale = MSUSY;


   if (errno == ERANGE) {
#ifdef ENABLE_VERBOSE
      ERROR("SplitMSSM_high_scale_constraint<Two_scale>: Overflow error"
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
