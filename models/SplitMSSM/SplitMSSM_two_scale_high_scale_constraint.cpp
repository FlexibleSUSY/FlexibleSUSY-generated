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

// File generated at Fri 20 Oct 2017 08:38:41

#include "SplitMSSM_two_scale_high_scale_constraint.hpp"
#include "SplitMSSM_two_scale_model.hpp"
#include "SplitMSSM_info.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#include <cmath>
#include <cerrno>
#include <cstring>

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
#define MZPole Electroweak_constants::MZ
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define MODEL model
#define MODELCLASSNAME SplitMSSM<Two_scale>

SplitMSSM_high_scale_constraint<Two_scale>::SplitMSSM_high_scale_constraint(
   SplitMSSM<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void SplitMSSM_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();



   update_scale();

   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto LambdaLoopOrder = INPUTPARAMETER(LambdaLoopOrder);
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
      )) + IF(IsCloseRel(msq2(2,2),msu2(2,2),0.01), (0.005208333333333333*Quad(Yu(
      2,2))*Sqr(g3)*(-18 - Quad(AtInput - Mu/TanBeta)/Abs(msq2(2,2)*msu2(2,2)) + (
      12*Sqr(AtInput - Mu/TanBeta))/Sqrt(Abs(msq2(2,2)*msu2(2,2))) + 6*Log(msq2(2,
      2)/Sqr(SCALE))*(4 + Quad(AtInput - Mu/TanBeta)/Abs(msq2(2,2)*msu2(2,2)) - (
      12*Sqr(AtInput - Mu/TanBeta))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))) - 36*Sqr(Log(
      msq2(2,2)/Sqr(SCALE)))))/Quad(3.141592653589793), (-0.03125*Quad(Yu(2,2))*
      Sqr(g3)*(3 + 4*Log(Sqrt(Abs(msq2(2,2)/msu2(2,2)))) - 4*(1 + 3*Log(Sqrt(Abs(
      msq2(2,2)/msu2(2,2)))))*Log(msq2(2,2)/Sqr(SCALE)) + 8*Sqr(Log(Sqrt(Abs(msq2(
      2,2)/msu2(2,2))))) + (Sqr(AtInput - Mu/TanBeta)*((12*Sqrt(Abs(msq2(2,2)/msu2
      (2,2)))*Log(Sqrt(Abs(msq2(2,2)/msu2(2,2))))*(-1 + 2*Log(msq2(2,2)/Sqr(SCALE)
      )))/(-1 + Abs(msq2(2,2)/msu2(2,2))) - (16*(-2 + Abs(msq2(2,2)/msu2(2,2)))*
      Sqrt(Abs(msq2(2,2)/msu2(2,2)))*Sqr(Log(Sqrt(Abs(msq2(2,2)/msu2(2,2))))))/Sqr
      (-1 + Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2))) + (Quad(
      AtInput - Mu/TanBeta)*((6*Abs(msq2(2,2)/msu2(2,2))*(5 + Abs(msq2(2,2)/msu2(2
      ,2)))*Log(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Cube(-1 + Abs(msq2(2,2)/msu2(2,2)
      )) - (10*Abs(msq2(2,2)/msu2(2,2)))/Sqr(-1 + Abs(msq2(2,2)/msu2(2,2))) + (12*
      Abs(msq2(2,2)/msu2(2,2))*(1 - ((1 + Abs(msq2(2,2)/msu2(2,2)))*Log(Sqrt(Abs(
      msq2(2,2)/msu2(2,2)))))/(-1 + Abs(msq2(2,2)/msu2(2,2))))*Log(msq2(2,2)/Sqr(
      SCALE)))/Sqr(-1 + Abs(msq2(2,2)/msu2(2,2))) + (4*Abs(msq2(2,2)/msu2(2,2))*(
      -5 - 4*Abs(msq2(2,2)/msu2(2,2)) + Sqr(Abs(msq2(2,2)/msu2(2,2))))*Sqr(Log(
      Sqrt(Abs(msq2(2,2)/msu2(2,2))))))/Quad(-1 + Abs(msq2(2,2)/msu2(2,2)))))/Abs(
      msq2(2,2)*msu2(2,2)) + 6*Sqr(Log(msq2(2,2)/Sqr(SCALE)))))/Quad(
      3.141592653589793))*UnitStep(-2 + LambdaLoopOrder) + UnitStep(-1 +
      LambdaLoopOrder)*(0.006332573977646111*(-0.09*Quad(g1) - 0.3*Sqr(g1)*Sqr(g2)
      - Quad(g2)*(0.75 - 0.16666666666666666*Sqr(Cos(2*ArcTan(TanBeta))))) +
      0.006332573977646111*(0.00020833333333333335*Log(Sqr(mAInput)/Sqr(SCALE))*(
      261*Quad(g1) + 1325*Quad(g2) + 630*Sqr(g1)*Sqr(g2) - 4*Cos(4*ArcTan(TanBeta)
      )*(9*Quad(g1) + 175*Quad(g2) + 90*Sqr(g1)*Sqr(g2)) - 9*Cos(8*ArcTan(TanBeta)
      )*Sqr(3*Sqr(g1) + 5*Sqr(g2))) + 0.0033333333333333335*(6*(Log(msd2(0,0)/Sqr(
      SCALE)) + Log(msd2(1,1)/Sqr(SCALE)) + Log(msd2(2,2)/Sqr(SCALE)))*Quad(g1) +
      18*(Log(mse2(0,0)/Sqr(SCALE)) + Log(mse2(1,1)/Sqr(SCALE)) + Log(mse2(2,2)
      /Sqr(SCALE)))*Quad(g1) + 24*(Log(msu2(0,0)/Sqr(SCALE)) + Log(msu2(1,1)/Sqr(
      SCALE)) + Log(msu2(2,2)/Sqr(SCALE)))*Quad(g1) + 3*(Log(msq2(0,0)/Sqr(SCALE))
      + Log(msq2(1,1)/Sqr(SCALE)) + Log(msq2(2,2)/Sqr(SCALE)))*(Quad(g1) + 25*
      Quad(g2)) + (Log(msl2(0,0)/Sqr(SCALE)) + Log(msl2(1,1)/Sqr(SCALE)) + Log(
      msl2(2,2)/Sqr(SCALE)))*(9*Quad(g1) + 25*Quad(g2)))*Sqr(Cos(2*ArcTan(TanBeta)
      )) - 0.1875*Sqr(0.6*Sqr(g1) + Sqr(g2))*Sqr(Sin(4*ArcTan(TanBeta))) + 3*Log(
      msu2(2,2)/Sqr(SCALE))*Sqr(Yu(2,2))*(0.4*Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr
      (Yu(2,2))) + 3*Log(msq2(2,2)/Sqr(SCALE))*Sqr(Yu(2,2))*(0.5*Cos(2*ArcTan(
      TanBeta))*(-0.2*Sqr(g1) + Sqr(g2)) + Sqr(Yu(2,2))) + (6*Quad(Yu(2,2))*Sqr(
      AtInput - Mu/TanBeta)*(TCF(1)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))) - (
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
      )*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))*UnitStep(-1 +
      LambdaLoopOrder))/Sqrt(1 + Sqr(TanBeta))));
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
      UnitStep(-1 + LambdaLoopOrder))/Sqrt(1 + Sqr(TanBeta))));
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
      LambdaLoopOrder))/Sqrt(1 + Sqr(TanBeta))));
   MODEL->set_g2d(Re(g2/Sqrt(1 + Sqr(TanBeta)) + (0.006332573977646111*g2*(
      0.16666666666666666*Log(msl2(0,0)/Sqr(SCALE))*Sqr(g2) + 0.16666666666666666*
      Log(msl2(1,1)/Sqr(SCALE))*Sqr(g2) + 0.16666666666666666*Log(msl2(2,2)/Sqr(
      SCALE))*Sqr(g2) + 0.5*Log(msq2(0,0)/Sqr(SCALE))*Sqr(g2) + 0.5*Log(msq2(1,1)
      /Sqr(SCALE))*Sqr(g2) + 0.5*Log(msq2(2,2)/Sqr(SCALE))*Sqr(g2) - Sqr(g2)*(
      0.6666666666666666 + (0.6875*Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.0375*Sqr(
      g1)*(-2 + (7*Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.008333333333333333*Log(
      Sqr(mAInput)/Sqr(SCALE))*(20*Sqr(g2) + (3*(-9*Sqr(g1) + 35*Sqr(g2))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta))))*UnitStep(-1 + LambdaLoopOrder))/Sqrt(1 + Sqr(
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
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::g3);
   }
   if (MaxAbsValue(Lambdax) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Lambdax, MaxAbsValue(Lambdax), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Lambdax);
   }
   if (MaxAbsValue(Yu(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu0_0, MaxAbsValue(Yu(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu0_0);
   }

   if (MaxAbsValue(Yu(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu0_1, MaxAbsValue(Yu(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu0_1);
   }

   if (MaxAbsValue(Yu(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu0_2, MaxAbsValue(Yu(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu0_2);
   }

   if (MaxAbsValue(Yu(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu1_0, MaxAbsValue(Yu(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu1_0);
   }

   if (MaxAbsValue(Yu(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu1_1, MaxAbsValue(Yu(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu1_1);
   }

   if (MaxAbsValue(Yu(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu1_2, MaxAbsValue(Yu(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu1_2);
   }

   if (MaxAbsValue(Yu(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu2_0, MaxAbsValue(Yu(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu2_0);
   }

   if (MaxAbsValue(Yu(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu2_1, MaxAbsValue(Yu(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu2_1);
   }

   if (MaxAbsValue(Yu(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yu2_2, MaxAbsValue(Yu(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yu2_2);
   }
   if (MaxAbsValue(Yd(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd0_0, MaxAbsValue(Yd(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd0_0);
   }

   if (MaxAbsValue(Yd(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd0_1, MaxAbsValue(Yd(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd0_1);
   }

   if (MaxAbsValue(Yd(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd0_2, MaxAbsValue(Yd(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd0_2);
   }

   if (MaxAbsValue(Yd(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd1_0, MaxAbsValue(Yd(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd1_0);
   }

   if (MaxAbsValue(Yd(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd1_1, MaxAbsValue(Yd(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd1_1);
   }

   if (MaxAbsValue(Yd(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd1_2, MaxAbsValue(Yd(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd1_2);
   }

   if (MaxAbsValue(Yd(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd2_0, MaxAbsValue(Yd(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd2_0);
   }

   if (MaxAbsValue(Yd(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd2_1, MaxAbsValue(Yd(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd2_1);
   }

   if (MaxAbsValue(Yd(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Yd2_2, MaxAbsValue(Yd(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Yd2_2);
   }
   if (MaxAbsValue(Ye(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye0_0, MaxAbsValue(Ye(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye0_0);
   }

   if (MaxAbsValue(Ye(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye0_1, MaxAbsValue(Ye(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye0_1);
   }

   if (MaxAbsValue(Ye(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye0_2, MaxAbsValue(Ye(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye0_2);
   }

   if (MaxAbsValue(Ye(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye1_0, MaxAbsValue(Ye(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye1_0);
   }

   if (MaxAbsValue(Ye(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye1_1, MaxAbsValue(Ye(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye1_1);
   }

   if (MaxAbsValue(Ye(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye1_2, MaxAbsValue(Ye(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye1_2);
   }

   if (MaxAbsValue(Ye(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye2_0, MaxAbsValue(Ye(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye2_0);
   }

   if (MaxAbsValue(Ye(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye2_1, MaxAbsValue(Ye(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye2_1);
   }

   if (MaxAbsValue(Ye(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::Ye2_2, MaxAbsValue(Ye(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::Ye2_2);
   }
   if (MaxAbsValue(gYd) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::gYd, MaxAbsValue(gYd), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::gYd);
   }
   if (MaxAbsValue(g2d) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::g2d, MaxAbsValue(g2d), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::g2d);
   }
   if (MaxAbsValue(gYu) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::gYu, MaxAbsValue(gYu), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::gYu);
   }
   if (MaxAbsValue(g2u) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(SplitMSSM_info::g2u, MaxAbsValue(g2u), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(SplitMSSM_info::g2u);
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

void SplitMSSM_high_scale_constraint<Two_scale>::set_model(Model* model_)
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
   model = nullptr;
}

void SplitMSSM_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   initial_scale_guess = MSUSY;

   scale = initial_scale_guess;
}

void SplitMSSM_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   scale = MSUSY;


}

void SplitMSSM_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("SplitMSSM_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
