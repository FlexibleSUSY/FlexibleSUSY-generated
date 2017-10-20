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

// File generated at Fri 20 Oct 2017 08:36:33

#include "HGTHDMIIMSSMBC_two_scale_high_scale_constraint.hpp"
#include "HGTHDMIIMSSMBC_two_scale_model.hpp"
#include "HGTHDMIIMSSMBC_info.hpp"
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
#define MODELCLASSNAME HGTHDMIIMSSMBC<Two_scale>

HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::HGTHDMIIMSSMBC_high_scale_constraint(
   HGTHDMIIMSSMBC<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();



   update_scale();

   const auto MSUSY = INPUTPARAMETER(MSUSY);
   const auto LambdaLoopOrder = INPUTPARAMETER(LambdaLoopOrder);
   const auto AbInput = INPUTPARAMETER(AbInput);
   const auto AtauInput = INPUTPARAMETER(AtauInput);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto g3 = MODELPARAMETER(g3);
   const auto Mu = MODELPARAMETER(Mu);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_g1d(Re(g2));
   MODEL->set_g1dp(Re(0.7745966692414834*g1));
   MODEL->set_g2u(Re(g2));
   MODEL->set_g2up(Re(0.7745966692414834*g1));
   MODEL->set_Lambda1(Re(0.5*(0.25*(0.6*Sqr(g1) + Sqr(g2)) - (
      0.000053468657576480914*Quad(Mu)*Quad(Yu(2,2))*Sqr(g3)*UnitStep(-2 +
      LambdaLoopOrder))/Quad(MSUSY) + ((-0.0031662869888230555*Quad(Mu)*Quad(Yu(2,
      2)))/Quad(MSUSY) + (0.037995443865876666*Quad(Yd(2,2))*Sqr(AbInput)*(1 - (
      0.08333333333333333*Sqr(AbInput))/Sqr(MSUSY)))/Sqr(MSUSY) + (
      0.012665147955292222*Quad(Ye(2,2))*Sqr(AtauInput)*(1 - (0.08333333333333333*
      Sqr(AtauInput))/Sqr(MSUSY)))/Sqr(MSUSY) + 0.0015831434944115277*(0.6*Sqr(g1)
      + Sqr(g2))*((-3*Sqr(AbInput)*Sqr(Yd(2,2)))/Sqr(MSUSY) - (Sqr(AtauInput)*Sqr
      (Ye(2,2)))/Sqr(MSUSY) + (3*Sqr(Mu)*Sqr(Yu(2,2)))/Sqr(MSUSY)) -
      0.0005277144981371759*(0.6*Sqr(g1) + Sqr(g2))*((3*Sqr(AbInput)*Sqr(Yd(2,2)))
      /Sqr(MSUSY) + (Sqr(AtauInput)*Sqr(Ye(2,2)))/Sqr(MSUSY) + (3*Sqr(Mu)*Sqr(Yu(2
      ,2)))/Sqr(MSUSY)))*UnitStep(-1 + LambdaLoopOrder))));
   MODEL->set_Lambda2(Re(0.5*(0.25*(0.6*Sqr(g1) + Sqr(g2)) +
      0.000641623890917771*((-2*AtInput)/MSUSY + (0.3333333333333333*Cube(AtInput)
      )/Cube(MSUSY) - (0.08333333333333333*Quad(AtInput))/Quad(MSUSY))*Quad(Yu(2,2
      ))*Sqr(g3)*UnitStep(-2 + LambdaLoopOrder) + ((-0.0031662869888230555*Quad(Mu
      )*Quad(Yd(2,2)))/Quad(MSUSY) - (0.0010554289962743518*Quad(Mu)*Quad(Ye(2,2))
      )/Quad(MSUSY) + (0.037995443865876666*Quad(Yu(2,2))*Sqr(AtInput)*(1 - (
      0.08333333333333333*Sqr(AtInput))/Sqr(MSUSY)))/Sqr(MSUSY) -
      0.0015831434944115277*(0.6*Sqr(g1) + Sqr(g2))*((-3*Sqr(Mu)*Sqr(Yd(2,2)))/Sqr
      (MSUSY) - (Sqr(Mu)*Sqr(Ye(2,2)))/Sqr(MSUSY) + (3*Sqr(AtInput)*Sqr(Yu(2,2)))
      /Sqr(MSUSY)) - 0.0005277144981371759*(0.6*Sqr(g1) + Sqr(g2))*((3*Sqr(Mu)*Sqr
      (Yd(2,2)))/Sqr(MSUSY) + (Sqr(Mu)*Sqr(Ye(2,2)))/Sqr(MSUSY) + (3*Sqr(AtInput)*
      Sqr(Yu(2,2)))/Sqr(MSUSY)))*UnitStep(-1 + LambdaLoopOrder))));
   MODEL->set_Lambda3(Re(0.25*(-0.6*Sqr(g1) + Sqr(g2)) + (
      0.00016040597272944275*AtInput*(1 - (0.5*AtInput)/MSUSY)*Quad(Yu(2,2))*Sqr(
      g3)*Sqr(Mu)*UnitStep(-2 + LambdaLoopOrder))/Cube(MSUSY) + ((
      0.0010554289962743518*(3*Quad(Yd(2,2))*(3 - Sqr(AbInput)/Sqr(MSUSY)) + Quad(
      Ye(2,2))*(3 - Sqr(AtauInput)/Sqr(MSUSY)) + 3*Quad(Yu(2,2))*(3 - Sqr(AtInput)
      /Sqr(MSUSY)))*Sqr(Mu))/Sqr(MSUSY) + 0.0031662869888230555*(3*Sqr(
      AbInput/MSUSY + AtInput/MSUSY) - (6*Sqr(Mu))/Sqr(MSUSY) - Sqr(-((AbInput*
      AtInput)/Sqr(MSUSY)) + Sqr(Mu)/Sqr(MSUSY)))*Sqr(Yd(2,2))*Sqr(Yu(2,2)) -
      0.0007915717472057639*(-0.6*Sqr(g1) + Sqr(g2))*(3*(Sqr(AbInput)/Sqr(MSUSY) -
      Sqr(Mu)/Sqr(MSUSY))*Sqr(Yd(2,2)) + (Sqr(AtauInput)/Sqr(MSUSY) - Sqr(Mu)/Sqr
      (MSUSY))*Sqr(Ye(2,2)) + 3*(Sqr(AtInput)/Sqr(MSUSY) - Sqr(Mu)/Sqr(MSUSY))*Sqr
      (Yu(2,2))) - 0.00026385724906858796*(-0.6*Sqr(g1) + Sqr(g2))*(3*(Sqr(AbInput
      )/Sqr(MSUSY) + Sqr(Mu)/Sqr(MSUSY))*Sqr(Yd(2,2)) + (Sqr(AtauInput)/Sqr(MSUSY)
      + Sqr(Mu)/Sqr(MSUSY))*Sqr(Ye(2,2)) + 3*(Sqr(AtInput)/Sqr(MSUSY) + Sqr(Mu)
      /Sqr(MSUSY))*Sqr(Yu(2,2))))*UnitStep(-1 + LambdaLoopOrder)));
   MODEL->set_Lambda4(Re(-0.5*Sqr(g2) + (0.00016040597272944275*AtInput*(1 - (
      0.5*AtInput)/MSUSY)*Quad(Yu(2,2))*Sqr(g3)*Sqr(Mu)*UnitStep(-2 +
      LambdaLoopOrder))/Cube(MSUSY) + ((0.0010554289962743518*(3*Quad(Yd(2,2))*(3
      - Sqr(AbInput)/Sqr(MSUSY)) + Quad(Ye(2,2))*(3 - Sqr(AtauInput)/Sqr(MSUSY)) +
      3*Quad(Yu(2,2))*(3 - Sqr(AtInput)/Sqr(MSUSY)))*Sqr(Mu))/Sqr(MSUSY) -
      0.0031662869888230555*(3*Sqr(AbInput/MSUSY + AtInput/MSUSY) - (6*Sqr(Mu))
      /Sqr(MSUSY) - Sqr(-((AbInput*AtInput)/Sqr(MSUSY)) + Sqr(Mu)/Sqr(MSUSY)))*Sqr
      (Yd(2,2))*Sqr(Yu(2,2)) + 0.0015831434944115277*Sqr(g2)*(3*(Sqr(AbInput)/Sqr(
      MSUSY) - Sqr(Mu)/Sqr(MSUSY))*Sqr(Yd(2,2)) + (Sqr(AtauInput)/Sqr(MSUSY) - Sqr
      (Mu)/Sqr(MSUSY))*Sqr(Ye(2,2)) + 3*(Sqr(AtInput)/Sqr(MSUSY) - Sqr(Mu)/Sqr(
      MSUSY))*Sqr(Yu(2,2))) + 0.0005277144981371759*Sqr(g2)*(3*(Sqr(AbInput)/Sqr(
      MSUSY) + Sqr(Mu)/Sqr(MSUSY))*Sqr(Yd(2,2)) + (Sqr(AtauInput)/Sqr(MSUSY) + Sqr
      (Mu)/Sqr(MSUSY))*Sqr(Ye(2,2)) + 3*(Sqr(AtInput)/Sqr(MSUSY) + Sqr(Mu)/Sqr(
      MSUSY))*Sqr(Yu(2,2))))*UnitStep(-1 + LambdaLoopOrder)));
   MODEL->set_Lambda5(Re((0.00016040597272944275*AtInput*(1 - (0.5*AtInput)
      /MSUSY)*Quad(Yu(2,2))*Sqr(g3)*Sqr(Mu)*UnitStep(-2 + LambdaLoopOrder))/Cube(
      MSUSY) - (0.0010554289962743518*((3*Quad(Yd(2,2))*Sqr(AbInput))/Sqr(MSUSY) +
      (Quad(Ye(2,2))*Sqr(AtauInput))/Sqr(MSUSY) + (3*Quad(Yu(2,2))*Sqr(AtInput))
      /Sqr(MSUSY))*Sqr(Mu)*UnitStep(-1 + LambdaLoopOrder))/Sqr(MSUSY)));
   MODEL->set_Lambda6(Re((0.000053468657576480914*(-1 + AtInput/MSUSY)*Cube(Mu)
      *Quad(Yu(2,2))*Sqr(g3)*UnitStep(-2 + LambdaLoopOrder))/Cube(MSUSY) + ((
      0.0010554289962743518*Mu*((3*AbInput*Quad(Yd(2,2))*(-6 + Sqr(AbInput)/Sqr(
      MSUSY)))/MSUSY + (AtauInput*Quad(Ye(2,2))*(-6 + Sqr(AtauInput)/Sqr(MSUSY)))
      /MSUSY + (3*AtInput*Quad(Yu(2,2))*Sqr(Mu))/Cube(MSUSY)))/MSUSY + (
      0.0007915717472057639*Mu*(0.6*Sqr(g1) + Sqr(g2))*((3*AbInput*Sqr(Yd(2,2)))
      /MSUSY + (AtauInput*Sqr(Ye(2,2)))/MSUSY - (3*AtInput*Sqr(Yu(2,2)))/MSUSY))
      /MSUSY)*UnitStep(-1 + LambdaLoopOrder)));
   MODEL->set_Lambda7(Re((0.00016040597272944275*Mu*Quad(Yu(2,2))*Sqr(g3)*(2 +
      (0.3333333333333333*Cube(AtInput))/Cube(MSUSY) - Sqr(AtInput)/Sqr(MSUSY))*
      UnitStep(-2 + LambdaLoopOrder))/MSUSY + ((0.0010554289962743518*Mu*((3*
      AtInput*Quad(Yu(2,2))*(-6 + Sqr(AtInput)/Sqr(MSUSY)))/MSUSY + (3*AbInput*
      Quad(Yd(2,2))*Sqr(Mu))/Cube(MSUSY) + (AtauInput*Quad(Ye(2,2))*Sqr(Mu))/Cube(
      MSUSY)))/MSUSY - (0.0007915717472057639*Mu*(0.6*Sqr(g1) + Sqr(g2))*((3*
      AbInput*Sqr(Yd(2,2)))/MSUSY + (AtauInput*Sqr(Ye(2,2)))/MSUSY - (3*AtInput*
      Sqr(Yu(2,2)))/MSUSY))/MSUSY)*UnitStep(-1 + LambdaLoopOrder)));


   check_non_perturbative();
}

bool HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Lambda6 = MODELPARAMETER(Lambda6);
   const auto Lambda5 = MODELPARAMETER(Lambda5);
   const auto Lambda7 = MODELPARAMETER(Lambda7);
   const auto Lambda1 = MODELPARAMETER(Lambda1);
   const auto Lambda4 = MODELPARAMETER(Lambda4);
   const auto Lambda3 = MODELPARAMETER(Lambda3);
   const auto Lambda2 = MODELPARAMETER(Lambda2);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);
   const auto g1dp = MODELPARAMETER(g1dp);
   const auto g1d = MODELPARAMETER(g1d);
   const auto g2up = MODELPARAMETER(g2up);
   const auto g2u = MODELPARAMETER(g2u);

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g3);
   }
   if (MaxAbsValue(Lambda6) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda6, MaxAbsValue(Lambda6), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda6);
   }
   if (MaxAbsValue(Lambda5) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda5, MaxAbsValue(Lambda5), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda5);
   }
   if (MaxAbsValue(Lambda7) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda7, MaxAbsValue(Lambda7), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda7);
   }
   if (MaxAbsValue(Lambda1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda1, MaxAbsValue(Lambda1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda1);
   }
   if (MaxAbsValue(Lambda4) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda4, MaxAbsValue(Lambda4), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda4);
   }
   if (MaxAbsValue(Lambda3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda3, MaxAbsValue(Lambda3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda3);
   }
   if (MaxAbsValue(Lambda2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda2, MaxAbsValue(Lambda2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Lambda2);
   }
   if (MaxAbsValue(Yu(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu0_0, MaxAbsValue(Yu(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu0_0);
   }

   if (MaxAbsValue(Yu(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu0_1, MaxAbsValue(Yu(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu0_1);
   }

   if (MaxAbsValue(Yu(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu0_2, MaxAbsValue(Yu(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu0_2);
   }

   if (MaxAbsValue(Yu(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu1_0, MaxAbsValue(Yu(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu1_0);
   }

   if (MaxAbsValue(Yu(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu1_1, MaxAbsValue(Yu(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu1_1);
   }

   if (MaxAbsValue(Yu(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu1_2, MaxAbsValue(Yu(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu1_2);
   }

   if (MaxAbsValue(Yu(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu2_0, MaxAbsValue(Yu(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu2_0);
   }

   if (MaxAbsValue(Yu(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu2_1, MaxAbsValue(Yu(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu2_1);
   }

   if (MaxAbsValue(Yu(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu2_2, MaxAbsValue(Yu(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yu2_2);
   }
   if (MaxAbsValue(Yd(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd0_0, MaxAbsValue(Yd(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd0_0);
   }

   if (MaxAbsValue(Yd(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd0_1, MaxAbsValue(Yd(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd0_1);
   }

   if (MaxAbsValue(Yd(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd0_2, MaxAbsValue(Yd(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd0_2);
   }

   if (MaxAbsValue(Yd(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd1_0, MaxAbsValue(Yd(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd1_0);
   }

   if (MaxAbsValue(Yd(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd1_1, MaxAbsValue(Yd(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd1_1);
   }

   if (MaxAbsValue(Yd(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd1_2, MaxAbsValue(Yd(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd1_2);
   }

   if (MaxAbsValue(Yd(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd2_0, MaxAbsValue(Yd(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd2_0);
   }

   if (MaxAbsValue(Yd(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd2_1, MaxAbsValue(Yd(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd2_1);
   }

   if (MaxAbsValue(Yd(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd2_2, MaxAbsValue(Yd(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Yd2_2);
   }
   if (MaxAbsValue(Ye(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye0_0, MaxAbsValue(Ye(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye0_0);
   }

   if (MaxAbsValue(Ye(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye0_1, MaxAbsValue(Ye(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye0_1);
   }

   if (MaxAbsValue(Ye(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye0_2, MaxAbsValue(Ye(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye0_2);
   }

   if (MaxAbsValue(Ye(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye1_0, MaxAbsValue(Ye(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye1_0);
   }

   if (MaxAbsValue(Ye(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye1_1, MaxAbsValue(Ye(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye1_1);
   }

   if (MaxAbsValue(Ye(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye1_2, MaxAbsValue(Ye(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye1_2);
   }

   if (MaxAbsValue(Ye(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye2_0, MaxAbsValue(Ye(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye2_0);
   }

   if (MaxAbsValue(Ye(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye2_1, MaxAbsValue(Ye(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye2_1);
   }

   if (MaxAbsValue(Ye(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye2_2, MaxAbsValue(Ye(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::Ye2_2);
   }
   if (MaxAbsValue(g1dp) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g1dp, MaxAbsValue(g1dp), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g1dp);
   }
   if (MaxAbsValue(g1d) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g1d, MaxAbsValue(g1d), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g1d);
   }
   if (MaxAbsValue(g2up) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g2up, MaxAbsValue(g2up), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g2up);
   }
   if (MaxAbsValue(g2u) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g2u, MaxAbsValue(g2u), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HGTHDMIIMSSMBC_info::g2u);
   }


   return problem;
}

double HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const HGTHDMIIMSSMBC_input_parameters& HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

HGTHDMIIMSSMBC<Two_scale>* HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<HGTHDMIIMSSMBC<Two_scale>*>(model_);
}

void HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
}

void HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   initial_scale_guess = MSUSY;

   scale = initial_scale_guess;
}

void HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   scale = MSUSY;


}

void HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("HGTHDMIIMSSMBC_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
