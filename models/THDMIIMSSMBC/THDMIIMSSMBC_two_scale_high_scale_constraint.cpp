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

// File generated at Sun 4 Aug 2019 19:02:03

#include "THDMIIMSSMBC_two_scale_high_scale_constraint.hpp"
#include "THDMIIMSSMBC_two_scale_model.hpp"
#include "THDMIIMSSMBC_info.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "error.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#ifdef ENABLE_HIMALAYA
#include "HierarchyCalculator.hpp"
#include "version.hpp"
#endif

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
#define MODELCLASSNAME THDMIIMSSMBC<Two_scale>

#if defined(ENABLE_HIMALAYA) && Himalaya_VERSION_MAJOR >= 2
#define FSHimalayaMh23L [&] () {                                        \
      MODEL->calculate_DRbar_masses();                                  \
                                                                        \
      himalaya::Parameters pars;                                        \
      const auto g1 = MODELPARAMETER(g1); \
      const auto g2 = MODELPARAMETER(g2); \
      const auto g3 = MODELPARAMETER(g3); \
      const auto v2 = MODELPARAMETER(v2); \
      const auto v1 = MODELPARAMETER(v1); \
      const auto Yu = MODELPARAMETER(Yu); \
      const auto Yd = MODELPARAMETER(Yd); \
      const auto Ye = MODELPARAMETER(Ye); \
      const auto Vu = MODELPARAMETER(Vu); \
      const auto Ud = MODELPARAMETER(Ud); \
      const auto Uu = MODELPARAMETER(Uu); \
      const auto Ve = MODELPARAMETER(Ve); \
      const auto Ue = MODELPARAMETER(Ue); \
      const auto MAh = MODELPARAMETER(MAh); \
       \
       \
      pars.scale = MODELPARAMETER(scale); \
      pars.mu = Re(Mu); \
      pars.g1 = Re(g1); \
      pars.g2 = Re(g2); \
      pars.g3 = Re(g3); \
      pars.vd = Re(v1); \
      pars.vu = Re(v2); \
      pars.mq2 = Re(Vu); \
      pars.md2 = Re(Ud); \
      pars.mu2 = Re(Uu); \
      pars.ml2 = Re(Ve); \
      pars.me2 = Re(Ue); \
      pars.Au(2,2) = Re(TrilinearUp); \
      pars.Ad(2,2) = Re(TrilinearDown); \
      pars.Ae(2,2) = Re(TrilinearLepton); \
      pars.Yu = Re(Yu); \
      pars.Yd = Re(Yd); \
      pars.Ye = Re(Ye); \
      pars.M1 = 0; \
      pars.M2 = 0; \
      pars.MG = MGluino; \
      pars.MA = MAh; \
       \
      const double msbar_scheme = 1; \
      const double lambda_3L_eft = 1; \
      const double lambda_3L_uncertainty = 0; \
       \
                                                                        \
      double lambda_3L = 0.;                                            \
                                                                        \
      try {                                                             \
         const bool verbose = false;                                    \
         himalaya::HierarchyCalculator hc(pars, verbose);               \
                                                                        \
         const auto ho = hc.calculateDMh3L(false);                      \
                                                                        \
         lambda_3L =                                                    \
            lambda_3L_eft * (                                           \
               ho.getDLambda(3)                                         \
               + msbar_scheme*ho.getDLambdaDRbarPrimeToMSbarShift(3)    \
               + lambda_3L_uncertainty*ho.getDLambdaUncertainty(3)      \
            );                                                          \
                                                                        \
         VERBOSE_MSG("Himalaya top (hierarchy, Dlambda_3L) = ("         \
                     << ho.getSuitableHierarchy() << ", "               \
                     << lambda_3L <<")");                               \
      } catch (const std::exception& e) {                               \
         model->get_problems().flag_bad_mass(THDMIIMSSMBC_info::hh); \
         WARNING(e.what());                                             \
         VERBOSE_MSG(pars);                                             \
      }                                                                 \
                                                                        \
      return lambda_3L;                                                 \
   }()
#else
#define FSHimalayaMh23L [] () {                                         \
      throw HimalayaError("The 3-loop corrections to lambda "           \
                          "require Himalaya 2.0.0 (or higher), but "    \
                          "FlexibleSUSY has not been configured with "  \
                          "this Himalaya version!");                    \
      return 0.;                                                        \
   }()
#endif

THDMIIMSSMBC_high_scale_constraint<Two_scale>::THDMIIMSSMBC_high_scale_constraint(
   THDMIIMSSMBC<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void THDMIIMSSMBC_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   update_scale();

   const auto MSUSY = INPUTPARAMETER(MSUSY);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto LambdaLoopOrder = INPUTPARAMETER(LambdaLoopOrder);
   const auto AbInput = INPUTPARAMETER(AbInput);
   const auto AtauInput = INPUTPARAMETER(AtauInput);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_Lambda1(Re(0.5*(0.25*(0.6*Sqr(g1) + Sqr(g2)) - (
      0.000053468657576480914*Quad(MuInput)*Quad(Yu(2,2))*Sqr(g3)*UnitStep(-2 +
      LambdaLoopOrder))/Quad(MSUSY) + ((-0.0031662869888230555*Quad(MuInput)*Quad(
      Yu(2,2)))/Quad(MSUSY) + (0.037995443865876666*Quad(Yd(2,2))*Sqr(AbInput)*(1
      - (0.08333333333333333*Sqr(AbInput))/Sqr(MSUSY)))/Sqr(MSUSY) + (
      0.012665147955292222*Quad(Ye(2,2))*Sqr(AtauInput)*(1 - (0.08333333333333333*
      Sqr(AtauInput))/Sqr(MSUSY)))/Sqr(MSUSY) + 0.0015831434944115277*(0.6*Sqr(g1)
      + Sqr(g2))*((-3*Sqr(AbInput)*Sqr(Yd(2,2)))/Sqr(MSUSY) - (Sqr(AtauInput)*Sqr(
      Ye(2,2)))/Sqr(MSUSY) + (3*Sqr(MuInput)*Sqr(Yu(2,2)))/Sqr(MSUSY)) -
      0.0005277144981371759*(0.6*Sqr(g1) + Sqr(g2))*((3*Sqr(AbInput)*Sqr(Yd(2,2)))
      /Sqr(MSUSY) + (Sqr(AtauInput)*Sqr(Ye(2,2)))/Sqr(MSUSY) + (3*Sqr(MuInput)*Sqr
      (Yu(2,2)))/Sqr(MSUSY)))*UnitStep(-1 + LambdaLoopOrder))));
   MODEL->set_Lambda2(Re(0.5*(0.25*(0.6*Sqr(g1) + Sqr(g2)) + 0.000641623890917771*
      ((-2*AtInput)/MSUSY + (0.3333333333333333*Cube(AtInput))/Cube(MSUSY) - (
      0.08333333333333333*Quad(AtInput))/Quad(MSUSY))*Quad(Yu(2,2))*Sqr(g3)*
      UnitStep(-2 + LambdaLoopOrder) + ((-0.0031662869888230555*Quad(MuInput)*Quad
      (Yd(2,2)))/Quad(MSUSY) - (0.0010554289962743518*Quad(MuInput)*Quad(Ye(2,2)))
      /Quad(MSUSY) + (0.037995443865876666*Quad(Yu(2,2))*Sqr(AtInput)*(1 - (
      0.08333333333333333*Sqr(AtInput))/Sqr(MSUSY)))/Sqr(MSUSY) -
      0.0015831434944115277*(0.6*Sqr(g1) + Sqr(g2))*((-3*Sqr(MuInput)*Sqr(Yd(2,2))
      )/Sqr(MSUSY) - (Sqr(MuInput)*Sqr(Ye(2,2)))/Sqr(MSUSY) + (3*Sqr(AtInput)*Sqr(
      Yu(2,2)))/Sqr(MSUSY)) - 0.0005277144981371759*(0.6*Sqr(g1) + Sqr(g2))*((3*
      Sqr(MuInput)*Sqr(Yd(2,2)))/Sqr(MSUSY) + (Sqr(MuInput)*Sqr(Ye(2,2)))/Sqr(
      MSUSY) + (3*Sqr(AtInput)*Sqr(Yu(2,2)))/Sqr(MSUSY)))*UnitStep(-1 +
      LambdaLoopOrder))));
   MODEL->set_Lambda3(Re(0.25*(-0.6*Sqr(g1) + Sqr(g2)) + (0.00008020298636472138*
      AtInput*(1 - (0.5*AtInput)/MSUSY)*Quad(Yu(2,2))*Sqr(g3)*Sqr(MuInput)*
      UnitStep(-2 + LambdaLoopOrder))/Cube(MSUSY) + ((0.0010554289962743518*(3*
      Quad(Yd(2,2))*(3 - Sqr(AbInput)/Sqr(MSUSY)) + Quad(Ye(2,2))*(3 - Sqr(
      AtauInput)/Sqr(MSUSY)) + 3*Quad(Yu(2,2))*(3 - Sqr(AtInput)/Sqr(MSUSY)))*Sqr(
      MuInput))/Sqr(MSUSY) + 0.0031662869888230555*(3*Sqr(AbInput/MSUSY + AtInput/
      MSUSY) - (6*Sqr(MuInput))/Sqr(MSUSY) - Sqr(-((AbInput*AtInput)/Sqr(MSUSY)) +
      Sqr(MuInput)/Sqr(MSUSY)))*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 0.0007915717472057639*
      (-0.6*Sqr(g1) + Sqr(g2))*(3*(Sqr(AbInput)/Sqr(MSUSY) - Sqr(MuInput)/Sqr(
      MSUSY))*Sqr(Yd(2,2)) + (Sqr(AtauInput)/Sqr(MSUSY) - Sqr(MuInput)/Sqr(MSUSY))
      *Sqr(Ye(2,2)) + 3*(Sqr(AtInput)/Sqr(MSUSY) - Sqr(MuInput)/Sqr(MSUSY))*Sqr(Yu
      (2,2))) - 0.00026385724906858796*(-0.6*Sqr(g1) + Sqr(g2))*(3*(Sqr(AbInput)/
      Sqr(MSUSY) + Sqr(MuInput)/Sqr(MSUSY))*Sqr(Yd(2,2)) + (Sqr(AtauInput)/Sqr(
      MSUSY) + Sqr(MuInput)/Sqr(MSUSY))*Sqr(Ye(2,2)) + 3*(Sqr(AtInput)/Sqr(MSUSY)
      + Sqr(MuInput)/Sqr(MSUSY))*Sqr(Yu(2,2))))*UnitStep(-1 + LambdaLoopOrder)));
   MODEL->set_Lambda4(Re(-0.5*Sqr(g2) + (0.00008020298636472138*AtInput*(1 - (0.5*
      AtInput)/MSUSY)*Quad(Yu(2,2))*Sqr(g3)*Sqr(MuInput)*UnitStep(-2 +
      LambdaLoopOrder))/Cube(MSUSY) + ((0.0010554289962743518*(3*Quad(Yd(2,2))*(3
      - Sqr(AbInput)/Sqr(MSUSY)) + Quad(Ye(2,2))*(3 - Sqr(AtauInput)/Sqr(MSUSY)) +
      3*Quad(Yu(2,2))*(3 - Sqr(AtInput)/Sqr(MSUSY)))*Sqr(MuInput))/Sqr(MSUSY) -
      0.0031662869888230555*(3*Sqr(AbInput/MSUSY + AtInput/MSUSY) - (6*Sqr(MuInput
      ))/Sqr(MSUSY) - Sqr(-((AbInput*AtInput)/Sqr(MSUSY)) + Sqr(MuInput)/Sqr(MSUSY
      )))*Sqr(Yd(2,2))*Sqr(Yu(2,2)) + 0.0015831434944115277*Sqr(g2)*(3*(Sqr(
      AbInput)/Sqr(MSUSY) - Sqr(MuInput)/Sqr(MSUSY))*Sqr(Yd(2,2)) + (Sqr(AtauInput
      )/Sqr(MSUSY) - Sqr(MuInput)/Sqr(MSUSY))*Sqr(Ye(2,2)) + 3*(Sqr(AtInput)/Sqr(
      MSUSY) - Sqr(MuInput)/Sqr(MSUSY))*Sqr(Yu(2,2))) + 0.0005277144981371759*Sqr(
      g2)*(3*(Sqr(AbInput)/Sqr(MSUSY) + Sqr(MuInput)/Sqr(MSUSY))*Sqr(Yd(2,2)) + (
      Sqr(AtauInput)/Sqr(MSUSY) + Sqr(MuInput)/Sqr(MSUSY))*Sqr(Ye(2,2)) + 3*(Sqr(
      AtInput)/Sqr(MSUSY) + Sqr(MuInput)/Sqr(MSUSY))*Sqr(Yu(2,2))))*UnitStep(-1 +
      LambdaLoopOrder)));
   MODEL->set_Lambda5(Re((-0.0010554289962743518*((3*Quad(Yd(2,2))*Sqr(AbInput))/
      Sqr(MSUSY) + (Quad(Ye(2,2))*Sqr(AtauInput))/Sqr(MSUSY) + (3*Quad(Yu(2,2))*
      Sqr(AtInput))/Sqr(MSUSY))*Sqr(MuInput)*UnitStep(-1 + LambdaLoopOrder))/Sqr(
      MSUSY)));
   MODEL->set_Lambda6(Re((0.000053468657576480914*(-1 + AtInput/MSUSY)*Cube(
      MuInput)*Quad(Yu(2,2))*Sqr(g3)*UnitStep(-2 + LambdaLoopOrder))/Cube(MSUSY) +
      ((0.0010554289962743518*MuInput*((3*AbInput*Quad(Yd(2,2))*(-6 + Sqr(AbInput)
      /Sqr(MSUSY)))/MSUSY + (AtauInput*Quad(Ye(2,2))*(-6 + Sqr(AtauInput)/Sqr(
      MSUSY)))/MSUSY + (3*AtInput*Quad(Yu(2,2))*Sqr(MuInput))/Cube(MSUSY)))/MSUSY
      + (0.0007915717472057639*MuInput*(0.6*Sqr(g1) + Sqr(g2))*((3*AbInput*Sqr(Yd(
      2,2)))/MSUSY + (AtauInput*Sqr(Ye(2,2)))/MSUSY - (3*AtInput*Sqr(Yu(2,2)))/
      MSUSY))/MSUSY)*UnitStep(-1 + LambdaLoopOrder)));
   MODEL->set_Lambda7(Re((0.00016040597272944275*MuInput*Quad(Yu(2,2))*Sqr(g3)*(2
      + (0.3333333333333333*Cube(AtInput))/Cube(MSUSY) - Sqr(AtInput)/Sqr(MSUSY))*
      UnitStep(-2 + LambdaLoopOrder))/MSUSY + ((0.0010554289962743518*MuInput*((3*
      AtInput*Quad(Yu(2,2))*(-6 + Sqr(AtInput)/Sqr(MSUSY)))/MSUSY + (3*AbInput*
      Quad(Yd(2,2))*Sqr(MuInput))/Cube(MSUSY) + (AtauInput*Quad(Ye(2,2))*Sqr(
      MuInput))/Cube(MSUSY)))/MSUSY - (0.0007915717472057639*MuInput*(0.6*Sqr(g1)
      + Sqr(g2))*((3*AbInput*Sqr(Yd(2,2)))/MSUSY + (AtauInput*Sqr(Ye(2,2)))/MSUSY
      - (3*AtInput*Sqr(Yu(2,2)))/MSUSY))/MSUSY)*UnitStep(-1 + LambdaLoopOrder)));


   check_non_perturbative();
}

bool THDMIIMSSMBC_high_scale_constraint<Two_scale>::check_non_perturbative()
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

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::g3);
   }
   if (MaxAbsValue(Lambda6) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda6, MaxAbsValue(Lambda6), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda6);
   }
   if (MaxAbsValue(Lambda5) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda5, MaxAbsValue(Lambda5), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda5);
   }
   if (MaxAbsValue(Lambda7) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda7, MaxAbsValue(Lambda7), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda7);
   }
   if (MaxAbsValue(Lambda1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda1, MaxAbsValue(Lambda1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda1);
   }
   if (MaxAbsValue(Lambda4) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda4, MaxAbsValue(Lambda4), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda4);
   }
   if (MaxAbsValue(Lambda3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda3, MaxAbsValue(Lambda3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda3);
   }
   if (MaxAbsValue(Lambda2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda2, MaxAbsValue(Lambda2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Lambda2);
   }
   if (MaxAbsValue(Yu(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu0_0, MaxAbsValue(Yu(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu0_0);
   }

   if (MaxAbsValue(Yu(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu0_1, MaxAbsValue(Yu(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu0_1);
   }

   if (MaxAbsValue(Yu(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu0_2, MaxAbsValue(Yu(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu0_2);
   }

   if (MaxAbsValue(Yu(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu1_0, MaxAbsValue(Yu(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu1_0);
   }

   if (MaxAbsValue(Yu(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu1_1, MaxAbsValue(Yu(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu1_1);
   }

   if (MaxAbsValue(Yu(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu1_2, MaxAbsValue(Yu(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu1_2);
   }

   if (MaxAbsValue(Yu(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu2_0, MaxAbsValue(Yu(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu2_0);
   }

   if (MaxAbsValue(Yu(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu2_1, MaxAbsValue(Yu(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu2_1);
   }

   if (MaxAbsValue(Yu(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu2_2, MaxAbsValue(Yu(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yu2_2);
   }
   if (MaxAbsValue(Yd(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd0_0, MaxAbsValue(Yd(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd0_0);
   }

   if (MaxAbsValue(Yd(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd0_1, MaxAbsValue(Yd(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd0_1);
   }

   if (MaxAbsValue(Yd(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd0_2, MaxAbsValue(Yd(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd0_2);
   }

   if (MaxAbsValue(Yd(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd1_0, MaxAbsValue(Yd(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd1_0);
   }

   if (MaxAbsValue(Yd(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd1_1, MaxAbsValue(Yd(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd1_1);
   }

   if (MaxAbsValue(Yd(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd1_2, MaxAbsValue(Yd(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd1_2);
   }

   if (MaxAbsValue(Yd(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd2_0, MaxAbsValue(Yd(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd2_0);
   }

   if (MaxAbsValue(Yd(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd2_1, MaxAbsValue(Yd(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd2_1);
   }

   if (MaxAbsValue(Yd(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd2_2, MaxAbsValue(Yd(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Yd2_2);
   }
   if (MaxAbsValue(Ye(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye0_0, MaxAbsValue(Ye(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye0_0);
   }

   if (MaxAbsValue(Ye(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye0_1, MaxAbsValue(Ye(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye0_1);
   }

   if (MaxAbsValue(Ye(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye0_2, MaxAbsValue(Ye(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye0_2);
   }

   if (MaxAbsValue(Ye(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye1_0, MaxAbsValue(Ye(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye1_0);
   }

   if (MaxAbsValue(Ye(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye1_1, MaxAbsValue(Ye(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye1_1);
   }

   if (MaxAbsValue(Ye(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye1_2, MaxAbsValue(Ye(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye1_2);
   }

   if (MaxAbsValue(Ye(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye2_0, MaxAbsValue(Ye(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye2_0);
   }

   if (MaxAbsValue(Ye(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye2_1, MaxAbsValue(Ye(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye2_1);
   }

   if (MaxAbsValue(Ye(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye2_2, MaxAbsValue(Ye(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(THDMIIMSSMBC_info::Ye2_2);
   }


   return problem;
}

double THDMIIMSSMBC_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double THDMIIMSSMBC_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const THDMIIMSSMBC_input_parameters& THDMIIMSSMBC_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

THDMIIMSSMBC<Two_scale>* THDMIIMSSMBC_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void THDMIIMSSMBC_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<THDMIIMSSMBC<Two_scale>*>(model_);
}

void THDMIIMSSMBC_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void THDMIIMSSMBC_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
}

void THDMIIMSSMBC_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   initial_scale_guess = MSUSY;

   scale = initial_scale_guess;
}

void THDMIIMSSMBC_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);

   scale = MSUSY;


}

void THDMIIMSSMBC_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("THDMIIMSSMBC_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
