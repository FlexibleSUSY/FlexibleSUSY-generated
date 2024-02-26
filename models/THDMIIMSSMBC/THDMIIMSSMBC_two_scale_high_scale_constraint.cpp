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


#include "THDMIIMSSMBC_two_scale_high_scale_constraint.hpp"
#include "THDMIIMSSMBC_two_scale_model.hpp"
#include "THDMIIMSSMBC_info.hpp"
#include "config.h"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "error.hpp"
#include "minimizer.hpp"
#include "raii.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"

#ifdef ENABLE_HIMALAYA
#include "himalaya/HierarchyCalculator.hpp"
#include "himalaya/version.hpp"
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
   const auto M2Input = INPUTPARAMETER(M2Input);
   const auto M1Input = INPUTPARAMETER(M1Input);
   const auto mslInput = INPUTPARAMETER(mslInput);
   const auto msqInput = INPUTPARAMETER(msqInput);
   const auto msdInput = INPUTPARAMETER(msdInput);
   const auto mseInput = INPUTPARAMETER(mseInput);
   const auto msuInput = INPUTPARAMETER(msuInput);
   const auto AdInput = INPUTPARAMETER(AdInput);
   const auto AeInput = INPUTPARAMETER(AeInput);
   const auto AuInput = INPUTPARAMETER(AuInput);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_Lambda1(Re(0.5*(0.25*(0.6*Sqr(g1) + Sqr(g2)) - (
      0.000053468657576480914*Quad(MuInput)*Quad(Yu(2,2))*Sqr(g3)*UnitStep(-2 +
      LambdaLoopOrder))/Quad(MSUSY) + UnitStep(-1 + LambdaLoopOrder)*(0.5*(-
      0.00037995443865876665*(2*Log(Sqr(Abs(MuInput))/Sqr(SCALE)) +
      0.3333333333333333*Log(Sqr(msdInput(0))/Sqr(SCALE)) + 0.3333333333333333*Log
      (Sqr(msdInput(1))/Sqr(SCALE)) + 0.3333333333333333*Log(Sqr(msdInput(2))/Sqr(
      SCALE)) + Log(Sqr(mseInput(0))/Sqr(SCALE)) + Log(Sqr(mseInput(1))/Sqr(SCALE)
      ) + Log(Sqr(mseInput(2))/Sqr(SCALE)) + 0.5*Log(Sqr(mslInput(0))/Sqr(SCALE))
      + 0.5*Log(Sqr(mslInput(1))/Sqr(SCALE)) + 0.5*Log(Sqr(mslInput(2))/Sqr(SCALE)
      ) + 0.16666666666666666*Log(Sqr(msqInput(0))/Sqr(SCALE)) +
      0.16666666666666666*Log(Sqr(msqInput(1))/Sqr(SCALE)) + 0.16666666666666666*
      Log(Sqr(msqInput(2))/Sqr(SCALE)) + 1.3333333333333333*Log(Sqr(msuInput(0))/
      Sqr(SCALE)) + 1.3333333333333333*Log(Sqr(msuInput(1))/Sqr(SCALE)) +
      1.3333333333333333*Log(Sqr(msuInput(2))/Sqr(SCALE)))*Quad(g1) -
      0.0005277144981371759*(-4 + 8*Log(Sqr(M2Input)/Sqr(SCALE)) + 4*Log(Sqr(Abs(
      MuInput))/Sqr(SCALE)) + Log(Sqr(mslInput(0))/Sqr(SCALE)) + Log(Sqr(mslInput(
      1))/Sqr(SCALE)) + Log(Sqr(mslInput(2))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(0))/
      Sqr(SCALE)) + 3*Log(Sqr(msqInput(1))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(2))/
      Sqr(SCALE)))*Quad(g2)) + Re(0.0031662869888230555*(3*Sqr(AdInput(0,0))*Sqr(
      Yd(0,0))*TCDB0(msdInput(0),msqInput(0)) + 3*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*
      TCDB0(msdInput(0),msqInput(1)) + 3*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*TCDB0(
      msdInput(0),msqInput(2)) + 3*Sqr(AdInput(1,0))*Sqr(Yd(1,0))*TCDB0(msdInput(1
      ),msqInput(0)) + 3*Sqr(AdInput(1,1))*Sqr(Yd(1,1))*TCDB0(msdInput(1),msqInput
      (1)) + 3*Sqr(AdInput(1,2))*Sqr(Yd(1,2))*TCDB0(msdInput(1),msqInput(2)) + 3*
      Sqr(AdInput(2,0))*Sqr(Yd(2,0))*TCDB0(msdInput(2),msqInput(0)) + 3*Sqr(
      AdInput(2,1))*Sqr(Yd(2,1))*TCDB0(msdInput(2),msqInput(1)) + 3*Sqr(AdInput(2,
      2))*Sqr(Yd(2,2))*TCDB0(msdInput(2),msqInput(2)) + Sqr(AeInput(0,0))*Sqr(Ye(0
      ,0))*TCDB0(mseInput(0),mslInput(0)) + Sqr(AeInput(0,1))*Sqr(Ye(0,1))*TCDB0(
      mseInput(0),mslInput(1)) + Sqr(AeInput(0,2))*Sqr(Ye(0,2))*TCDB0(mseInput(0),
      mslInput(2)) + Sqr(AeInput(1,0))*Sqr(Ye(1,0))*TCDB0(mseInput(1),mslInput(0))
      + Sqr(AeInput(1,1))*Sqr(Ye(1,1))*TCDB0(mseInput(1),mslInput(1)) + Sqr(
      AeInput(1,2))*Sqr(Ye(1,2))*TCDB0(mseInput(1),mslInput(2)) + Sqr(AeInput(2,0)
      )*Sqr(Ye(2,0))*TCDB0(mseInput(2),mslInput(0)) + Sqr(AeInput(2,1))*Sqr(Ye(2,1
      ))*TCDB0(mseInput(2),mslInput(1)) + Sqr(AeInput(2,2))*Sqr(Ye(2,2))*TCDB0(
      mseInput(2),mslInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*TCDB0(msuInput(0
      ),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*TCDB0(msuInput(0),msqInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*TCDB0(msuInput(0),msqInput(2)) + 3*
      Sqr(Abs(MuInput))*Sqr(Yu(1,0))*TCDB0(msuInput(1),msqInput(0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(1,1))*TCDB0(msuInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yu(1,2))*TCDB0(msuInput(1),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,0
      ))*TCDB0(msuInput(2),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*TCDB0(
      msuInput(2),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,2))*TCDB0(msuInput(2
      ),msqInput(2))) - 0.0007915717472057639*(0.6*Sqr(g1)*TCW(Abs(M1Input),Abs(
      MuInput),SCALE) + 3*Sqr(g2)*TCW(Abs(M2Input),Abs(MuInput),SCALE)))*(0.6*Sqr(
      g1) + Sqr(g2)) + 0.006332573977646111*((-0.03*Quad(g1) + 0.6*Sqr(g1)*(Sqr(Yd
      (0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))) - 3*Sqr(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) +
      Sqr(Yd(0,2))))*TCB0(msdInput(0),msdInput(0),SCALE) - 3*Sqr(Yd(0,0)*Yd(1,0) +
      Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2))*TCB0(msdInput(0),msdInput(1),SCALE) - 3*
      Sqr(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2))*TCB0(msdInput(0),
      msdInput(2),SCALE) - 3*Sqr(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,
      2))*TCB0(msdInput(1),msdInput(0),SCALE) + (-0.03*Quad(g1) + 0.6*Sqr(g1)*(Sqr
      (Yd(1,0)) + Sqr(Yd(1,1)) + Sqr(Yd(1,2))) - 3*Sqr(Sqr(Yd(1,0)) + Sqr(Yd(1,1))
      + Sqr(Yd(1,2))))*TCB0(msdInput(1),msdInput(1),SCALE) - 3*Sqr(Yd(1,0)*Yd(2,0)
      + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2))*TCB0(msdInput(1),msdInput(2),SCALE) - 3
      *Sqr(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2))*TCB0(msdInput(2),
      msdInput(0),SCALE) - 3*Sqr(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,
      2))*TCB0(msdInput(2),msdInput(1),SCALE) + (-0.03*Quad(g1) + 0.6*Sqr(g1)*(Sqr
      (Yd(2,0)) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))) - 3*Sqr(Sqr(Yd(2,0)) + Sqr(Yd(2,1))
      + Sqr(Yd(2,2))))*TCB0(msdInput(2),msdInput(2),SCALE) + (-0.09*Quad(g1) + 0.6
      *Sqr(g1)*(Sqr(Ye(0,0)) + Sqr(Ye(0,1)) + Sqr(Ye(0,2))) - Sqr(Sqr(Ye(0,0)) +
      Sqr(Ye(0,1)) + Sqr(Ye(0,2))))*TCB0(mseInput(0),mseInput(0),SCALE) - Sqr(Ye(0
      ,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye(1,2))*TCB0(mseInput(0),mseInput(1
      ),SCALE) - Sqr(Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2))*TCB0(
      mseInput(0),mseInput(2),SCALE) - Sqr(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(
      0,2)*Ye(1,2))*TCB0(mseInput(1),mseInput(0),SCALE) + (-0.09*Quad(g1) + 0.6*
      Sqr(g1)*(Sqr(Ye(1,0)) + Sqr(Ye(1,1)) + Sqr(Ye(1,2))) - Sqr(Sqr(Ye(1,0)) +
      Sqr(Ye(1,1)) + Sqr(Ye(1,2))))*TCB0(mseInput(1),mseInput(1),SCALE) - Sqr(Ye(1
      ,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2)*Ye(2,2))*TCB0(mseInput(1),mseInput(2
      ),SCALE) - Sqr(Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2))*TCB0(
      mseInput(2),mseInput(0),SCALE) - Sqr(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(
      1,2)*Ye(2,2))*TCB0(mseInput(2),mseInput(1),SCALE) + (-0.09*Quad(g1) + 0.6*
      Sqr(g1)*(Sqr(Ye(2,0)) + Sqr(Ye(2,1)) + Sqr(Ye(2,2))) - Sqr(Sqr(Ye(2,0)) +
      Sqr(Ye(2,1)) + Sqr(Ye(2,2))))*TCB0(mseInput(2),mseInput(2),SCALE) + (0.125*(
      -0.36*Quad(g1) - Quad(g2)) + 0.5*(-0.6*Sqr(g1) + Sqr(g2))*(Sqr(Ye(0,0)) +
      Sqr(Ye(1,0)) + Sqr(Ye(2,0))) - Sqr(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0)
      )))*TCB0(mslInput(0),mslInput(0),SCALE) - Sqr(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1
      ,1) + Ye(2,0)*Ye(2,1))*TCB0(mslInput(0),mslInput(1),SCALE) - Sqr(Ye(0,0)*Ye(
      0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2))*TCB0(mslInput(0),mslInput(2),SCALE
      ) - Sqr(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1))*TCB0(mslInput(1
      ),mslInput(0),SCALE) + (0.125*(-0.36*Quad(g1) - Quad(g2)) + 0.5*(-0.6*Sqr(g1
      ) + Sqr(g2))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))) - Sqr(Sqr(Ye(0,1))
      + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCB0(mslInput(1),mslInput(1),SCALE) - Sqr(Ye
      (0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2))*TCB0(mslInput(1),mslInput
      (2),SCALE) - Sqr(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2))*TCB0(
      mslInput(2),mslInput(0),SCALE) - Sqr(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(
      2,1)*Ye(2,2))*TCB0(mslInput(2),mslInput(1),SCALE) + (0.125*(-0.36*Quad(g1) -
      Quad(g2)) + 0.5*(-0.6*Sqr(g1) + Sqr(g2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(
      Ye(2,2))) - Sqr(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCB0(mslInput(2
      ),mslInput(2),SCALE) + (0.041666666666666664*(-0.36*Quad(g1) - 9*Quad(g2)) +
      0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))) -
      3*Sqr(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))))*TCB0(msqInput(0),msqInput
      (0),SCALE) - 3*Sqr(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*TCB0
      (msqInput(0),msqInput(1),SCALE) - 3*Sqr(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) +
      Yd(2,0)*Yd(2,2))*TCB0(msqInput(0),msqInput(2),SCALE) - 3*Sqr(Yd(0,0)*Yd(0,1)
      + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*TCB0(msqInput(1),msqInput(0),SCALE) + (
      0.041666666666666664*(-0.36*Quad(g1) - 9*Quad(g2)) + 0.5*(0.6*Sqr(g1) + 3*
      Sqr(g2))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) - 3*Sqr(Sqr(Yd(0,1)) +
      Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*TCB0(msqInput(1),msqInput(1),SCALE) - 3*Sqr(Yd
      (0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*TCB0(msqInput(1),msqInput
      (2),SCALE) - 3*Sqr(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*TCB0
      (msqInput(2),msqInput(0),SCALE) - 3*Sqr(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) +
      Yd(2,1)*Yd(2,2))*TCB0(msqInput(2),msqInput(1),SCALE) + (0.041666666666666664
      *(-0.36*Quad(g1) - 9*Quad(g2)) + 0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*(Sqr(Yd(0,2))
      + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) - 3*Sqr(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(
      2,2))))*TCB0(msqInput(2),msqInput(2),SCALE) - 0.12*Quad(g1)*TCB0(msuInput(0)
      ,msuInput(0),SCALE) - 0.12*Quad(g1)*TCB0(msuInput(1),msuInput(1),SCALE) -
      0.12*Quad(g1)*TCB0(msuInput(2),msuInput(2),SCALE) + (0.6*Sqr(g1)*Sqr(AdInput
      (0,0))*Sqr(Yd(0,0)) - 6*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*(Sqr(Yd(0,0)) + Sqr(
      Yd(0,1)) + Sqr(Yd(0,2))))*TCC0(msdInput(0),msdInput(0),msqInput(0)) + (0.6*
      Sqr(g1)*Sqr(AdInput(0,1))*Sqr(Yd(0,1)) - 6*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*(
      Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))*TCC0(msdInput(0),msdInput(0),
      msqInput(1)) + (0.6*Sqr(g1)*Sqr(AdInput(0,2))*Sqr(Yd(0,2)) - 6*Sqr(AdInput(0
      ,2))*Sqr(Yd(0,2))*(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))*TCC0(
      msdInput(0),msdInput(0),msqInput(2)) + (0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(
      AdInput(0,0))*Sqr(Yd(0,0)) - 6*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*(Sqr(Yd(0,0))
      + Sqr(Yd(1,0)) + Sqr(Yd(2,0))))*TCC0(msdInput(0),msqInput(0),msqInput(0)) +
      (0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AdInput(0,1))*Sqr(Yd(0,1)) - 6*Sqr(
      AdInput(0,1))*Sqr(Yd(0,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*
      TCC0(msdInput(0),msqInput(1),msqInput(1)) + (0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*
      Sqr(AdInput(0,2))*Sqr(Yd(0,2)) - 6*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*(Sqr(Yd(0,
      2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))))*TCC0(msdInput(0),msqInput(2),msqInput(2)
      ) + (0.6*Sqr(g1)*Sqr(AdInput(1,0))*Sqr(Yd(1,0)) - 6*Sqr(AdInput(1,0))*Sqr(Yd
      (1,0))*(Sqr(Yd(1,0)) + Sqr(Yd(1,1)) + Sqr(Yd(1,2))))*TCC0(msdInput(1),
      msdInput(1),msqInput(0)) + (0.6*Sqr(g1)*Sqr(AdInput(1,1))*Sqr(Yd(1,1)) - 6*
      Sqr(AdInput(1,1))*Sqr(Yd(1,1))*(Sqr(Yd(1,0)) + Sqr(Yd(1,1)) + Sqr(Yd(1,2))))
      *TCC0(msdInput(1),msdInput(1),msqInput(1)) + (0.6*Sqr(g1)*Sqr(AdInput(1,2))*
      Sqr(Yd(1,2)) - 6*Sqr(AdInput(1,2))*Sqr(Yd(1,2))*(Sqr(Yd(1,0)) + Sqr(Yd(1,1))
      + Sqr(Yd(1,2))))*TCC0(msdInput(1),msdInput(1),msqInput(2)) + (0.5*(0.6*Sqr(
      g1) + 3*Sqr(g2))*Sqr(AdInput(1,0))*Sqr(Yd(1,0)) - 6*Sqr(AdInput(1,0))*Sqr(Yd
      (1,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))))*TCC0(msdInput(1),
      msqInput(0),msqInput(0)) + (0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AdInput(1,1))*
      Sqr(Yd(1,1)) - 6*Sqr(AdInput(1,1))*Sqr(Yd(1,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1))
      + Sqr(Yd(2,1))))*TCC0(msdInput(1),msqInput(1),msqInput(1)) + (0.5*(0.6*Sqr(
      g1) + 3*Sqr(g2))*Sqr(AdInput(1,2))*Sqr(Yd(1,2)) - 6*Sqr(AdInput(1,2))*Sqr(Yd
      (1,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))))*TCC0(msdInput(1),
      msqInput(2),msqInput(2)) + (0.6*Sqr(g1)*Sqr(AdInput(2,0))*Sqr(Yd(2,0)) - 6*
      Sqr(AdInput(2,0))*Sqr(Yd(2,0))*(Sqr(Yd(2,0)) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))))
      *TCC0(msdInput(2),msdInput(2),msqInput(0)) + (0.6*Sqr(g1)*Sqr(AdInput(2,1))*
      Sqr(Yd(2,1)) - 6*Sqr(AdInput(2,1))*Sqr(Yd(2,1))*(Sqr(Yd(2,0)) + Sqr(Yd(2,1))
      + Sqr(Yd(2,2))))*TCC0(msdInput(2),msdInput(2),msqInput(1)) + (0.6*Sqr(g1)*
      Sqr(AdInput(2,2))*Sqr(Yd(2,2)) - 6*Sqr(AdInput(2,2))*Sqr(Yd(2,2))*(Sqr(Yd(2,
      0)) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))))*TCC0(msdInput(2),msdInput(2),msqInput(2)
      ) + (0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AdInput(2,0))*Sqr(Yd(2,0)) - 6*Sqr(
      AdInput(2,0))*Sqr(Yd(2,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))))*
      TCC0(msdInput(2),msqInput(0),msqInput(0)) + (0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*
      Sqr(AdInput(2,1))*Sqr(Yd(2,1)) - 6*Sqr(AdInput(2,1))*Sqr(Yd(2,1))*(Sqr(Yd(0,
      1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*TCC0(msdInput(2),msqInput(1),msqInput(1)
      ) + (0.5*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AdInput(2,2))*Sqr(Yd(2,2)) - 6*Sqr(
      AdInput(2,2))*Sqr(Yd(2,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))))*
      TCC0(msdInput(2),msqInput(2),msqInput(2)) + (0.6*Sqr(g1)*Sqr(AeInput(0,0))*
      Sqr(Ye(0,0)) - 2*Sqr(AeInput(0,0))*Sqr(Ye(0,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0))
      + Sqr(Ye(2,0))))*TCC0(mseInput(0),mseInput(0),mslInput(0)) + (0.6*Sqr(g1)*
      Sqr(AeInput(0,1))*Sqr(Ye(0,1)) - 2*Sqr(AeInput(0,1))*Sqr(Ye(0,1))*(Sqr(Ye(0,
      0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(0),mseInput(0),mslInput(1)
      ) + (0.6*Sqr(g1)*Sqr(AeInput(0,2))*Sqr(Ye(0,2)) - 2*Sqr(AeInput(0,2))*Sqr(Ye
      (0,2))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(0),
      mseInput(0),mslInput(2)) + (0.5*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(AeInput(0,0))*
      Sqr(Ye(0,0)) - 2*Sqr(AeInput(0,0))*Sqr(Ye(0,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0))
      + Sqr(Ye(2,0))))*TCC0(mseInput(0),mslInput(0),mslInput(0)) + (0.5*(-0.6*Sqr(
      g1) + Sqr(g2))*Sqr(AeInput(0,1))*Sqr(Ye(0,1)) - 2*Sqr(AeInput(0,1))*Sqr(Ye(0
      ,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(0),mslInput
      (1),mslInput(1)) + (0.5*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(AeInput(0,2))*Sqr(Ye(0,
      2)) - 2*Sqr(AeInput(0,2))*Sqr(Ye(0,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye
      (2,2))))*TCC0(mseInput(0),mslInput(2),mslInput(2)) + (0.6*Sqr(g1)*Sqr(
      AeInput(1,0))*Sqr(Ye(1,0)) - 2*Sqr(AeInput(1,0))*Sqr(Ye(1,0))*(Sqr(Ye(0,1))
      + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(1),mseInput(1),mslInput(0)) +
      (0.6*Sqr(g1)*Sqr(AeInput(1,1))*Sqr(Ye(1,1)) - 2*Sqr(AeInput(1,1))*Sqr(Ye(1,1
      ))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(1),mseInput(1
      ),mslInput(1)) + (0.6*Sqr(g1)*Sqr(AeInput(1,2))*Sqr(Ye(1,2)) - 2*Sqr(AeInput
      (1,2))*Sqr(Ye(1,2))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(
      mseInput(1),mseInput(1),mslInput(2)) + (0.5*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(
      AeInput(1,0))*Sqr(Ye(1,0)) - 2*Sqr(AeInput(1,0))*Sqr(Ye(1,0))*(Sqr(Ye(0,0))
      + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(1),mslInput(0),mslInput(0)) +
      (0.5*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(AeInput(1,1))*Sqr(Ye(1,1)) - 2*Sqr(AeInput
      (1,1))*Sqr(Ye(1,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(
      mseInput(1),mslInput(1),mslInput(1)) + (0.5*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(
      AeInput(1,2))*Sqr(Ye(1,2)) - 2*Sqr(AeInput(1,2))*Sqr(Ye(1,2))*(Sqr(Ye(0,2))
      + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(1),mslInput(2),mslInput(2)) +
      (0.6*Sqr(g1)*Sqr(AeInput(2,0))*Sqr(Ye(2,0)) - 2*Sqr(AeInput(2,0))*Sqr(Ye(2,0
      ))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(2),mseInput(2
      ),mslInput(0)) + (0.6*Sqr(g1)*Sqr(AeInput(2,1))*Sqr(Ye(2,1)) - 2*Sqr(AeInput
      (2,1))*Sqr(Ye(2,1))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(
      mseInput(2),mseInput(2),mslInput(1)) + (0.6*Sqr(g1)*Sqr(AeInput(2,2))*Sqr(Ye
      (2,2)) - 2*Sqr(AeInput(2,2))*Sqr(Ye(2,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr
      (Ye(2,2))))*TCC0(mseInput(2),mseInput(2),mslInput(2)) + (0.5*(-0.6*Sqr(g1) +
      Sqr(g2))*Sqr(AeInput(2,0))*Sqr(Ye(2,0)) - 2*Sqr(AeInput(2,0))*Sqr(Ye(2,0))*(
      Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(2),mslInput(0),
      mslInput(0)) + (0.5*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(AeInput(2,1))*Sqr(Ye(2,1))
      - 2*Sqr(AeInput(2,1))*Sqr(Ye(2,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1
      ))))*TCC0(mseInput(2),mslInput(1),mslInput(1)) + (0.5*(-0.6*Sqr(g1) + Sqr(g2
      ))*Sqr(AeInput(2,2))*Sqr(Ye(2,2)) - 2*Sqr(AeInput(2,2))*Sqr(Ye(2,2))*(Sqr(Ye
      (0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(2),mslInput(2),mslInput
      (2)) + 0.5*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*TCC0(
      msqInput(0),msqInput(0),msuInput(0)) + 0.5*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs
      (MuInput))*Sqr(Yu(1,0))*TCC0(msqInput(0),msqInput(0),msuInput(1)) + 0.5*(0.6
      *Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*TCC0(msqInput(0),
      msqInput(0),msuInput(2)) - 1.2*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*TCC0(
      msqInput(0),msuInput(0),msuInput(0)) - 1.2*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(
      1,0))*TCC0(msqInput(0),msuInput(1),msuInput(1)) - 1.2*Sqr(g1)*Sqr(Abs(
      MuInput))*Sqr(Yu(2,0))*TCC0(msqInput(0),msuInput(2),msuInput(2)) + 0.5*(0.6*
      Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*TCC0(msqInput(1),
      msqInput(1),msuInput(0)) + 0.5*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*
      Sqr(Yu(1,1))*TCC0(msqInput(1),msqInput(1),msuInput(1)) + 0.5*(0.6*Sqr(g1) -
      3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*TCC0(msqInput(1),msqInput(1),
      msuInput(2)) - 1.2*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*TCC0(msqInput(1),
      msuInput(0),msuInput(0)) - 1.2*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(1,1))*TCC0(
      msqInput(1),msuInput(1),msuInput(1)) - 1.2*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(
      2,1))*TCC0(msqInput(1),msuInput(2),msuInput(2)) + 0.5*(0.6*Sqr(g1) - 3*Sqr(
      g2))*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*TCC0(msqInput(2),msqInput(2),msuInput(0)
      ) + 0.5*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(1,2))*TCC0(
      msqInput(2),msqInput(2),msuInput(1)) + 0.5*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs
      (MuInput))*Sqr(Yu(2,2))*TCC0(msqInput(2),msqInput(2),msuInput(2)) - 1.2*Sqr(
      g1)*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*TCC0(msqInput(2),msuInput(0),msuInput(0))
      - 1.2*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(1,2))*TCC0(msqInput(2),msuInput(1),
      msuInput(1)) - 1.2*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(2,2))*TCC0(msqInput(2),
      msuInput(2),msuInput(2)) - 3*Quad(AdInput(0,0))*Quad(Yd(0,0))*TCD0(msdInput(
      0),msdInput(0),msqInput(0),msqInput(0)) - 3*Sqr(AdInput(0,0))*Sqr(AdInput(0,
      1))*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(0),
      msqInput(1)) - 3*Sqr(AdInput(0,0))*Sqr(AdInput(0,2))*Sqr(Yd(0,0))*Sqr(Yd(0,2
      ))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(2)) - 3*Sqr(AdInput(0,0
      ))*Sqr(AdInput(0,1))*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),
      msqInput(1),msqInput(0)) - 3*Quad(AdInput(0,1))*Quad(Yd(0,1))*TCD0(msdInput(
      0),msdInput(0),msqInput(1),msqInput(1)) - 3*Sqr(AdInput(0,1))*Sqr(AdInput(0,
      2))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(1),
      msqInput(2)) - 3*Sqr(AdInput(0,0))*Sqr(AdInput(0,2))*Sqr(Yd(0,0))*Sqr(Yd(0,2
      ))*TCD0(msdInput(0),msdInput(0),msqInput(2),msqInput(0)) - 3*Sqr(AdInput(0,1
      ))*Sqr(AdInput(0,2))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),
      msqInput(2),msqInput(1)) - 3*Quad(AdInput(0,2))*Quad(Yd(0,2))*TCD0(msdInput(
      0),msdInput(0),msqInput(2),msqInput(2)) - 3*Sqr(AdInput(0,0))*Sqr(AdInput(1,
      0))*Sqr(Yd(0,0))*Sqr(Yd(1,0))*TCD0(msdInput(0),msdInput(1),msqInput(0),
      msqInput(0)) - 3*Sqr(AdInput(0,1))*Sqr(AdInput(1,1))*Sqr(Yd(0,1))*Sqr(Yd(1,1
      ))*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(1)) - 3*Sqr(AdInput(0,2
      ))*Sqr(AdInput(1,2))*Sqr(Yd(0,2))*Sqr(Yd(1,2))*TCD0(msdInput(0),msdInput(1),
      msqInput(2),msqInput(2)) - 3*Sqr(AdInput(0,0))*Sqr(AdInput(2,0))*Sqr(Yd(0,0)
      )*Sqr(Yd(2,0))*TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(0)) - 3*Sqr
      (AdInput(0,1))*Sqr(AdInput(2,1))*Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(1)) - 3*Sqr(AdInput(0,2))*Sqr(AdInput(2,2))
      *Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput
      (2)) - 3*Sqr(AdInput(0,0))*Sqr(AdInput(1,0))*Sqr(Yd(0,0))*Sqr(Yd(1,0))*TCD0(
      msdInput(1),msdInput(0),msqInput(0),msqInput(0)) - 3*Sqr(AdInput(0,1))*Sqr(
      AdInput(1,1))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(0),
      msqInput(1),msqInput(1)) - 3*Sqr(AdInput(0,2))*Sqr(AdInput(1,2))*Sqr(Yd(0,2)
      )*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(2)) - 3*
      Quad(AdInput(1,0))*Quad(Yd(1,0))*TCD0(msdInput(1),msdInput(1),msqInput(0),
      msqInput(0)) - 3*Sqr(AdInput(1,0))*Sqr(AdInput(1,1))*Sqr(Yd(1,0))*Sqr(Yd(1,1
      ))*TCD0(msdInput(1),msdInput(1),msqInput(0),msqInput(1)) - 3*Sqr(AdInput(1,0
      ))*Sqr(AdInput(1,2))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),
      msqInput(0),msqInput(2)) - 3*Sqr(AdInput(1,0))*Sqr(AdInput(1,1))*Sqr(Yd(1,0)
      )*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(0)) - 3*
      Quad(AdInput(1,1))*Quad(Yd(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(1),
      msqInput(1)) - 3*Sqr(AdInput(1,1))*Sqr(AdInput(1,2))*Sqr(Yd(1,1))*Sqr(Yd(1,2
      ))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(2)) - 3*Sqr(AdInput(1,0
      ))*Sqr(AdInput(1,2))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),
      msqInput(2),msqInput(0)) - 3*Sqr(AdInput(1,1))*Sqr(AdInput(1,2))*Sqr(Yd(1,1)
      )*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2),msqInput(1)) - 3*
      Quad(AdInput(1,2))*Quad(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2),
      msqInput(2)) - 3*Sqr(AdInput(1,0))*Sqr(AdInput(2,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0
      ))*TCD0(msdInput(1),msdInput(2),msqInput(0),msqInput(0)) - 3*Sqr(AdInput(1,1
      ))*Sqr(AdInput(2,1))*Sqr(Yd(1,1))*Sqr(Yd(2,1))*TCD0(msdInput(1),msdInput(2),
      msqInput(1),msqInput(1)) - 3*Sqr(AdInput(1,2))*Sqr(AdInput(2,2))*Sqr(Yd(1,2)
      )*Sqr(Yd(2,2))*TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput(2)) - 3*Sqr
      (AdInput(0,0))*Sqr(AdInput(2,0))*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(2),
      msdInput(0),msqInput(0),msqInput(0)) - 3*Sqr(AdInput(0,1))*Sqr(AdInput(2,1))
      *Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput
      (1)) - 3*Sqr(AdInput(0,2))*Sqr(AdInput(2,2))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(0),msqInput(2),msqInput(2)) - 3*Sqr(AdInput(1,0))*Sqr(
      AdInput(2,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(2),msdInput(1),
      msqInput(0),msqInput(0)) - 3*Sqr(AdInput(1,1))*Sqr(AdInput(2,1))*Sqr(Yd(1,1)
      )*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(1)) - 3*Sqr
      (AdInput(1,2))*Sqr(AdInput(2,2))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(msdInput(2),
      msdInput(1),msqInput(2),msqInput(2)) - 3*Quad(AdInput(2,0))*Quad(Yd(2,0))*
      TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(0)) - 3*Sqr(AdInput(2,0))*
      Sqr(AdInput(2,1))*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(2),
      msqInput(0),msqInput(1)) - 3*Sqr(AdInput(2,0))*Sqr(AdInput(2,2))*Sqr(Yd(2,0)
      )*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(2)) - 3*Sqr
      (AdInput(2,0))*Sqr(AdInput(2,1))*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(2),msqInput(1),msqInput(0)) - 3*Quad(AdInput(2,1))*Quad(Yd(2,1))*
      TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(1)) - 3*Sqr(AdInput(2,1))*
      Sqr(AdInput(2,2))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(1),msqInput(2)) - 3*Sqr(AdInput(2,0))*Sqr(AdInput(2,2))*Sqr(Yd(2,0)
      )*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(0)) - 3*Sqr
      (AdInput(2,1))*Sqr(AdInput(2,2))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),
      msdInput(2),msqInput(2),msqInput(1)) - 3*Quad(AdInput(2,2))*Quad(Yd(2,2))*
      TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(2)) - Quad(AeInput(0,0))*
      Quad(Ye(0,0))*TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(0)) - Sqr(
      AeInput(0,0))*Sqr(AeInput(0,1))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(mseInput(0),
      mseInput(0),mslInput(0),mslInput(1)) - Sqr(AeInput(0,0))*Sqr(AeInput(0,2))*
      Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(
      2)) - Sqr(AeInput(0,0))*Sqr(AeInput(0,1))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(
      mseInput(0),mseInput(0),mslInput(1),mslInput(0)) - Quad(AeInput(0,1))*Quad(
      Ye(0,1))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(1)) - Sqr(AeInput
      (0,1))*Sqr(AeInput(0,2))*Sqr(Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput
      (0),mslInput(1),mslInput(2)) - Sqr(AeInput(0,0))*Sqr(AeInput(0,2))*Sqr(Ye(0,
      0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(0)) - Sqr
      (AeInput(0,1))*Sqr(AeInput(0,2))*Sqr(Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),
      mseInput(0),mslInput(2),mslInput(1)) - Quad(AeInput(0,2))*Quad(Ye(0,2))*TCD0
      (mseInput(0),mseInput(0),mslInput(2),mslInput(2)) - Sqr(AeInput(0,0))*Sqr(
      AeInput(1,0))*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(0),mseInput(1),
      mslInput(0),mslInput(0)) - Sqr(AeInput(0,1))*Sqr(AeInput(1,1))*Sqr(Ye(0,1))*
      Sqr(Ye(1,1))*TCD0(mseInput(0),mseInput(1),mslInput(1),mslInput(1)) - Sqr(
      AeInput(0,2))*Sqr(AeInput(1,2))*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(0),
      mseInput(1),mslInput(2),mslInput(2)) - Sqr(AeInput(0,0))*Sqr(AeInput(2,0))*
      Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(
      0)) - Sqr(AeInput(0,1))*Sqr(AeInput(2,1))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(1)) - Sqr(AeInput(0,2))*Sqr(
      AeInput(2,2))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(2)) - Sqr(AeInput(0,0))*Sqr(AeInput(1,0))*Sqr(Ye(0,0))*
      Sqr(Ye(1,0))*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(0)) - Sqr(
      AeInput(0,1))*Sqr(AeInput(1,1))*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0(mseInput(1),
      mseInput(0),mslInput(1),mslInput(1)) - Sqr(AeInput(0,2))*Sqr(AeInput(1,2))*
      Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(0),mslInput(2),mslInput(
      2)) - Quad(AeInput(1,0))*Quad(Ye(1,0))*TCD0(mseInput(1),mseInput(1),mslInput
      (0),mslInput(0)) - Sqr(AeInput(1,0))*Sqr(AeInput(1,1))*Sqr(Ye(1,0))*Sqr(Ye(1
      ,1))*TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(1)) - Sqr(AeInput(1,0
      ))*Sqr(AeInput(1,2))*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(0),mslInput(2)) - Sqr(AeInput(1,0))*Sqr(AeInput(1,1))*Sqr(Ye(1,0))*
      Sqr(Ye(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(0)) - Quad(
      AeInput(1,1))*Quad(Ye(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1),
      mslInput(1)) - Sqr(AeInput(1,1))*Sqr(AeInput(1,2))*Sqr(Ye(1,1))*Sqr(Ye(1,2))
      *TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(2)) - Sqr(AeInput(1,0))*
      Sqr(AeInput(1,2))*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(2),mslInput(0)) - Sqr(AeInput(1,1))*Sqr(AeInput(1,2))*Sqr(Ye(1,1))*
      Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(2),mslInput(1)) - Quad(
      AeInput(1,2))*Quad(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(2),
      mslInput(2)) - Sqr(AeInput(1,0))*Sqr(AeInput(2,0))*Sqr(Ye(1,0))*Sqr(Ye(2,0))
      *TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(0)) - Sqr(AeInput(1,1))*
      Sqr(AeInput(2,1))*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(mseInput(1),mseInput(2),
      mslInput(1),mslInput(1)) - Sqr(AeInput(1,2))*Sqr(AeInput(2,2))*Sqr(Ye(1,2))*
      Sqr(Ye(2,2))*TCD0(mseInput(1),mseInput(2),mslInput(2),mslInput(2)) - Sqr(
      AeInput(0,0))*Sqr(AeInput(2,0))*Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),
      mseInput(0),mslInput(0),mslInput(0)) - Sqr(AeInput(0,1))*Sqr(AeInput(2,1))*
      Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(0),mslInput(1),mslInput(
      1)) - Sqr(AeInput(0,2))*Sqr(AeInput(2,2))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(2),mseInput(0),mslInput(2),mslInput(2)) - Sqr(AeInput(1,0))*Sqr(
      AeInput(2,0))*Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(1),
      mslInput(0),mslInput(0)) - Sqr(AeInput(1,1))*Sqr(AeInput(2,1))*Sqr(Ye(1,1))*
      Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(1),mslInput(1),mslInput(1)) - Sqr(
      AeInput(1,2))*Sqr(AeInput(2,2))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(mseInput(2),
      mseInput(1),mslInput(2),mslInput(2)) - Quad(AeInput(2,0))*Quad(Ye(2,0))*TCD0
      (mseInput(2),mseInput(2),mslInput(0),mslInput(0)) - Sqr(AeInput(2,0))*Sqr(
      AeInput(2,1))*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),
      mslInput(0),mslInput(1)) - Sqr(AeInput(2,0))*Sqr(AeInput(2,2))*Sqr(Ye(2,0))*
      Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(0),mslInput(2)) - Sqr(
      AeInput(2,0))*Sqr(AeInput(2,1))*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),
      mseInput(2),mslInput(1),mslInput(0)) - Quad(AeInput(2,1))*Quad(Ye(2,1))*TCD0
      (mseInput(2),mseInput(2),mslInput(1),mslInput(1)) - Sqr(AeInput(2,1))*Sqr(
      AeInput(2,2))*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),
      mslInput(1),mslInput(2)) - Sqr(AeInput(2,0))*Sqr(AeInput(2,2))*Sqr(Ye(2,0))*
      Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(0)) - Sqr(
      AeInput(2,1))*Sqr(AeInput(2,2))*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),
      mseInput(2),mslInput(2),mslInput(1)) - Quad(AeInput(2,2))*Quad(Ye(2,2))*TCD0
      (mseInput(2),mseInput(2),mslInput(2),mslInput(2)) - 3*Quad(Abs(MuInput))*
      Quad(Yu(0,0))*TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(0)) - 3*Quad
      (Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),
      msuInput(0),msuInput(1)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(2,0))*
      TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(2)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),msuInput
      (0)) - 3*Quad(Abs(MuInput))*Quad(Yu(1,0))*TCD0(msqInput(0),msqInput(0),
      msuInput(1),msuInput(1)) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*
      TCD0(msqInput(0),msqInput(0),msuInput(1),msuInput(2)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(2),msuInput
      (0)) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(2),msuInput(1)) - 3*Quad(Abs(MuInput))*Quad(Yu(2,0))*
      TCD0(msqInput(0),msqInput(0),msuInput(2),msuInput(2)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput
      (0)) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(msqInput(0),
      msqInput(1),msuInput(1),msuInput(1)) - 3*Quad(Abs(MuInput))*Sqr(Yu(2,0))*Sqr
      (Yu(2,1))*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(2)) - 3*Quad(Abs
      (MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(0),msqInput(2),msuInput(0
      ),msuInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(
      msqInput(0),msqInput(2),msuInput(1),msuInput(1)) - 3*Quad(Abs(MuInput))*Sqr(
      Yu(2,0))*Sqr(Yu(2,2))*TCD0(msqInput(0),msqInput(2),msuInput(2),msuInput(2))
      - 3*Quad(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(1),msqInput(0
      ),msuInput(0),msuInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*
      TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(1)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(0),msuInput(2),msuInput
      (2)) - 3*Quad(Abs(MuInput))*Quad(Yu(0,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(0),msuInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(1,1))*
      TCD0(msqInput(1),msqInput(1),msuInput(0),msuInput(1)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(0),msuInput
      (2)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(1,1))*TCD0(msqInput(1),
      msqInput(1),msuInput(1),msuInput(0)) - 3*Quad(Abs(MuInput))*Quad(Yu(1,1))*
      TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(1)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(1,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput
      (2)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),
      msqInput(1),msuInput(2),msuInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,1))*Sqr
      (Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(2),msuInput(1)) - 3*Quad(Abs
      (MuInput))*Quad(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(2),msuInput(2
      )) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(1),
      msqInput(2),msuInput(0),msuInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,1))*Sqr
      (Yu(1,2))*TCD0(msqInput(1),msqInput(2),msuInput(1),msuInput(1)) - 3*Quad(Abs
      (MuInput))*Sqr(Yu(2,1))*Sqr(Yu(2,2))*TCD0(msqInput(1),msqInput(2),msuInput(2
      ),msuInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(
      msqInput(2),msqInput(0),msuInput(0),msuInput(0)) - 3*Quad(Abs(MuInput))*Sqr(
      Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(0),msuInput(1),msuInput(1))
      - 3*Quad(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(0
      ),msuInput(2),msuInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*
      TCD0(msqInput(2),msqInput(1),msuInput(0),msuInput(0)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(1),msuInput(1),msuInput
      (1)) - 3*Quad(Abs(MuInput))*Sqr(Yu(2,1))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(1),msuInput(2),msuInput(2)) - 3*Quad(Abs(MuInput))*Quad(Yu(0,2))*
      TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput(0)) - 3*Quad(Abs(MuInput))
      *Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput
      (1)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(0),msuInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,2))*Sqr
      (Yu(1,2))*TCD0(msqInput(2),msqInput(2),msuInput(1),msuInput(0)) - 3*Quad(Abs
      (MuInput))*Quad(Yu(1,2))*TCD0(msqInput(2),msqInput(2),msuInput(1),msuInput(1
      )) - 3*Quad(Abs(MuInput))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(1),msuInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yu(0,2))*Sqr
      (Yu(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(0)) - 3*Quad(Abs
      (MuInput))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2
      ),msuInput(1)) - 3*Quad(Abs(MuInput))*Quad(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(2),msuInput(2)) + 0.36*Quad(g1)*(-0.25 + 0.5*Sqr(Abs(
      M1Input))*TCD2t(M1Input,M1Input,Abs(MuInput),Abs(MuInput)) + 0.5*TCD4t(
      M1Input,M1Input,Abs(MuInput),Abs(MuInput),SCALE)) + 0.6*Sqr(g1)*Sqr(g2)*(-
      0.5 + M1Input*M2Input*TCD2t(M1Input,M2Input,Abs(MuInput),Abs(MuInput)) +
      TCD4t(M1Input,M2Input,Abs(MuInput),Abs(MuInput),SCALE)) + Quad(g2)*(-0.75 +
      0.5*Sqr(Abs(M2Input))*TCD2t(M2Input,M2Input,Abs(MuInput),Abs(MuInput)) + 2.5
      *TCD4t(M2Input,M2Input,Abs(MuInput),Abs(MuInput),SCALE)) - 3*AdInput(0,0)*
      AdInput(0,1)*AdInput(1,0)*AdInput(1,1)*TCD0(msdInput(0),msdInput(1),msqInput
      (0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*AdInput(0,0)*AdInput(0,
      1)*AdInput(1,0)*AdInput(1,1)*TCD0(msdInput(0),msdInput(1),msqInput(1),
      msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*AdInput(0,0)*AdInput(0,1)*
      AdInput(1,0)*AdInput(1,1)*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(
      1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*AdInput(0,0)*AdInput(0,1)*AdInput(1,
      0)*AdInput(1,1)*TCD0(msdInput(1),msdInput(0),msqInput(1),msqInput(0))*Yd(0,0
      )*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*AdInput(0,0)*AdInput(0,2)*AdInput(1,0)*AdInput
      (1,2)*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*
      Yd(1,0)*Yd(1,2) - 3*AdInput(0,0)*AdInput(0,2)*AdInput(1,0)*AdInput(1,2)*TCD0
      (msdInput(0),msdInput(1),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd
      (1,2) - 3*AdInput(0,0)*AdInput(0,2)*AdInput(1,0)*AdInput(1,2)*TCD0(msdInput(
      1),msdInput(0),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*
      AdInput(0,0)*AdInput(0,2)*AdInput(1,0)*AdInput(1,2)*TCD0(msdInput(1),
      msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(1,1)*AdInput(1,2)*TCD0(msdInput(0),
      msdInput(1),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(1,1)*AdInput(1,2)*TCD0(msdInput(0),
      msdInput(1),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(1,1)*AdInput(1,2)*TCD0(msdInput(1),
      msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(1,1)*AdInput(1,2)*TCD0(msdInput(1),
      msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 6*
      AdInput(0,0)*AdInput(1,0)*TCC0(msdInput(0),msdInput(1),msqInput(0))*Yd(0,0)*
      Yd(1,0)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 6*AdInput(0,
      0)*AdInput(1,0)*TCC0(msdInput(1),msdInput(0),msqInput(0))*Yd(0,0)*Yd(1,0)*(
      Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 6*AdInput(0,1)*
      AdInput(1,1)*TCC0(msdInput(0),msdInput(1),msqInput(1))*Yd(0,1)*Yd(1,1)*(Yd(0
      ,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 6*AdInput(0,1)*AdInput(1,
      1)*TCC0(msdInput(1),msdInput(0),msqInput(1))*Yd(0,1)*Yd(1,1)*(Yd(0,0)*Yd(1,0
      ) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 6*AdInput(0,2)*AdInput(1,2)*TCC0(
      msdInput(0),msdInput(1),msqInput(2))*Yd(0,2)*Yd(1,2)*(Yd(0,0)*Yd(1,0) + Yd(0
      ,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 6*AdInput(0,2)*AdInput(1,2)*TCC0(msdInput(1
      ),msdInput(0),msqInput(2))*Yd(0,2)*Yd(1,2)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1
      ) + Yd(0,2)*Yd(1,2)) - 3*AdInput(0,0)*AdInput(0,1)*AdInput(2,0)*AdInput(2,1)
      *TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,
      0)*Yd(2,1) - 3*AdInput(0,0)*AdInput(0,1)*AdInput(2,0)*AdInput(2,1)*TCD0(
      msdInput(0),msdInput(2),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(
      2,1) - 3*AdInput(0,0)*AdInput(0,1)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput(2
      ),msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*
      AdInput(0,0)*AdInput(0,1)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput(2),
      msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*
      AdInput(1,0)*AdInput(1,1)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput(1),
      msdInput(2),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*
      AdInput(1,0)*AdInput(1,1)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput(1),
      msdInput(2),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*
      AdInput(1,0)*AdInput(1,1)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput(2),
      msdInput(1),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*
      AdInput(1,0)*AdInput(1,1)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput(2),
      msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 6*
      AdInput(0,0)*AdInput(0,1)*TCC0(msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*
      Yd(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 6*AdInput(0,
      0)*AdInput(0,1)*TCC0(msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*(
      Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 6*AdInput(1,0)*
      AdInput(1,1)*TCC0(msdInput(1),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*(Yd(0
      ,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 6*AdInput(1,0)*AdInput(1,
      1)*TCC0(msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,1
      ) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 6*AdInput(2,0)*AdInput(2,1)*TCC0(
      msdInput(2),msqInput(0),msqInput(1))*Yd(2,0)*Yd(2,1)*(Yd(0,0)*Yd(0,1) + Yd(1
      ,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 6*AdInput(2,0)*AdInput(2,1)*TCC0(msdInput(2
      ),msqInput(1),msqInput(0))*Yd(2,0)*Yd(2,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1
      ) + Yd(2,0)*Yd(2,1)) - 3*AdInput(0,0)*AdInput(0,2)*AdInput(2,0)*AdInput(2,2)
      *TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,
      0)*Yd(2,2) - 3*AdInput(0,0)*AdInput(0,2)*AdInput(2,0)*AdInput(2,2)*TCD0(
      msdInput(0),msdInput(2),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(
      2,2) - 3*AdInput(0,0)*AdInput(0,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2
      ),msdInput(0),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(0,0)*AdInput(0,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2),
      msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(1,0)*AdInput(1,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(1),
      msdInput(2),msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(1,0)*AdInput(1,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(1,0)*AdInput(1,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2),
      msdInput(1),msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(1,0)*AdInput(1,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2),
      msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(0),
      msdInput(2),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(2),
      msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(0,1)*AdInput(0,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(2),
      msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(1,1)*AdInput(1,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(1),
      msdInput(2),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(1,1)*AdInput(1,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(1,1)*AdInput(1,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(2),
      msdInput(1),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(1,1)*AdInput(1,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(2),
      msdInput(1),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 6*
      AdInput(0,0)*AdInput(2,0)*TCC0(msdInput(0),msdInput(2),msqInput(0))*Yd(0,0)*
      Yd(2,0)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 6*AdInput(0,
      0)*AdInput(2,0)*TCC0(msdInput(2),msdInput(0),msqInput(0))*Yd(0,0)*Yd(2,0)*(
      Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 6*AdInput(0,1)*
      AdInput(2,1)*TCC0(msdInput(0),msdInput(2),msqInput(1))*Yd(0,1)*Yd(2,1)*(Yd(0
      ,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 6*AdInput(0,1)*AdInput(2,
      1)*TCC0(msdInput(2),msdInput(0),msqInput(1))*Yd(0,1)*Yd(2,1)*(Yd(0,0)*Yd(2,0
      ) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 6*AdInput(0,2)*AdInput(2,2)*TCC0(
      msdInput(0),msdInput(2),msqInput(2))*Yd(0,2)*Yd(2,2)*(Yd(0,0)*Yd(2,0) + Yd(0
      ,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 6*AdInput(0,2)*AdInput(2,2)*TCC0(msdInput(2
      ),msdInput(0),msqInput(2))*Yd(0,2)*Yd(2,2)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1
      ) + Yd(0,2)*Yd(2,2)) - 6*AdInput(1,0)*AdInput(2,0)*TCC0(msdInput(1),msdInput
      (2),msqInput(0))*Yd(1,0)*Yd(2,0)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2
      )*Yd(2,2)) - 6*AdInput(1,0)*AdInput(2,0)*TCC0(msdInput(2),msdInput(1),
      msqInput(0))*Yd(1,0)*Yd(2,0)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd
      (2,2)) - 6*AdInput(1,1)*AdInput(2,1)*TCC0(msdInput(1),msdInput(2),msqInput(1
      ))*Yd(1,1)*Yd(2,1)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 6
      *AdInput(1,1)*AdInput(2,1)*TCC0(msdInput(2),msdInput(1),msqInput(1))*Yd(1,1)
      *Yd(2,1)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 6*AdInput(1
      ,2)*AdInput(2,2)*TCC0(msdInput(1),msdInput(2),msqInput(2))*Yd(1,2)*Yd(2,2)*(
      Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 6*AdInput(1,2)*
      AdInput(2,2)*TCC0(msdInput(2),msdInput(1),msqInput(2))*Yd(1,2)*Yd(2,2)*(Yd(1
      ,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 6*AdInput(0,0)*AdInput(0,
      2)*TCC0(msdInput(0),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2
      ) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 6*AdInput(0,0)*AdInput(0,2)*TCC0(
      msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1
      ,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 6*AdInput(1,0)*AdInput(1,2)*TCC0(msdInput(1
      ),msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2
      ) + Yd(2,0)*Yd(2,2)) - 6*AdInput(1,0)*AdInput(1,2)*TCC0(msdInput(1),msqInput
      (2),msqInput(0))*Yd(1,0)*Yd(1,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0
      )*Yd(2,2)) - 6*AdInput(2,0)*AdInput(2,2)*TCC0(msdInput(2),msqInput(0),
      msqInput(2))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd
      (2,2)) - 6*AdInput(2,0)*AdInput(2,2)*TCC0(msdInput(2),msqInput(2),msqInput(0
      ))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 6
      *AdInput(0,1)*AdInput(0,2)*TCC0(msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)
      *Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 6*AdInput(0
      ,1)*AdInput(0,2)*TCC0(msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*(
      Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 6*AdInput(1,1)*
      AdInput(1,2)*TCC0(msdInput(1),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*(Yd(0
      ,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 6*AdInput(1,1)*AdInput(1,
      2)*TCC0(msdInput(1),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*(Yd(0,1)*Yd(0,2
      ) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 6*AdInput(2,1)*AdInput(2,2)*TCC0(
      msdInput(2),msqInput(1),msqInput(2))*Yd(2,1)*Yd(2,2)*(Yd(0,1)*Yd(0,2) + Yd(1
      ,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 6*AdInput(2,1)*AdInput(2,2)*TCC0(msdInput(2
      ),msqInput(2),msqInput(1))*Yd(2,1)*Yd(2,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2
      ) + Yd(2,1)*Yd(2,2)) - AeInput(0,0)*AeInput(0,1)*AeInput(1,0)*AeInput(1,1)*
      TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1,0
      )*Ye(1,1) - AeInput(0,0)*AeInput(0,1)*AeInput(1,0)*AeInput(1,1)*TCD0(
      mseInput(0),mseInput(1),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(
      1,1) - AeInput(0,0)*AeInput(0,1)*AeInput(1,0)*AeInput(1,1)*TCD0(mseInput(1),
      mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) -
      AeInput(0,0)*AeInput(0,1)*AeInput(1,0)*AeInput(1,1)*TCD0(mseInput(1),
      mseInput(0),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) -
      AeInput(0,0)*AeInput(0,2)*AeInput(1,0)*AeInput(1,2)*TCD0(mseInput(0),
      mseInput(1),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) -
      AeInput(0,0)*AeInput(0,2)*AeInput(1,0)*AeInput(1,2)*TCD0(mseInput(0),
      mseInput(1),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) -
      AeInput(0,0)*AeInput(0,2)*AeInput(1,0)*AeInput(1,2)*TCD0(mseInput(1),
      mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) -
      AeInput(0,0)*AeInput(0,2)*AeInput(1,0)*AeInput(1,2)*TCD0(mseInput(1),
      mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(1,1)*AeInput(1,2)*TCD0(mseInput(0),
      mseInput(1),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(1,1)*AeInput(1,2)*TCD0(mseInput(0),
      mseInput(1),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(1,1)*AeInput(1,2)*TCD0(mseInput(1),
      mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(1,1)*AeInput(1,2)*TCD0(mseInput(1),
      mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) -
      AeInput(0,0)*AeInput(0,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(0),
      mseInput(2),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) -
      AeInput(0,0)*AeInput(0,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(0),
      mseInput(2),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) -
      AeInput(0,0)*AeInput(0,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(2),
      mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) -
      AeInput(0,0)*AeInput(0,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(2),
      mseInput(0),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) -
      AeInput(1,0)*AeInput(1,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(1),
      mseInput(2),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) -
      AeInput(1,0)*AeInput(1,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(1),
      mseInput(2),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) -
      AeInput(1,0)*AeInput(1,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(2),
      mseInput(1),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) -
      AeInput(1,0)*AeInput(1,1)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(2),
      mseInput(1),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) - 2*
      AeInput(0,0)*AeInput(0,1)*TCC0(mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*
      Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2*AeInput(0,
      0)*AeInput(0,1)*TCC0(mseInput(0),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*(
      Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2*AeInput(0,0)*
      AeInput(1,0)*TCC0(mseInput(0),mseInput(1),mslInput(0))*Ye(0,0)*Ye(1,0)*(Ye(0
      ,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2*AeInput(0,0)*AeInput(1,
      0)*TCC0(mseInput(1),mseInput(0),mslInput(0))*Ye(0,0)*Ye(1,0)*(Ye(0,0)*Ye(0,1
      ) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2*AeInput(0,1)*AeInput(1,1)*TCC0(
      mseInput(0),mseInput(1),mslInput(1))*Ye(0,1)*Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1
      ,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2*AeInput(0,1)*AeInput(1,1)*TCC0(mseInput(1
      ),mseInput(0),mslInput(1))*Ye(0,1)*Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1
      ) + Ye(2,0)*Ye(2,1)) - 2*AeInput(1,0)*AeInput(1,1)*TCC0(mseInput(1),mslInput
      (0),mslInput(1))*Ye(1,0)*Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0
      )*Ye(2,1)) - 2*AeInput(1,0)*AeInput(1,1)*TCC0(mseInput(1),mslInput(1),
      mslInput(0))*Ye(1,0)*Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye
      (2,1)) - 2*AeInput(0,2)*AeInput(1,2)*TCC0(mseInput(0),mseInput(1),mslInput(2
      ))*Ye(0,2)*Ye(1,2)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2
      *AeInput(0,2)*AeInput(1,2)*TCC0(mseInput(1),mseInput(0),mslInput(2))*Ye(0,2)
      *Ye(1,2)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2*AeInput(2
      ,0)*AeInput(2,1)*TCC0(mseInput(2),mslInput(0),mslInput(1))*Ye(2,0)*Ye(2,1)*(
      Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - 2*AeInput(2,0)*
      AeInput(2,1)*TCC0(mseInput(2),mslInput(1),mslInput(0))*Ye(2,0)*Ye(2,1)*(Ye(0
      ,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - AeInput(0,0)*AeInput(0,2)
      *AeInput(2,0)*AeInput(2,2)*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput
      (2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) - AeInput(0,0)*AeInput(0,2)*AeInput(2,0
      )*AeInput(2,2)*TCD0(mseInput(0),mseInput(2),mslInput(2),mslInput(0))*Ye(0,0)
      *Ye(0,2)*Ye(2,0)*Ye(2,2) - AeInput(0,0)*AeInput(0,2)*AeInput(2,0)*AeInput(2,
      2)*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(
      2,0)*Ye(2,2) - AeInput(0,0)*AeInput(0,2)*AeInput(2,0)*AeInput(2,2)*TCD0(
      mseInput(2),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(
      2,2) - AeInput(1,0)*AeInput(1,2)*AeInput(2,0)*AeInput(2,2)*TCD0(mseInput(1),
      mseInput(2),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) -
      AeInput(1,0)*AeInput(1,2)*AeInput(2,0)*AeInput(2,2)*TCD0(mseInput(1),
      mseInput(2),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) -
      AeInput(1,0)*AeInput(1,2)*AeInput(2,0)*AeInput(2,2)*TCD0(mseInput(2),
      mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) -
      AeInput(1,0)*AeInput(1,2)*AeInput(2,0)*AeInput(2,2)*TCD0(mseInput(2),
      mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(0),
      mseInput(2),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(0),
      mseInput(2),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2),
      mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) -
      AeInput(0,1)*AeInput(0,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2),
      mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) -
      AeInput(1,1)*AeInput(1,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(1),
      mseInput(2),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) -
      AeInput(1,1)*AeInput(1,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(1),
      mseInput(2),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) -
      AeInput(1,1)*AeInput(1,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2),
      mseInput(1),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) -
      AeInput(1,1)*AeInput(1,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2),
      mseInput(1),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) - 2*
      AeInput(0,0)*AeInput(0,2)*TCC0(mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*
      Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(0,
      0)*AeInput(0,2)*TCC0(mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*(
      Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(1,0)*
      AeInput(1,2)*TCC0(mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*(Ye(0
      ,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(1,0)*AeInput(1,
      2)*TCC0(mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*(Ye(0,0)*Ye(0,2
      ) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(0,0)*AeInput(2,0)*TCC0(
      mseInput(0),mseInput(2),mslInput(0))*Ye(0,0)*Ye(2,0)*(Ye(0,0)*Ye(0,2) + Ye(1
      ,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(0,0)*AeInput(2,0)*TCC0(mseInput(2
      ),mseInput(0),mslInput(0))*Ye(0,0)*Ye(2,0)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2
      ) + Ye(2,0)*Ye(2,2)) - 2*AeInput(0,1)*AeInput(2,1)*TCC0(mseInput(0),mseInput
      (2),mslInput(1))*Ye(0,1)*Ye(2,1)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0
      )*Ye(2,2)) - 2*AeInput(0,1)*AeInput(2,1)*TCC0(mseInput(2),mseInput(0),
      mslInput(1))*Ye(0,1)*Ye(2,1)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye
      (2,2)) - 2*AeInput(0,2)*AeInput(2,2)*TCC0(mseInput(0),mseInput(2),mslInput(2
      ))*Ye(0,2)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2
      *AeInput(0,2)*AeInput(2,2)*TCC0(mseInput(2),mseInput(0),mslInput(2))*Ye(0,2)
      *Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(2
      ,0)*AeInput(2,2)*TCC0(mseInput(2),mslInput(0),mslInput(2))*Ye(2,0)*Ye(2,2)*(
      Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(2,0)*
      AeInput(2,2)*TCC0(mseInput(2),mslInput(2),mslInput(0))*Ye(2,0)*Ye(2,2)*(Ye(0
      ,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - 2*AeInput(0,1)*AeInput(0,
      2)*TCC0(mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*(Ye(0,1)*Ye(0,2
      ) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2*AeInput(0,1)*AeInput(0,2)*TCC0(
      mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*(Ye(0,1)*Ye(0,2) + Ye(1
      ,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2*AeInput(1,1)*AeInput(1,2)*TCC0(mseInput(1
      ),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2
      ) + Ye(2,1)*Ye(2,2)) - 2*AeInput(1,1)*AeInput(1,2)*TCC0(mseInput(1),mslInput
      (2),mslInput(1))*Ye(1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1
      )*Ye(2,2)) - 2*AeInput(1,0)*AeInput(2,0)*TCC0(mseInput(1),mseInput(2),
      mslInput(0))*Ye(1,0)*Ye(2,0)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye
      (2,2)) - 2*AeInput(1,0)*AeInput(2,0)*TCC0(mseInput(2),mseInput(1),mslInput(0
      ))*Ye(1,0)*Ye(2,0)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2
      *AeInput(1,1)*AeInput(2,1)*TCC0(mseInput(1),mseInput(2),mslInput(1))*Ye(1,1)
      *Ye(2,1)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2*AeInput(1
      ,1)*AeInput(2,1)*TCC0(mseInput(2),mseInput(1),mslInput(1))*Ye(1,1)*Ye(2,1)*(
      Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2*AeInput(1,2)*
      AeInput(2,2)*TCC0(mseInput(1),mseInput(2),mslInput(2))*Ye(1,2)*Ye(2,2)*(Ye(0
      ,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2*AeInput(1,2)*AeInput(2,
      2)*TCC0(mseInput(2),mseInput(1),mslInput(2))*Ye(1,2)*Ye(2,2)*(Ye(0,1)*Ye(0,2
      ) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2*AeInput(2,1)*AeInput(2,2)*TCC0(
      mseInput(2),mslInput(1),mslInput(2))*Ye(2,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1
      ,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 2*AeInput(2,1)*AeInput(2,2)*TCC0(mseInput(2
      ),mslInput(2),mslInput(1))*Ye(2,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2
      ) + Ye(2,1)*Ye(2,2)) - 3*Quad(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),
      msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0
      ,1)*Yu(1,0)*Yu(1,1) - 3*Quad(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),
      msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0
      ,1)*Yu(1,0)*Yu(1,1) - 3*Quad(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),
      msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0
      ,2)*Yu(1,0)*Yu(1,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),
      msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0
      ,2)*Yu(1,0)*Yu(1,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),
      msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(1),msuInput(0))*Yu(0,1)*Yu(0
      ,2)*Yu(1,1)*Yu(1,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),
      msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(1),msuInput(0))*Yu(0,1)*Yu(0
      ,2)*Yu(1,1)*Yu(1,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),
      msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0
      ,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),
      msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0
      ,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),
      msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1
      ,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),
      msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1
      ,1)*Yu(2,0)*Yu(2,1) - 3*Quad(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),
      msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(2),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0
      ,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),
      msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0
      ,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),
      msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(2),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1
      ,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),
      msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1
      ,2)*Yu(2,0)*Yu(2,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),
      msuInput(0),msuInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(0))*Yu(0,1)*Yu(0
      ,2)*Yu(2,1)*Yu(2,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),
      msuInput(0),msuInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(0))*Yu(0,1)*Yu(0
      ,2)*Yu(2,1)*Yu(2,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),
      msuInput(1),msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(1))*Yu(1,1)*Yu(1
      ,2)*Yu(2,1)*Yu(2,2) - 3*Quad(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),
      msuInput(1),msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*Quad(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(1
      ,2)*Yu(2,1)*Yu(2,2))))));
   MODEL->set_Lambda2(Re(0.5*(0.25*(0.6*Sqr(g1) + Sqr(g2)) + 0.000641623890917771*
      ((-2*AuInput(2,2))/MSUSY + (0.3333333333333333*Cube(AuInput(2,2)))/Cube(
      MSUSY) - (0.08333333333333333*Quad(AuInput(2,2)))/Quad(MSUSY))*Quad(Yu(2,2))
      *Sqr(g3)*UnitStep(-2 + LambdaLoopOrder) + UnitStep(-1 + LambdaLoopOrder)*(
      0.5*(-0.00037995443865876665*(2*Log(Sqr(Abs(MuInput))/Sqr(SCALE)) +
      0.3333333333333333*Log(Sqr(msdInput(0))/Sqr(SCALE)) + 0.3333333333333333*Log
      (Sqr(msdInput(1))/Sqr(SCALE)) + 0.3333333333333333*Log(Sqr(msdInput(2))/Sqr(
      SCALE)) + Log(Sqr(mseInput(0))/Sqr(SCALE)) + Log(Sqr(mseInput(1))/Sqr(SCALE)
      ) + Log(Sqr(mseInput(2))/Sqr(SCALE)) + 0.5*Log(Sqr(mslInput(0))/Sqr(SCALE))
      + 0.5*Log(Sqr(mslInput(1))/Sqr(SCALE)) + 0.5*Log(Sqr(mslInput(2))/Sqr(SCALE)
      ) + 0.16666666666666666*Log(Sqr(msqInput(0))/Sqr(SCALE)) +
      0.16666666666666666*Log(Sqr(msqInput(1))/Sqr(SCALE)) + 0.16666666666666666*
      Log(Sqr(msqInput(2))/Sqr(SCALE)) + 1.3333333333333333*Log(Sqr(msuInput(0))/
      Sqr(SCALE)) + 1.3333333333333333*Log(Sqr(msuInput(1))/Sqr(SCALE)) +
      1.3333333333333333*Log(Sqr(msuInput(2))/Sqr(SCALE)))*Quad(g1) -
      0.0005277144981371759*(-4 + 8*Log(Sqr(M2Input)/Sqr(SCALE)) + 4*Log(Sqr(Abs(
      MuInput))/Sqr(SCALE)) + Log(Sqr(mslInput(0))/Sqr(SCALE)) + Log(Sqr(mslInput(
      1))/Sqr(SCALE)) + Log(Sqr(mslInput(2))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(0))/
      Sqr(SCALE)) + 3*Log(Sqr(msqInput(1))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(2))/
      Sqr(SCALE)))*Quad(g2)) + Re(0.0031662869888230555*(3*Sqr(Abs(MuInput))*Sqr(
      Yd(0,0))*TCDB0(msdInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*
      TCDB0(msdInput(0),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,2))*TCDB0(
      msdInput(0),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*TCDB0(msdInput(1
      ),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*TCDB0(msdInput(1),msqInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(1,2))*TCDB0(msdInput(1),msqInput(2)) + 3*
      Sqr(Abs(MuInput))*Sqr(Yd(2,0))*TCDB0(msdInput(2),msqInput(0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,1))*TCDB0(msdInput(2),msqInput(1)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yd(2,2))*TCDB0(msdInput(2),msqInput(2)) + Sqr(Abs(MuInput))*Sqr(Ye(0,0))
      *TCDB0(mseInput(0),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(Ye(0,1))*TCDB0(
      mseInput(0),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(Ye(0,2))*TCDB0(mseInput(0),
      mslInput(2)) + Sqr(Abs(MuInput))*Sqr(Ye(1,0))*TCDB0(mseInput(1),mslInput(0))
      + Sqr(Abs(MuInput))*Sqr(Ye(1,1))*TCDB0(mseInput(1),mslInput(1)) + Sqr(Abs(
      MuInput))*Sqr(Ye(1,2))*TCDB0(mseInput(1),mslInput(2)) + Sqr(Abs(MuInput))*
      Sqr(Ye(2,0))*TCDB0(mseInput(2),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(Ye(2,1))
      *TCDB0(mseInput(2),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(Ye(2,2))*TCDB0(
      mseInput(2),mslInput(2)) + 3*Sqr(AuInput(0,0))*Sqr(Yu(0,0))*TCDB0(msuInput(0
      ),msqInput(0)) + 3*Sqr(AuInput(0,1))*Sqr(Yu(0,1))*TCDB0(msuInput(0),msqInput
      (1)) + 3*Sqr(AuInput(0,2))*Sqr(Yu(0,2))*TCDB0(msuInput(0),msqInput(2)) + 3*
      Sqr(AuInput(1,0))*Sqr(Yu(1,0))*TCDB0(msuInput(1),msqInput(0)) + 3*Sqr(
      AuInput(1,1))*Sqr(Yu(1,1))*TCDB0(msuInput(1),msqInput(1)) + 3*Sqr(AuInput(1,
      2))*Sqr(Yu(1,2))*TCDB0(msuInput(1),msqInput(2)) + 3*Sqr(AuInput(2,0))*Sqr(Yu
      (2,0))*TCDB0(msuInput(2),msqInput(0)) + 3*Sqr(AuInput(2,1))*Sqr(Yu(2,1))*
      TCDB0(msuInput(2),msqInput(1)) + 3*Sqr(AuInput(2,2))*Sqr(Yu(2,2))*TCDB0(
      msuInput(2),msqInput(2))) - 0.0007915717472057639*(0.6*Sqr(g1)*TCW(Abs(
      M1Input),Abs(MuInput),SCALE) + 3*Sqr(g2)*TCW(Abs(M2Input),Abs(MuInput),SCALE
      )))*(0.6*Sqr(g1) + Sqr(g2)) + 0.006332573977646111*(-0.03*Quad(g1)*TCB0(
      msdInput(0),msdInput(0),SCALE) - 0.03*Quad(g1)*TCB0(msdInput(1),msdInput(1),
      SCALE) - 0.03*Quad(g1)*TCB0(msdInput(2),msdInput(2),SCALE) - 0.09*Quad(g1)*
      TCB0(mseInput(0),mseInput(0),SCALE) - 0.09*Quad(g1)*TCB0(mseInput(1),
      mseInput(1),SCALE) - 0.09*Quad(g1)*TCB0(mseInput(2),mseInput(2),SCALE) +
      0.125*(-0.36*Quad(g1) - Quad(g2))*TCB0(mslInput(0),mslInput(0),SCALE) +
      0.125*(-0.36*Quad(g1) - Quad(g2))*TCB0(mslInput(1),mslInput(1),SCALE) +
      0.125*(-0.36*Quad(g1) - Quad(g2))*TCB0(mslInput(2),mslInput(2),SCALE) + (
      0.041666666666666664*(-0.36*Quad(g1) - 9*Quad(g2)) + 0.5*(-0.6*Sqr(g1) + 3*
      Sqr(g2))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))) - 3*Sqr(Sqr(Yu(0,0)) +
      Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCB0(msqInput(0),msqInput(0),SCALE) - 3*Sqr(Yu
      (0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))*TCB0(msqInput(0),msqInput
      (1),SCALE) - 3*Sqr(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))*TCB0
      (msqInput(0),msqInput(2),SCALE) - 3*Sqr(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) +
      Yu(2,0)*Yu(2,1))*TCB0(msqInput(1),msqInput(0),SCALE) + (0.041666666666666664
      *(-0.36*Quad(g1) - 9*Quad(g2)) + 0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*(Sqr(Yu(0,1)
      ) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))) - 3*Sqr(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(
      Yu(2,1))))*TCB0(msqInput(1),msqInput(1),SCALE) - 3*Sqr(Yu(0,1)*Yu(0,2) + Yu(
      1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))*TCB0(msqInput(1),msqInput(2),SCALE) - 3*Sqr(
      Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))*TCB0(msqInput(2),
      msqInput(0),SCALE) - 3*Sqr(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,
      2))*TCB0(msqInput(2),msqInput(1),SCALE) + (0.041666666666666664*(-0.36*Quad(
      g1) - 9*Quad(g2)) + 0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,
      2)) + Sqr(Yu(2,2))) - 3*Sqr(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*
      TCB0(msqInput(2),msqInput(2),SCALE) + (-0.12*Quad(g1) + 1.2*Sqr(g1)*(Sqr(Yu(
      0,0)) + Sqr(Yu(0,1)) + Sqr(Yu(0,2))) - 3*Sqr(Sqr(Yu(0,0)) + Sqr(Yu(0,1)) +
      Sqr(Yu(0,2))))*TCB0(msuInput(0),msuInput(0),SCALE) - 3*Sqr(Yu(0,0)*Yu(1,0) +
      Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2))*TCB0(msuInput(0),msuInput(1),SCALE) - 3*
      Sqr(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2))*TCB0(msuInput(0),
      msuInput(2),SCALE) - 3*Sqr(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,
      2))*TCB0(msuInput(1),msuInput(0),SCALE) + (-0.12*Quad(g1) + 1.2*Sqr(g1)*(Sqr
      (Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))) - 3*Sqr(Sqr(Yu(1,0)) + Sqr(Yu(1,1))
      + Sqr(Yu(1,2))))*TCB0(msuInput(1),msuInput(1),SCALE) - 3*Sqr(Yu(1,0)*Yu(2,0)
      + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2))*TCB0(msuInput(1),msuInput(2),SCALE) - 3
      *Sqr(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2))*TCB0(msuInput(2),
      msuInput(0),SCALE) - 3*Sqr(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,
      2))*TCB0(msuInput(2),msuInput(1),SCALE) + (-0.12*Quad(g1) + 1.2*Sqr(g1)*(Sqr
      (Yu(2,0)) + Sqr(Yu(2,1)) + Sqr(Yu(2,2))) - 3*Sqr(Sqr(Yu(2,0)) + Sqr(Yu(2,1))
      + Sqr(Yu(2,2))))*TCB0(msuInput(2),msuInput(2),SCALE) - 0.6*Sqr(g1)*Sqr(Abs(
      MuInput))*Sqr(Yd(0,0))*TCC0(msdInput(0),msdInput(0),msqInput(0)) - 0.6*Sqr(
      g1)*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*TCC0(msdInput(0),msdInput(0),msqInput(1))
      - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(0,2))*TCC0(msdInput(0),msdInput(0),
      msqInput(2)) + 0.5*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(0,0))
      *TCC0(msdInput(0),msqInput(0),msqInput(0)) + 0.5*(-0.6*Sqr(g1) - 3*Sqr(g2))*
      Sqr(Abs(MuInput))*Sqr(Yd(0,1))*TCC0(msdInput(0),msqInput(1),msqInput(1)) +
      0.5*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(0,2))*TCC0(msdInput(
      0),msqInput(2),msqInput(2)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*
      TCC0(msdInput(1),msdInput(1),msqInput(0)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Yd(1,1))*TCC0(msdInput(1),msdInput(1),msqInput(1)) - 0.6*Sqr(g1)*Sqr(Abs
      (MuInput))*Sqr(Yd(1,2))*TCC0(msdInput(1),msdInput(1),msqInput(2)) + 0.5*(-
      0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*TCC0(msdInput(1),
      msqInput(0),msqInput(0)) + 0.5*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*
      Sqr(Yd(1,1))*TCC0(msdInput(1),msqInput(1),msqInput(1)) + 0.5*(-0.6*Sqr(g1) -
      3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(1,2))*TCC0(msdInput(1),msqInput(2),
      msqInput(2)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*TCC0(msdInput(2),
      msdInput(2),msqInput(0)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*TCC0(
      msdInput(2),msdInput(2),msqInput(1)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(
      2,2))*TCC0(msdInput(2),msdInput(2),msqInput(2)) + 0.5*(-0.6*Sqr(g1) - 3*Sqr(
      g2))*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*TCC0(msdInput(2),msqInput(0),msqInput(0)
      ) + 0.5*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*TCC0(
      msdInput(2),msqInput(1),msqInput(1)) + 0.5*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(
      Abs(MuInput))*Sqr(Yd(2,2))*TCC0(msdInput(2),msqInput(2),msqInput(2)) - 0.6*
      Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*TCC0(mseInput(0),mseInput(0),mslInput
      (0)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*TCC0(mseInput(0),mseInput(
      0),mslInput(1)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(0,2))*TCC0(mseInput(0
      ),mseInput(0),mslInput(2)) + 0.5*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Abs(MuInput))*
      Sqr(Ye(0,0))*TCC0(mseInput(0),mslInput(0),mslInput(0)) + 0.5*(0.6*Sqr(g1) -
      Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*TCC0(mseInput(0),mslInput(1),
      mslInput(1)) + 0.5*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(0,2))*
      TCC0(mseInput(0),mslInput(2),mslInput(2)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Ye(1,0))*TCC0(mseInput(1),mseInput(1),mslInput(0)) - 0.6*Sqr(g1)*Sqr(Abs
      (MuInput))*Sqr(Ye(1,1))*TCC0(mseInput(1),mseInput(1),mslInput(1)) - 0.6*Sqr(
      g1)*Sqr(Abs(MuInput))*Sqr(Ye(1,2))*TCC0(mseInput(1),mseInput(1),mslInput(2))
      + 0.5*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*TCC0(mseInput(1
      ),mslInput(0),mslInput(0)) + 0.5*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Abs(MuInput))*
      Sqr(Ye(1,1))*TCC0(mseInput(1),mslInput(1),mslInput(1)) + 0.5*(0.6*Sqr(g1) -
      Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(1,2))*TCC0(mseInput(1),mslInput(2),
      mslInput(2)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*TCC0(mseInput(2),
      mseInput(2),mslInput(0)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(2,1))*TCC0(
      mseInput(2),mseInput(2),mslInput(1)) - 0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(
      2,2))*TCC0(mseInput(2),mseInput(2),mslInput(2)) + 0.5*(0.6*Sqr(g1) - Sqr(g2)
      )*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*TCC0(mseInput(2),mslInput(0),mslInput(0)) +
      0.5*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(2,1))*TCC0(mseInput(2),
      mslInput(1),mslInput(1)) + 0.5*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Abs(MuInput))*Sqr
      (Ye(2,2))*TCC0(mseInput(2),mslInput(2),mslInput(2)) + (0.5*(-0.6*Sqr(g1) + 3
      *Sqr(g2))*Sqr(AuInput(0,0))*Sqr(Yu(0,0)) - 6*Sqr(AuInput(0,0))*Sqr(Yu(0,0))*
      (Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),
      msuInput(0)) + (0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AuInput(1,0))*Sqr(Yu(1,0)
      ) - 6*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2
      ,0))))*TCC0(msqInput(0),msqInput(0),msuInput(1)) + (0.5*(-0.6*Sqr(g1) + 3*
      Sqr(g2))*Sqr(AuInput(2,0))*Sqr(Yu(2,0)) - 6*Sqr(AuInput(2,0))*Sqr(Yu(2,0))*(
      Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),
      msuInput(2)) + (1.2*Sqr(g1)*Sqr(AuInput(0,0))*Sqr(Yu(0,0)) - 6*Sqr(AuInput(0
      ,0))*Sqr(Yu(0,0))*(Sqr(Yu(0,0)) + Sqr(Yu(0,1)) + Sqr(Yu(0,2))))*TCC0(
      msqInput(0),msuInput(0),msuInput(0)) + (1.2*Sqr(g1)*Sqr(AuInput(1,0))*Sqr(Yu
      (1,0)) - 6*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr
      (Yu(1,2))))*TCC0(msqInput(0),msuInput(1),msuInput(1)) + (1.2*Sqr(g1)*Sqr(
      AuInput(2,0))*Sqr(Yu(2,0)) - 6*Sqr(AuInput(2,0))*Sqr(Yu(2,0))*(Sqr(Yu(2,0))
      + Sqr(Yu(2,1)) + Sqr(Yu(2,2))))*TCC0(msqInput(0),msuInput(2),msuInput(2)) +
      (0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AuInput(0,1))*Sqr(Yu(0,1)) - 6*Sqr(
      AuInput(0,1))*Sqr(Yu(0,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*
      TCC0(msqInput(1),msqInput(1),msuInput(0)) + (0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*
      Sqr(AuInput(1,1))*Sqr(Yu(1,1)) - 6*Sqr(AuInput(1,1))*Sqr(Yu(1,1))*(Sqr(Yu(0,
      1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(msqInput(1),msqInput(1),msuInput(1)
      ) + (0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AuInput(2,1))*Sqr(Yu(2,1)) - 6*Sqr(
      AuInput(2,1))*Sqr(Yu(2,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*
      TCC0(msqInput(1),msqInput(1),msuInput(2)) + (1.2*Sqr(g1)*Sqr(AuInput(0,1))*
      Sqr(Yu(0,1)) - 6*Sqr(AuInput(0,1))*Sqr(Yu(0,1))*(Sqr(Yu(0,0)) + Sqr(Yu(0,1))
      + Sqr(Yu(0,2))))*TCC0(msqInput(1),msuInput(0),msuInput(0)) + (1.2*Sqr(g1)*
      Sqr(AuInput(1,1))*Sqr(Yu(1,1)) - 6*Sqr(AuInput(1,1))*Sqr(Yu(1,1))*(Sqr(Yu(1,
      0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))*TCC0(msqInput(1),msuInput(1),msuInput(1)
      ) + (1.2*Sqr(g1)*Sqr(AuInput(2,1))*Sqr(Yu(2,1)) - 6*Sqr(AuInput(2,1))*Sqr(Yu
      (2,1))*(Sqr(Yu(2,0)) + Sqr(Yu(2,1)) + Sqr(Yu(2,2))))*TCC0(msqInput(1),
      msuInput(2),msuInput(2)) + (0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AuInput(0,2))
      *Sqr(Yu(0,2)) - 6*Sqr(AuInput(0,2))*Sqr(Yu(0,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)
      ) + Sqr(Yu(2,2))))*TCC0(msqInput(2),msqInput(2),msuInput(0)) + (0.5*(-0.6*
      Sqr(g1) + 3*Sqr(g2))*Sqr(AuInput(1,2))*Sqr(Yu(1,2)) - 6*Sqr(AuInput(1,2))*
      Sqr(Yu(1,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(msqInput(2),
      msqInput(2),msuInput(1)) + (0.5*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(AuInput(2,2))
      *Sqr(Yu(2,2)) - 6*Sqr(AuInput(2,2))*Sqr(Yu(2,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)
      ) + Sqr(Yu(2,2))))*TCC0(msqInput(2),msqInput(2),msuInput(2)) + (1.2*Sqr(g1)*
      Sqr(AuInput(0,2))*Sqr(Yu(0,2)) - 6*Sqr(AuInput(0,2))*Sqr(Yu(0,2))*(Sqr(Yu(0,
      0)) + Sqr(Yu(0,1)) + Sqr(Yu(0,2))))*TCC0(msqInput(2),msuInput(0),msuInput(0)
      ) + (1.2*Sqr(g1)*Sqr(AuInput(1,2))*Sqr(Yu(1,2)) - 6*Sqr(AuInput(1,2))*Sqr(Yu
      (1,2))*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))*TCC0(msqInput(2),
      msuInput(1),msuInput(1)) + (1.2*Sqr(g1)*Sqr(AuInput(2,2))*Sqr(Yu(2,2)) - 6*
      Sqr(AuInput(2,2))*Sqr(Yu(2,2))*(Sqr(Yu(2,0)) + Sqr(Yu(2,1)) + Sqr(Yu(2,2))))
      *TCC0(msqInput(2),msuInput(2),msuInput(2)) - 3*Quad(Abs(MuInput))*Quad(Yd(0,
      0))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(0)) - 3*Quad(Abs(
      MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(0)
      ,msqInput(1)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,2))*TCD0(msdInput
      (0),msdInput(0),msqInput(0),msqInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,0))
      *Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1),msqInput(0)) - 3*Quad
      (Abs(MuInput))*Quad(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1),
      msqInput(1)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(
      0),msdInput(0),msqInput(1),msqInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,0))*
      Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),msqInput(0)) - 3*Quad(
      Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),
      msqInput(2),msqInput(1)) - 3*Quad(Abs(MuInput))*Quad(Yd(0,2))*TCD0(msdInput(
      0),msdInput(0),msqInput(2),msqInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,0))*
      Sqr(Yd(1,0))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(0)) - 3*Quad(
      Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(0),msdInput(1),
      msqInput(1),msqInput(1)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,2))*Sqr(Yd(1,2))*
      TCD0(msdInput(0),msdInput(1),msqInput(2),msqInput(2)) - 3*Quad(Abs(MuInput))
      *Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput
      (0)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(1)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,2))*Sqr
      (Yd(2,2))*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput(2)) - 3*Quad(Abs
      (MuInput))*Sqr(Yd(0,0))*Sqr(Yd(1,0))*TCD0(msdInput(1),msdInput(0),msqInput(0
      ),msqInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(
      msdInput(1),msdInput(0),msqInput(1),msqInput(1)) - 3*Quad(Abs(MuInput))*Sqr(
      Yd(0,2))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(2))
      - 3*Quad(Abs(MuInput))*Quad(Yd(1,0))*TCD0(msdInput(1),msdInput(1),msqInput(0
      ),msqInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(
      msdInput(1),msdInput(1),msqInput(0),msqInput(1)) - 3*Quad(Abs(MuInput))*Sqr(
      Yd(1,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(0),msqInput(2))
      - 3*Quad(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1
      ),msqInput(1),msqInput(0)) - 3*Quad(Abs(MuInput))*Quad(Yd(1,1))*TCD0(
      msdInput(1),msdInput(1),msqInput(1),msqInput(1)) - 3*Quad(Abs(MuInput))*Sqr(
      Yd(1,1))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(2))
      - 3*Quad(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1
      ),msqInput(2),msqInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yd(1,1))*Sqr(Yd(1,2))*
      TCD0(msdInput(1),msdInput(1),msqInput(2),msqInput(1)) - 3*Quad(Abs(MuInput))
      *Quad(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2),msqInput(2)) - 3*
      Quad(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(1),msdInput(2),
      msqInput(0),msqInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yd(1,1))*Sqr(Yd(2,1))*
      TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(1)) - 3*Quad(Abs(MuInput))
      *Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput
      (2)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(2),
      msdInput(0),msqInput(0),msqInput(0)) - 3*Quad(Abs(MuInput))*Sqr(Yd(0,1))*Sqr
      (Yd(2,1))*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(1)) - 3*Quad(Abs
      (MuInput))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(0),msqInput(2
      ),msqInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(
      msdInput(2),msdInput(1),msqInput(0),msqInput(0)) - 3*Quad(Abs(MuInput))*Sqr(
      Yd(1,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(1))
      - 3*Quad(Abs(MuInput))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(1
      ),msqInput(2),msqInput(2)) - 3*Quad(Abs(MuInput))*Quad(Yd(2,0))*TCD0(
      msdInput(2),msdInput(2),msqInput(0),msqInput(0)) - 3*Quad(Abs(MuInput))*Sqr(
      Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(1))
      - 3*Quad(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2
      ),msqInput(0),msqInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,1))*
      TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(0)) - 3*Quad(Abs(MuInput))
      *Quad(Yd(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(1)) - 3*
      Quad(Abs(MuInput))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(1),msqInput(2)) - 3*Quad(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,2))*
      TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(0)) - 3*Quad(Abs(MuInput))
      *Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput
      (1)) - 3*Quad(Abs(MuInput))*Quad(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(2),msqInput(2)) - Quad(Abs(MuInput))*Quad(Ye(0,0))*TCD0(mseInput(0)
      ,mseInput(0),mslInput(0),mslInput(0)) - Quad(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(
      Ye(0,1))*TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(1)) - Quad(Abs(
      MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(0)
      ,mslInput(2)) - Quad(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(mseInput(0
      ),mseInput(0),mslInput(1),mslInput(0)) - Quad(Abs(MuInput))*Quad(Ye(0,1))*
      TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*
      Sqr(Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(
      2)) - Quad(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput
      (0),mslInput(2),mslInput(0)) - Quad(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(0,2))*
      TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(1)) - Quad(Abs(MuInput))*
      Quad(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(2)) - Quad(
      Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(0),mseInput(1),
      mslInput(0),mslInput(0)) - Quad(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0
      (mseInput(0),mseInput(1),mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*Sqr(
      Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(0),mseInput(1),mslInput(2),mslInput(2))
      - Quad(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(0),mseInput(2),
      mslInput(0),mslInput(0)) - Quad(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0
      (mseInput(0),mseInput(2),mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*Sqr(
      Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(0),mseInput(2),mslInput(2),mslInput(2))
      - Quad(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(1),mseInput(0),
      mslInput(0),mslInput(0)) - Quad(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0
      (mseInput(1),mseInput(0),mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*Sqr(
      Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(0),mslInput(2),mslInput(2))
      - Quad(Abs(MuInput))*Quad(Ye(1,0))*TCD0(mseInput(1),mseInput(1),mslInput(0),
      mslInput(0)) - Quad(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1)
      ,mseInput(1),mslInput(0),mslInput(1)) - Quad(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(
      Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(2)) - Quad(Abs(
      MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1)
      ,mslInput(0)) - Quad(Abs(MuInput))*Quad(Ye(1,1))*TCD0(mseInput(1),mseInput(1
      ),mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*Sqr(Ye(1,1))*Sqr(Ye(1,2))*
      TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(2)) - Quad(Abs(MuInput))*
      Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(2),mslInput(
      0)) - Quad(Abs(MuInput))*Sqr(Ye(1,1))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput
      (1),mslInput(2),mslInput(1)) - Quad(Abs(MuInput))*Quad(Ye(1,2))*TCD0(
      mseInput(1),mseInput(1),mslInput(2),mslInput(2)) - Quad(Abs(MuInput))*Sqr(Ye
      (1,0))*Sqr(Ye(2,0))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(0)) -
      Quad(Abs(MuInput))*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(mseInput(1),mseInput(2),
      mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0
      (mseInput(1),mseInput(2),mslInput(2),mslInput(2)) - Quad(Abs(MuInput))*Sqr(
      Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(0))
      - Quad(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(0),
      mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0
      (mseInput(2),mseInput(0),mslInput(2),mslInput(2)) - Quad(Abs(MuInput))*Sqr(
      Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(1),mslInput(0),mslInput(0))
      - Quad(Abs(MuInput))*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(1),
      mslInput(1),mslInput(1)) - Quad(Abs(MuInput))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0
      (mseInput(2),mseInput(1),mslInput(2),mslInput(2)) - Quad(Abs(MuInput))*Quad(
      Ye(2,0))*TCD0(mseInput(2),mseInput(2),mslInput(0),mslInput(0)) - Quad(Abs(
      MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(0)
      ,mslInput(1)) - Quad(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0(mseInput(2
      ),mseInput(2),mslInput(0),mslInput(2)) - Quad(Abs(MuInput))*Sqr(Ye(2,0))*Sqr
      (Ye(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(0)) - Quad(Abs(
      MuInput))*Quad(Ye(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(1)
      ) - Quad(Abs(MuInput))*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2
      ),mslInput(1),mslInput(2)) - Quad(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*
      TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(0)) - Quad(Abs(MuInput))*
      Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(
      1)) - Quad(Abs(MuInput))*Quad(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput
      (2),mslInput(2)) - 3*Quad(AuInput(0,0))*Quad(Yu(0,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(0),msuInput(0)) - 3*Sqr(AuInput(0,0))*Sqr(AuInput(1,0))
      *Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput
      (1)) - 3*Sqr(AuInput(0,0))*Sqr(AuInput(2,0))*Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(
      msqInput(0),msqInput(0),msuInput(0),msuInput(2)) - 3*Sqr(AuInput(0,0))*Sqr(
      AuInput(1,0))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),
      msuInput(1),msuInput(0)) - 3*Quad(AuInput(1,0))*Quad(Yu(1,0))*TCD0(msqInput(
      0),msqInput(0),msuInput(1),msuInput(1)) - 3*Sqr(AuInput(1,0))*Sqr(AuInput(2,
      0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),
      msuInput(2)) - 3*Sqr(AuInput(0,0))*Sqr(AuInput(2,0))*Sqr(Yu(0,0))*Sqr(Yu(2,0
      ))*TCD0(msqInput(0),msqInput(0),msuInput(2),msuInput(0)) - 3*Sqr(AuInput(1,0
      ))*Sqr(AuInput(2,0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),
      msuInput(2),msuInput(1)) - 3*Quad(AuInput(2,0))*Quad(Yu(2,0))*TCD0(msqInput(
      0),msqInput(0),msuInput(2),msuInput(2)) - 3*Sqr(AuInput(0,0))*Sqr(AuInput(0,
      1))*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(0),msqInput(1),msuInput(0),
      msuInput(0)) - 3*Sqr(AuInput(1,0))*Sqr(AuInput(1,1))*Sqr(Yu(1,0))*Sqr(Yu(1,1
      ))*TCD0(msqInput(0),msqInput(1),msuInput(1),msuInput(1)) - 3*Sqr(AuInput(2,0
      ))*Sqr(AuInput(2,1))*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(0),msqInput(1),
      msuInput(2),msuInput(2)) - 3*Sqr(AuInput(0,0))*Sqr(AuInput(0,2))*Sqr(Yu(0,0)
      )*Sqr(Yu(0,2))*TCD0(msqInput(0),msqInput(2),msuInput(0),msuInput(0)) - 3*Sqr
      (AuInput(1,0))*Sqr(AuInput(1,2))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(0),
      msqInput(2),msuInput(1),msuInput(1)) - 3*Sqr(AuInput(2,0))*Sqr(AuInput(2,2))
      *Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(msqInput(0),msqInput(2),msuInput(2),msuInput
      (2)) - 3*Sqr(AuInput(0,0))*Sqr(AuInput(0,1))*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(
      msqInput(1),msqInput(0),msuInput(0),msuInput(0)) - 3*Sqr(AuInput(1,0))*Sqr(
      AuInput(1,1))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(0),
      msuInput(1),msuInput(1)) - 3*Sqr(AuInput(2,0))*Sqr(AuInput(2,1))*Sqr(Yu(2,0)
      )*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(0),msuInput(2),msuInput(2)) - 3*
      Quad(AuInput(0,1))*Quad(Yu(0,1))*TCD0(msqInput(1),msqInput(1),msuInput(0),
      msuInput(0)) - 3*Sqr(AuInput(0,1))*Sqr(AuInput(1,1))*Sqr(Yu(0,1))*Sqr(Yu(1,1
      ))*TCD0(msqInput(1),msqInput(1),msuInput(0),msuInput(1)) - 3*Sqr(AuInput(0,1
      ))*Sqr(AuInput(2,1))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(0),msuInput(2)) - 3*Sqr(AuInput(0,1))*Sqr(AuInput(1,1))*Sqr(Yu(0,1)
      )*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(0)) - 3*
      Quad(AuInput(1,1))*Quad(Yu(1,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),
      msuInput(1)) - 3*Sqr(AuInput(1,1))*Sqr(AuInput(2,1))*Sqr(Yu(1,1))*Sqr(Yu(2,1
      ))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(2)) - 3*Sqr(AuInput(0,1
      ))*Sqr(AuInput(2,1))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(2),msuInput(0)) - 3*Sqr(AuInput(1,1))*Sqr(AuInput(2,1))*Sqr(Yu(1,1)
      )*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(2),msuInput(1)) - 3*
      Quad(AuInput(2,1))*Quad(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(2),
      msuInput(2)) - 3*Sqr(AuInput(0,1))*Sqr(AuInput(0,2))*Sqr(Yu(0,1))*Sqr(Yu(0,2
      ))*TCD0(msqInput(1),msqInput(2),msuInput(0),msuInput(0)) - 3*Sqr(AuInput(1,1
      ))*Sqr(AuInput(1,2))*Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(1),msqInput(2),
      msuInput(1),msuInput(1)) - 3*Sqr(AuInput(2,1))*Sqr(AuInput(2,2))*Sqr(Yu(2,1)
      )*Sqr(Yu(2,2))*TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(2)) - 3*Sqr
      (AuInput(0,0))*Sqr(AuInput(0,2))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(2),
      msqInput(0),msuInput(0),msuInput(0)) - 3*Sqr(AuInput(1,0))*Sqr(AuInput(1,2))
      *Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(0),msuInput(1),msuInput
      (1)) - 3*Sqr(AuInput(2,0))*Sqr(AuInput(2,2))*Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(
      msqInput(2),msqInput(0),msuInput(2),msuInput(2)) - 3*Sqr(AuInput(0,1))*Sqr(
      AuInput(0,2))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(2),msqInput(1),
      msuInput(0),msuInput(0)) - 3*Sqr(AuInput(1,1))*Sqr(AuInput(1,2))*Sqr(Yu(1,1)
      )*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(1),msuInput(1),msuInput(1)) - 3*Sqr
      (AuInput(2,1))*Sqr(AuInput(2,2))*Sqr(Yu(2,1))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(1),msuInput(2),msuInput(2)) - 3*Quad(AuInput(0,2))*Quad(Yu(0,2))*
      TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput(0)) - 3*Sqr(AuInput(0,2))*
      Sqr(AuInput(1,2))*Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),
      msuInput(0),msuInput(1)) - 3*Sqr(AuInput(0,2))*Sqr(AuInput(2,2))*Sqr(Yu(0,2)
      )*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput(2)) - 3*Sqr
      (AuInput(0,2))*Sqr(AuInput(1,2))*Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(1),msuInput(0)) - 3*Quad(AuInput(1,2))*Quad(Yu(1,2))*
      TCD0(msqInput(2),msqInput(2),msuInput(1),msuInput(1)) - 3*Sqr(AuInput(1,2))*
      Sqr(AuInput(2,2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),
      msuInput(1),msuInput(2)) - 3*Sqr(AuInput(0,2))*Sqr(AuInput(2,2))*Sqr(Yu(0,2)
      )*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(0)) - 3*Sqr
      (AuInput(1,2))*Sqr(AuInput(2,2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(2),msuInput(1)) - 3*Quad(AuInput(2,2))*Quad(Yu(2,2))*
      TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(2)) + 0.36*Quad(g1)*(-0.25
       + 0.5*Sqr(Abs(M1Input))*TCD2t(M1Input,M1Input,Abs(MuInput),Abs(MuInput)) +
      0.5*TCD4t(M1Input,M1Input,Abs(MuInput),Abs(MuInput),SCALE)) + 0.6*Sqr(g1)*
      Sqr(g2)*(-0.5 + M1Input*M2Input*TCD2t(M1Input,M2Input,Abs(MuInput),Abs(
      MuInput)) + TCD4t(M1Input,M2Input,Abs(MuInput),Abs(MuInput),SCALE)) + Quad(
      g2)*(-0.75 + 0.5*Sqr(Abs(M2Input))*TCD2t(M2Input,M2Input,Abs(MuInput),Abs(
      MuInput)) + 2.5*TCD4t(M2Input,M2Input,Abs(MuInput),Abs(MuInput),SCALE)) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(1))*Yd(
      0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*Quad(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(1),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(1))*Yd(
      0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*Quad(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(2))*Yd(
      0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(1),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(2))*Yd(
      0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(2))*Yd(
      0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(1),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(1),msqInput(2))*Yd(
      0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(1))*Yd(
      0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*Quad(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(1))*Yd(
      0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*Quad(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(0),msqInput(1))*Yd(
      1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*Quad(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(2),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(1))*Yd(
      1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*Quad(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(2))*Yd(
      0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(2),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(2))*Yd(
      0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(0),msqInput(2))*Yd(
      1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(2))*Yd(
      1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(1),msqInput(2))*Yd(
      0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(2),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(2))*Yd(
      0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(2))*Yd(
      1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*
      Quad(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(2))*Yd(
      1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*Quad(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(1),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - Quad(
      Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,0)*
      Ye(0,1)*Ye(1,0)*Ye(1,1) - Quad(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),
      mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) - Quad(Abs(MuInput)
      )*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1
      ,0)*Ye(1,1) - Quad(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) - Quad(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(
      1,2) - Quad(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(2),mslInput(
      0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) - Quad(Abs(MuInput))*TCD0(mseInput(1),
      mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) - Quad(
      Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*
      Ye(0,2)*Ye(1,0)*Ye(1,2) - Quad(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),
      mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) - Quad(Abs(MuInput)
      )*TCD0(mseInput(0),mseInput(1),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1
      ,1)*Ye(1,2) - Quad(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(1),
      mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) - Quad(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(
      1,2) - Quad(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(
      1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) - Quad(Abs(MuInput))*TCD0(mseInput(0),
      mseInput(2),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) - Quad(
      Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*
      Ye(0,1)*Ye(2,0)*Ye(2,1) - Quad(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),
      mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) - Quad(Abs(MuInput)
      )*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2
      ,0)*Ye(2,1) - Quad(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(1),
      mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) - Quad(Abs(MuInput))*TCD0(
      mseInput(2),mseInput(1),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(
      2,1) - Quad(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),mslInput(
      0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) - Quad(Abs(MuInput))*TCD0(mseInput(0),
      mseInput(2),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) - Quad(
      Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(2),mslInput(0))*Ye(0,0)*
      Ye(0,2)*Ye(2,0)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),
      mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) - Quad(Abs(MuInput)
      )*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2
      ,0)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(0),
      mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(
      2,2) - Quad(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(0),mslInput(
      2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(mseInput(2),
      mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) - Quad(
      Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(1),mslInput(2))*Ye(0,1)*
      Ye(0,2)*Ye(2,1)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) - Quad(Abs(MuInput)
      )*TCD0(mseInput(2),mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2
      ,1)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(2),
      mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(
      2,2) - Quad(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(2),mslInput(
      1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) - Quad(Abs(MuInput))*TCD0(mseInput(2),
      mseInput(1),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) - Quad(
      Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(1))*Ye(1,1)*
      Ye(1,2)*Ye(2,1)*Ye(2,2) - 3*AuInput(0,0)*AuInput(0,1)*AuInput(1,0)*AuInput(1
      ,1)*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,1)*Yu
      (1,0)*Yu(1,1) - 3*AuInput(0,0)*AuInput(0,1)*AuInput(1,0)*AuInput(1,1)*TCD0(
      msqInput(0),msqInput(1),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(
      1,1) - 3*AuInput(0,0)*AuInput(0,1)*AuInput(1,0)*AuInput(1,1)*TCD0(msqInput(1
      ),msqInput(0),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*
      AuInput(0,0)*AuInput(0,1)*AuInput(1,0)*AuInput(1,1)*TCD0(msqInput(1),
      msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*
      AuInput(0,0)*AuInput(0,2)*AuInput(1,0)*AuInput(1,2)*TCD0(msqInput(0),
      msqInput(2),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*
      AuInput(0,0)*AuInput(0,2)*AuInput(1,0)*AuInput(1,2)*TCD0(msqInput(0),
      msqInput(2),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*
      AuInput(0,0)*AuInput(0,2)*AuInput(1,0)*AuInput(1,2)*TCD0(msqInput(2),
      msqInput(0),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*
      AuInput(0,0)*AuInput(0,2)*AuInput(1,0)*AuInput(1,2)*TCD0(msqInput(2),
      msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(1,1)*AuInput(1,2)*TCD0(msqInput(1),
      msqInput(2),msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(1,1)*AuInput(1,2)*TCD0(msqInput(1),
      msqInput(2),msuInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(1,1)*AuInput(1,2)*TCD0(msqInput(2),
      msqInput(1),msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(1,1)*AuInput(1,2)*TCD0(msqInput(2),
      msqInput(1),msuInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 6*
      AuInput(0,0)*AuInput(1,0)*TCC0(msqInput(0),msuInput(0),msuInput(1))*Yu(0,0)*
      Yu(1,0)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 6*AuInput(0,
      0)*AuInput(1,0)*TCC0(msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(1,0)*(
      Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 6*AuInput(0,1)*
      AuInput(1,1)*TCC0(msqInput(1),msuInput(0),msuInput(1))*Yu(0,1)*Yu(1,1)*(Yu(0
      ,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 6*AuInput(0,1)*AuInput(1,
      1)*TCC0(msqInput(1),msuInput(1),msuInput(0))*Yu(0,1)*Yu(1,1)*(Yu(0,0)*Yu(1,0
      ) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 6*AuInput(0,2)*AuInput(1,2)*TCC0(
      msqInput(2),msuInput(0),msuInput(1))*Yu(0,2)*Yu(1,2)*(Yu(0,0)*Yu(1,0) + Yu(0
      ,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 6*AuInput(0,2)*AuInput(1,2)*TCC0(msqInput(2
      ),msuInput(1),msuInput(0))*Yu(0,2)*Yu(1,2)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1
      ) + Yu(0,2)*Yu(1,2)) - 3*AuInput(0,0)*AuInput(0,1)*AuInput(2,0)*AuInput(2,1)
      *TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,
      0)*Yu(2,1) - 3*AuInput(0,0)*AuInput(0,1)*AuInput(2,0)*AuInput(2,1)*TCD0(
      msqInput(0),msqInput(1),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(
      2,1) - 3*AuInput(0,0)*AuInput(0,1)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput(1
      ),msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) - 3*
      AuInput(0,0)*AuInput(0,1)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput(1),
      msqInput(0),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) - 3*
      AuInput(1,0)*AuInput(1,1)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput(0),
      msqInput(1),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*
      AuInput(1,0)*AuInput(1,1)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput(0),
      msqInput(1),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*
      AuInput(1,0)*AuInput(1,1)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput(1),
      msqInput(0),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*
      AuInput(1,0)*AuInput(1,1)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput(1),
      msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 6*
      AuInput(0,0)*AuInput(0,1)*TCC0(msqInput(0),msqInput(1),msuInput(0))*Yu(0,0)*
      Yu(0,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 6*AuInput(0,
      0)*AuInput(0,1)*TCC0(msqInput(1),msqInput(0),msuInput(0))*Yu(0,0)*Yu(0,1)*(
      Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 6*AuInput(1,0)*
      AuInput(1,1)*TCC0(msqInput(0),msqInput(1),msuInput(1))*Yu(1,0)*Yu(1,1)*(Yu(0
      ,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 6*AuInput(1,0)*AuInput(1,
      1)*TCC0(msqInput(1),msqInput(0),msuInput(1))*Yu(1,0)*Yu(1,1)*(Yu(0,0)*Yu(0,1
      ) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 6*AuInput(2,0)*AuInput(2,1)*TCC0(
      msqInput(0),msqInput(1),msuInput(2))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1
      ,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 6*AuInput(2,0)*AuInput(2,1)*TCC0(msqInput(1
      ),msqInput(0),msuInput(2))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1
      ) + Yu(2,0)*Yu(2,1)) - 3*AuInput(0,0)*AuInput(0,2)*AuInput(2,0)*AuInput(2,2)
      *TCD0(msqInput(0),msqInput(2),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,
      0)*Yu(2,2) - 3*AuInput(0,0)*AuInput(0,2)*AuInput(2,0)*AuInput(2,2)*TCD0(
      msqInput(0),msqInput(2),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(
      2,2) - 3*AuInput(0,0)*AuInput(0,2)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(2
      ),msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) - 3*
      AuInput(0,0)*AuInput(0,2)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(2),
      msqInput(0),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) - 3*
      AuInput(1,0)*AuInput(1,2)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(0),
      msqInput(2),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*
      AuInput(1,0)*AuInput(1,2)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(0),
      msqInput(2),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*
      AuInput(1,0)*AuInput(1,2)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(2),
      msqInput(0),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*
      AuInput(1,0)*AuInput(1,2)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(2),
      msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(1),
      msqInput(2),msuInput(0),msuInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(1),
      msqInput(2),msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(2),
      msqInput(1),msuInput(0),msuInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(0,1)*AuInput(0,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(2),
      msqInput(1),msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(1,1)*AuInput(1,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(1),
      msqInput(2),msuInput(1),msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(1,1)*AuInput(1,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(1),
      msqInput(2),msuInput(2),msuInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(1,1)*AuInput(1,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(2),
      msqInput(1),msuInput(1),msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(1,1)*AuInput(1,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(2),
      msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 6*
      AuInput(0,0)*AuInput(2,0)*TCC0(msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*
      Yu(2,0)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 6*AuInput(0,
      0)*AuInput(2,0)*TCC0(msqInput(0),msuInput(2),msuInput(0))*Yu(0,0)*Yu(2,0)*(
      Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 6*AuInput(0,1)*
      AuInput(2,1)*TCC0(msqInput(1),msuInput(0),msuInput(2))*Yu(0,1)*Yu(2,1)*(Yu(0
      ,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 6*AuInput(0,1)*AuInput(2,
      1)*TCC0(msqInput(1),msuInput(2),msuInput(0))*Yu(0,1)*Yu(2,1)*(Yu(0,0)*Yu(2,0
      ) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 6*AuInput(0,2)*AuInput(2,2)*TCC0(
      msqInput(2),msuInput(0),msuInput(2))*Yu(0,2)*Yu(2,2)*(Yu(0,0)*Yu(2,0) + Yu(0
      ,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 6*AuInput(0,2)*AuInput(2,2)*TCC0(msqInput(2
      ),msuInput(2),msuInput(0))*Yu(0,2)*Yu(2,2)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1
      ) + Yu(0,2)*Yu(2,2)) - 6*AuInput(1,0)*AuInput(2,0)*TCC0(msqInput(0),msuInput
      (1),msuInput(2))*Yu(1,0)*Yu(2,0)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2
      )*Yu(2,2)) - 6*AuInput(1,0)*AuInput(2,0)*TCC0(msqInput(0),msuInput(2),
      msuInput(1))*Yu(1,0)*Yu(2,0)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu
      (2,2)) - 6*AuInput(1,1)*AuInput(2,1)*TCC0(msqInput(1),msuInput(1),msuInput(2
      ))*Yu(1,1)*Yu(2,1)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 6
      *AuInput(1,1)*AuInput(2,1)*TCC0(msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)
      *Yu(2,1)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 6*AuInput(1
      ,2)*AuInput(2,2)*TCC0(msqInput(2),msuInput(1),msuInput(2))*Yu(1,2)*Yu(2,2)*(
      Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 6*AuInput(1,2)*
      AuInput(2,2)*TCC0(msqInput(2),msuInput(2),msuInput(1))*Yu(1,2)*Yu(2,2)*(Yu(1
      ,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 6*AuInput(0,0)*AuInput(0,
      2)*TCC0(msqInput(0),msqInput(2),msuInput(0))*Yu(0,0)*Yu(0,2)*(Yu(0,0)*Yu(0,2
      ) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) - 6*AuInput(0,0)*AuInput(0,2)*TCC0(
      msqInput(2),msqInput(0),msuInput(0))*Yu(0,0)*Yu(0,2)*(Yu(0,0)*Yu(0,2) + Yu(1
      ,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) - 6*AuInput(1,0)*AuInput(1,2)*TCC0(msqInput(0
      ),msqInput(2),msuInput(1))*Yu(1,0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2
      ) + Yu(2,0)*Yu(2,2)) - 6*AuInput(1,0)*AuInput(1,2)*TCC0(msqInput(2),msqInput
      (0),msuInput(1))*Yu(1,0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0
      )*Yu(2,2)) - 6*AuInput(2,0)*AuInput(2,2)*TCC0(msqInput(0),msqInput(2),
      msuInput(2))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu
      (2,2)) - 6*AuInput(2,0)*AuInput(2,2)*TCC0(msqInput(2),msqInput(0),msuInput(2
      ))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) - 6
      *AuInput(0,1)*AuInput(0,2)*TCC0(msqInput(1),msqInput(2),msuInput(0))*Yu(0,1)
      *Yu(0,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 6*AuInput(0
      ,1)*AuInput(0,2)*TCC0(msqInput(2),msqInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*(
      Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 6*AuInput(1,1)*
      AuInput(1,2)*TCC0(msqInput(1),msqInput(2),msuInput(1))*Yu(1,1)*Yu(1,2)*(Yu(0
      ,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 6*AuInput(1,1)*AuInput(1,
      2)*TCC0(msqInput(2),msqInput(1),msuInput(1))*Yu(1,1)*Yu(1,2)*(Yu(0,1)*Yu(0,2
      ) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 6*AuInput(2,1)*AuInput(2,2)*TCC0(
      msqInput(1),msqInput(2),msuInput(2))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1
      ,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 6*AuInput(2,1)*AuInput(2,2)*TCC0(msqInput(2
      ),msqInput(1),msuInput(2))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2
      ) + Yu(2,1)*Yu(2,2)))))));
   MODEL->set_Lambda3(Re(0.25*(-0.6*Sqr(g1) - Sqr(g2)) + 0.5*Sqr(g2) + (
      0.00008020298636472138*AuInput(2,2)*(1 - (0.5*AuInput(2,2))/MSUSY)*Quad(Yu(2
      ,2))*Sqr(g3)*Sqr(MuInput)*UnitStep(-2 + LambdaLoopOrder))/Cube(MSUSY) +
      UnitStep(-1 + LambdaLoopOrder)*(-0.0005277144981371759*(-4 + 8*Log(Sqr(
      M2Input)/Sqr(SCALE)) + 4*Log(Sqr(Abs(MuInput))/Sqr(SCALE)) + Log(Sqr(
      mslInput(0))/Sqr(SCALE)) + Log(Sqr(mslInput(1))/Sqr(SCALE)) + Log(Sqr(
      mslInput(2))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(0))/Sqr(SCALE)) + 3*Log(Sqr(
      msqInput(1))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(2))/Sqr(SCALE)))*Quad(g2) +
      0.5*(0.00037995443865876665*(2*Log(Sqr(Abs(MuInput))/Sqr(SCALE)) +
      0.3333333333333333*Log(Sqr(msdInput(0))/Sqr(SCALE)) + 0.3333333333333333*Log
      (Sqr(msdInput(1))/Sqr(SCALE)) + 0.3333333333333333*Log(Sqr(msdInput(2))/Sqr(
      SCALE)) + Log(Sqr(mseInput(0))/Sqr(SCALE)) + Log(Sqr(mseInput(1))/Sqr(SCALE)
      ) + Log(Sqr(mseInput(2))/Sqr(SCALE)) + 0.5*Log(Sqr(mslInput(0))/Sqr(SCALE))
      + 0.5*Log(Sqr(mslInput(1))/Sqr(SCALE)) + 0.5*Log(Sqr(mslInput(2))/Sqr(SCALE)
      ) + 0.16666666666666666*Log(Sqr(msqInput(0))/Sqr(SCALE)) +
      0.16666666666666666*Log(Sqr(msqInput(1))/Sqr(SCALE)) + 0.16666666666666666*
      Log(Sqr(msqInput(2))/Sqr(SCALE)) + 1.3333333333333333*Log(Sqr(msuInput(0))/
      Sqr(SCALE)) + 1.3333333333333333*Log(Sqr(msuInput(1))/Sqr(SCALE)) +
      1.3333333333333333*Log(Sqr(msuInput(2))/Sqr(SCALE)))*Quad(g1) +
      0.0005277144981371759*(-4 + 8*Log(Sqr(M2Input)/Sqr(SCALE)) + 4*Log(Sqr(Abs(
      MuInput))/Sqr(SCALE)) + Log(Sqr(mslInput(0))/Sqr(SCALE)) + Log(Sqr(mslInput(
      1))/Sqr(SCALE)) + Log(Sqr(mslInput(2))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(0))/
      Sqr(SCALE)) + 3*Log(Sqr(msqInput(1))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(2))/
      Sqr(SCALE)))*Quad(g2)) + 0.5*(Re(0.0031662869888230555*(3*Sqr(AdInput(0,0))*
      Sqr(Yd(0,0))*TCDB0(msdInput(0),msqInput(0)) + 3*Sqr(AdInput(0,1))*Sqr(Yd(0,1
      ))*TCDB0(msdInput(0),msqInput(1)) + 3*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*TCDB0(
      msdInput(0),msqInput(2)) + 3*Sqr(AdInput(1,0))*Sqr(Yd(1,0))*TCDB0(msdInput(1
      ),msqInput(0)) + 3*Sqr(AdInput(1,1))*Sqr(Yd(1,1))*TCDB0(msdInput(1),msqInput
      (1)) + 3*Sqr(AdInput(1,2))*Sqr(Yd(1,2))*TCDB0(msdInput(1),msqInput(2)) + 3*
      Sqr(AdInput(2,0))*Sqr(Yd(2,0))*TCDB0(msdInput(2),msqInput(0)) + 3*Sqr(
      AdInput(2,1))*Sqr(Yd(2,1))*TCDB0(msdInput(2),msqInput(1)) + 3*Sqr(AdInput(2,
      2))*Sqr(Yd(2,2))*TCDB0(msdInput(2),msqInput(2)) + Sqr(AeInput(0,0))*Sqr(Ye(0
      ,0))*TCDB0(mseInput(0),mslInput(0)) + Sqr(AeInput(0,1))*Sqr(Ye(0,1))*TCDB0(
      mseInput(0),mslInput(1)) + Sqr(AeInput(0,2))*Sqr(Ye(0,2))*TCDB0(mseInput(0),
      mslInput(2)) + Sqr(AeInput(1,0))*Sqr(Ye(1,0))*TCDB0(mseInput(1),mslInput(0))
      + Sqr(AeInput(1,1))*Sqr(Ye(1,1))*TCDB0(mseInput(1),mslInput(1)) + Sqr(
      AeInput(1,2))*Sqr(Ye(1,2))*TCDB0(mseInput(1),mslInput(2)) + Sqr(AeInput(2,0)
      )*Sqr(Ye(2,0))*TCDB0(mseInput(2),mslInput(0)) + Sqr(AeInput(2,1))*Sqr(Ye(2,1
      ))*TCDB0(mseInput(2),mslInput(1)) + Sqr(AeInput(2,2))*Sqr(Ye(2,2))*TCDB0(
      mseInput(2),mslInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*TCDB0(msuInput(0
      ),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*TCDB0(msuInput(0),msqInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*TCDB0(msuInput(0),msqInput(2)) + 3*
      Sqr(Abs(MuInput))*Sqr(Yu(1,0))*TCDB0(msuInput(1),msqInput(0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(1,1))*TCDB0(msuInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yu(1,2))*TCDB0(msuInput(1),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,0
      ))*TCDB0(msuInput(2),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*TCDB0(
      msuInput(2),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,2))*TCDB0(msuInput(2
      ),msqInput(2))) - 0.0007915717472057639*(0.6*Sqr(g1)*TCW(Abs(M1Input),Abs(
      MuInput),SCALE) + 3*Sqr(g2)*TCW(Abs(M2Input),Abs(MuInput),SCALE))) + Re(
      0.0031662869888230555*(3*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*TCDB0(msdInput(0),
      msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*TCDB0(msdInput(0),msqInput(1
      )) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,2))*TCDB0(msdInput(0),msqInput(2)) + 3*Sqr
      (Abs(MuInput))*Sqr(Yd(1,0))*TCDB0(msdInput(1),msqInput(0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(1,1))*TCDB0(msdInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yd(1,2))*TCDB0(msdInput(1),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,0
      ))*TCDB0(msdInput(2),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*TCDB0(
      msdInput(2),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,2))*TCDB0(msdInput(2
      ),msqInput(2)) + Sqr(Abs(MuInput))*Sqr(Ye(0,0))*TCDB0(mseInput(0),mslInput(0
      )) + Sqr(Abs(MuInput))*Sqr(Ye(0,1))*TCDB0(mseInput(0),mslInput(1)) + Sqr(Abs
      (MuInput))*Sqr(Ye(0,2))*TCDB0(mseInput(0),mslInput(2)) + Sqr(Abs(MuInput))*
      Sqr(Ye(1,0))*TCDB0(mseInput(1),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(Ye(1,1))
      *TCDB0(mseInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(Ye(1,2))*TCDB0(
      mseInput(1),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(Ye(2,0))*TCDB0(mseInput(2),
      mslInput(0)) + Sqr(Abs(MuInput))*Sqr(Ye(2,1))*TCDB0(mseInput(2),mslInput(1))
      + Sqr(Abs(MuInput))*Sqr(Ye(2,2))*TCDB0(mseInput(2),mslInput(2)) + 3*Sqr(
      AuInput(0,0))*Sqr(Yu(0,0))*TCDB0(msuInput(0),msqInput(0)) + 3*Sqr(AuInput(0,
      1))*Sqr(Yu(0,1))*TCDB0(msuInput(0),msqInput(1)) + 3*Sqr(AuInput(0,2))*Sqr(Yu
      (0,2))*TCDB0(msuInput(0),msqInput(2)) + 3*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*
      TCDB0(msuInput(1),msqInput(0)) + 3*Sqr(AuInput(1,1))*Sqr(Yu(1,1))*TCDB0(
      msuInput(1),msqInput(1)) + 3*Sqr(AuInput(1,2))*Sqr(Yu(1,2))*TCDB0(msuInput(1
      ),msqInput(2)) + 3*Sqr(AuInput(2,0))*Sqr(Yu(2,0))*TCDB0(msuInput(2),msqInput
      (0)) + 3*Sqr(AuInput(2,1))*Sqr(Yu(2,1))*TCDB0(msuInput(2),msqInput(1)) + 3*
      Sqr(AuInput(2,2))*Sqr(Yu(2,2))*TCDB0(msuInput(2),msqInput(2))) -
      0.0007915717472057639*(0.6*Sqr(g1)*TCW(Abs(M1Input),Abs(MuInput),SCALE) + 3*
      Sqr(g2)*TCW(Abs(M2Input),Abs(MuInput),SCALE))))*(-0.6*Sqr(g1) - Sqr(g2)) + (
      Re(0.0031662869888230555*(3*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*TCDB0(msdInput(0)
      ,msqInput(0)) + 3*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*TCDB0(msdInput(0),msqInput(
      1)) + 3*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*TCDB0(msdInput(0),msqInput(2)) + 3*
      Sqr(AdInput(1,0))*Sqr(Yd(1,0))*TCDB0(msdInput(1),msqInput(0)) + 3*Sqr(
      AdInput(1,1))*Sqr(Yd(1,1))*TCDB0(msdInput(1),msqInput(1)) + 3*Sqr(AdInput(1,
      2))*Sqr(Yd(1,2))*TCDB0(msdInput(1),msqInput(2)) + 3*Sqr(AdInput(2,0))*Sqr(Yd
      (2,0))*TCDB0(msdInput(2),msqInput(0)) + 3*Sqr(AdInput(2,1))*Sqr(Yd(2,1))*
      TCDB0(msdInput(2),msqInput(1)) + 3*Sqr(AdInput(2,2))*Sqr(Yd(2,2))*TCDB0(
      msdInput(2),msqInput(2)) + Sqr(AeInput(0,0))*Sqr(Ye(0,0))*TCDB0(mseInput(0),
      mslInput(0)) + Sqr(AeInput(0,1))*Sqr(Ye(0,1))*TCDB0(mseInput(0),mslInput(1))
      + Sqr(AeInput(0,2))*Sqr(Ye(0,2))*TCDB0(mseInput(0),mslInput(2)) + Sqr(
      AeInput(1,0))*Sqr(Ye(1,0))*TCDB0(mseInput(1),mslInput(0)) + Sqr(AeInput(1,1)
      )*Sqr(Ye(1,1))*TCDB0(mseInput(1),mslInput(1)) + Sqr(AeInput(1,2))*Sqr(Ye(1,2
      ))*TCDB0(mseInput(1),mslInput(2)) + Sqr(AeInput(2,0))*Sqr(Ye(2,0))*TCDB0(
      mseInput(2),mslInput(0)) + Sqr(AeInput(2,1))*Sqr(Ye(2,1))*TCDB0(mseInput(2),
      mslInput(1)) + Sqr(AeInput(2,2))*Sqr(Ye(2,2))*TCDB0(mseInput(2),mslInput(2))
      + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*TCDB0(msuInput(0),msqInput(0)) + 3*Sqr(
      Abs(MuInput))*Sqr(Yu(0,1))*TCDB0(msuInput(0),msqInput(1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(0,2))*TCDB0(msuInput(0),msqInput(2)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yu(1,0))*TCDB0(msuInput(1),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(1,1
      ))*TCDB0(msuInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(1,2))*TCDB0(
      msuInput(1),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*TCDB0(msuInput(2
      ),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*TCDB0(msuInput(2),msqInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,2))*TCDB0(msuInput(2),msqInput(2))) -
      0.0007915717472057639*(0.6*Sqr(g1)*TCW(Abs(M1Input),Abs(MuInput),SCALE) + 3*
      Sqr(g2)*TCW(Abs(M2Input),Abs(MuInput),SCALE))) + Re(0.0031662869888230555*(3
      *Sqr(Abs(MuInput))*Sqr(Yd(0,0))*TCDB0(msdInput(0),msqInput(0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(0,1))*TCDB0(msdInput(0),msqInput(1)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yd(0,2))*TCDB0(msdInput(0),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(1,0
      ))*TCDB0(msdInput(1),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*TCDB0(
      msdInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(1,2))*TCDB0(msdInput(1
      ),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*TCDB0(msdInput(2),msqInput
      (0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*TCDB0(msdInput(2),msqInput(1)) + 3*
      Sqr(Abs(MuInput))*Sqr(Yd(2,2))*TCDB0(msdInput(2),msqInput(2)) + Sqr(Abs(
      MuInput))*Sqr(Ye(0,0))*TCDB0(mseInput(0),mslInput(0)) + Sqr(Abs(MuInput))*
      Sqr(Ye(0,1))*TCDB0(mseInput(0),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(Ye(0,2))
      *TCDB0(mseInput(0),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(Ye(1,0))*TCDB0(
      mseInput(1),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(Ye(1,1))*TCDB0(mseInput(1),
      mslInput(1)) + Sqr(Abs(MuInput))*Sqr(Ye(1,2))*TCDB0(mseInput(1),mslInput(2))
      + Sqr(Abs(MuInput))*Sqr(Ye(2,0))*TCDB0(mseInput(2),mslInput(0)) + Sqr(Abs(
      MuInput))*Sqr(Ye(2,1))*TCDB0(mseInput(2),mslInput(1)) + Sqr(Abs(MuInput))*
      Sqr(Ye(2,2))*TCDB0(mseInput(2),mslInput(2)) + 3*Sqr(AuInput(0,0))*Sqr(Yu(0,0
      ))*TCDB0(msuInput(0),msqInput(0)) + 3*Sqr(AuInput(0,1))*Sqr(Yu(0,1))*TCDB0(
      msuInput(0),msqInput(1)) + 3*Sqr(AuInput(0,2))*Sqr(Yu(0,2))*TCDB0(msuInput(0
      ),msqInput(2)) + 3*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*TCDB0(msuInput(1),msqInput
      (0)) + 3*Sqr(AuInput(1,1))*Sqr(Yu(1,1))*TCDB0(msuInput(1),msqInput(1)) + 3*
      Sqr(AuInput(1,2))*Sqr(Yu(1,2))*TCDB0(msuInput(1),msqInput(2)) + 3*Sqr(
      AuInput(2,0))*Sqr(Yu(2,0))*TCDB0(msuInput(2),msqInput(0)) + 3*Sqr(AuInput(2,
      1))*Sqr(Yu(2,1))*TCDB0(msuInput(2),msqInput(1)) + 3*Sqr(AuInput(2,2))*Sqr(Yu
      (2,2))*TCDB0(msuInput(2),msqInput(2))) - 0.0007915717472057639*(0.6*Sqr(g1)*
      TCW(Abs(M1Input),Abs(MuInput),SCALE) + 3*Sqr(g2)*TCW(Abs(M2Input),Abs(
      MuInput),SCALE))))*Sqr(g2) + 0.006332573977646111*((0.03*Quad(g1) - 0.3*Sqr(
      g1)*(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))*TCB0(msdInput(0),msdInput(
      0),SCALE) + (0.03*Quad(g1) - 0.3*Sqr(g1)*(Sqr(Yd(1,0)) + Sqr(Yd(1,1)) + Sqr(
      Yd(1,2))))*TCB0(msdInput(1),msdInput(1),SCALE) + (0.03*Quad(g1) - 0.3*Sqr(g1
      )*(Sqr(Yd(2,0)) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))))*TCB0(msdInput(2),msdInput(2)
      ,SCALE) + (0.09*Quad(g1) - 0.3*Sqr(g1)*(Sqr(Ye(0,0)) + Sqr(Ye(0,1)) + Sqr(Ye
      (0,2))))*TCB0(mseInput(0),mseInput(0),SCALE) + (0.09*Quad(g1) - 0.3*Sqr(g1)*
      (Sqr(Ye(1,0)) + Sqr(Ye(1,1)) + Sqr(Ye(1,2))))*TCB0(mseInput(1),mseInput(1),
      SCALE) + (0.09*Quad(g1) - 0.3*Sqr(g1)*(Sqr(Ye(2,0)) + Sqr(Ye(2,1)) + Sqr(Ye(
      2,2))))*TCB0(mseInput(2),mseInput(2),SCALE) + (0.125*(0.36*Quad(g1) + Quad(
      g2)) + 0.25*(0.6*Sqr(g1) - Sqr(g2))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,
      0))))*TCB0(mslInput(0),mslInput(0),SCALE) + (0.125*(0.36*Quad(g1) + Quad(g2)
      ) + 0.25*(0.6*Sqr(g1) - Sqr(g2))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))
      ))*TCB0(mslInput(1),mslInput(1),SCALE) + (0.125*(0.36*Quad(g1) + Quad(g2)) +
      0.25*(0.6*Sqr(g1) - Sqr(g2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*
      TCB0(mslInput(2),mslInput(2),SCALE) + (0.041666666666666664*(0.36*Quad(g1) +
      9*Quad(g2)) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) +
      Sqr(Yd(2,0))) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0))
      + Sqr(Yu(2,0))))*TCB0(msqInput(0),msqInput(0),SCALE) + (0.041666666666666664
      *(0.36*Quad(g1) + 9*Quad(g2)) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*(Sqr(Yd(0,1)
      ) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*(Sqr(Yu(0,
      1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCB0(msqInput(1),msqInput(1),SCALE) + (
      0.041666666666666664*(0.36*Quad(g1) + 9*Quad(g2)) + 0.25*(-0.6*Sqr(g1) - 3*
      Sqr(g2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) + 0.25*(0.6*Sqr(g1) -
      3*Sqr(g2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCB0(msqInput(2),
      msqInput(2),SCALE) + (0.12*Quad(g1) - 0.6*Sqr(g1)*(Sqr(Yu(0,0)) + Sqr(Yu(0,1
      )) + Sqr(Yu(0,2))))*TCB0(msuInput(0),msuInput(0),SCALE) + (0.12*Quad(g1) -
      0.6*Sqr(g1)*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))*TCB0(msuInput(1),
      msuInput(1),SCALE) + (0.12*Quad(g1) - 0.6*Sqr(g1)*(Sqr(Yu(2,0)) + Sqr(Yu(2,1
      )) + Sqr(Yu(2,2))))*TCB0(msuInput(2),msuInput(2),SCALE) + (0.3*Sqr(g1)*Sqr(
      Abs(MuInput))*Sqr(Yd(0,0)) - 0.3*Sqr(g1)*Sqr(AdInput(0,0))*Sqr(Yd(0,0)) - 3*
      Sqr(Abs(MuInput))*Sqr(Yd(0,0))*(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))
      *TCC0(msdInput(0),msdInput(0),msqInput(0)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Yd(0,1)) - 0.3*Sqr(g1)*Sqr(AdInput(0,1))*Sqr(Yd(0,1)) - 3*Sqr(Abs(
      MuInput))*Sqr(Yd(0,1))*(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))*TCC0(
      msdInput(0),msdInput(0),msqInput(1)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd
      (0,2)) - 0.3*Sqr(g1)*Sqr(AdInput(0,2))*Sqr(Yd(0,2)) - 3*Sqr(Abs(MuInput))*
      Sqr(Yd(0,2))*(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))*TCC0(msdInput(0),
      msdInput(0),msqInput(2)) + (0.25*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))
      *Sqr(Yd(0,0)) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AdInput(0,0))*Sqr(Yd(0,0
      )) - 3*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(
      2,0))))*TCC0(msdInput(0),msqInput(0),msqInput(0)) + (0.25*(0.6*Sqr(g1) + 3*
      Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(0,1)) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*
      Sqr(AdInput(0,1))*Sqr(Yd(0,1)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*(Sqr(Yd(0,
      1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*TCC0(msdInput(0),msqInput(1),msqInput(1)
      ) + (0.25*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(0,2)) + 0.25*(-
      0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AdInput(0,2))*Sqr(Yd(0,2)) - 3*Sqr(Abs(MuInput)
      )*Sqr(Yd(0,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))))*TCC0(msdInput(0
      ),msqInput(2),msqInput(2)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(1,0)) -
      0.3*Sqr(g1)*Sqr(AdInput(1,0))*Sqr(Yd(1,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(1,0)
      )*(Sqr(Yd(1,0)) + Sqr(Yd(1,1)) + Sqr(Yd(1,2))))*TCC0(msdInput(1),msdInput(1)
      ,msqInput(0)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(1,1)) - 0.3*Sqr(g1)*
      Sqr(AdInput(1,1))*Sqr(Yd(1,1)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*(Sqr(Yd(1,
      0)) + Sqr(Yd(1,1)) + Sqr(Yd(1,2))))*TCC0(msdInput(1),msdInput(1),msqInput(1)
      ) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(1,2)) - 0.3*Sqr(g1)*Sqr(AdInput(1,
      2))*Sqr(Yd(1,2)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(1,2))*(Sqr(Yd(1,0)) + Sqr(Yd(1
      ,1)) + Sqr(Yd(1,2))))*TCC0(msdInput(1),msdInput(1),msqInput(2)) + (0.25*(0.6
      *Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(1,0)) + 0.25*(-0.6*Sqr(g1) -
      3*Sqr(g2))*Sqr(AdInput(1,0))*Sqr(Yd(1,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(1,0))
      *(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))))*TCC0(msdInput(1),msqInput(0),
      msqInput(0)) + (0.25*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(1,1)
      ) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AdInput(1,1))*Sqr(Yd(1,1)) - 3*Sqr(
      Abs(MuInput))*Sqr(Yd(1,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*
      TCC0(msdInput(1),msqInput(1),msqInput(1)) + (0.25*(0.6*Sqr(g1) + 3*Sqr(g2))*
      Sqr(Abs(MuInput))*Sqr(Yd(1,2)) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AdInput
      (1,2))*Sqr(Yd(1,2)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(1,2))*(Sqr(Yd(0,2)) + Sqr(
      Yd(1,2)) + Sqr(Yd(2,2))))*TCC0(msdInput(1),msqInput(2),msqInput(2)) + (0.3*
      Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yd(2,0)) - 0.3*Sqr(g1)*Sqr(AdInput(2,0))*Sqr(
      Yd(2,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*(Sqr(Yd(2,0)) + Sqr(Yd(2,1)) +
      Sqr(Yd(2,2))))*TCC0(msdInput(2),msdInput(2),msqInput(0)) + (0.3*Sqr(g1)*Sqr(
      Abs(MuInput))*Sqr(Yd(2,1)) - 0.3*Sqr(g1)*Sqr(AdInput(2,1))*Sqr(Yd(2,1)) - 3*
      Sqr(Abs(MuInput))*Sqr(Yd(2,1))*(Sqr(Yd(2,0)) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))))
      *TCC0(msdInput(2),msdInput(2),msqInput(1)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,2)) - 0.3*Sqr(g1)*Sqr(AdInput(2,2))*Sqr(Yd(2,2)) - 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,2))*(Sqr(Yd(2,0)) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))))*TCC0(
      msdInput(2),msdInput(2),msqInput(2)) + (0.25*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(
      Abs(MuInput))*Sqr(Yd(2,0)) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AdInput(2,0
      ))*Sqr(Yd(2,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,
      0)) + Sqr(Yd(2,0))))*TCC0(msdInput(2),msqInput(0),msqInput(0)) + (0.25*(0.6*
      Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(2,1)) + 0.25*(-0.6*Sqr(g1) - 3
      *Sqr(g2))*Sqr(AdInput(2,1))*Sqr(Yd(2,1)) - 3*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*
      (Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*TCC0(msdInput(2),msqInput(1),
      msqInput(1)) + (0.25*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yd(2,2)
      ) + 0.25*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AdInput(2,2))*Sqr(Yd(2,2)) - 3*Sqr(
      Abs(MuInput))*Sqr(Yd(2,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))))*
      TCC0(msdInput(2),msqInput(2),msqInput(2)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Ye(0,0)) - 0.3*Sqr(g1)*Sqr(AeInput(0,0))*Sqr(Ye(0,0)) - Sqr(Abs(MuInput)
      )*Sqr(Ye(0,0))*(Sqr(Ye(0,0)) + Sqr(Ye(0,1)) + Sqr(Ye(0,2))))*TCC0(mseInput(0
      ),mseInput(0),mslInput(0)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(0,1)) -
      0.3*Sqr(g1)*Sqr(AeInput(0,1))*Sqr(Ye(0,1)) - Sqr(Abs(MuInput))*Sqr(Ye(0,1))*
      (Sqr(Ye(0,0)) + Sqr(Ye(0,1)) + Sqr(Ye(0,2))))*TCC0(mseInput(0),mseInput(0),
      mslInput(1)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(0,2)) - 0.3*Sqr(g1)*Sqr
      (AeInput(0,2))*Sqr(Ye(0,2)) - Sqr(Abs(MuInput))*Sqr(Ye(0,2))*(Sqr(Ye(0,0)) +
      Sqr(Ye(0,1)) + Sqr(Ye(0,2))))*TCC0(mseInput(0),mseInput(0),mslInput(2)) + (
      0.25*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(0,0)) + 0.25*(0.6*Sqr
      (g1) - Sqr(g2))*Sqr(AeInput(0,0))*Sqr(Ye(0,0)) - Sqr(Abs(MuInput))*Sqr(Ye(0,
      0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(0),mslInput(
      0),mslInput(0)) + (0.25*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(0,
      1)) + 0.25*(0.6*Sqr(g1) - Sqr(g2))*Sqr(AeInput(0,1))*Sqr(Ye(0,1)) - Sqr(Abs(
      MuInput))*Sqr(Ye(0,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(
      mseInput(0),mslInput(1),mslInput(1)) + (0.25*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(
      Abs(MuInput))*Sqr(Ye(0,2)) + 0.25*(0.6*Sqr(g1) - Sqr(g2))*Sqr(AeInput(0,2))*
      Sqr(Ye(0,2)) - Sqr(Abs(MuInput))*Sqr(Ye(0,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) +
      Sqr(Ye(2,2))))*TCC0(mseInput(0),mslInput(2),mslInput(2)) + (0.3*Sqr(g1)*Sqr(
      Abs(MuInput))*Sqr(Ye(1,0)) - 0.3*Sqr(g1)*Sqr(AeInput(1,0))*Sqr(Ye(1,0)) -
      Sqr(Abs(MuInput))*Sqr(Ye(1,0))*(Sqr(Ye(1,0)) + Sqr(Ye(1,1)) + Sqr(Ye(1,2))))
      *TCC0(mseInput(1),mseInput(1),mslInput(0)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Ye(1,1)) - 0.3*Sqr(g1)*Sqr(AeInput(1,1))*Sqr(Ye(1,1)) - Sqr(Abs(MuInput)
      )*Sqr(Ye(1,1))*(Sqr(Ye(1,0)) + Sqr(Ye(1,1)) + Sqr(Ye(1,2))))*TCC0(mseInput(1
      ),mseInput(1),mslInput(1)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(1,2)) -
      0.3*Sqr(g1)*Sqr(AeInput(1,2))*Sqr(Ye(1,2)) - Sqr(Abs(MuInput))*Sqr(Ye(1,2))*
      (Sqr(Ye(1,0)) + Sqr(Ye(1,1)) + Sqr(Ye(1,2))))*TCC0(mseInput(1),mseInput(1),
      mslInput(2)) + (0.25*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(1,0))
      + 0.25*(0.6*Sqr(g1) - Sqr(g2))*Sqr(AeInput(1,0))*Sqr(Ye(1,0)) - Sqr(Abs(
      MuInput))*Sqr(Ye(1,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(
      mseInput(1),mslInput(0),mslInput(0)) + (0.25*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(
      Abs(MuInput))*Sqr(Ye(1,1)) + 0.25*(0.6*Sqr(g1) - Sqr(g2))*Sqr(AeInput(1,1))*
      Sqr(Ye(1,1)) - Sqr(Abs(MuInput))*Sqr(Ye(1,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) +
      Sqr(Ye(2,1))))*TCC0(mseInput(1),mslInput(1),mslInput(1)) + (0.25*(-0.6*Sqr(
      g1) + Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(1,2)) + 0.25*(0.6*Sqr(g1) - Sqr(g2))
      *Sqr(AeInput(1,2))*Sqr(Ye(1,2)) - Sqr(Abs(MuInput))*Sqr(Ye(1,2))*(Sqr(Ye(0,2
      )) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(1),mslInput(2),mslInput(2))
      + (0.3*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Ye(2,0)) - 0.3*Sqr(g1)*Sqr(AeInput(2,0)
      )*Sqr(Ye(2,0)) - Sqr(Abs(MuInput))*Sqr(Ye(2,0))*(Sqr(Ye(2,0)) + Sqr(Ye(2,1))
      + Sqr(Ye(2,2))))*TCC0(mseInput(2),mseInput(2),mslInput(0)) + (0.3*Sqr(g1)*
      Sqr(Abs(MuInput))*Sqr(Ye(2,1)) - 0.3*Sqr(g1)*Sqr(AeInput(2,1))*Sqr(Ye(2,1))
      - Sqr(Abs(MuInput))*Sqr(Ye(2,1))*(Sqr(Ye(2,0)) + Sqr(Ye(2,1)) + Sqr(Ye(2,2))
      ))*TCC0(mseInput(2),mseInput(2),mslInput(1)) + (0.3*Sqr(g1)*Sqr(Abs(MuInput)
      )*Sqr(Ye(2,2)) - 0.3*Sqr(g1)*Sqr(AeInput(2,2))*Sqr(Ye(2,2)) - Sqr(Abs(
      MuInput))*Sqr(Ye(2,2))*(Sqr(Ye(2,0)) + Sqr(Ye(2,1)) + Sqr(Ye(2,2))))*TCC0(
      mseInput(2),mseInput(2),mslInput(2)) + (0.25*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(
      Abs(MuInput))*Sqr(Ye(2,0)) + 0.25*(0.6*Sqr(g1) - Sqr(g2))*Sqr(AeInput(2,0))*
      Sqr(Ye(2,0)) - Sqr(Abs(MuInput))*Sqr(Ye(2,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) +
      Sqr(Ye(2,0))))*TCC0(mseInput(2),mslInput(0),mslInput(0)) + (0.25*(-0.6*Sqr(
      g1) + Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(2,1)) + 0.25*(0.6*Sqr(g1) - Sqr(g2))
      *Sqr(AeInput(2,1))*Sqr(Ye(2,1)) - Sqr(Abs(MuInput))*Sqr(Ye(2,1))*(Sqr(Ye(0,1
      )) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(2),mslInput(1),mslInput(1))
      + (0.25*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Ye(2,2)) + 0.25*(0.6*
      Sqr(g1) - Sqr(g2))*Sqr(AeInput(2,2))*Sqr(Ye(2,2)) - Sqr(Abs(MuInput))*Sqr(Ye
      (2,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(2),
      mslInput(2),mslInput(2)) + (0.25*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput)
      )*Sqr(Yu(0,0)) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AuInput(0,0))*Sqr(Yu(0,0
      )) - 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(
      2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(0)) + (0.25*(-0.6*Sqr(g1) + 3*
      Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(1,0)) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr
      (AuInput(1,0))*Sqr(Yu(1,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*(Sqr(Yu(0,0))
      + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(1)) +
      (0.25*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(2,0)) + 0.25*(0.6*
      Sqr(g1) - 3*Sqr(g2))*Sqr(AuInput(2,0))*Sqr(Yu(2,0)) - 3*Sqr(Abs(MuInput))*
      Sqr(Yu(2,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),
      msqInput(0),msuInput(2)) + (0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(0,0)) - 0.6
      *Sqr(g1)*Sqr(AuInput(0,0))*Sqr(Yu(0,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*(
      Sqr(Yu(0,0)) + Sqr(Yu(0,1)) + Sqr(Yu(0,2))))*TCC0(msqInput(0),msuInput(0),
      msuInput(0)) + (0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(1,0)) - 0.6*Sqr(g1)*Sqr
      (AuInput(1,0))*Sqr(Yu(1,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*(Sqr(Yu(1,0))
      + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))*TCC0(msqInput(0),msuInput(1),msuInput(1)) +
      (0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(2,0)) - 0.6*Sqr(g1)*Sqr(AuInput(2,0))*
      Sqr(Yu(2,0)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*(Sqr(Yu(2,0)) + Sqr(Yu(2,1))
      + Sqr(Yu(2,2))))*TCC0(msqInput(0),msuInput(2),msuInput(2)) + (0.25*(-0.6*Sqr
      (g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(0,1)) + 0.25*(0.6*Sqr(g1) - 3*Sqr
      (g2))*Sqr(AuInput(0,1))*Sqr(Yu(0,1)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*(Sqr
      (Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(msqInput(1),msqInput(1),
      msuInput(0)) + (0.25*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(1,1
      )) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AuInput(1,1))*Sqr(Yu(1,1)) - 3*Sqr(
      Abs(MuInput))*Sqr(Yu(1,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*
      TCC0(msqInput(1),msqInput(1),msuInput(1)) + (0.25*(-0.6*Sqr(g1) + 3*Sqr(g2))
      *Sqr(Abs(MuInput))*Sqr(Yu(2,1)) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AuInput
      (2,1))*Sqr(Yu(2,1)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*(Sqr(Yu(0,1)) + Sqr(
      Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(msqInput(1),msqInput(1),msuInput(2)) + (0.6*
      Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(0,1)) - 0.6*Sqr(g1)*Sqr(AuInput(0,1))*Sqr(
      Yu(0,1)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*(Sqr(Yu(0,0)) + Sqr(Yu(0,1)) +
      Sqr(Yu(0,2))))*TCC0(msqInput(1),msuInput(0),msuInput(0)) + (0.6*Sqr(g1)*Sqr(
      Abs(MuInput))*Sqr(Yu(1,1)) - 0.6*Sqr(g1)*Sqr(AuInput(1,1))*Sqr(Yu(1,1)) - 3*
      Sqr(Abs(MuInput))*Sqr(Yu(1,1))*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))
      *TCC0(msqInput(1),msuInput(1),msuInput(1)) + (0.6*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Yu(2,1)) - 0.6*Sqr(g1)*Sqr(AuInput(2,1))*Sqr(Yu(2,1)) - 3*Sqr(Abs(
      MuInput))*Sqr(Yu(2,1))*(Sqr(Yu(2,0)) + Sqr(Yu(2,1)) + Sqr(Yu(2,2))))*TCC0(
      msqInput(1),msuInput(2),msuInput(2)) + (0.25*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(
      Abs(MuInput))*Sqr(Yu(0,2)) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AuInput(0,2)
      )*Sqr(Yu(0,2)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2
      )) + Sqr(Yu(2,2))))*TCC0(msqInput(2),msqInput(2),msuInput(0)) + (0.25*(-0.6*
      Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(1,2)) + 0.25*(0.6*Sqr(g1) - 3*
      Sqr(g2))*Sqr(AuInput(1,2))*Sqr(Yu(1,2)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(1,2))*(
      Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(msqInput(2),msqInput(2),
      msuInput(1)) + (0.25*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Abs(MuInput))*Sqr(Yu(2,2
      )) + 0.25*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(AuInput(2,2))*Sqr(Yu(2,2)) - 3*Sqr(
      Abs(MuInput))*Sqr(Yu(2,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*
      TCC0(msqInput(2),msqInput(2),msuInput(2)) + (0.6*Sqr(g1)*Sqr(Abs(MuInput))*
      Sqr(Yu(0,2)) - 0.6*Sqr(g1)*Sqr(AuInput(0,2))*Sqr(Yu(0,2)) - 3*Sqr(Abs(
      MuInput))*Sqr(Yu(0,2))*(Sqr(Yu(0,0)) + Sqr(Yu(0,1)) + Sqr(Yu(0,2))))*TCC0(
      msqInput(2),msuInput(0),msuInput(0)) + (0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu
      (1,2)) - 0.6*Sqr(g1)*Sqr(AuInput(1,2))*Sqr(Yu(1,2)) - 3*Sqr(Abs(MuInput))*
      Sqr(Yu(1,2))*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))*TCC0(msqInput(2),
      msuInput(1),msuInput(1)) + (0.6*Sqr(g1)*Sqr(Abs(MuInput))*Sqr(Yu(2,2)) - 0.6
      *Sqr(g1)*Sqr(AuInput(2,2))*Sqr(Yu(2,2)) - 3*Sqr(Abs(MuInput))*Sqr(Yu(2,2))*(
      Sqr(Yu(2,0)) + Sqr(Yu(2,1)) + Sqr(Yu(2,2))))*TCC0(msqInput(2),msuInput(2),
      msuInput(2)) - 6*Quad(Yd(0,0))*Sqr(Abs(MuInput))*Sqr(AdInput(0,0))*TCD0(
      msdInput(0),msdInput(0),msqInput(0),msqInput(0)) - 6*Quad(Yd(0,1))*Sqr(Abs(
      MuInput))*Sqr(AdInput(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1),
      msqInput(1)) - 6*Quad(Yd(0,2))*Sqr(Abs(MuInput))*Sqr(AdInput(0,2))*TCD0(
      msdInput(0),msdInput(0),msqInput(2),msqInput(2)) - 6*Quad(Yd(1,0))*Sqr(Abs(
      MuInput))*Sqr(AdInput(1,0))*TCD0(msdInput(1),msdInput(1),msqInput(0),
      msqInput(0)) - 6*Quad(Yd(1,1))*Sqr(Abs(MuInput))*Sqr(AdInput(1,1))*TCD0(
      msdInput(1),msdInput(1),msqInput(1),msqInput(1)) - 6*Quad(Yd(1,2))*Sqr(Abs(
      MuInput))*Sqr(AdInput(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2),
      msqInput(2)) - 6*Quad(Yd(2,0))*Sqr(Abs(MuInput))*Sqr(AdInput(2,0))*TCD0(
      msdInput(2),msdInput(2),msqInput(0),msqInput(0)) - 6*Quad(Yd(2,1))*Sqr(Abs(
      MuInput))*Sqr(AdInput(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(1),
      msqInput(1)) - 6*Quad(Yd(2,2))*Sqr(Abs(MuInput))*Sqr(AdInput(2,2))*TCD0(
      msdInput(2),msdInput(2),msqInput(2),msqInput(2)) - 2*Quad(Ye(0,0))*Sqr(Abs(
      MuInput))*Sqr(AeInput(0,0))*TCD0(mseInput(0),mseInput(0),mslInput(0),
      mslInput(0)) - 2*Quad(Ye(0,1))*Sqr(Abs(MuInput))*Sqr(AeInput(0,1))*TCD0(
      mseInput(0),mseInput(0),mslInput(1),mslInput(1)) - 2*Quad(Ye(0,2))*Sqr(Abs(
      MuInput))*Sqr(AeInput(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(2),
      mslInput(2)) - 2*Quad(Ye(1,0))*Sqr(Abs(MuInput))*Sqr(AeInput(1,0))*TCD0(
      mseInput(1),mseInput(1),mslInput(0),mslInput(0)) - 2*Quad(Ye(1,1))*Sqr(Abs(
      MuInput))*Sqr(AeInput(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1),
      mslInput(1)) - 2*Quad(Ye(1,2))*Sqr(Abs(MuInput))*Sqr(AeInput(1,2))*TCD0(
      mseInput(1),mseInput(1),mslInput(2),mslInput(2)) - 2*Quad(Ye(2,0))*Sqr(Abs(
      MuInput))*Sqr(AeInput(2,0))*TCD0(mseInput(2),mseInput(2),mslInput(0),
      mslInput(0)) - 2*Quad(Ye(2,1))*Sqr(Abs(MuInput))*Sqr(AeInput(2,1))*TCD0(
      mseInput(2),mseInput(2),mslInput(1),mslInput(1)) - 2*Quad(Ye(2,2))*Sqr(Abs(
      MuInput))*Sqr(AeInput(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2),
      mslInput(2)) - 6*Quad(Yu(0,0))*Sqr(Abs(MuInput))*Sqr(AuInput(0,0))*TCD0(
      msqInput(0),msqInput(0),msuInput(0),msuInput(0)) - 6*Quad(Yu(1,0))*Sqr(Abs(
      MuInput))*Sqr(AuInput(1,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),
      msuInput(1)) - 6*Quad(Yu(2,0))*Sqr(Abs(MuInput))*Sqr(AuInput(2,0))*TCD0(
      msqInput(0),msqInput(0),msuInput(2),msuInput(2)) - 6*Quad(Yu(0,1))*Sqr(Abs(
      MuInput))*Sqr(AuInput(0,1))*TCD0(msqInput(1),msqInput(1),msuInput(0),
      msuInput(0)) - 6*Quad(Yu(1,1))*Sqr(Abs(MuInput))*Sqr(AuInput(1,1))*TCD0(
      msqInput(1),msqInput(1),msuInput(1),msuInput(1)) - 6*Quad(Yu(2,1))*Sqr(Abs(
      MuInput))*Sqr(AuInput(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(2),
      msuInput(2)) - 6*Quad(Yu(0,2))*Sqr(Abs(MuInput))*Sqr(AuInput(0,2))*TCD0(
      msqInput(2),msqInput(2),msuInput(0),msuInput(0)) - 6*Quad(Yu(1,2))*Sqr(Abs(
      MuInput))*Sqr(AuInput(1,2))*TCD0(msqInput(2),msqInput(2),msuInput(1),
      msuInput(1)) - 6*Quad(Yu(2,2))*Sqr(Abs(MuInput))*Sqr(AuInput(2,2))*TCD0(
      msqInput(2),msqInput(2),msuInput(2),msuInput(2)) + 0.36*Quad(g1)*(-0.25 + (
      0.5*Sqr(Abs(M1Input)) + Sqr(Abs(MuInput)))*TCD2t(M1Input,M1Input,Abs(MuInput
      ),Abs(MuInput)) + 0.5*TCD4t(M1Input,M1Input,Abs(MuInput),Abs(MuInput),SCALE)
      ) + 0.6*Sqr(g1)*Sqr(g2)*(-0.5 + (M1Input*M2Input + 2*Sqr(Abs(MuInput)))*
      TCD2t(M1Input,M2Input,Abs(MuInput),Abs(MuInput)) + TCD4t(M1Input,M2Input,Abs
      (MuInput),Abs(MuInput),SCALE)) + Quad(g2)*(-0.75 + (2.5*Sqr(Abs(M2Input)) +
      3*Sqr(Abs(MuInput)))*TCD2t(M2Input,M2Input,Abs(MuInput),Abs(MuInput)) + 0.5*
      TCD4t(M2Input,M2Input,Abs(MuInput),Abs(MuInput),SCALE)) - 3*AdInput(0,1)*Sqr
      (Abs(MuInput))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)
      *Yd(0,1)*(AdInput(0,0)*Yd(0,0)*Yd(0,1) + AdInput(0,1)*Yd(0,0)*Yd(0,1)) - 3*
      AdInput(0,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(0),msqInput(1),
      msqInput(0))*Yd(0,0)*Yd(0,1)*(AdInput(0,0)*Yd(0,0)*Yd(0,1) + AdInput(0,1)*Yd
      (0,0)*Yd(0,1)) - 3*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      0),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*(AdInput(0,0)*Yd(0,0)*Yd(0,2) +
      AdInput(0,2)*Yd(0,0)*Yd(0,2)) - 3*AdInput(0,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*(AdInput(0,
      0)*Yd(0,0)*Yd(0,2) + AdInput(0,2)*Yd(0,0)*Yd(0,2)) - 3*AdInput(0,2)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0
      ,2)*(AdInput(0,1)*Yd(0,1)*Yd(0,2) + AdInput(0,2)*Yd(0,1)*Yd(0,2)) - 3*
      AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(0),msqInput(2),
      msqInput(1))*Yd(0,1)*Yd(0,2)*(AdInput(0,1)*Yd(0,1)*Yd(0,2) + AdInput(0,2)*Yd
      (0,1)*Yd(0,2)) - 3*AdInput(0,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      1),msqInput(0),msqInput(0))*Yd(0,0)*Yd(1,0)*(AdInput(0,0)*Yd(0,0)*Yd(1,0) +
      AdInput(1,0)*Yd(0,0)*Yd(1,0)) - 3*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(0),msqInput(0),msqInput(0))*Yd(0,0)*Yd(1,0)*(AdInput(0,
      0)*Yd(0,0)*Yd(1,0) + AdInput(1,0)*Yd(0,0)*Yd(1,0)) - 3*AdInput(0,0)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(0))*Yd(0,0)*(
      AdInput(0,1)*Yd(0,1)*Yd(1,0) + AdInput(1,0)*Yd(0,1)*Yd(1,0))*Yd(1,1) - 3*
      AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),
      msqInput(1))*Yd(0,0)*(AdInput(0,1)*Yd(0,1)*Yd(1,0) + AdInput(1,0)*Yd(0,1)*Yd
      (1,0))*Yd(1,1) - 3*AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      1),msqInput(0),msqInput(1))*Yd(0,1)*Yd(1,0)*(AdInput(0,0)*Yd(0,0)*Yd(1,1) +
      AdInput(1,1)*Yd(0,0)*Yd(1,1)) - 3*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(0),msqInput(1),msqInput(0))*Yd(0,1)*Yd(1,0)*(AdInput(0,
      0)*Yd(0,0)*Yd(1,1) + AdInput(1,1)*Yd(0,0)*Yd(1,1)) - 3*AdInput(0,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(1))*Yd(0,1)*Yd(1
      ,1)*(AdInput(0,1)*Yd(0,1)*Yd(1,1) + AdInput(1,1)*Yd(0,1)*Yd(1,1)) - 3*
      AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(1),
      msqInput(1))*Yd(0,1)*Yd(1,1)*(AdInput(0,1)*Yd(0,1)*Yd(1,1) + AdInput(1,1)*Yd
      (0,1)*Yd(1,1)) - 3*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(
      1),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*(AdInput(1,0)*Yd(1,0)*Yd(1,1) +
      AdInput(1,1)*Yd(1,0)*Yd(1,1)) - 3*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*(AdInput(1,
      0)*Yd(1,0)*Yd(1,1) + AdInput(1,1)*Yd(1,0)*Yd(1,1)) - 3*AdInput(0,0)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(2),msqInput(0))*Yd(0,0)*(
      AdInput(0,2)*Yd(0,2)*Yd(1,0) + AdInput(1,0)*Yd(0,2)*Yd(1,0))*Yd(1,2) - 3*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),
      msqInput(2))*Yd(0,0)*(AdInput(0,2)*Yd(0,2)*Yd(1,0) + AdInput(1,0)*Yd(0,2)*Yd
      (1,0))*Yd(1,2) - 3*AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      1),msqInput(2),msqInput(1))*Yd(0,1)*(AdInput(0,2)*Yd(0,2)*Yd(1,1) + AdInput(
      1,1)*Yd(0,2)*Yd(1,1))*Yd(1,2) - 3*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*(AdInput(0,2)*Yd(0,
      2)*Yd(1,1) + AdInput(1,1)*Yd(0,2)*Yd(1,1))*Yd(1,2) - 3*AdInput(0,2)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(2))*Yd(0,2)*Yd(1
      ,0)*(AdInput(0,0)*Yd(0,0)*Yd(1,2) + AdInput(1,2)*Yd(0,0)*Yd(1,2)) - 3*
      AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(2),
      msqInput(0))*Yd(0,2)*Yd(1,0)*(AdInput(0,0)*Yd(0,0)*Yd(1,2) + AdInput(1,2)*Yd
      (0,0)*Yd(1,2)) - 3*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      1),msqInput(1),msqInput(2))*Yd(0,2)*Yd(1,1)*(AdInput(0,1)*Yd(0,1)*Yd(1,2) +
      AdInput(1,2)*Yd(0,1)*Yd(1,2)) - 3*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(0),msqInput(2),msqInput(1))*Yd(0,2)*Yd(1,1)*(AdInput(0,
      1)*Yd(0,1)*Yd(1,2) + AdInput(1,2)*Yd(0,1)*Yd(1,2)) - 3*Sqr(Abs(MuInput))*
      TCC0(msdInput(0),msdInput(1),msqInput(0))*Yd(0,0)*Yd(1,0)*(Yd(0,0)*Yd(1,0) +
      Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),
      msdInput(0),msqInput(0))*Yd(0,0)*Yd(1,0)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1)
      + Yd(0,2)*Yd(1,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(0),msdInput(1),
      msqInput(1))*Yd(0,1)*Yd(1,1)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd
      (1,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),msdInput(0),msqInput(1))*Yd(0,
      1)*Yd(1,1)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 3*Sqr(Abs
      (MuInput))*TCC0(msdInput(0),msdInput(1),msqInput(2))*Yd(0,2)*Yd(1,2)*(Yd(0,0
      )*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 3*Sqr(Abs(MuInput))*TCC0(
      msdInput(1),msdInput(0),msqInput(2))*Yd(0,2)*Yd(1,2)*(Yd(0,0)*Yd(1,0) + Yd(0
      ,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) - 3*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(1),msqInput(2),msqInput(2))*Yd(0,2)*Yd(1,2)*(AdInput(0,
      2)*Yd(0,2)*Yd(1,2) + AdInput(1,2)*Yd(0,2)*Yd(1,2)) - 3*AdInput(1,2)*Sqr(Abs(
      MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(2))*Yd(0,2)*Yd(1
      ,2)*(AdInput(0,2)*Yd(0,2)*Yd(1,2) + AdInput(1,2)*Yd(0,2)*Yd(1,2)) - 3*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(1),msqInput(0),
      msqInput(2))*Yd(1,0)*Yd(1,2)*(AdInput(1,0)*Yd(1,0)*Yd(1,2) + AdInput(1,2)*Yd
      (1,0)*Yd(1,2)) - 3*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(
      1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*(AdInput(1,0)*Yd(1,0)*Yd(1,2) +
      AdInput(1,2)*Yd(1,0)*Yd(1,2)) - 3*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(1),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*(AdInput(1,
      1)*Yd(1,1)*Yd(1,2) + AdInput(1,2)*Yd(1,1)*Yd(1,2)) - 3*AdInput(1,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(1),msdInput(1),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1
      ,2)*(AdInput(1,1)*Yd(1,1)*Yd(1,2) + AdInput(1,2)*Yd(1,1)*Yd(1,2)) - 3*
      AdInput(0,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(0),
      msqInput(0))*Yd(0,0)*Yd(2,0)*(AdInput(0,0)*Yd(0,0)*Yd(2,0) + AdInput(2,0)*Yd
      (0,0)*Yd(2,0)) - 3*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(
      0),msqInput(0),msqInput(0))*Yd(0,0)*Yd(2,0)*(AdInput(0,0)*Yd(0,0)*Yd(2,0) +
      AdInput(2,0)*Yd(0,0)*Yd(2,0)) - 3*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(2),msqInput(0),msqInput(0))*Yd(1,0)*Yd(2,0)*(AdInput(1,
      0)*Yd(1,0)*Yd(2,0) + AdInput(2,0)*Yd(1,0)*Yd(2,0)) - 3*AdInput(2,0)*Sqr(Abs(
      MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(0))*Yd(1,0)*Yd(2
      ,0)*(AdInput(1,0)*Yd(1,0)*Yd(2,0) + AdInput(2,0)*Yd(1,0)*Yd(2,0)) - 3*
      AdInput(0,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(1),
      msqInput(0))*Yd(0,0)*(AdInput(0,1)*Yd(0,1)*Yd(2,0) + AdInput(2,0)*Yd(0,1)*Yd
      (2,0))*Yd(2,1) - 3*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(
      0),msqInput(0),msqInput(1))*Yd(0,0)*(AdInput(0,1)*Yd(0,1)*Yd(2,0) + AdInput(
      2,0)*Yd(0,1)*Yd(2,0))*Yd(2,1) - 3*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(2),msqInput(1),msqInput(0))*Yd(1,0)*(AdInput(1,1)*Yd(1,
      1)*Yd(2,0) + AdInput(2,0)*Yd(1,1)*Yd(2,0))*Yd(2,1) - 3*AdInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(1))*Yd(1,0)*(
      AdInput(1,1)*Yd(1,1)*Yd(2,0) + AdInput(2,0)*Yd(1,1)*Yd(2,0))*Yd(2,1) - 3*
      AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(0),
      msqInput(1))*Yd(0,1)*Yd(2,0)*(AdInput(0,0)*Yd(0,0)*Yd(2,1) + AdInput(2,1)*Yd
      (0,0)*Yd(2,1)) - 3*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(
      0),msqInput(1),msqInput(0))*Yd(0,1)*Yd(2,0)*(AdInput(0,0)*Yd(0,0)*Yd(2,1) +
      AdInput(2,1)*Yd(0,0)*Yd(2,1)) - 3*AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(2),msqInput(1),msqInput(1))*Yd(0,1)*Yd(2,1)*(AdInput(0,
      1)*Yd(0,1)*Yd(2,1) + AdInput(2,1)*Yd(0,1)*Yd(2,1)) - 3*AdInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(1))*Yd(0,1)*Yd(2
      ,1)*(AdInput(0,1)*Yd(0,1)*Yd(2,1) + AdInput(2,1)*Yd(0,1)*Yd(2,1)) - 3*
      AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(0),
      msqInput(1))*Yd(1,1)*Yd(2,0)*(AdInput(1,0)*Yd(1,0)*Yd(2,1) + AdInput(2,1)*Yd
      (1,0)*Yd(2,1)) - 3*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(
      1),msqInput(1),msqInput(0))*Yd(1,1)*Yd(2,0)*(AdInput(1,0)*Yd(1,0)*Yd(2,1) +
      AdInput(2,1)*Yd(1,0)*Yd(2,1)) - 3*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(2),msqInput(1),msqInput(1))*Yd(1,1)*Yd(2,1)*(AdInput(1,
      1)*Yd(1,1)*Yd(2,1) + AdInput(2,1)*Yd(1,1)*Yd(2,1)) - 3*AdInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(1))*Yd(1,1)*Yd(2
      ,1)*(AdInput(1,1)*Yd(1,1)*Yd(2,1) + AdInput(2,1)*Yd(1,1)*Yd(2,1)) - 3*Sqr(
      Abs(MuInput))*TCC0(msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*(Yd(
      0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*Sqr(Abs(MuInput))*TCC0
      (msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*(Yd(0,0)*Yd(0,1) + Yd(
      1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),
      msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1)
      + Yd(2,0)*Yd(2,1)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),msqInput(1),
      msqInput(0))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd
      (2,1)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(2),msqInput(0),msqInput(1))*Yd(2,
      0)*Yd(2,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*Sqr(Abs
      (MuInput))*TCC0(msdInput(2),msqInput(1),msqInput(0))*Yd(2,0)*Yd(2,1)*(Yd(0,0
      )*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*AdInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(1))*Yd(2,0)*Yd(2
      ,1)*(AdInput(2,0)*Yd(2,0)*Yd(2,1) + AdInput(2,1)*Yd(2,0)*Yd(2,1)) - 3*
      AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(2),msqInput(1),
      msqInput(0))*Yd(2,0)*Yd(2,1)*(AdInput(2,0)*Yd(2,0)*Yd(2,1) + AdInput(2,1)*Yd
      (2,0)*Yd(2,1)) - 3*AdInput(0,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      2),msqInput(2),msqInput(0))*Yd(0,0)*(AdInput(0,2)*Yd(0,2)*Yd(2,0) + AdInput(
      2,0)*Yd(0,2)*Yd(2,0))*Yd(2,2) - 3*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(2),msdInput(0),msqInput(0),msqInput(2))*Yd(0,0)*(AdInput(0,2)*Yd(0,
      2)*Yd(2,0) + AdInput(2,0)*Yd(0,2)*Yd(2,0))*Yd(2,2) - 3*AdInput(1,0)*Sqr(Abs(
      MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput(0))*Yd(1,0)*(
      AdInput(1,2)*Yd(1,2)*Yd(2,0) + AdInput(2,0)*Yd(1,2)*Yd(2,0))*Yd(2,2) - 3*
      AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),
      msqInput(2))*Yd(1,0)*(AdInput(1,2)*Yd(1,2)*Yd(2,0) + AdInput(2,0)*Yd(1,2)*Yd
      (2,0))*Yd(2,2) - 3*AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      2),msqInput(2),msqInput(1))*Yd(0,1)*(AdInput(0,2)*Yd(0,2)*Yd(2,1) + AdInput(
      2,1)*Yd(0,2)*Yd(2,1))*Yd(2,2) - 3*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(2),msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*(AdInput(0,2)*Yd(0,
      2)*Yd(2,1) + AdInput(2,1)*Yd(0,2)*Yd(2,1))*Yd(2,2) - 3*AdInput(1,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput(1))*Yd(1,1)*(
      AdInput(1,2)*Yd(1,2)*Yd(2,1) + AdInput(2,1)*Yd(1,2)*Yd(2,1))*Yd(2,2) - 3*
      AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(1),
      msqInput(2))*Yd(1,1)*(AdInput(1,2)*Yd(1,2)*Yd(2,1) + AdInput(2,1)*Yd(1,2)*Yd
      (2,1))*Yd(2,2) - 3*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(
      2),msqInput(0),msqInput(2))*Yd(0,2)*Yd(2,0)*(AdInput(0,0)*Yd(0,0)*Yd(2,2) +
      AdInput(2,2)*Yd(0,0)*Yd(2,2)) - 3*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(2),msdInput(0),msqInput(2),msqInput(0))*Yd(0,2)*Yd(2,0)*(AdInput(0,
      0)*Yd(0,0)*Yd(2,2) + AdInput(2,2)*Yd(0,0)*Yd(2,2)) - 3*AdInput(0,2)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(1),msqInput(2))*Yd(0,2)*Yd(2
      ,1)*(AdInput(0,1)*Yd(0,1)*Yd(2,2) + AdInput(2,2)*Yd(0,1)*Yd(2,2)) - 3*
      AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(2),
      msqInput(1))*Yd(0,2)*Yd(2,1)*(AdInput(0,1)*Yd(0,1)*Yd(2,2) + AdInput(2,2)*Yd
      (0,1)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(0),msdInput(2),msqInput(0
      ))*Yd(0,0)*Yd(2,0)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 3
      *Sqr(Abs(MuInput))*TCC0(msdInput(2),msdInput(0),msqInput(0))*Yd(0,0)*Yd(2,0)
      *(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 3*Sqr(Abs(MuInput))
      *TCC0(msdInput(0),msdInput(2),msqInput(1))*Yd(0,1)*Yd(2,1)*(Yd(0,0)*Yd(2,0)
      + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(2),
      msdInput(0),msqInput(1))*Yd(0,1)*Yd(2,1)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1)
      + Yd(0,2)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(0),msdInput(2),
      msqInput(2))*Yd(0,2)*Yd(2,2)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(2),msdInput(0),msqInput(2))*Yd(0,
      2)*Yd(2,2)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) - 3*AdInput
      (0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput(2)
      )*Yd(0,2)*Yd(2,2)*(AdInput(0,2)*Yd(0,2)*Yd(2,2) + AdInput(2,2)*Yd(0,2)*Yd(2,
      2)) - 3*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput
      (2),msqInput(2))*Yd(0,2)*Yd(2,2)*(AdInput(0,2)*Yd(0,2)*Yd(2,2) + AdInput(2,2
      )*Yd(0,2)*Yd(2,2)) - 3*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(2),msqInput(0),msqInput(2))*Yd(1,2)*Yd(2,0)*(AdInput(1,0)*Yd(1,0)*
      Yd(2,2) + AdInput(2,2)*Yd(1,0)*Yd(2,2)) - 3*AdInput(2,0)*Sqr(Abs(MuInput))*
      TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(0))*Yd(1,2)*Yd(2,0)*(
      AdInput(1,0)*Yd(1,0)*Yd(2,2) + AdInput(2,2)*Yd(1,0)*Yd(2,2)) - 3*AdInput(1,2
      )*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(2))*Yd
      (1,2)*Yd(2,1)*(AdInput(1,1)*Yd(1,1)*Yd(2,2) + AdInput(2,2)*Yd(1,1)*Yd(2,2))
      - 3*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(2),
      msqInput(1))*Yd(1,2)*Yd(2,1)*(AdInput(1,1)*Yd(1,1)*Yd(2,2) + AdInput(2,2)*Yd
      (1,1)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),msdInput(2),msqInput(0
      ))*Yd(1,0)*Yd(2,0)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 3
      *Sqr(Abs(MuInput))*TCC0(msdInput(2),msdInput(1),msqInput(0))*Yd(1,0)*Yd(2,0)
      *(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 3*Sqr(Abs(MuInput))
      *TCC0(msdInput(1),msdInput(2),msqInput(1))*Yd(1,1)*Yd(2,1)*(Yd(1,0)*Yd(2,0)
      + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(2),
      msdInput(1),msqInput(1))*Yd(1,1)*Yd(2,1)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1)
      + Yd(1,2)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),msdInput(2),
      msqInput(2))*Yd(1,2)*Yd(2,2)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(2),msdInput(1),msqInput(2))*Yd(1,
      2)*Yd(2,2)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) - 3*AdInput
      (1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput(2)
      )*Yd(1,2)*Yd(2,2)*(AdInput(1,2)*Yd(1,2)*Yd(2,2) + AdInput(2,2)*Yd(1,2)*Yd(2,
      2)) - 3*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput
      (2),msqInput(2))*Yd(1,2)*Yd(2,2)*(AdInput(1,2)*Yd(1,2)*Yd(2,2) + AdInput(2,2
      )*Yd(1,2)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(0),msqInput(0),
      msqInput(2))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(0),msqInput(2),msqInput(0))*Yd(0,
      0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*Sqr(Abs
      (MuInput))*TCC0(msdInput(1),msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*(Yd(0,0
      )*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(
      msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*(Yd(0,0)*Yd(0,2) + Yd(1
      ,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(2),
      msqInput(0),msqInput(2))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2)
      + Yd(2,0)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(2),msqInput(2),
      msqInput(0))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd
      (2,2)) - 3*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(2),
      msqInput(0),msqInput(2))*Yd(2,0)*Yd(2,2)*(AdInput(2,0)*Yd(2,0)*Yd(2,2) +
      AdInput(2,2)*Yd(2,0)*Yd(2,2)) - 3*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(
      msdInput(2),msdInput(2),msqInput(2),msqInput(0))*Yd(2,0)*Yd(2,2)*(AdInput(2,
      0)*Yd(2,0)*Yd(2,2) + AdInput(2,2)*Yd(2,0)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*
      TCC0(msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) +
      Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(0),
      msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2)
      + Yd(2,1)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),msqInput(1),
      msqInput(2))*Yd(1,1)*Yd(1,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msdInput(1),msqInput(2),msqInput(1))*Yd(1,
      1)*Yd(1,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*Sqr(Abs
      (MuInput))*TCC0(msdInput(2),msqInput(1),msqInput(2))*Yd(2,1)*Yd(2,2)*(Yd(0,1
      )*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(
      msdInput(2),msqInput(2),msqInput(1))*Yd(2,1)*Yd(2,2)*(Yd(0,1)*Yd(0,2) + Yd(1
      ,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(2),msdInput(2),msqInput(1),msqInput(2))*Yd(2,1)*Yd(2,2)*(AdInput(2,
      1)*Yd(2,1)*Yd(2,2) + AdInput(2,2)*Yd(2,1)*Yd(2,2)) - 3*AdInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(1))*Yd(2,1)*Yd(2
      ,2)*(AdInput(2,1)*Yd(2,1)*Yd(2,2) + AdInput(2,2)*Yd(2,1)*Yd(2,2)) - AeInput(
      0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(1))
      *Ye(0,0)*Ye(0,1)*(AeInput(0,0)*Ye(0,0)*Ye(0,1) + AeInput(0,1)*Ye(0,0)*Ye(0,1
      )) - AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(0),mslInput(1)
      ,mslInput(0))*Ye(0,0)*Ye(0,1)*(AeInput(0,0)*Ye(0,0)*Ye(0,1) + AeInput(0,1)*
      Ye(0,0)*Ye(0,1)) - AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(
      0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*(AeInput(0,0)*Ye(0,0)*Ye(0,2) +
      AeInput(0,2)*Ye(0,0)*Ye(0,2)) - AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput
      (0),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*(AeInput(0,0)*Ye(0,
      0)*Ye(0,2) + AeInput(0,2)*Ye(0,0)*Ye(0,2)) - AeInput(0,2)*Sqr(Abs(MuInput))*
      TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*(
      AeInput(0,1)*Ye(0,1)*Ye(0,2) + AeInput(0,2)*Ye(0,1)*Ye(0,2)) - AeInput(0,1)*
      Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(1))*Ye(0
      ,1)*Ye(0,2)*(AeInput(0,1)*Ye(0,1)*Ye(0,2) + AeInput(0,2)*Ye(0,1)*Ye(0,2)) -
      AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),
      mslInput(0))*Ye(0,0)*Ye(1,0)*(AeInput(0,0)*Ye(0,0)*Ye(1,0) + AeInput(1,0)*Ye
      (0,0)*Ye(1,0)) - AeInput(1,0)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0)
      ,mslInput(0),mslInput(0))*Ye(0,0)*Ye(1,0)*(AeInput(0,0)*Ye(0,0)*Ye(1,0) +
      AeInput(1,0)*Ye(0,0)*Ye(1,0)) - AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput
      (0),mseInput(1),mslInput(1),mslInput(0))*Ye(0,0)*(AeInput(0,1)*Ye(0,1)*Ye(1,
      0) + AeInput(1,0)*Ye(0,1)*Ye(1,0))*Ye(1,1) - AeInput(1,1)*Sqr(Abs(MuInput))*
      TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*(AeInput(0,1)*
      Ye(0,1)*Ye(1,0) + AeInput(1,0)*Ye(0,1)*Ye(1,0))*Ye(1,1) - AeInput(0,1)*Sqr(
      Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,1)*
      Ye(1,0)*(AeInput(0,0)*Ye(0,0)*Ye(1,1) + AeInput(1,1)*Ye(0,0)*Ye(1,1)) -
      AeInput(1,0)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,1)*Ye(1,0)*(AeInput(0,0)*Ye(0,0)*Ye(1,1) + AeInput(1,1)*Ye
      (0,0)*Ye(1,1)) - AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1)
      ,mslInput(1),mslInput(1))*Ye(0,1)*Ye(1,1)*(AeInput(0,1)*Ye(0,1)*Ye(1,1) +
      AeInput(1,1)*Ye(0,1)*Ye(1,1)) - AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput
      (1),mseInput(0),mslInput(1),mslInput(1))*Ye(0,1)*Ye(1,1)*(AeInput(0,1)*Ye(0,
      1)*Ye(1,1) + AeInput(1,1)*Ye(0,1)*Ye(1,1)) - AeInput(1,1)*Sqr(Abs(MuInput))*
      TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*(
      AeInput(1,0)*Ye(1,0)*Ye(1,1) + AeInput(1,1)*Ye(1,0)*Ye(1,1)) - AeInput(1,0)*
      Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(0))*Ye(1
      ,0)*Ye(1,1)*(AeInput(1,0)*Ye(1,0)*Ye(1,1) + AeInput(1,1)*Ye(1,0)*Ye(1,1)) -
      AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(2),
      mslInput(0))*Ye(0,0)*(AeInput(0,2)*Ye(0,2)*Ye(1,0) + AeInput(1,0)*Ye(0,2)*Ye
      (1,0))*Ye(1,2) - AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0)
      ,mslInput(0),mslInput(2))*Ye(0,0)*(AeInput(0,2)*Ye(0,2)*Ye(1,0) + AeInput(1,
      0)*Ye(0,2)*Ye(1,0))*Ye(1,2) - AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(0
      ),mseInput(1),mslInput(2),mslInput(1))*Ye(0,1)*(AeInput(0,2)*Ye(0,2)*Ye(1,1)
      + AeInput(1,1)*Ye(0,2)*Ye(1,1))*Ye(1,2) - AeInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(mseInput(1),mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*(AeInput(0,2)*
      Ye(0,2)*Ye(1,1) + AeInput(1,1)*Ye(0,2)*Ye(1,1))*Ye(1,2) - AeInput(0,2)*Sqr(
      Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(2))*Ye(0,2)*
      Ye(1,0)*(AeInput(0,0)*Ye(0,0)*Ye(1,2) + AeInput(1,2)*Ye(0,0)*Ye(1,2)) -
      AeInput(1,0)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(2),
      mslInput(0))*Ye(0,2)*Ye(1,0)*(AeInput(0,0)*Ye(0,0)*Ye(1,2) + AeInput(1,2)*Ye
      (0,0)*Ye(1,2)) - AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1)
      ,mslInput(1),mslInput(2))*Ye(0,2)*Ye(1,1)*(AeInput(0,1)*Ye(0,1)*Ye(1,2) +
      AeInput(1,2)*Ye(0,1)*Ye(1,2)) - AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput
      (1),mseInput(0),mslInput(2),mslInput(1))*Ye(0,2)*Ye(1,1)*(AeInput(0,1)*Ye(0,
      1)*Ye(1,2) + AeInput(1,2)*Ye(0,1)*Ye(1,2)) - Sqr(Abs(MuInput))*TCC0(mseInput
      (0),mseInput(1),mslInput(0))*Ye(0,0)*Ye(1,0)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1
      ,1) + Ye(0,2)*Ye(1,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(1),mseInput(0),
      mslInput(0))*Ye(0,0)*Ye(1,0)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye
      (1,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(0),mseInput(1),mslInput(1))*Ye(0,1)
      *Ye(1,1)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye(1,2)) - Sqr(Abs(
      MuInput))*TCC0(mseInput(1),mseInput(0),mslInput(1))*Ye(0,1)*Ye(1,1)*(Ye(0,0)
      *Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye(1,2)) - Sqr(Abs(MuInput))*TCC0(
      mseInput(0),mseInput(1),mslInput(2))*Ye(0,2)*Ye(1,2)*(Ye(0,0)*Ye(1,0) + Ye(0
      ,1)*Ye(1,1) + Ye(0,2)*Ye(1,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(1),mseInput
      (0),mslInput(2))*Ye(0,2)*Ye(1,2)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2
      )*Ye(1,2)) - AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),
      mslInput(2),mslInput(2))*Ye(0,2)*Ye(1,2)*(AeInput(0,2)*Ye(0,2)*Ye(1,2) +
      AeInput(1,2)*Ye(0,2)*Ye(1,2)) - AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput
      (1),mseInput(0),mslInput(2),mslInput(2))*Ye(0,2)*Ye(1,2)*(AeInput(0,2)*Ye(0,
      2)*Ye(1,2) + AeInput(1,2)*Ye(0,2)*Ye(1,2)) - AeInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*(
      AeInput(1,0)*Ye(1,0)*Ye(1,2) + AeInput(1,2)*Ye(1,0)*Ye(1,2)) - AeInput(1,0)*
      Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(1),mslInput(2),mslInput(0))*Ye(1
      ,0)*Ye(1,2)*(AeInput(1,0)*Ye(1,0)*Ye(1,2) + AeInput(1,2)*Ye(1,0)*Ye(1,2)) -
      AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(1),mslInput(1),
      mslInput(2))*Ye(1,1)*Ye(1,2)*(AeInput(1,1)*Ye(1,1)*Ye(1,2) + AeInput(1,2)*Ye
      (1,1)*Ye(1,2)) - AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(1)
      ,mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*(AeInput(1,1)*Ye(1,1)*Ye(1,2) +
      AeInput(1,2)*Ye(1,1)*Ye(1,2)) - AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput
      (0),mseInput(2),mslInput(0),mslInput(0))*Ye(0,0)*Ye(2,0)*(AeInput(0,0)*Ye(0,
      0)*Ye(2,0) + AeInput(2,0)*Ye(0,0)*Ye(2,0)) - AeInput(2,0)*Sqr(Abs(MuInput))*
      TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(0))*Ye(0,0)*Ye(2,0)*(
      AeInput(0,0)*Ye(0,0)*Ye(2,0) + AeInput(2,0)*Ye(0,0)*Ye(2,0)) - AeInput(1,0)*
      Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(0))*Ye(1
      ,0)*Ye(2,0)*(AeInput(1,0)*Ye(1,0)*Ye(2,0) + AeInput(2,0)*Ye(1,0)*Ye(2,0)) -
      AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(0),
      mslInput(0))*Ye(1,0)*Ye(2,0)*(AeInput(1,0)*Ye(1,0)*Ye(2,0) + AeInput(2,0)*Ye
      (1,0)*Ye(2,0)) - AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2)
      ,mslInput(1),mslInput(0))*Ye(0,0)*(AeInput(0,1)*Ye(0,1)*Ye(2,0) + AeInput(2,
      0)*Ye(0,1)*Ye(2,0))*Ye(2,1) - AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2
      ),mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*(AeInput(0,1)*Ye(0,1)*Ye(2,0)
      + AeInput(2,0)*Ye(0,1)*Ye(2,0))*Ye(2,1) - AeInput(1,0)*Sqr(Abs(MuInput))*
      TCD0(mseInput(1),mseInput(2),mslInput(1),mslInput(0))*Ye(1,0)*(AeInput(1,1)*
      Ye(1,1)*Ye(2,0) + AeInput(2,0)*Ye(1,1)*Ye(2,0))*Ye(2,1) - AeInput(2,1)*Sqr(
      Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(0),mslInput(1))*Ye(1,0)*
      (AeInput(1,1)*Ye(1,1)*Ye(2,0) + AeInput(2,0)*Ye(1,1)*Ye(2,0))*Ye(2,1) -
      AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(0),
      mslInput(1))*Ye(0,1)*Ye(2,0)*(AeInput(0,0)*Ye(0,0)*Ye(2,1) + AeInput(2,1)*Ye
      (0,0)*Ye(2,1)) - AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0)
      ,mslInput(1),mslInput(0))*Ye(0,1)*Ye(2,0)*(AeInput(0,0)*Ye(0,0)*Ye(2,1) +
      AeInput(2,1)*Ye(0,0)*Ye(2,1)) - AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput
      (0),mseInput(2),mslInput(1),mslInput(1))*Ye(0,1)*Ye(2,1)*(AeInput(0,1)*Ye(0,
      1)*Ye(2,1) + AeInput(2,1)*Ye(0,1)*Ye(2,1)) - AeInput(2,1)*Sqr(Abs(MuInput))*
      TCD0(mseInput(2),mseInput(0),mslInput(1),mslInput(1))*Ye(0,1)*Ye(2,1)*(
      AeInput(0,1)*Ye(0,1)*Ye(2,1) + AeInput(2,1)*Ye(0,1)*Ye(2,1)) - AeInput(1,1)*
      Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(1))*Ye(1
      ,1)*Ye(2,0)*(AeInput(1,0)*Ye(1,0)*Ye(2,1) + AeInput(2,1)*Ye(1,0)*Ye(2,1)) -
      AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),
      mslInput(0))*Ye(1,1)*Ye(2,0)*(AeInput(1,0)*Ye(1,0)*Ye(2,1) + AeInput(2,1)*Ye
      (1,0)*Ye(2,1)) - AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2)
      ,mslInput(1),mslInput(1))*Ye(1,1)*Ye(2,1)*(AeInput(1,1)*Ye(1,1)*Ye(2,1) +
      AeInput(2,1)*Ye(1,1)*Ye(2,1)) - AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput
      (2),mseInput(1),mslInput(1),mslInput(1))*Ye(1,1)*Ye(2,1)*(AeInput(1,1)*Ye(1,
      1)*Ye(2,1) + AeInput(2,1)*Ye(1,1)*Ye(2,1)) - Sqr(Abs(MuInput))*TCC0(mseInput
      (0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1
      ,1) + Ye(2,0)*Ye(2,1)) - Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye
      (2,1)) - Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput(0),mslInput(1))*Ye(1,0)
      *Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - Sqr(Abs(
      MuInput))*TCC0(mseInput(1),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*(Ye(0,0)
      *Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - Sqr(Abs(MuInput))*TCC0(
      mseInput(2),mslInput(0),mslInput(1))*Ye(2,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1
      ,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) - Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput
      (1),mslInput(0))*Ye(2,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0
      )*Ye(2,1)) - AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(2),
      mslInput(0),mslInput(1))*Ye(2,0)*Ye(2,1)*(AeInput(2,0)*Ye(2,0)*Ye(2,1) +
      AeInput(2,1)*Ye(2,0)*Ye(2,1)) - AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput
      (2),mseInput(2),mslInput(1),mslInput(0))*Ye(2,0)*Ye(2,1)*(AeInput(2,0)*Ye(2,
      0)*Ye(2,1) + AeInput(2,1)*Ye(2,0)*Ye(2,1)) - AeInput(0,0)*Sqr(Abs(MuInput))*
      TCD0(mseInput(0),mseInput(2),mslInput(2),mslInput(0))*Ye(0,0)*(AeInput(0,2)*
      Ye(0,2)*Ye(2,0) + AeInput(2,0)*Ye(0,2)*Ye(2,0))*Ye(2,2) - AeInput(2,2)*Sqr(
      Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*
      (AeInput(0,2)*Ye(0,2)*Ye(2,0) + AeInput(2,0)*Ye(0,2)*Ye(2,0))*Ye(2,2) -
      AeInput(1,0)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(2),
      mslInput(0))*Ye(1,0)*(AeInput(1,2)*Ye(1,2)*Ye(2,0) + AeInput(2,0)*Ye(1,2)*Ye
      (2,0))*Ye(2,2) - AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1)
      ,mslInput(0),mslInput(2))*Ye(1,0)*(AeInput(1,2)*Ye(1,2)*Ye(2,0) + AeInput(2,
      0)*Ye(1,2)*Ye(2,0))*Ye(2,2) - AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(0
      ),mseInput(2),mslInput(2),mslInput(1))*Ye(0,1)*(AeInput(0,2)*Ye(0,2)*Ye(2,1)
      + AeInput(2,1)*Ye(0,2)*Ye(2,1))*Ye(2,2) - AeInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(mseInput(2),mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*(AeInput(0,2)*
      Ye(0,2)*Ye(2,1) + AeInput(2,1)*Ye(0,2)*Ye(2,1))*Ye(2,2) - AeInput(1,1)*Sqr(
      Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(2),mslInput(1))*Ye(1,1)*
      (AeInput(1,2)*Ye(1,2)*Ye(2,1) + AeInput(2,1)*Ye(1,2)*Ye(2,1))*Ye(2,2) -
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),
      mslInput(2))*Ye(1,1)*(AeInput(1,2)*Ye(1,2)*Ye(2,1) + AeInput(2,1)*Ye(1,2)*Ye
      (2,1))*Ye(2,2) - AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2)
      ,mslInput(0),mslInput(2))*Ye(0,2)*Ye(2,0)*(AeInput(0,0)*Ye(0,0)*Ye(2,2) +
      AeInput(2,2)*Ye(0,0)*Ye(2,2)) - AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput
      (2),mseInput(0),mslInput(2),mslInput(0))*Ye(0,2)*Ye(2,0)*(AeInput(0,0)*Ye(0,
      0)*Ye(2,2) + AeInput(2,2)*Ye(0,0)*Ye(2,2)) - AeInput(0,2)*Sqr(Abs(MuInput))*
      TCD0(mseInput(0),mseInput(2),mslInput(1),mslInput(2))*Ye(0,2)*Ye(2,1)*(
      AeInput(0,1)*Ye(0,1)*Ye(2,2) + AeInput(2,2)*Ye(0,1)*Ye(2,2)) - AeInput(2,1)*
      Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(1))*Ye(0
      ,2)*Ye(2,1)*(AeInput(0,1)*Ye(0,1)*Ye(2,2) + AeInput(2,2)*Ye(0,1)*Ye(2,2)) -
      Sqr(Abs(MuInput))*TCC0(mseInput(0),mseInput(2),mslInput(0))*Ye(0,0)*Ye(2,0)*
      (Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) - Sqr(Abs(MuInput))*
      TCC0(mseInput(2),mseInput(0),mslInput(0))*Ye(0,0)*Ye(2,0)*(Ye(0,0)*Ye(2,0) +
      Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(0),
      mseInput(2),mslInput(1))*Ye(0,1)*Ye(2,1)*(Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1)
      + Ye(0,2)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(2),mseInput(0),mslInput
      (1))*Ye(0,1)*Ye(2,1)*(Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) -
      Sqr(Abs(MuInput))*TCC0(mseInput(0),mseInput(2),mslInput(2))*Ye(0,2)*Ye(2,2)*
      (Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) - Sqr(Abs(MuInput))*
      TCC0(mseInput(2),mseInput(0),mslInput(2))*Ye(0,2)*Ye(2,2)*(Ye(0,0)*Ye(2,0) +
      Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) - AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(2),mslInput(2),mslInput(2))*Ye(0,2)*Ye(2,2)*(AeInput(0,
      2)*Ye(0,2)*Ye(2,2) + AeInput(2,2)*Ye(0,2)*Ye(2,2)) - AeInput(2,2)*Sqr(Abs(
      MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(2))*Ye(0,2)*Ye(2
      ,2)*(AeInput(0,2)*Ye(0,2)*Ye(2,2) + AeInput(2,2)*Ye(0,2)*Ye(2,2)) - AeInput(
      1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(2))
      *Ye(1,2)*Ye(2,0)*(AeInput(1,0)*Ye(1,0)*Ye(2,2) + AeInput(2,2)*Ye(1,0)*Ye(2,2
      )) - AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(2)
      ,mslInput(0))*Ye(1,2)*Ye(2,0)*(AeInput(1,0)*Ye(1,0)*Ye(2,2) + AeInput(2,2)*
      Ye(1,0)*Ye(2,2)) - AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(
      2),mslInput(1),mslInput(2))*Ye(1,2)*Ye(2,1)*(AeInput(1,1)*Ye(1,1)*Ye(2,2) +
      AeInput(2,2)*Ye(1,1)*Ye(2,2)) - AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput
      (2),mseInput(1),mslInput(2),mslInput(1))*Ye(1,2)*Ye(2,1)*(AeInput(1,1)*Ye(1,
      1)*Ye(2,2) + AeInput(2,2)*Ye(1,1)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput
      (1),mseInput(2),mslInput(0))*Ye(1,0)*Ye(2,0)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2
      ,1) + Ye(1,2)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(2),mseInput(1),
      mslInput(0))*Ye(1,0)*Ye(2,0)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2)*Ye
      (2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(1),mseInput(2),mslInput(1))*Ye(1,1)
      *Ye(2,1)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2)*Ye(2,2)) - Sqr(Abs(
      MuInput))*TCC0(mseInput(2),mseInput(1),mslInput(1))*Ye(1,1)*Ye(2,1)*(Ye(1,0)
      *Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(
      mseInput(1),mseInput(2),mslInput(2))*Ye(1,2)*Ye(2,2)*(Ye(1,0)*Ye(2,0) + Ye(1
      ,1)*Ye(2,1) + Ye(1,2)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(2),mseInput
      (1),mslInput(2))*Ye(1,2)*Ye(2,2)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2
      )*Ye(2,2)) - AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),
      mslInput(2),mslInput(2))*Ye(1,2)*Ye(2,2)*(AeInput(1,2)*Ye(1,2)*Ye(2,2) +
      AeInput(2,2)*Ye(1,2)*Ye(2,2)) - AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput
      (2),mseInput(1),mslInput(2),mslInput(2))*Ye(1,2)*Ye(2,2)*(AeInput(1,2)*Ye(1,
      2)*Ye(2,2) + AeInput(2,2)*Ye(1,2)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput
      (0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1
      ,2) + Ye(2,0)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput(2),
      mslInput(0))*Ye(0,0)*Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye
      (2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)
      *Ye(1,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - Sqr(Abs(
      MuInput))*TCC0(mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*(Ye(0,0)
      *Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(
      mseInput(2),mslInput(0),mslInput(2))*Ye(2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1
      ,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput
      (2),mslInput(0))*Ye(2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0
      )*Ye(2,2)) - AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(2),
      mslInput(0),mslInput(2))*Ye(2,0)*Ye(2,2)*(AeInput(2,0)*Ye(2,0)*Ye(2,2) +
      AeInput(2,2)*Ye(2,0)*Ye(2,2)) - AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput
      (2),mseInput(2),mslInput(2),mslInput(0))*Ye(2,0)*Ye(2,2)*(AeInput(2,0)*Ye(2,
      0)*Ye(2,2) + AeInput(2,2)*Ye(2,0)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput
      (0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1
      ,2) + Ye(2,1)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput(2),
      mslInput(1))*Ye(0,1)*Ye(0,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye
      (2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput(1),mslInput(2))*Ye(1,1)
      *Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - Sqr(Abs(
      MuInput))*TCC0(mseInput(1),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*(Ye(0,1)
      *Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(
      mseInput(2),mslInput(1),mslInput(2))*Ye(2,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1
      ,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput
      (2),mslInput(1))*Ye(2,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1
      )*Ye(2,2)) - AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(2),
      mslInput(1),mslInput(2))*Ye(2,1)*Ye(2,2)*(AeInput(2,1)*Ye(2,1)*Ye(2,2) +
      AeInput(2,2)*Ye(2,1)*Ye(2,2)) - AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput
      (2),mseInput(2),mslInput(2),mslInput(1))*Ye(2,1)*Ye(2,2)*(AeInput(2,1)*Ye(2,
      1)*Ye(2,2) + AeInput(2,2)*Ye(2,1)*Ye(2,2)) - 3*AuInput(0,1)*Sqr(Abs(MuInput)
      )*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(0))*Yu(0,0)*Yu(0,1)*(
      AuInput(0,0)*Yu(0,0)*Yu(0,1) + AuInput(0,1)*Yu(0,0)*Yu(0,1)) - 3*AuInput(0,0
      )*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(0),msuInput(0))*Yu
      (0,0)*Yu(0,1)*(AuInput(0,0)*Yu(0,0)*Yu(0,1) + AuInput(0,1)*Yu(0,0)*Yu(0,1))
      - 3*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msuInput(0),
      msuInput(0))*Yu(0,0)*Yu(0,2)*(AuInput(0,0)*Yu(0,0)*Yu(0,2) + AuInput(0,2)*Yu
      (0,0)*Yu(0,2)) - 3*AuInput(0,0)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(
      0),msuInput(0),msuInput(0))*Yu(0,0)*Yu(0,2)*(AuInput(0,0)*Yu(0,0)*Yu(0,2) +
      AuInput(0,2)*Yu(0,0)*Yu(0,2)) - 3*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(2),msuInput(0),msuInput(0))*Yu(0,1)*Yu(0,2)*(AuInput(0,
      1)*Yu(0,1)*Yu(0,2) + AuInput(0,2)*Yu(0,1)*Yu(0,2)) - 3*AuInput(0,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(0),msuInput(0))*Yu(0,1)*Yu(0
      ,2)*(AuInput(0,1)*Yu(0,1)*Yu(0,2) + AuInput(0,2)*Yu(0,1)*Yu(0,2)) - 3*
      AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(0),msuInput(0),
      msuInput(1))*Yu(0,0)*Yu(1,0)*(AuInput(0,0)*Yu(0,0)*Yu(1,0) + AuInput(1,0)*Yu
      (0,0)*Yu(1,0)) - 3*AuInput(0,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(
      0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(1,0)*(AuInput(0,0)*Yu(0,0)*Yu(1,0) +
      AuInput(1,0)*Yu(0,0)*Yu(1,0)) - 3*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(0),msqInput(1),msuInput(0),msuInput(1))*Yu(0,0)*(AuInput(0,1)*Yu(0,
      1)*Yu(1,0) + AuInput(1,0)*Yu(0,1)*Yu(1,0))*Yu(1,1) - 3*AuInput(0,0)*Sqr(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*(
      AuInput(0,1)*Yu(0,1)*Yu(1,0) + AuInput(1,0)*Yu(0,1)*Yu(1,0))*Yu(1,1) - 3*
      AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(1),
      msuInput(0))*Yu(0,1)*Yu(1,0)*(AuInput(0,0)*Yu(0,0)*Yu(1,1) + AuInput(1,1)*Yu
      (0,0)*Yu(1,1)) - 3*AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(
      0),msuInput(0),msuInput(1))*Yu(0,1)*Yu(1,0)*(AuInput(0,0)*Yu(0,0)*Yu(1,1) +
      AuInput(1,1)*Yu(0,0)*Yu(1,1)) - 3*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(1),msuInput(0),msuInput(1))*Yu(0,1)*Yu(1,1)*(AuInput(0,
      1)*Yu(0,1)*Yu(1,1) + AuInput(1,1)*Yu(0,1)*Yu(1,1)) - 3*AuInput(0,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(0))*Yu(0,1)*Yu(1
      ,1)*(AuInput(0,1)*Yu(0,1)*Yu(1,1) + AuInput(1,1)*Yu(0,1)*Yu(1,1)) - 3*
      AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(1),
      msuInput(1))*Yu(1,0)*Yu(1,1)*(AuInput(1,0)*Yu(1,0)*Yu(1,1) + AuInput(1,1)*Yu
      (1,0)*Yu(1,1)) - 3*AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(
      0),msuInput(1),msuInput(1))*Yu(1,0)*Yu(1,1)*(AuInput(1,0)*Yu(1,0)*Yu(1,1) +
      AuInput(1,1)*Yu(1,0)*Yu(1,1)) - 3*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(0),msqInput(2),msuInput(0),msuInput(1))*Yu(0,0)*(AuInput(0,2)*Yu(0,
      2)*Yu(1,0) + AuInput(1,0)*Yu(0,2)*Yu(1,0))*Yu(1,2) - 3*AuInput(0,0)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*(
      AuInput(0,2)*Yu(0,2)*Yu(1,0) + AuInput(1,0)*Yu(0,2)*Yu(1,0))*Yu(1,2) - 3*
      AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(0),
      msuInput(1))*Yu(0,1)*(AuInput(0,2)*Yu(0,2)*Yu(1,1) + AuInput(1,1)*Yu(0,2)*Yu
      (1,1))*Yu(1,2) - 3*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(
      1),msuInput(1),msuInput(0))*Yu(0,1)*(AuInput(0,2)*Yu(0,2)*Yu(1,1) + AuInput(
      1,1)*Yu(0,2)*Yu(1,1))*Yu(1,2) - 3*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(0),msqInput(2),msuInput(1),msuInput(0))*Yu(0,2)*Yu(1,0)*(AuInput(0,
      0)*Yu(0,0)*Yu(1,2) + AuInput(1,2)*Yu(0,0)*Yu(1,2)) - 3*AuInput(1,0)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(0),msuInput(1))*Yu(0,2)*Yu(1
      ,0)*(AuInput(0,0)*Yu(0,0)*Yu(1,2) + AuInput(1,2)*Yu(0,0)*Yu(1,2)) - 3*
      AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(1),
      msuInput(0))*Yu(0,2)*Yu(1,1)*(AuInput(0,1)*Yu(0,1)*Yu(1,2) + AuInput(1,2)*Yu
      (0,1)*Yu(1,2)) - 3*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(
      1),msuInput(0),msuInput(1))*Yu(0,2)*Yu(1,1)*(AuInput(0,1)*Yu(0,1)*Yu(1,2) +
      AuInput(1,2)*Yu(0,1)*Yu(1,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(0),
      msuInput(0),msuInput(1))*Yu(0,0)*Yu(1,0)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1)
      + Yu(0,2)*Yu(1,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(0),msuInput(1),
      msuInput(0))*Yu(0,0)*Yu(1,0)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu
      (1,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(1),msuInput(0),msuInput(1))*Yu(0,
      1)*Yu(1,1)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 3*Sqr(Abs
      (MuInput))*TCC0(msqInput(1),msuInput(1),msuInput(0))*Yu(0,1)*Yu(1,1)*(Yu(0,0
      )*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 3*Sqr(Abs(MuInput))*TCC0(
      msqInput(2),msuInput(0),msuInput(1))*Yu(0,2)*Yu(1,2)*(Yu(0,0)*Yu(1,0) + Yu(0
      ,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),
      msuInput(1),msuInput(0))*Yu(0,2)*Yu(1,2)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1)
      + Yu(0,2)*Yu(1,2)) - 3*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),
      msqInput(2),msuInput(0),msuInput(1))*Yu(0,2)*Yu(1,2)*(AuInput(0,2)*Yu(0,2)*
      Yu(1,2) + AuInput(1,2)*Yu(0,2)*Yu(1,2)) - 3*AuInput(0,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(2),msqInput(2),msuInput(1),msuInput(0))*Yu(0,2)*Yu(1,2)*(
      AuInput(0,2)*Yu(0,2)*Yu(1,2) + AuInput(1,2)*Yu(0,2)*Yu(1,2)) - 3*AuInput(1,2
      )*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(1))*Yu
      (1,0)*Yu(1,2)*(AuInput(1,0)*Yu(1,0)*Yu(1,2) + AuInput(1,2)*Yu(1,0)*Yu(1,2))
      - 3*AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(1),
      msuInput(1))*Yu(1,0)*Yu(1,2)*(AuInput(1,0)*Yu(1,0)*Yu(1,2) + AuInput(1,2)*Yu
      (1,0)*Yu(1,2)) - 3*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(
      2),msuInput(1),msuInput(1))*Yu(1,1)*Yu(1,2)*(AuInput(1,1)*Yu(1,1)*Yu(1,2) +
      AuInput(1,2)*Yu(1,1)*Yu(1,2)) - 3*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(1),msuInput(1),msuInput(1))*Yu(1,1)*Yu(1,2)*(AuInput(1,
      1)*Yu(1,1)*Yu(1,2) + AuInput(1,2)*Yu(1,1)*Yu(1,2)) - 3*AuInput(2,0)*Sqr(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*Yu(2
      ,0)*(AuInput(0,0)*Yu(0,0)*Yu(2,0) + AuInput(2,0)*Yu(0,0)*Yu(2,0)) - 3*
      AuInput(0,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(0),msuInput(2),
      msuInput(0))*Yu(0,0)*Yu(2,0)*(AuInput(0,0)*Yu(0,0)*Yu(2,0) + AuInput(2,0)*Yu
      (0,0)*Yu(2,0)) - 3*AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(
      0),msuInput(1),msuInput(2))*Yu(1,0)*Yu(2,0)*(AuInput(1,0)*Yu(1,0)*Yu(2,0) +
      AuInput(2,0)*Yu(1,0)*Yu(2,0)) - 3*AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(
      msqInput(0),msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(2,0)*(AuInput(1,
      0)*Yu(1,0)*Yu(2,0) + AuInput(2,0)*Yu(1,0)*Yu(2,0)) - 3*AuInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(2))*Yu(0,0)*(
      AuInput(0,1)*Yu(0,1)*Yu(2,0) + AuInput(2,0)*Yu(0,1)*Yu(2,0))*Yu(2,1) - 3*
      AuInput(0,0)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(2),
      msuInput(0))*Yu(0,0)*(AuInput(0,1)*Yu(0,1)*Yu(2,0) + AuInput(2,0)*Yu(0,1)*Yu
      (2,0))*Yu(2,1) - 3*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(
      1),msuInput(1),msuInput(2))*Yu(1,0)*(AuInput(1,1)*Yu(1,1)*Yu(2,0) + AuInput(
      2,0)*Yu(1,1)*Yu(2,0))*Yu(2,1) - 3*AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*(AuInput(1,1)*Yu(1,
      1)*Yu(2,0) + AuInput(2,0)*Yu(1,1)*Yu(2,0))*Yu(2,1) - 3*AuInput(0,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(0))*Yu(0,1)*Yu(2
      ,0)*(AuInput(0,0)*Yu(0,0)*Yu(2,1) + AuInput(2,1)*Yu(0,0)*Yu(2,1)) - 3*
      AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(0),
      msuInput(2))*Yu(0,1)*Yu(2,0)*(AuInput(0,0)*Yu(0,0)*Yu(2,1) + AuInput(2,1)*Yu
      (0,0)*Yu(2,1)) - 3*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(
      1),msuInput(0),msuInput(2))*Yu(0,1)*Yu(2,1)*(AuInput(0,1)*Yu(0,1)*Yu(2,1) +
      AuInput(2,1)*Yu(0,1)*Yu(2,1)) - 3*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(1),msuInput(2),msuInput(0))*Yu(0,1)*Yu(2,1)*(AuInput(0,
      1)*Yu(0,1)*Yu(2,1) + AuInput(2,1)*Yu(0,1)*Yu(2,1)) - 3*AuInput(1,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(2
      ,0)*(AuInput(1,0)*Yu(1,0)*Yu(2,1) + AuInput(2,1)*Yu(1,0)*Yu(2,1)) - 3*
      AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(1),
      msuInput(2))*Yu(1,1)*Yu(2,0)*(AuInput(1,0)*Yu(1,0)*Yu(2,1) + AuInput(2,1)*Yu
      (1,0)*Yu(2,1)) - 3*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(
      1),msuInput(1),msuInput(2))*Yu(1,1)*Yu(2,1)*(AuInput(1,1)*Yu(1,1)*Yu(2,1) +
      AuInput(2,1)*Yu(1,1)*Yu(2,1)) - 3*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(2,1)*(AuInput(1,
      1)*Yu(1,1)*Yu(2,1) + AuInput(2,1)*Yu(1,1)*Yu(2,1)) - 3*Sqr(Abs(MuInput))*
      TCC0(msqInput(0),msqInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*(Yu(0,0)*Yu(0,1) +
      Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(1),
      msqInput(0),msuInput(0))*Yu(0,0)*Yu(0,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1)
      + Yu(2,0)*Yu(2,1)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(0),msqInput(1),
      msuInput(1))*Yu(1,0)*Yu(1,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu
      (2,1)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(1),msqInput(0),msuInput(1))*Yu(1,
      0)*Yu(1,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 3*Sqr(Abs
      (MuInput))*TCC0(msqInput(0),msqInput(1),msuInput(2))*Yu(2,0)*Yu(2,1)*(Yu(0,0
      )*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 3*Sqr(Abs(MuInput))*TCC0(
      msqInput(1),msqInput(0),msuInput(2))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1
      ,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) - 3*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(0),msqInput(1),msuInput(2),msuInput(2))*Yu(2,0)*Yu(2,1)*(AuInput(2,
      0)*Yu(2,0)*Yu(2,1) + AuInput(2,1)*Yu(2,0)*Yu(2,1)) - 3*AuInput(2,0)*Sqr(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(0),msuInput(2),msuInput(2))*Yu(2,0)*Yu(2
      ,1)*(AuInput(2,0)*Yu(2,0)*Yu(2,1) + AuInput(2,1)*Yu(2,0)*Yu(2,1)) - 3*
      AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msuInput(0),
      msuInput(2))*Yu(0,0)*(AuInput(0,2)*Yu(0,2)*Yu(2,0) + AuInput(2,0)*Yu(0,2)*Yu
      (2,0))*Yu(2,2) - 3*AuInput(0,0)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(
      0),msuInput(2),msuInput(0))*Yu(0,0)*(AuInput(0,2)*Yu(0,2)*Yu(2,0) + AuInput(
      2,0)*Yu(0,2)*Yu(2,0))*Yu(2,2) - 3*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(0),msqInput(2),msuInput(1),msuInput(2))*Yu(1,0)*(AuInput(1,2)*Yu(1,
      2)*Yu(2,0) + AuInput(2,0)*Yu(1,2)*Yu(2,0))*Yu(2,2) - 3*AuInput(1,0)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*(
      AuInput(1,2)*Yu(1,2)*Yu(2,0) + AuInput(2,0)*Yu(1,2)*Yu(2,0))*Yu(2,2) - 3*
      AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(0),
      msuInput(2))*Yu(0,1)*(AuInput(0,2)*Yu(0,2)*Yu(2,1) + AuInput(2,1)*Yu(0,2)*Yu
      (2,1))*Yu(2,2) - 3*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(
      1),msuInput(2),msuInput(0))*Yu(0,1)*(AuInput(0,2)*Yu(0,2)*Yu(2,1) + AuInput(
      2,1)*Yu(0,2)*Yu(2,1))*Yu(2,2) - 3*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(2),msuInput(1),msuInput(2))*Yu(1,1)*(AuInput(1,2)*Yu(1,
      2)*Yu(2,1) + AuInput(2,1)*Yu(1,2)*Yu(2,1))*Yu(2,2) - 3*AuInput(1,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*(
      AuInput(1,2)*Yu(1,2)*Yu(2,1) + AuInput(2,1)*Yu(1,2)*Yu(2,1))*Yu(2,2) - 3*
      AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msuInput(2),
      msuInput(0))*Yu(0,2)*Yu(2,0)*(AuInput(0,0)*Yu(0,0)*Yu(2,2) + AuInput(2,2)*Yu
      (0,0)*Yu(2,2)) - 3*AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(
      0),msuInput(0),msuInput(2))*Yu(0,2)*Yu(2,0)*(AuInput(0,0)*Yu(0,0)*Yu(2,2) +
      AuInput(2,2)*Yu(0,0)*Yu(2,2)) - 3*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(2),msuInput(2),msuInput(0))*Yu(0,2)*Yu(2,1)*(AuInput(0,
      1)*Yu(0,1)*Yu(2,2) + AuInput(2,2)*Yu(0,1)*Yu(2,2)) - 3*AuInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(0),msuInput(2))*Yu(0,2)*Yu(2
      ,1)*(AuInput(0,1)*Yu(0,1)*Yu(2,2) + AuInput(2,2)*Yu(0,1)*Yu(2,2)) - 3*Sqr(
      Abs(MuInput))*TCC0(msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*Yu(2,0)*(Yu(
      0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0
      (msqInput(0),msuInput(2),msuInput(0))*Yu(0,0)*Yu(2,0)*(Yu(0,0)*Yu(2,0) + Yu(
      0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(1),
      msuInput(0),msuInput(2))*Yu(0,1)*Yu(2,1)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1)
      + Yu(0,2)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(1),msuInput(2),
      msuInput(0))*Yu(0,1)*Yu(2,1)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),msuInput(0),msuInput(2))*Yu(0,
      2)*Yu(2,2)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 3*Sqr(Abs
      (MuInput))*TCC0(msqInput(2),msuInput(2),msuInput(0))*Yu(0,2)*Yu(2,2)*(Yu(0,0
      )*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) - 3*AuInput(2,2)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput(2))*Yu(0,2)*Yu(2
      ,2)*(AuInput(0,2)*Yu(0,2)*Yu(2,2) + AuInput(2,2)*Yu(0,2)*Yu(2,2)) - 3*
      AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(2),msuInput(2),
      msuInput(0))*Yu(0,2)*Yu(2,2)*(AuInput(0,2)*Yu(0,2)*Yu(2,2) + AuInput(2,2)*Yu
      (0,2)*Yu(2,2)) - 3*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(
      2),msuInput(2),msuInput(1))*Yu(1,2)*Yu(2,0)*(AuInput(1,0)*Yu(1,0)*Yu(2,2) +
      AuInput(2,2)*Yu(1,0)*Yu(2,2)) - 3*AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(0),msuInput(1),msuInput(2))*Yu(1,2)*Yu(2,0)*(AuInput(1,
      0)*Yu(1,0)*Yu(2,2) + AuInput(2,2)*Yu(1,0)*Yu(2,2)) - 3*AuInput(1,2)*Sqr(Abs(
      MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(1))*Yu(1,2)*Yu(2
      ,1)*(AuInput(1,1)*Yu(1,1)*Yu(2,2) + AuInput(2,2)*Yu(1,1)*Yu(2,2)) - 3*
      AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(1),
      msuInput(2))*Yu(1,2)*Yu(2,1)*(AuInput(1,1)*Yu(1,1)*Yu(2,2) + AuInput(2,2)*Yu
      (1,1)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(0),msuInput(1),msuInput(2
      ))*Yu(1,0)*Yu(2,0)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 3
      *Sqr(Abs(MuInput))*TCC0(msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(2,0)
      *(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 3*Sqr(Abs(MuInput))
      *TCC0(msqInput(1),msuInput(1),msuInput(2))*Yu(1,1)*Yu(2,1)*(Yu(1,0)*Yu(2,0)
      + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(1),
      msuInput(2),msuInput(1))*Yu(1,1)*Yu(2,1)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1)
      + Yu(1,2)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),msuInput(1),
      msuInput(2))*Yu(1,2)*Yu(2,2)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),msuInput(2),msuInput(1))*Yu(1,
      2)*Yu(2,2)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) - 3*AuInput
      (2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(2),msuInput(1),msuInput(2)
      )*Yu(1,2)*Yu(2,2)*(AuInput(1,2)*Yu(1,2)*Yu(2,2) + AuInput(2,2)*Yu(1,2)*Yu(2,
      2)) - 3*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(2),msuInput
      (2),msuInput(1))*Yu(1,2)*Yu(2,2)*(AuInput(1,2)*Yu(1,2)*Yu(2,2) + AuInput(2,2
      )*Yu(1,2)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(0),msqInput(2),
      msuInput(0))*Yu(0,0)*Yu(0,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),msqInput(0),msuInput(0))*Yu(0,
      0)*Yu(0,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) - 3*Sqr(Abs
      (MuInput))*TCC0(msqInput(0),msqInput(2),msuInput(1))*Yu(1,0)*Yu(1,2)*(Yu(0,0
      )*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(
      msqInput(2),msqInput(0),msuInput(1))*Yu(1,0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1
      ,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(0),
      msqInput(2),msuInput(2))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2)
      + Yu(2,0)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),msqInput(0),
      msuInput(2))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu
      (2,2)) - 3*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),
      msuInput(2),msuInput(2))*Yu(2,0)*Yu(2,2)*(AuInput(2,0)*Yu(2,0)*Yu(2,2) +
      AuInput(2,2)*Yu(2,0)*Yu(2,2)) - 3*AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(0),msuInput(2),msuInput(2))*Yu(2,0)*Yu(2,2)*(AuInput(2,
      0)*Yu(2,0)*Yu(2,2) + AuInput(2,2)*Yu(2,0)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*
      TCC0(msqInput(1),msqInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2) +
      Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),
      msqInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2)
      + Yu(2,1)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(1),msqInput(2),
      msuInput(1))*Yu(1,1)*Yu(1,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu
      (2,2)) - 3*Sqr(Abs(MuInput))*TCC0(msqInput(2),msqInput(1),msuInput(1))*Yu(1,
      1)*Yu(1,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 3*Sqr(Abs
      (MuInput))*TCC0(msqInput(1),msqInput(2),msuInput(2))*Yu(2,1)*Yu(2,2)*(Yu(0,1
      )*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 3*Sqr(Abs(MuInput))*TCC0(
      msqInput(2),msqInput(1),msuInput(2))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1
      ,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) - 3*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(2),msuInput(2),msuInput(2))*Yu(2,1)*Yu(2,2)*(AuInput(2,
      1)*Yu(2,1)*Yu(2,2) + AuInput(2,2)*Yu(2,1)*Yu(2,2)) - 3*AuInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(2))*Yu(2,1)*Yu(2
      ,2)*(AuInput(2,1)*Yu(2,1)*Yu(2,2) + AuInput(2,2)*Yu(2,1)*Yu(2,2))) +
      0.006332573977646111*(-3*Sqr(Yd(0,0)*Yu(0,0) + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(
      0,2))*TCB0(msdInput(0),msuInput(0),SCALE) - 3*Sqr(Yd(0,0)*Yu(1,0) + Yd(0,1)*
      Yu(1,1) + Yd(0,2)*Yu(1,2))*TCB0(msdInput(0),msuInput(1),SCALE) - 3*Sqr(Yd(0,
      0)*Yu(2,0) + Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu(2,2))*TCB0(msdInput(0),msuInput(2)
      ,SCALE) - 3*Sqr(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2))*TCB0(
      msdInput(1),msuInput(0),SCALE) - 3*Sqr(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) +
      Yd(1,2)*Yu(1,2))*TCB0(msdInput(1),msuInput(1),SCALE) - 3*Sqr(Yd(1,0)*Yu(2,0)
      + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2))*TCB0(msdInput(1),msuInput(2),SCALE) - 3
      *Sqr(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu(0,2))*TCB0(msdInput(2),
      msuInput(0),SCALE) - 3*Sqr(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu(1,
      2))*TCB0(msdInput(2),msuInput(1),SCALE) - 3*Sqr(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu
      (2,1) + Yd(2,2)*Yu(2,2))*TCB0(msdInput(2),msuInput(2),SCALE) + (-0.25*Quad(
      g2) + 0.5*Sqr(g2)*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCB0(
      mslInput(0),mslInput(0),SCALE) + (-0.25*Quad(g2) + 0.5*Sqr(g2)*(Sqr(Ye(0,1))
      + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCB0(mslInput(1),mslInput(1),SCALE) + (-0.25
      *Quad(g2) + 0.5*Sqr(g2)*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCB0(
      mslInput(2),mslInput(2),SCALE) + (-0.75*Quad(g2) + 1.5*Sqr(g2)*(Sqr(Yd(0,0))
      + Sqr(Yd(1,0)) + Sqr(Yd(2,0))) + 1.5*Sqr(g2)*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) +
      Sqr(Yu(2,0))) - 3*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0)))*(Sqr(Yu(0,0))
      + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCB0(msqInput(0),msqInput(0),SCALE) + (-0.75
      *Quad(g2) + 1.5*Sqr(g2)*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) + 1.5*
      Sqr(g2)*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))) - 3*(Sqr(Yd(0,1)) + Sqr
      (Yd(1,1)) + Sqr(Yd(2,1)))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCB0
      (msqInput(1),msqInput(1),SCALE) + (-0.75*Quad(g2) + 1.5*Sqr(g2)*(Sqr(Yd(0,2)
      ) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) + 1.5*Sqr(g2)*(Sqr(Yu(0,2)) + Sqr(Yu(1,2))
      + Sqr(Yu(2,2))) - 3*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*(Sqr(Yu(0,2
      )) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCB0(msqInput(2),msqInput(2),SCALE) + (-
      1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Yd(0,0)) + 1.5*Sqr(g2)*Sqr(AdInput(0,0))*
      Sqr(Yd(0,0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0))
      + Sqr(Yd(2,0))) - 3*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,
      0)) + Sqr(Yu(2,0))))*TCC0(msdInput(0),msqInput(0),msqInput(0)) + (-1.5*Sqr(
      g2)*Sqr(Abs(MuInput))*Sqr(Yd(0,1)) + 1.5*Sqr(g2)*Sqr(AdInput(0,1))*Sqr(Yd(0,
      1)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd
      (2,1))) - 3*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) +
      Sqr(Yu(2,1))))*TCC0(msdInput(0),msqInput(1),msqInput(1)) + (-1.5*Sqr(g2)*Sqr
      (Abs(MuInput))*Sqr(Yd(0,2)) + 1.5*Sqr(g2)*Sqr(AdInput(0,2))*Sqr(Yd(0,2)) + 3
      *Sqr(Abs(MuInput))*Sqr(Yd(0,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))
      - 3*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2
      ))))*TCC0(msdInput(0),msqInput(2),msqInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(
      MuInput))*Sqr(Yd(1,0)) + 1.5*Sqr(g2)*Sqr(AdInput(1,0))*Sqr(Yd(1,0)) + 3*Sqr(
      Abs(MuInput))*Sqr(Yd(1,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))) - 3*
      Sqr(AdInput(1,0))*Sqr(Yd(1,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))
      *TCC0(msdInput(1),msqInput(0),msqInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))
      *Sqr(Yd(1,1)) + 1.5*Sqr(g2)*Sqr(AdInput(1,1))*Sqr(Yd(1,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(1,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) - 3*Sqr(
      AdInput(1,1))*Sqr(Yd(1,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*
      TCC0(msdInput(1),msqInput(1),msqInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(1,2)) + 1.5*Sqr(g2)*Sqr(AdInput(1,2))*Sqr(Yd(1,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(1,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) - 3*Sqr(
      AdInput(1,2))*Sqr(Yd(1,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*
      TCC0(msdInput(1),msqInput(2),msqInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,0)) + 1.5*Sqr(g2)*Sqr(AdInput(2,0))*Sqr(Yd(2,0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))) - 3*Sqr(
      AdInput(2,0))*Sqr(Yd(2,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*
      TCC0(msdInput(2),msqInput(0),msqInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,1)) + 1.5*Sqr(g2)*Sqr(AdInput(2,1))*Sqr(Yd(2,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) - 3*Sqr(
      AdInput(2,1))*Sqr(Yd(2,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*
      TCC0(msdInput(2),msqInput(1),msqInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,2)) + 1.5*Sqr(g2)*Sqr(AdInput(2,2))*Sqr(Yd(2,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) - 3*Sqr(
      AdInput(2,2))*Sqr(Yd(2,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*
      TCC0(msdInput(2),msqInput(2),msqInput(2)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Ye(0,0)) + 0.5*Sqr(g2)*Sqr(AeInput(0,0))*Sqr(Ye(0,0)) + Sqr(Abs(MuInput)
      )*Sqr(Ye(0,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(0
      ),mslInput(0),mslInput(0)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(0,1)) +
      0.5*Sqr(g2)*Sqr(AeInput(0,1))*Sqr(Ye(0,1)) + Sqr(Abs(MuInput))*Sqr(Ye(0,1))*
      (Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(0),mslInput(1),
      mslInput(1)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(0,2)) + 0.5*Sqr(g2)*
      Sqr(AeInput(0,2))*Sqr(Ye(0,2)) + Sqr(Abs(MuInput))*Sqr(Ye(0,2))*(Sqr(Ye(0,2)
      ) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(0),mslInput(2),mslInput(2))
      + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(1,0)) + 0.5*Sqr(g2)*Sqr(AeInput(1,0
      ))*Sqr(Ye(1,0)) + Sqr(Abs(MuInput))*Sqr(Ye(1,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)
      ) + Sqr(Ye(2,0))))*TCC0(mseInput(1),mslInput(0),mslInput(0)) + (-0.5*Sqr(g2)
      *Sqr(Abs(MuInput))*Sqr(Ye(1,1)) + 0.5*Sqr(g2)*Sqr(AeInput(1,1))*Sqr(Ye(1,1))
      + Sqr(Abs(MuInput))*Sqr(Ye(1,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))
      ))*TCC0(mseInput(1),mslInput(1),mslInput(1)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput
      ))*Sqr(Ye(1,2)) + 0.5*Sqr(g2)*Sqr(AeInput(1,2))*Sqr(Ye(1,2)) + Sqr(Abs(
      MuInput))*Sqr(Ye(1,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(
      mseInput(1),mslInput(2),mslInput(2)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Ye(2,0)) + 0.5*Sqr(g2)*Sqr(AeInput(2,0))*Sqr(Ye(2,0)) + Sqr(Abs(MuInput))*
      Sqr(Ye(2,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(2),
      mslInput(0),mslInput(0)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(2,1)) +
      0.5*Sqr(g2)*Sqr(AeInput(2,1))*Sqr(Ye(2,1)) + Sqr(Abs(MuInput))*Sqr(Ye(2,1))*
      (Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(2),mslInput(1),
      mslInput(1)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(2,2)) + 0.5*Sqr(g2)*
      Sqr(AeInput(2,2))*Sqr(Ye(2,2)) + Sqr(Abs(MuInput))*Sqr(Ye(2,2))*(Sqr(Ye(0,2)
      ) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(2),mslInput(2),mslInput(2))
      + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Yu(0,0)) + 1.5*Sqr(g2)*Sqr(AuInput(0,0
      ))*Sqr(Yu(0,0)) - 3*Sqr(AuInput(0,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(
      2,0)))*Sqr(Yu(0,0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*(Sqr(Yu(0,0)) + Sqr(
      Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(0)) + (-1.5*
      Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Yu(1,0)) + 1.5*Sqr(g2)*Sqr(AuInput(1,0))*Sqr(
      Yu(1,0)) - 3*Sqr(AuInput(1,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0)))*
      Sqr(Yu(1,0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0))
      + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(1)) + (-1.5*Sqr(g2)*
      Sqr(Abs(MuInput))*Sqr(Yu(2,0)) + 1.5*Sqr(g2)*Sqr(AuInput(2,0))*Sqr(Yu(2,0))
      - 3*Sqr(AuInput(2,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0)))*Sqr(Yu(2,
      0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu
      (2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(
      MuInput))*Sqr(Yu(0,1)) + 1.5*Sqr(g2)*Sqr(AuInput(0,1))*Sqr(Yu(0,1)) - 3*Sqr(
      AuInput(0,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1)))*Sqr(Yu(0,1)) + 3*
      Sqr(Abs(MuInput))*Sqr(Yu(0,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))
      *TCC0(msqInput(1),msqInput(1),msuInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))
      *Sqr(Yu(1,1)) + 1.5*Sqr(g2)*Sqr(AuInput(1,1))*Sqr(Yu(1,1)) - 3*Sqr(AuInput(1
      ,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1)))*Sqr(Yu(1,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(1,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(
      msqInput(1),msqInput(1),msuInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(2,1)) + 1.5*Sqr(g2)*Sqr(AuInput(2,1))*Sqr(Yu(2,1)) - 3*Sqr(AuInput(2,1))*
      (Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1)))*Sqr(Yu(2,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(2,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(
      msqInput(1),msqInput(1),msuInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(0,2)) + 1.5*Sqr(g2)*Sqr(AuInput(0,2))*Sqr(Yu(0,2)) - 3*Sqr(AuInput(0,2))*
      (Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*Sqr(Yu(0,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(0,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(
      msqInput(2),msqInput(2),msuInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(1,2)) + 1.5*Sqr(g2)*Sqr(AuInput(1,2))*Sqr(Yu(1,2)) - 3*Sqr(AuInput(1,2))*
      (Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*Sqr(Yu(1,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(1,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(
      msqInput(2),msqInput(2),msuInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(2,2)) + 1.5*Sqr(g2)*Sqr(AuInput(2,2))*Sqr(Yu(2,2)) - 3*Sqr(AuInput(2,2))*
      (Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*Sqr(Yu(2,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(2,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(
      msqInput(2),msqInput(2),msuInput(2)) + 3*Quad(Yd(0,0))*Sqr(Abs(MuInput))*Sqr
      (AdInput(0,0))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(0)) + 3*
      AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(
      msdInput(0),msdInput(0),msqInput(0),msqInput(1)) + 3*AdInput(0,0)*AdInput(0,
      2)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),
      msqInput(0),msqInput(2)) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*Sqr
      (Yd(0,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1),msqInput(0))
      + 3*Quad(Yd(0,1))*Sqr(Abs(MuInput))*Sqr(AdInput(0,1))*TCD0(msdInput(0),
      msdInput(0),msqInput(1),msqInput(1)) + 3*AdInput(0,1)*AdInput(0,2)*Sqr(Abs(
      MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(1)
      ,msqInput(2)) + 3*AdInput(0,0)*AdInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*
      Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),msqInput(0)) + 3*
      AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(
      msdInput(0),msdInput(0),msqInput(2),msqInput(1)) + 3*Quad(Yd(0,2))*Sqr(Abs(
      MuInput))*Sqr(AdInput(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),
      msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*Sqr(Yd(1,0
      ))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(0),msdInput(1),
      msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,2))*Sqr(Yd(0,2)
      )*Sqr(Yd(1,2))*TCD0(msdInput(0),msdInput(1),msqInput(2),msqInput(2)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(0),
      msdInput(2),msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,1))
      *Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(0),msdInput(2),msqInput(1),msqInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(0),msdInput(2),msqInput(2),msqInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))
      *Yd(0,0)*Yu(0,0)) + AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0))*TCD0(msdInput
      (0),msqInput(0),msqInput(0),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,0)
      *Yu(1,0)) + AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0))*TCD0(msdInput(0),
      msqInput(0),msqInput(0),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,0)*Yu(
      2,0)) + AdInput(0,0)*AuInput(2,0)*Yd(0,0)*Yu(2,0))*TCD0(msdInput(0),msqInput
      (0),msqInput(0),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,1)*Yu(0,1)) +
      AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1))*TCD0(msdInput(0),msqInput(1),
      msqInput(1),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)) +
      AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1))*TCD0(msdInput(0),msqInput(1),
      msqInput(1),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)) +
      AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu(2,1))*TCD0(msdInput(0),msqInput(1),
      msqInput(1),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,2)*Yu(0,2)) +
      AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2))*TCD0(msdInput(0),msqInput(2),
      msqInput(2),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) +
      AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2))*TCD0(msdInput(0),msqInput(2),
      msqInput(2),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) +
      AdInput(0,2)*AuInput(2,2)*Yd(0,2)*Yu(2,2))*TCD0(msdInput(0),msqInput(2),
      msqInput(2),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(1,0))*Sqr(Yd(0,0)
      )*Sqr(Yd(1,0))*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(0)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(1,1))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(1),
      msdInput(0),msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(1,2))
      *Sqr(Yd(0,2))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput
      (2)) + 3*Quad(Yd(1,0))*Sqr(Abs(MuInput))*Sqr(AdInput(1,0))*TCD0(msdInput(1),
      msdInput(1),msqInput(0),msqInput(0)) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(
      MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(0)
      ,msqInput(1)) + 3*AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*
      Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(0),msqInput(2)) + 3*
      AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(
      msdInput(1),msdInput(1),msqInput(1),msqInput(0)) + 3*Quad(Yd(1,1))*Sqr(Abs(
      MuInput))*Sqr(AdInput(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(1),
      msqInput(1)) + 3*AdInput(1,1)*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*
      Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(2)) + 3*
      AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(
      msdInput(1),msdInput(1),msqInput(2),msqInput(0)) + 3*AdInput(1,1)*AdInput(1,
      2)*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),
      msqInput(2),msqInput(1)) + 3*Quad(Yd(1,2))*Sqr(Abs(MuInput))*Sqr(AdInput(1,2
      ))*TCD0(msdInput(1),msdInput(1),msqInput(2),msqInput(2)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AdInput(1,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(1),msdInput(2),
      msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(1,1))*Sqr(Yd(1,1)
      )*Sqr(Yd(2,1))*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(1)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(1,2))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,0)*Yu(
      0,0)) + AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0))*TCD0(msdInput(1),msqInput
      (0),msqInput(0),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,0)*Yu(1,0)) +
      AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0))*TCD0(msdInput(1),msqInput(0),
      msqInput(0),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0)) +
      AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0))*TCD0(msdInput(1),msqInput(0),
      msqInput(0),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,1)*Yu(0,1)) +
      AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1))*TCD0(msdInput(1),msqInput(1),
      msqInput(1),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,1)*Yu(1,1)) +
      AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1))*TCD0(msdInput(1),msqInput(1),
      msqInput(1),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,1)*Yu(2,1)) +
      AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1))*TCD0(msdInput(1),msqInput(1),
      msqInput(1),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,2)*Yu(0,2)) +
      AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2))*TCD0(msdInput(1),msqInput(2),
      msqInput(2),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,2)*Yu(1,2)) +
      AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2))*TCD0(msdInput(1),msqInput(2),
      msqInput(2),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,2)*Yu(2,2)) +
      AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2))*TCD0(msdInput(1),msqInput(2),
      msqInput(2),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,0))*Sqr(Yd(0,0)
      )*Sqr(Yd(2,0))*TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(0)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(2,1))*Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(0),msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,2))
      *Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput
      (2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(
      msdInput(2),msdInput(1),msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(
      AdInput(2,1))*Sqr(Yd(1,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(1),
      msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,2))*Sqr(Yd(1,2)
      )*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(2)) + 3*
      Quad(Yd(2,0))*Sqr(Abs(MuInput))*Sqr(AdInput(2,0))*TCD0(msdInput(2),msdInput(
      2),msqInput(0),msqInput(0)) + 3*AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(
      1)) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,2)
      )*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(2)) + 3*AdInput(2,0)*
      AdInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(2),msqInput(1),msqInput(0)) + 3*Quad(Yd(2,1))*Sqr(Abs(MuInput))*Sqr
      (AdInput(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(1)) + 3*
      AdInput(2,1)*AdInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(2),msqInput(1),msqInput(2)) + 3*AdInput(2,0)*AdInput(2,
      2)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(2),msqInput(0)) + 3*AdInput(2,1)*AdInput(2,2)*Sqr(Abs(MuInput))*Sqr
      (Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(1))
      + 3*Quad(Yd(2,2))*Sqr(Abs(MuInput))*Sqr(AdInput(2,2))*TCD0(msdInput(2),
      msdInput(2),msqInput(2),msqInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,0)*Yu(
      0,0)) + AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0))*TCD0(msdInput(2),msqInput
      (0),msqInput(0),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,0)*Yu(1,0)) +
      AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0))*TCD0(msdInput(2),msqInput(0),
      msqInput(0),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,0)*Yu(2,0)) +
      AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0))*TCD0(msdInput(2),msqInput(0),
      msqInput(0),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,1)*Yu(0,1)) +
      AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1))*TCD0(msdInput(2),msqInput(1),
      msqInput(1),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,1)*Yu(1,1)) +
      AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1))*TCD0(msdInput(2),msqInput(1),
      msqInput(1),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,1)*Yu(2,1)) +
      AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1))*TCD0(msdInput(2),msqInput(1),
      msqInput(1),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,2)*Yu(0,2)) +
      AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2))*TCD0(msdInput(2),msqInput(2),
      msqInput(2),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,2)*Yu(1,2)) +
      AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2))*TCD0(msdInput(2),msqInput(2),
      msqInput(2),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,2)*Yu(2,2)) +
      AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2))*TCD0(msdInput(2),msqInput(2),
      msqInput(2),msuInput(2)) + Quad(Ye(0,0))*Sqr(Abs(MuInput))*Sqr(AeInput(0,0))
      *TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(0)) + AeInput(0,0)*
      AeInput(0,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(mseInput(0),
      mseInput(0),mslInput(0),mslInput(1)) + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(
      MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(0)
      ,mslInput(2)) + AeInput(0,0)*AeInput(0,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr
      (Ye(0,1))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(0)) + Quad(Ye(0,
      1))*Sqr(Abs(MuInput))*Sqr(AeInput(0,1))*TCD0(mseInput(0),mseInput(0),
      mslInput(1),mslInput(1)) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*Sqr(
      Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(2))
      + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0
      (mseInput(0),mseInput(0),mslInput(2),mslInput(0)) + AeInput(0,1)*AeInput(0,2
      )*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),
      mslInput(2),mslInput(1)) + Quad(Ye(0,2))*Sqr(Abs(MuInput))*Sqr(AeInput(0,2))
      *TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*
      Sqr(AeInput(0,0))*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(0),mseInput(1),
      mslInput(0),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(AeInput(0,1))*Sqr(Ye(0,1))*
      Sqr(Ye(1,1))*TCD0(mseInput(0),mseInput(1),mslInput(1),mslInput(1)) + Sqr(Abs
      (MuInput))*Sqr(AeInput(0,2))*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(0),
      mseInput(1),mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(0,0))*
      Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(
      0)) + Sqr(Abs(MuInput))*Sqr(AeInput(0,1))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(
      AeInput(0,2))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,0))*Sqr(Ye(0,0))*
      Sqr(Ye(1,0))*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(0)) + Sqr(Abs
      (MuInput))*Sqr(AeInput(1,1))*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0(mseInput(1),
      mseInput(0),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,2))*
      Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(0),mslInput(2),mslInput(
      2)) + Quad(Ye(1,0))*Sqr(Abs(MuInput))*Sqr(AeInput(1,0))*TCD0(mseInput(1),
      mseInput(1),mslInput(0),mslInput(0)) + AeInput(1,0)*AeInput(1,1)*Sqr(Abs(
      MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(0)
      ,mslInput(1)) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr
      (Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(2)) + AeInput(1,
      0)*AeInput(1,1)*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1)
      ,mseInput(1),mslInput(1),mslInput(0)) + Quad(Ye(1,1))*Sqr(Abs(MuInput))*Sqr(
      AeInput(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(1)) +
      AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(Ye(1,1))*Sqr(Ye(1,2))*TCD0(
      mseInput(1),mseInput(1),mslInput(1),mslInput(2)) + AeInput(1,0)*AeInput(1,2)
      *Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(2),mslInput(0)) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(
      Ye(1,1))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(2),mslInput(1))
      + Quad(Ye(1,2))*Sqr(Abs(MuInput))*Sqr(AeInput(1,2))*TCD0(mseInput(1),
      mseInput(1),mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,0))*
      Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(
      0)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,1))*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(
      AeInput(1,2))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(mseInput(1),mseInput(2),
      mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,0))*Sqr(Ye(0,0))*
      Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(0)) + Sqr(Abs
      (MuInput))*Sqr(AeInput(2,1))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(mseInput(2),
      mseInput(0),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,2))*
      Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(
      2)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,0))*Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(
      mseInput(2),mseInput(1),mslInput(0),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(
      AeInput(2,1))*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(1),
      mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,2))*Sqr(Ye(1,2))*
      Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(2)) + Quad(Ye
      (2,0))*Sqr(Abs(MuInput))*Sqr(AeInput(2,0))*TCD0(mseInput(2),mseInput(2),
      mslInput(0),mslInput(0)) + AeInput(2,0)*AeInput(2,1)*Sqr(Abs(MuInput))*Sqr(
      Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(0),mslInput(1))
      + AeInput(2,0)*AeInput(2,2)*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0
      (mseInput(2),mseInput(2),mslInput(0),mslInput(2)) + AeInput(2,0)*AeInput(2,1
      )*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),
      mslInput(1),mslInput(0)) + Quad(Ye(2,1))*Sqr(Abs(MuInput))*Sqr(AeInput(2,1))
      *TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(1)) + AeInput(2,1)*
      AeInput(2,2)*Sqr(Abs(MuInput))*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),
      mseInput(2),mslInput(1),mslInput(2)) + AeInput(2,0)*AeInput(2,2)*Sqr(Abs(
      MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2)
      ,mslInput(0)) + AeInput(2,1)*AeInput(2,2)*Sqr(Abs(MuInput))*Sqr(Ye(2,1))*Sqr
      (Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(1)) + Quad(Ye(2,
      2))*Sqr(Abs(MuInput))*Sqr(AeInput(2,2))*TCD0(mseInput(2),mseInput(2),
      mslInput(2),mslInput(2)) + 3*Quad(Yu(0,0))*Sqr(Abs(MuInput))*Sqr(AuInput(0,0
      ))*TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(0)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AuInput(1,0))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),
      msuInput(0),msuInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(2,0))*Sqr(Yu(0,0)
      )*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(2)) + 3*Sqr
      (Abs(MuInput))*Sqr(AuInput(0,0))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(1),msuInput(0)) + 3*Quad(Yu(1,0))*Sqr(Abs(MuInput))*Sqr
      (AuInput(1,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),msuInput(1)) + 3*Sqr
      (Abs(MuInput))*Sqr(AuInput(2,0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(1),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(0,0))
      *Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(2),msuInput
      (0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(
      msqInput(0),msqInput(0),msuInput(2),msuInput(1)) + 3*Quad(Yu(2,0))*Sqr(Abs(
      MuInput))*Sqr(AuInput(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(2),
      msuInput(2)) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*
      Sqr(Yu(0,1))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(0)) + 3*
      AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(
      msqInput(0),msqInput(1),msuInput(1),msuInput(1)) + 3*AuInput(2,0)*AuInput(2,
      1)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(0),msqInput(1),
      msuInput(2),msuInput(2)) + 3*AuInput(0,0)*AuInput(0,2)*Sqr(Abs(MuInput))*Sqr
      (Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(0),msqInput(2),msuInput(0),msuInput(0))
      + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*
      TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(1)) + 3*AuInput(2,0)*
      AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(msqInput(0),
      msqInput(2),msuInput(2),msuInput(2)) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(
      MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(1),msqInput(0),msuInput(0)
      ,msuInput(0)) + 3*AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*
      Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(1)) + 3*
      AuInput(2,0)*AuInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(
      msqInput(1),msqInput(0),msuInput(2),msuInput(2)) + 3*Quad(Yu(0,1))*Sqr(Abs(
      MuInput))*Sqr(AuInput(0,1))*TCD0(msqInput(1),msqInput(1),msuInput(0),
      msuInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,1))*Sqr(Yu(0,1))*Sqr(Yu(1,1
      ))*TCD0(msqInput(1),msqInput(1),msuInput(0),msuInput(1)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AuInput(2,1))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(0),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(0,1))*Sqr(Yu(0,1)
      )*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(0)) + 3*
      Quad(Yu(1,1))*Sqr(Abs(MuInput))*Sqr(AuInput(1,1))*TCD0(msqInput(1),msqInput(
      1),msuInput(1),msuInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(2,1))*Sqr(Yu(1
      ,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(2)) + 3*
      Sqr(Abs(MuInput))*Sqr(AuInput(0,1))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(
      1),msqInput(1),msuInput(2),msuInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,
      1))*Sqr(Yu(1,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(2),
      msuInput(1)) + 3*Quad(Yu(2,1))*Sqr(Abs(MuInput))*Sqr(AuInput(2,1))*TCD0(
      msqInput(1),msqInput(1),msuInput(2),msuInput(2)) + 3*AuInput(0,1)*AuInput(0,
      2)*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(1),msqInput(2),
      msuInput(0),msuInput(0)) + 3*AuInput(1,1)*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr
      (Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(1),msqInput(2),msuInput(1),msuInput(1))
      + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*Sqr(Yu(2,2))*
      TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(2)) + 3*AuInput(0,0)*
      AuInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(2),
      msqInput(0),msuInput(0),msuInput(0)) + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(
      MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(0),msuInput(1)
      ,msuInput(1)) + 3*AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*
      Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(2)) + 3*
      AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(
      msqInput(2),msqInput(1),msuInput(0),msuInput(0)) + 3*AuInput(1,1)*AuInput(1,
      2)*Sqr(Abs(MuInput))*Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(1),
      msuInput(1),msuInput(1)) + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr
      (Yu(2,1))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(2))
      + 3*Quad(Yu(0,2))*Sqr(Abs(MuInput))*Sqr(AuInput(0,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(0),msuInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,2))
      *Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(2,2))*Sqr(Yu(0,2))*Sqr(Yu(2,2))*TCD0(
      msqInput(2),msqInput(2),msuInput(0),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(
      AuInput(0,2))*Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),
      msuInput(1),msuInput(0)) + 3*Quad(Yu(1,2))*Sqr(Abs(MuInput))*Sqr(AuInput(1,2
      ))*TCD0(msqInput(2),msqInput(2),msuInput(1),msuInput(1)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AuInput(2,2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),
      msuInput(1),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(0,2))*Sqr(Yu(0,2)
      )*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(0)) + 3*Sqr
      (Abs(MuInput))*Sqr(AuInput(1,2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(2),msuInput(1)) + 3*Quad(Yu(2,2))*Sqr(Abs(MuInput))*Sqr
      (AuInput(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(2)) - 0.36*
      Quad(g1)*Sqr(Abs(MuInput))*TCD2t(M1Input,M1Input,Abs(MuInput),Abs(MuInput))
      + 0.6*Sqr(g1)*Sqr(g2)*(1 + (-2*M1Input*M2Input + 2*Sqr(Abs(MuInput)))*TCD2t(
      M1Input,M2Input,Abs(MuInput),Abs(MuInput)) - 2*TCD4t(M1Input,M2Input,Abs(
      MuInput),Abs(MuInput),SCALE)) + Quad(g2)*((-2*Sqr(Abs(M2Input)) - 3*Sqr(Abs(
      MuInput)))*TCD2t(M2Input,M2Input,Abs(MuInput),Abs(MuInput)) + 2*TCD4t(
      M2Input,M2Input,Abs(MuInput),Abs(MuInput),SCALE)) + 3*AdInput(0,0)*AdInput(0
      ,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(1))*
      Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0
      ,1)*Yd(1,0)*Yd(1,1) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(
      1,1) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*
      AdInput(0,0)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),
      msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*AdInput(0,0)*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(2),
      msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*AdInput(1,0)*AdInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(2))*Yd(0
      ,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0
      )*Yd(1,2) + 3*AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(1),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*
      AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),
      msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*AdInput(1,1)*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(1),
      msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*AdInput(1,1)*AdInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(1))*Yd(0
      ,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*
      TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0
      )*Yd(2,1) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*
      AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),
      msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(2,0)*
      AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),
      msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(1,0)*AdInput(1,1)*
      Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(0),msqInput(1))*Yd(1
      ,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*
      TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0
      )*Yd(2,1) + 3*AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(1),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*
      AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),
      msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(0,0)*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(0),
      msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(0,0)*AdInput(0,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput(0))*Yd(0
      ,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0
      )*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*
      AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),
      msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(1,0)*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(2),
      msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(2))*Yd(1
      ,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0
      )*Yd(2,2) + 3*AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*
      AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),
      msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(2,1)*
      AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),
      msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(2,1)*AdInput(2,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput(1))*Yd(0
      ,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(1,1)*AdInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1
      )*Yd(2,2) + 3*AdInput(1,1)*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*
      AdInput(2,1)*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),
      msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(2,1)*
      AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(2),
      msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + AeInput(0,0)*AeInput(0,1)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,0)
      *Ye(0,1)*Ye(1,0)*Ye(1,1) + AeInput(0,0)*AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(
      1,1) + AeInput(1,0)*AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + AeInput(1,0)*
      AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + AeInput(0,0)*AeInput(0,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(2))*Ye(0,0)
      *Ye(0,2)*Ye(1,0)*Ye(1,2) + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(
      1,2) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + AeInput(1,0)*
      AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(2),
      mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + AeInput(0,1)*AeInput(0,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(1),mslInput(2))*Ye(0,1)
      *Ye(0,2)*Ye(1,1)*Ye(1,2) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(
      1,2) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + AeInput(1,1)*
      AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(2),
      mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + AeInput(0,0)*AeInput(0,1)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(1))*Ye(0,0)
      *Ye(0,1)*Ye(2,0)*Ye(2,1) + AeInput(0,0)*AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(
      2,1) + AeInput(2,0)*AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput
      (0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + AeInput(2,0)*
      AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + AeInput(1,0)*AeInput(1,1)*Sqr
      (Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(1))*Ye(1,0)
      *Ye(1,1)*Ye(2,0)*Ye(2,1) + AeInput(1,0)*AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(
      2,1) + AeInput(2,0)*AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput
      (1),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + AeInput(2,0)*
      AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),
      mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + Sqr(Abs(MuInput))*TCC0(
      mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1
      ,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput
      (1),mslInput(0))*Ye(0,0)*Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0
      )*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput(0),mslInput(1))*Ye(
      1,0)*Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs
      (MuInput))*TCC0(mseInput(1),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*(Ye(0,0
      )*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(2),mslInput(0),mslInput(1))*Ye(2,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1
      ,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput
      (1),mslInput(0))*Ye(2,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0
      )*Ye(2,1)) + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),
      mseInput(2),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) +
      AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(0),
      mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)
      *Ye(0,2)*Ye(2,0)*Ye(2,2) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(
      2,2) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (2),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(0),
      mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)
      *Ye(1,2)*Ye(2,0)*Ye(2,2) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(
      2,2) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput
      (2),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(1),
      mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)
      *Ye(0,2)*Ye(2,1)*Ye(2,2) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(
      2,2) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (2),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),
      mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(1))*Ye(1,1)
      *Ye(1,2)*Ye(2,1)*Ye(2,2) + Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput(0),
      mslInput(2))*Ye(0,0)*Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye
      (2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)
      *Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs(
      MuInput))*TCC0(mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*(Ye(0,0)
      *Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*(Ye(0,0)*Ye(0,2) + Ye(1
      ,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput
      (0),mslInput(2))*Ye(2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0
      )*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput(2),mslInput(0))*Ye(
      2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs
      (MuInput))*TCC0(mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*(Ye(0,1
      )*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*(Ye(0,1)*Ye(0,2) + Ye(1
      ,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput
      (1),mslInput(2))*Ye(1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1
      )*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput(2),mslInput(1))*Ye(
      1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs
      (MuInput))*TCC0(mseInput(2),mslInput(1),mslInput(2))*Ye(2,1)*Ye(2,2)*(Ye(0,1
      )*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(2),mslInput(2),mslInput(1))*Ye(2,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1
      ,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 3*TCD0(msdInput(0),msqInput(0),msqInput(1),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(0,1)*Yu(0,1)) - 3*TCD0(msdInput(0),msqInput(1),msqInput(0),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(0,1)*Yu(0,1)) - 3*TCD0(msdInput(1),msqInput(0),msqInput(1),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(0,1)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(0),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(0,1)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(1),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(0,1)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(0),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(0,1)) - 3*TCD0(msdInput(0),msqInput(0),msqInput(2),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(0),msqInput(2),msqInput(0),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(0),msqInput(1),msqInput(2),
      msuInput(0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(0,1)*Yu(0,1))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(0),msqInput(2),msqInput(1),
      msuInput(0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(0,1)*Yu(0,1))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(0),msqInput(2),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(0),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(2),
      msuInput(0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(0,1))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(1),
      msuInput(0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(0,1))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(2),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(0),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(2),
      msuInput(0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(0,1))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(1),
      msuInput(0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(0,1))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) + TCC0(msdInput(0),msqInput(0),msuInput(0))*(-6*
      AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0)*(Yd(0,0)*Yu(0,0) + Yd(0,1)*Yu(0,1)
      + Yd(0,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(0,0)*Yu(0,0)*(Yd(0,0)*Yu(0,0) +
      Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2))) + TCC0(msdInput(0),msqInput(1),msuInput(
      0))*(-6*AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1)*(Yd(0,0)*Yu(0,0) + Yd(0,1)
      *Yu(0,1) + Yd(0,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(0,1)*Yu(0,1)*(Yd(0,0)*
      Yu(0,0) + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2))) + TCC0(msdInput(0),msqInput(2)
      ,msuInput(0))*(-6*AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2)*(Yd(0,0)*Yu(0,0)
      + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(0,2)*Yu(0,2)*(
      Yd(0,0)*Yu(0,0) + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2))) + TCC0(msdInput(1),
      msqInput(0),msuInput(0))*(-6*AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0)*(Yd(1
      ,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(1,
      0)*Yu(0,0)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2))) + TCC0(
      msdInput(1),msqInput(1),msuInput(0))*(-6*AdInput(1,1)*AuInput(0,1)*Yd(1,1)*
      Yu(0,1)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2)) + 6*Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(0,1)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,
      2))) + TCC0(msdInput(1),msqInput(2),msuInput(0))*(-6*AdInput(1,2)*AuInput(0,
      2)*Yd(1,2)*Yu(0,2)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2)) + 6
      *Sqr(Abs(MuInput))*Yd(1,2)*Yu(0,2)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1
      ,2)*Yu(0,2))) + TCC0(msdInput(2),msqInput(0),msuInput(0))*(-6*AdInput(2,0)*
      AuInput(0,0)*Yd(2,0)*Yu(0,0)*(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu
      (0,2)) + 6*Sqr(Abs(MuInput))*Yd(2,0)*Yu(0,0)*(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0
      ,1) + Yd(2,2)*Yu(0,2))) + TCC0(msdInput(2),msqInput(1),msuInput(0))*(-6*
      AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1)*(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0,1)
      + Yd(2,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(2,1)*Yu(0,1)*(Yd(2,0)*Yu(0,0) +
      Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu(0,2))) + TCC0(msdInput(2),msqInput(2),msuInput(
      0))*(-6*AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2)*(Yd(2,0)*Yu(0,0) + Yd(2,1)
      *Yu(0,1) + Yd(2,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(2,2)*Yu(0,2)*(Yd(2,0)*
      Yu(0,0) + Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu(0,2))) + 3*AuInput(1,0)*AuInput(1,1)*
      Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(1))*Yu(0
      ,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(1),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0
      )*Yu(1,1) + 3*AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),
      msqInput(0),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*
      AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),
      msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*TCD0(msdInput(0
      ),msqInput(0),msqInput(1),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu
      (1,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1
      )*Yu(1,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)) - 3*TCD0(msdInput(0),msqInput
      (1),msqInput(0),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0) -
      Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1
      ) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)) - 3*TCD0(msdInput(1),msqInput(0),
      msqInput(1),msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(
      Abs(MuInput))*Yd(1,0)*Yu(1,0))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) -
      Sqr(Abs(MuInput))*Yd(1,1)*Yu(1,1)) - 3*TCD0(msdInput(1),msqInput(1),msqInput
      (0),msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(Abs(
      MuInput))*Yd(1,0)*Yu(1,0))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) - Sqr(
      Abs(MuInput))*Yd(1,1)*Yu(1,1)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(1),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(1,1)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(0),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(1,1)) + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(2),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0
      )*Yu(1,2) + 3*AuInput(0,0)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),
      msqInput(2),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*
      AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),
      msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*AuInput(0,0)*
      AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(1),
      msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*AuInput(1,1)*AuInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(0),msuInput(1))*Yu(0
      ,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(1),msqInput(2),msuInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1
      )*Yu(1,2) + 3*AuInput(1,1)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),
      msqInput(1),msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*
      AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),
      msuInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*TCD0(msdInput(0
      ),msqInput(0),msqInput(2),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu
      (1,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2
      )*Yu(1,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(0),msqInput
      (2),msqInput(0),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0) -
      Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2
      ) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(0),msqInput(1),
      msqInput(2),msuInput(1))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1) - Sqr(
      Abs(MuInput))*Yd(0,1)*Yu(1,1))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2) -
      Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(0),msqInput(2),msqInput
      (1),msuInput(1))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1) - Sqr(Abs(
      MuInput))*Yd(0,1)*Yu(1,1))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2) - Sqr(
      Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(0),msqInput(2),
      msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(1,0))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(0),
      msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(1,0))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(2),
      msuInput(1))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(1,1))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(1),
      msuInput(1))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(1,1))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(2),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(0),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(2),
      msuInput(1))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(1,1))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(1),
      msuInput(1))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(1,1))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) + TCC0(msdInput(0),msqInput(0),msuInput(1))*(-6*
      AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0)*(Yd(0,0)*Yu(1,0) + Yd(0,1)*Yu(1,1)
      + Yd(0,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0)*(Yd(0,0)*Yu(1,0) +
      Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2))) + TCC0(msdInput(0),msqInput(1),msuInput(
      1))*(-6*AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1)*(Yd(0,0)*Yu(1,0) + Yd(0,1)
      *Yu(1,1) + Yd(0,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)*(Yd(0,0)*
      Yu(1,0) + Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2))) + TCC0(msdInput(0),msqInput(2)
      ,msuInput(1))*(-6*AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2)*(Yd(0,0)*Yu(1,0)
      + Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)*(
      Yd(0,0)*Yu(1,0) + Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2))) + TCC0(msdInput(1),
      msqInput(0),msuInput(1))*(-6*AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0)*(Yd(1
      ,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(1,
      0)*Yu(1,0)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2))) + TCC0(
      msdInput(1),msqInput(1),msuInput(1))*(-6*AdInput(1,1)*AuInput(1,1)*Yd(1,1)*
      Yu(1,1)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2)) + 6*Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(1,1)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,
      2))) + TCC0(msdInput(1),msqInput(2),msuInput(1))*(-6*AdInput(1,2)*AuInput(1,
      2)*Yd(1,2)*Yu(1,2)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2)) + 6
      *Sqr(Abs(MuInput))*Yd(1,2)*Yu(1,2)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1
      ,2)*Yu(1,2))) + TCC0(msdInput(2),msqInput(0),msuInput(1))*(-6*AdInput(2,0)*
      AuInput(1,0)*Yd(2,0)*Yu(1,0)*(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu
      (1,2)) + 6*Sqr(Abs(MuInput))*Yd(2,0)*Yu(1,0)*(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1
      ,1) + Yd(2,2)*Yu(1,2))) + TCC0(msdInput(2),msqInput(1),msuInput(1))*(-6*
      AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1)*(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1,1)
      + Yd(2,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(2,1)*Yu(1,1)*(Yd(2,0)*Yu(1,0) +
      Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu(1,2))) + TCC0(msdInput(2),msqInput(2),msuInput(
      1))*(-6*AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2)*(Yd(2,0)*Yu(1,0) + Yd(2,1)
      *Yu(1,1) + Yd(2,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(2,2)*Yu(1,2)*(Yd(2,0)*
      Yu(1,0) + Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu(1,2))) + 3*AuInput(2,0)*AuInput(2,1)*
      Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(2))*Yu(0
      ,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0
      )*Yu(2,1) + 3*AuInput(2,0)*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),
      msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*
      AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),
      msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(2,0)*
      AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(1),
      msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(1,0)*AuInput(1,1)*
      Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(1))*Yu(1
      ,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(2,0)*AuInput(2,1)*Sqr(Abs(MuInput))*
      TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0
      )*Yu(2,1) + 3*AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),
      msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*
      TCD0(msdInput(0),msqInput(0),msqInput(1),msuInput(2))*(AdInput(0,0)*AuInput(
      2,0)*Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,1)*
      AuInput(2,1)*Yd(0,1)*Yu(2,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)) - 3*TCD0(
      msdInput(0),msqInput(1),msqInput(0),msuInput(2))*(AdInput(0,0)*AuInput(2,0)*
      Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,1)*AuInput(2
      ,1)*Yd(0,1)*Yu(2,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)) - 3*TCD0(msdInput(1
      ),msqInput(0),msqInput(1),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu
      (2,0) - Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1
      )*Yu(2,1) - Sqr(Abs(MuInput))*Yd(1,1)*Yu(2,1)) - 3*TCD0(msdInput(1),msqInput
      (1),msqInput(0),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0) -
      Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1
      ) - Sqr(Abs(MuInput))*Yd(1,1)*Yu(2,1)) - 3*TCD0(msdInput(2),msqInput(0),
      msqInput(1),msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(
      Abs(MuInput))*Yd(2,0)*Yu(2,0))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) -
      Sqr(Abs(MuInput))*Yd(2,1)*Yu(2,1)) - 3*TCD0(msdInput(2),msqInput(1),msqInput
      (0),msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(Abs(
      MuInput))*Yd(2,0)*Yu(2,0))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) - Sqr(
      Abs(MuInput))*Yd(2,1)*Yu(2,1)) - 3*TCB0(msqInput(0),msqInput(1),SCALE)*(Yd(0
      ,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*(Yu(0,0)*Yu(0,1) + Yu(1,0)*
      Yu(1,1) + Yu(2,0)*Yu(2,1)) - 3*TCB0(msqInput(1),msqInput(0),SCALE)*(Yd(0,0)*
      Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1
      ,1) + Yu(2,0)*Yu(2,1)) + TCC0(msdInput(0),msqInput(0),msqInput(1))*(3*Sqr(
      Abs(MuInput))*Yd(0,0)*Yd(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*
      Yd(2,1)) - 3*AdInput(0,0)*AdInput(0,1)*Yd(0,0)*Yd(0,1)*(Yu(0,0)*Yu(0,1) + Yu
      (1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msdInput(0),msqInput(1),msqInput(0)
      )*(3*Sqr(Abs(MuInput))*Yd(0,0)*Yd(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) +
      Yd(2,0)*Yd(2,1)) - 3*AdInput(0,0)*AdInput(0,1)*Yd(0,0)*Yd(0,1)*(Yu(0,0)*Yu(0
      ,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msdInput(1),msqInput(0),
      msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)
      *Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*AdInput(1,0)*AdInput(1,1)*Yd(1,0)*Yd(1,1)*(
      Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msdInput(1),
      msqInput(1),msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,
      1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*AdInput(1,0)*AdInput(1,1)*Yd(1,0
      )*Yd(1,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(
      msdInput(2),msqInput(0),msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(2,0)*Yd(2,1)*(
      Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*AdInput(2,0)*
      AdInput(2,1)*Yd(2,0)*Yd(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu
      (2,1))) + TCC0(msdInput(2),msqInput(1),msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(
      2,0)*Yd(2,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*
      AdInput(2,0)*AdInput(2,1)*Yd(2,0)*Yd(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1)
      + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(0),msqInput(1),msuInput(0))*(-3*AuInput(
      0,0)*AuInput(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(0
      ,0)*Yu(0,1) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)
      *Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(1),msqInput(0),msuInput(0))*(-3
      *AuInput(0,0)*AuInput(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2
      ,1))*Yu(0,0)*Yu(0,1) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,1)*(Yu(0,0)*Yu(0,1)
      + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(0),msqInput(1),
      msuInput(1))*(-3*AuInput(1,0)*AuInput(1,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1
      ) + Yd(2,0)*Yd(2,1))*Yu(1,0)*Yu(1,1) + 3*Sqr(Abs(MuInput))*Yu(1,0)*Yu(1,1)*(
      Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(1),
      msqInput(0),msuInput(1))*(-3*AuInput(1,0)*AuInput(1,1)*(Yd(0,0)*Yd(0,1) + Yd
      (1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(1,0)*Yu(1,1) + 3*Sqr(Abs(MuInput))*Yu(1,
      0)*Yu(1,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(
      msqInput(0),msqInput(1),msuInput(2))*(-3*AuInput(2,0)*AuInput(2,1)*(Yd(0,0)*
      Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(2,0)*Yu(2,1) + 3*Sqr(Abs(
      MuInput))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,
      1))) + TCC0(msqInput(1),msqInput(0),msuInput(2))*(-3*AuInput(2,0)*AuInput(2,
      1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(2,0)*Yu(2,1) + 3
      *Sqr(Abs(MuInput))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2
      ,0)*Yu(2,1))) + 3*AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(
      0),msqInput(2),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*
      AuInput(0,0)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),
      msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(2,0)*
      AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(0),
      msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(0,0)*AuInput(0,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(0))*Yu(0
      ,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0
      )*Yu(2,2) + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),
      msqInput(2),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*
      AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),
      msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(1,0)*
      AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(2),
      msuInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(2,1)*AuInput(2,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(0),msuInput(2))*Yu(0
      ,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1
      )*Yu(2,2) + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),
      msqInput(1),msuInput(0),msuInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*
      AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),
      msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(2,1)*
      AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(1),
      msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(1,1)*AuInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(1))*Yu(1
      ,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(2),msqInput(1),msuInput(1),msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1
      )*Yu(2,2) + 3*AuInput(1,1)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),
      msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*
      TCD0(msdInput(0),msqInput(0),msqInput(2),msuInput(2))*(AdInput(0,0)*AuInput(
      2,0)*Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,2)*
      AuInput(2,2)*Yd(0,2)*Yu(2,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(
      msdInput(0),msqInput(2),msqInput(0),msuInput(2))*(AdInput(0,0)*AuInput(2,0)*
      Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,2)*AuInput(2
      ,2)*Yd(0,2)*Yu(2,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(msdInput(0
      ),msqInput(1),msqInput(2),msuInput(2))*(AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu
      (2,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1))*(AdInput(0,2)*AuInput(2,2)*Yd(0,2
      )*Yu(2,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(msdInput(0),msqInput
      (2),msqInput(1),msuInput(2))*(AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu(2,1) -
      Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1))*(AdInput(0,2)*AuInput(2,2)*Yd(0,2)*Yu(2,2
      ) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(0),
      msqInput(2),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0) - Sqr(
      Abs(MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) -
      Sqr(Abs(MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput
      (0),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0) - Sqr(Abs(
      MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) - Sqr(
      Abs(MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(2),
      msuInput(2))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(2,1))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(1),
      msuInput(2))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(2,1))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(2),
      msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(2,0))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(0),
      msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(2,0))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(2),
      msuInput(2))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(2,1))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(1),
      msuInput(2))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(2,1))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCB0(msqInput(0),msqInput(2),SCALE)*(Yd(0,0)*
      Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1
      ,2) + Yu(2,0)*Yu(2,2)) - 3*TCB0(msqInput(2),msqInput(0),SCALE)*(Yd(0,0)*Yd(0
      ,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2)
      + Yu(2,0)*Yu(2,2)) - 3*TCB0(msqInput(1),msqInput(2),SCALE)*(Yd(0,1)*Yd(0,2)
      + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu
      (2,1)*Yu(2,2)) - 3*TCB0(msqInput(2),msqInput(1),SCALE)*(Yd(0,1)*Yd(0,2) + Yd
      (1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1
      )*Yu(2,2)) + TCC0(msdInput(0),msqInput(0),msuInput(2))*(-6*AdInput(0,0)*
      AuInput(2,0)*Yd(0,0)*Yu(2,0)*(Yd(0,0)*Yu(2,0) + Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu
      (2,2)) + 6*Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0)*(Yd(0,0)*Yu(2,0) + Yd(0,1)*Yu(2
      ,1) + Yd(0,2)*Yu(2,2))) + TCC0(msdInput(0),msqInput(1),msuInput(2))*(-6*
      AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu(2,1)*(Yd(0,0)*Yu(2,0) + Yd(0,1)*Yu(2,1)
      + Yd(0,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)*(Yd(0,0)*Yu(2,0) +
      Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu(2,2))) + TCC0(msdInput(0),msqInput(2),msuInput(
      2))*(-6*AdInput(0,2)*AuInput(2,2)*Yd(0,2)*Yu(2,2)*(Yd(0,0)*Yu(2,0) + Yd(0,1)
      *Yu(2,1) + Yd(0,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)*(Yd(0,0)*
      Yu(2,0) + Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu(2,2))) + TCC0(msdInput(1),msqInput(0)
      ,msuInput(2))*(-6*AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0)*(Yd(1,0)*Yu(2,0)
      + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0)*(
      Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2))) + TCC0(msdInput(1),
      msqInput(1),msuInput(2))*(-6*AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1)*(Yd(1
      ,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(1,
      1)*Yu(2,1)*(Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2))) + TCC0(
      msdInput(1),msqInput(2),msuInput(2))*(-6*AdInput(1,2)*AuInput(2,2)*Yd(1,2)*
      Yu(2,2)*(Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2)) + 6*Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(2,2)*(Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,
      2))) + TCC0(msdInput(2),msqInput(0),msuInput(2))*(-6*AdInput(2,0)*AuInput(2,
      0)*Yd(2,0)*Yu(2,0)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1) + Yd(2,2)*Yu(2,2)) + 6
      *Sqr(Abs(MuInput))*Yd(2,0)*Yu(2,0)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1) + Yd(2
      ,2)*Yu(2,2))) + TCC0(msdInput(2),msqInput(1),msuInput(2))*(-6*AdInput(2,1)*
      AuInput(2,1)*Yd(2,1)*Yu(2,1)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1) + Yd(2,2)*Yu
      (2,2)) + 6*Sqr(Abs(MuInput))*Yd(2,1)*Yu(2,1)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2
      ,1) + Yd(2,2)*Yu(2,2))) + TCC0(msdInput(2),msqInput(2),msuInput(2))*(-6*
      AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1)
      + Yd(2,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(2,2)*Yu(2,2)*(Yd(2,0)*Yu(2,0) +
      Yd(2,1)*Yu(2,1) + Yd(2,2)*Yu(2,2))) + TCC0(msdInput(0),msqInput(0),msqInput(
      2))*(3*Sqr(Abs(MuInput))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2)
      + Yd(2,0)*Yd(2,2)) - 3*AdInput(0,0)*AdInput(0,2)*Yd(0,0)*Yd(0,2)*(Yu(0,0)*Yu
      (0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(0),msqInput(2),
      msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)
      *Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*AdInput(0,0)*AdInput(0,2)*Yd(0,0)*Yd(0,2)*(
      Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(1),
      msqInput(0),msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,2)*(Yd(0,0)*Yd(0,
      2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*AdInput(1,0)*AdInput(1,2)*Yd(1,0
      )*Yd(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(
      msdInput(1),msqInput(2),msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,2)*(
      Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*AdInput(1,0)*
      AdInput(1,2)*Yd(1,0)*Yd(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu
      (2,2))) + TCC0(msdInput(2),msqInput(0),msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(
      2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*
      AdInput(2,0)*AdInput(2,2)*Yd(2,0)*Yd(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2)
      + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(2),msqInput(2),msqInput(0))*(3*Sqr(Abs(
      MuInput))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,
      2)) - 3*AdInput(2,0)*AdInput(2,2)*Yd(2,0)*Yd(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)
      *Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msqInput(0),msqInput(2),msuInput(0))*(-3
      *AuInput(0,0)*AuInput(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2
      ,2))*Yu(0,0)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,2)*(Yu(0,0)*Yu(0,2)
      + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msqInput(2),msqInput(0),
      msuInput(0))*(-3*AuInput(0,0)*AuInput(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2
      ) + Yd(2,0)*Yd(2,2))*Yu(0,0)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,2)*(
      Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msqInput(0),
      msqInput(2),msuInput(1))*(-3*AuInput(1,0)*AuInput(1,2)*(Yd(0,0)*Yd(0,2) + Yd
      (1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(1,0)*Yu(1,2) + 3*Sqr(Abs(MuInput))*Yu(1,
      0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(
      msqInput(2),msqInput(0),msuInput(1))*(-3*AuInput(1,0)*AuInput(1,2)*(Yd(0,0)*
      Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(1,0)*Yu(1,2) + 3*Sqr(Abs(
      MuInput))*Yu(1,0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,
      2))) + TCC0(msqInput(0),msqInput(2),msuInput(2))*(-3*AuInput(2,0)*AuInput(2,
      2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(2,0)*Yu(2,2) + 3
      *Sqr(Abs(MuInput))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2
      ,0)*Yu(2,2))) + TCC0(msqInput(2),msqInput(0),msuInput(2))*(-3*AuInput(2,0)*
      AuInput(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(2,0)*
      Yu(2,2) + 3*Sqr(Abs(MuInput))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(
      1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(0),msqInput(1),msqInput(2))*(3*Sqr(
      Abs(MuInput))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*
      Yd(2,2)) - 3*AdInput(0,1)*AdInput(0,2)*Yd(0,1)*Yd(0,2)*(Yu(0,1)*Yu(0,2) + Yu
      (1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msdInput(0),msqInput(2),msqInput(1)
      )*(3*Sqr(Abs(MuInput))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) +
      Yd(2,1)*Yd(2,2)) - 3*AdInput(0,1)*AdInput(0,2)*Yd(0,1)*Yd(0,2)*(Yu(0,1)*Yu(0
      ,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msdInput(1),msqInput(1),
      msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(1,1)*Yd(1,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)
      *Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*AdInput(1,1)*AdInput(1,2)*Yd(1,1)*Yd(1,2)*(
      Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msdInput(1),
      msqInput(2),msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(1,1)*Yd(1,2)*(Yd(0,1)*Yd(0,
      2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*AdInput(1,1)*AdInput(1,2)*Yd(1,1
      )*Yd(1,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(
      msdInput(2),msqInput(1),msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(2,1)*Yd(2,2)*(
      Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*AdInput(2,1)*
      AdInput(2,2)*Yd(2,1)*Yd(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu
      (2,2))) + TCC0(msdInput(2),msqInput(2),msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(
      2,1)*Yd(2,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*
      AdInput(2,1)*AdInput(2,2)*Yd(2,1)*Yd(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2)
      + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(1),msqInput(2),msuInput(0))*(-3*AuInput(
      0,1)*AuInput(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(0
      ,1)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)
      *Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(2),msqInput(1),msuInput(0))*(-3
      *AuInput(0,1)*AuInput(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2
      ,2))*Yu(0,1)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2)
      + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(1),msqInput(2),
      msuInput(1))*(-3*AuInput(1,1)*AuInput(1,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2
      ) + Yd(2,1)*Yd(2,2))*Yu(1,1)*Yu(1,2) + 3*Sqr(Abs(MuInput))*Yu(1,1)*Yu(1,2)*(
      Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(2),
      msqInput(1),msuInput(1))*(-3*AuInput(1,1)*AuInput(1,2)*(Yd(0,1)*Yd(0,2) + Yd
      (1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(1,1)*Yu(1,2) + 3*Sqr(Abs(MuInput))*Yu(1,
      1)*Yu(1,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(
      msqInput(1),msqInput(2),msuInput(2))*(-3*AuInput(2,1)*AuInput(2,2)*(Yd(0,1)*
      Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(2,1)*Yu(2,2) + 3*Sqr(Abs(
      MuInput))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,
      2))) + TCC0(msqInput(2),msqInput(1),msuInput(2))*(-3*AuInput(2,1)*AuInput(2,
      2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(2,1)*Yu(2,2) + 3
      *Sqr(Abs(MuInput))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2
      ,1)*Yu(2,2)))))));
   MODEL->set_Lambda4(Re(-0.5*Sqr(g2) + (0.00008020298636472138*AuInput(2,2)*(1 -
      (0.5*AuInput(2,2))/MSUSY)*Quad(Yu(2,2))*Sqr(g3)*Sqr(MuInput)*UnitStep(-2 +
      LambdaLoopOrder))/Cube(MSUSY) + UnitStep(-1 + LambdaLoopOrder)*(
      0.0005277144981371759*(-4 + 8*Log(Sqr(M2Input)/Sqr(SCALE)) + 4*Log(Sqr(Abs(
      MuInput))/Sqr(SCALE)) + Log(Sqr(mslInput(0))/Sqr(SCALE)) + Log(Sqr(mslInput(
      1))/Sqr(SCALE)) + Log(Sqr(mslInput(2))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(0))/
      Sqr(SCALE)) + 3*Log(Sqr(msqInput(1))/Sqr(SCALE)) + 3*Log(Sqr(msqInput(2))/
      Sqr(SCALE)))*Quad(g2) - (Re(0.0031662869888230555*(3*Sqr(AdInput(0,0))*Sqr(
      Yd(0,0))*TCDB0(msdInput(0),msqInput(0)) + 3*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*
      TCDB0(msdInput(0),msqInput(1)) + 3*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*TCDB0(
      msdInput(0),msqInput(2)) + 3*Sqr(AdInput(1,0))*Sqr(Yd(1,0))*TCDB0(msdInput(1
      ),msqInput(0)) + 3*Sqr(AdInput(1,1))*Sqr(Yd(1,1))*TCDB0(msdInput(1),msqInput
      (1)) + 3*Sqr(AdInput(1,2))*Sqr(Yd(1,2))*TCDB0(msdInput(1),msqInput(2)) + 3*
      Sqr(AdInput(2,0))*Sqr(Yd(2,0))*TCDB0(msdInput(2),msqInput(0)) + 3*Sqr(
      AdInput(2,1))*Sqr(Yd(2,1))*TCDB0(msdInput(2),msqInput(1)) + 3*Sqr(AdInput(2,
      2))*Sqr(Yd(2,2))*TCDB0(msdInput(2),msqInput(2)) + Sqr(AeInput(0,0))*Sqr(Ye(0
      ,0))*TCDB0(mseInput(0),mslInput(0)) + Sqr(AeInput(0,1))*Sqr(Ye(0,1))*TCDB0(
      mseInput(0),mslInput(1)) + Sqr(AeInput(0,2))*Sqr(Ye(0,2))*TCDB0(mseInput(0),
      mslInput(2)) + Sqr(AeInput(1,0))*Sqr(Ye(1,0))*TCDB0(mseInput(1),mslInput(0))
      + Sqr(AeInput(1,1))*Sqr(Ye(1,1))*TCDB0(mseInput(1),mslInput(1)) + Sqr(
      AeInput(1,2))*Sqr(Ye(1,2))*TCDB0(mseInput(1),mslInput(2)) + Sqr(AeInput(2,0)
      )*Sqr(Ye(2,0))*TCDB0(mseInput(2),mslInput(0)) + Sqr(AeInput(2,1))*Sqr(Ye(2,1
      ))*TCDB0(mseInput(2),mslInput(1)) + Sqr(AeInput(2,2))*Sqr(Ye(2,2))*TCDB0(
      mseInput(2),mslInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*TCDB0(msuInput(0
      ),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*TCDB0(msuInput(0),msqInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*TCDB0(msuInput(0),msqInput(2)) + 3*
      Sqr(Abs(MuInput))*Sqr(Yu(1,0))*TCDB0(msuInput(1),msqInput(0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(1,1))*TCDB0(msuInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yu(1,2))*TCDB0(msuInput(1),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,0
      ))*TCDB0(msuInput(2),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*TCDB0(
      msuInput(2),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,2))*TCDB0(msuInput(2
      ),msqInput(2))) - 0.0007915717472057639*(0.6*Sqr(g1)*TCW(Abs(M1Input),Abs(
      MuInput),SCALE) + 3*Sqr(g2)*TCW(Abs(M2Input),Abs(MuInput),SCALE))) + Re(
      0.0031662869888230555*(3*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*TCDB0(msdInput(0),
      msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*TCDB0(msdInput(0),msqInput(1
      )) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,2))*TCDB0(msdInput(0),msqInput(2)) + 3*Sqr
      (Abs(MuInput))*Sqr(Yd(1,0))*TCDB0(msdInput(1),msqInput(0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(1,1))*TCDB0(msdInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*
      Sqr(Yd(1,2))*TCDB0(msdInput(1),msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,0
      ))*TCDB0(msdInput(2),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*TCDB0(
      msdInput(2),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(2,2))*TCDB0(msdInput(2
      ),msqInput(2)) + Sqr(Abs(MuInput))*Sqr(Ye(0,0))*TCDB0(mseInput(0),mslInput(0
      )) + Sqr(Abs(MuInput))*Sqr(Ye(0,1))*TCDB0(mseInput(0),mslInput(1)) + Sqr(Abs
      (MuInput))*Sqr(Ye(0,2))*TCDB0(mseInput(0),mslInput(2)) + Sqr(Abs(MuInput))*
      Sqr(Ye(1,0))*TCDB0(mseInput(1),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(Ye(1,1))
      *TCDB0(mseInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(Ye(1,2))*TCDB0(
      mseInput(1),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(Ye(2,0))*TCDB0(mseInput(2),
      mslInput(0)) + Sqr(Abs(MuInput))*Sqr(Ye(2,1))*TCDB0(mseInput(2),mslInput(1))
      + Sqr(Abs(MuInput))*Sqr(Ye(2,2))*TCDB0(mseInput(2),mslInput(2)) + 3*Sqr(
      AuInput(0,0))*Sqr(Yu(0,0))*TCDB0(msuInput(0),msqInput(0)) + 3*Sqr(AuInput(0,
      1))*Sqr(Yu(0,1))*TCDB0(msuInput(0),msqInput(1)) + 3*Sqr(AuInput(0,2))*Sqr(Yu
      (0,2))*TCDB0(msuInput(0),msqInput(2)) + 3*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*
      TCDB0(msuInput(1),msqInput(0)) + 3*Sqr(AuInput(1,1))*Sqr(Yu(1,1))*TCDB0(
      msuInput(1),msqInput(1)) + 3*Sqr(AuInput(1,2))*Sqr(Yu(1,2))*TCDB0(msuInput(1
      ),msqInput(2)) + 3*Sqr(AuInput(2,0))*Sqr(Yu(2,0))*TCDB0(msuInput(2),msqInput
      (0)) + 3*Sqr(AuInput(2,1))*Sqr(Yu(2,1))*TCDB0(msuInput(2),msqInput(1)) + 3*
      Sqr(AuInput(2,2))*Sqr(Yu(2,2))*TCDB0(msuInput(2),msqInput(2))) -
      0.0007915717472057639*(0.6*Sqr(g1)*TCW(Abs(M1Input),Abs(MuInput),SCALE) + 3*
      Sqr(g2)*TCW(Abs(M2Input),Abs(MuInput),SCALE))))*Sqr(g2) -
      0.006332573977646111*(-3*Sqr(Yd(0,0)*Yu(0,0) + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(
      0,2))*TCB0(msdInput(0),msuInput(0),SCALE) - 3*Sqr(Yd(0,0)*Yu(1,0) + Yd(0,1)*
      Yu(1,1) + Yd(0,2)*Yu(1,2))*TCB0(msdInput(0),msuInput(1),SCALE) - 3*Sqr(Yd(0,
      0)*Yu(2,0) + Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu(2,2))*TCB0(msdInput(0),msuInput(2)
      ,SCALE) - 3*Sqr(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2))*TCB0(
      msdInput(1),msuInput(0),SCALE) - 3*Sqr(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) +
      Yd(1,2)*Yu(1,2))*TCB0(msdInput(1),msuInput(1),SCALE) - 3*Sqr(Yd(1,0)*Yu(2,0)
      + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2))*TCB0(msdInput(1),msuInput(2),SCALE) - 3
      *Sqr(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu(0,2))*TCB0(msdInput(2),
      msuInput(0),SCALE) - 3*Sqr(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu(1,
      2))*TCB0(msdInput(2),msuInput(1),SCALE) - 3*Sqr(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu
      (2,1) + Yd(2,2)*Yu(2,2))*TCB0(msdInput(2),msuInput(2),SCALE) + (-0.25*Quad(
      g2) + 0.5*Sqr(g2)*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCB0(
      mslInput(0),mslInput(0),SCALE) + (-0.25*Quad(g2) + 0.5*Sqr(g2)*(Sqr(Ye(0,1))
      + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCB0(mslInput(1),mslInput(1),SCALE) + (-0.25
      *Quad(g2) + 0.5*Sqr(g2)*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCB0(
      mslInput(2),mslInput(2),SCALE) + (-0.75*Quad(g2) + 1.5*Sqr(g2)*(Sqr(Yd(0,0))
      + Sqr(Yd(1,0)) + Sqr(Yd(2,0))) + 1.5*Sqr(g2)*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) +
      Sqr(Yu(2,0))) - 3*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0)))*(Sqr(Yu(0,0))
      + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCB0(msqInput(0),msqInput(0),SCALE) + (-0.75
      *Quad(g2) + 1.5*Sqr(g2)*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) + 1.5*
      Sqr(g2)*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))) - 3*(Sqr(Yd(0,1)) + Sqr
      (Yd(1,1)) + Sqr(Yd(2,1)))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCB0
      (msqInput(1),msqInput(1),SCALE) + (-0.75*Quad(g2) + 1.5*Sqr(g2)*(Sqr(Yd(0,2)
      ) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) + 1.5*Sqr(g2)*(Sqr(Yu(0,2)) + Sqr(Yu(1,2))
      + Sqr(Yu(2,2))) - 3*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*(Sqr(Yu(0,2
      )) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCB0(msqInput(2),msqInput(2),SCALE) + (-
      1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Yd(0,0)) + 1.5*Sqr(g2)*Sqr(AdInput(0,0))*
      Sqr(Yd(0,0)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0))
      + Sqr(Yd(2,0))) - 3*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,
      0)) + Sqr(Yu(2,0))))*TCC0(msdInput(0),msqInput(0),msqInput(0)) + (-1.5*Sqr(
      g2)*Sqr(Abs(MuInput))*Sqr(Yd(0,1)) + 1.5*Sqr(g2)*Sqr(AdInput(0,1))*Sqr(Yd(0,
      1)) + 3*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd
      (2,1))) - 3*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) +
      Sqr(Yu(2,1))))*TCC0(msdInput(0),msqInput(1),msqInput(1)) + (-1.5*Sqr(g2)*Sqr
      (Abs(MuInput))*Sqr(Yd(0,2)) + 1.5*Sqr(g2)*Sqr(AdInput(0,2))*Sqr(Yd(0,2)) + 3
      *Sqr(Abs(MuInput))*Sqr(Yd(0,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))
      - 3*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2
      ))))*TCC0(msdInput(0),msqInput(2),msqInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(
      MuInput))*Sqr(Yd(1,0)) + 1.5*Sqr(g2)*Sqr(AdInput(1,0))*Sqr(Yd(1,0)) + 3*Sqr(
      Abs(MuInput))*Sqr(Yd(1,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))) - 3*
      Sqr(AdInput(1,0))*Sqr(Yd(1,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))
      *TCC0(msdInput(1),msqInput(0),msqInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))
      *Sqr(Yd(1,1)) + 1.5*Sqr(g2)*Sqr(AdInput(1,1))*Sqr(Yd(1,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(1,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) - 3*Sqr(
      AdInput(1,1))*Sqr(Yd(1,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*
      TCC0(msdInput(1),msqInput(1),msqInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(1,2)) + 1.5*Sqr(g2)*Sqr(AdInput(1,2))*Sqr(Yd(1,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(1,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) - 3*Sqr(
      AdInput(1,2))*Sqr(Yd(1,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*
      TCC0(msdInput(1),msqInput(2),msqInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,0)) + 1.5*Sqr(g2)*Sqr(AdInput(2,0))*Sqr(Yd(2,0)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))) - 3*Sqr(
      AdInput(2,0))*Sqr(Yd(2,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*
      TCC0(msdInput(2),msqInput(0),msqInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,1)) + 1.5*Sqr(g2)*Sqr(AdInput(2,1))*Sqr(Yd(2,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))) - 3*Sqr(
      AdInput(2,1))*Sqr(Yd(2,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*
      TCC0(msdInput(2),msqInput(1),msqInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,2)) + 1.5*Sqr(g2)*Sqr(AdInput(2,2))*Sqr(Yd(2,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yd(2,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))) - 3*Sqr(
      AdInput(2,2))*Sqr(Yd(2,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*
      TCC0(msdInput(2),msqInput(2),msqInput(2)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*
      Sqr(Ye(0,0)) + 0.5*Sqr(g2)*Sqr(AeInput(0,0))*Sqr(Ye(0,0)) + Sqr(Abs(MuInput)
      )*Sqr(Ye(0,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(0
      ),mslInput(0),mslInput(0)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(0,1)) +
      0.5*Sqr(g2)*Sqr(AeInput(0,1))*Sqr(Ye(0,1)) + Sqr(Abs(MuInput))*Sqr(Ye(0,1))*
      (Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(0),mslInput(1),
      mslInput(1)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(0,2)) + 0.5*Sqr(g2)*
      Sqr(AeInput(0,2))*Sqr(Ye(0,2)) + Sqr(Abs(MuInput))*Sqr(Ye(0,2))*(Sqr(Ye(0,2)
      ) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(0),mslInput(2),mslInput(2))
      + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(1,0)) + 0.5*Sqr(g2)*Sqr(AeInput(1,0
      ))*Sqr(Ye(1,0)) + Sqr(Abs(MuInput))*Sqr(Ye(1,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)
      ) + Sqr(Ye(2,0))))*TCC0(mseInput(1),mslInput(0),mslInput(0)) + (-0.5*Sqr(g2)
      *Sqr(Abs(MuInput))*Sqr(Ye(1,1)) + 0.5*Sqr(g2)*Sqr(AeInput(1,1))*Sqr(Ye(1,1))
      + Sqr(Abs(MuInput))*Sqr(Ye(1,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))
      ))*TCC0(mseInput(1),mslInput(1),mslInput(1)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput
      ))*Sqr(Ye(1,2)) + 0.5*Sqr(g2)*Sqr(AeInput(1,2))*Sqr(Ye(1,2)) + Sqr(Abs(
      MuInput))*Sqr(Ye(1,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(
      mseInput(1),mslInput(2),mslInput(2)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Ye(2,0)) + 0.5*Sqr(g2)*Sqr(AeInput(2,0))*Sqr(Ye(2,0)) + Sqr(Abs(MuInput))*
      Sqr(Ye(2,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(2),
      mslInput(0),mslInput(0)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(2,1)) +
      0.5*Sqr(g2)*Sqr(AeInput(2,1))*Sqr(Ye(2,1)) + Sqr(Abs(MuInput))*Sqr(Ye(2,1))*
      (Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(2),mslInput(1),
      mslInput(1)) + (-0.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Ye(2,2)) + 0.5*Sqr(g2)*
      Sqr(AeInput(2,2))*Sqr(Ye(2,2)) + Sqr(Abs(MuInput))*Sqr(Ye(2,2))*(Sqr(Ye(0,2)
      ) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(2),mslInput(2),mslInput(2))
      + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Yu(0,0)) + 1.5*Sqr(g2)*Sqr(AuInput(0,0
      ))*Sqr(Yu(0,0)) - 3*Sqr(AuInput(0,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(
      2,0)))*Sqr(Yu(0,0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*(Sqr(Yu(0,0)) + Sqr(
      Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(0)) + (-1.5*
      Sqr(g2)*Sqr(Abs(MuInput))*Sqr(Yu(1,0)) + 1.5*Sqr(g2)*Sqr(AuInput(1,0))*Sqr(
      Yu(1,0)) - 3*Sqr(AuInput(1,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0)))*
      Sqr(Yu(1,0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0))
      + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(1)) + (-1.5*Sqr(g2)*
      Sqr(Abs(MuInput))*Sqr(Yu(2,0)) + 1.5*Sqr(g2)*Sqr(AuInput(2,0))*Sqr(Yu(2,0))
      - 3*Sqr(AuInput(2,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0)))*Sqr(Yu(2,
      0)) + 3*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu
      (2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(
      MuInput))*Sqr(Yu(0,1)) + 1.5*Sqr(g2)*Sqr(AuInput(0,1))*Sqr(Yu(0,1)) - 3*Sqr(
      AuInput(0,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1)))*Sqr(Yu(0,1)) + 3*
      Sqr(Abs(MuInput))*Sqr(Yu(0,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))
      *TCC0(msqInput(1),msqInput(1),msuInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))
      *Sqr(Yu(1,1)) + 1.5*Sqr(g2)*Sqr(AuInput(1,1))*Sqr(Yu(1,1)) - 3*Sqr(AuInput(1
      ,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1)))*Sqr(Yu(1,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(1,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(
      msqInput(1),msqInput(1),msuInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(2,1)) + 1.5*Sqr(g2)*Sqr(AuInput(2,1))*Sqr(Yu(2,1)) - 3*Sqr(AuInput(2,1))*
      (Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1)))*Sqr(Yu(2,1)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(2,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(
      msqInput(1),msqInput(1),msuInput(2)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(0,2)) + 1.5*Sqr(g2)*Sqr(AuInput(0,2))*Sqr(Yu(0,2)) - 3*Sqr(AuInput(0,2))*
      (Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*Sqr(Yu(0,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(0,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(
      msqInput(2),msqInput(2),msuInput(0)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(1,2)) + 1.5*Sqr(g2)*Sqr(AuInput(1,2))*Sqr(Yu(1,2)) - 3*Sqr(AuInput(1,2))*
      (Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*Sqr(Yu(1,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(1,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(
      msqInput(2),msqInput(2),msuInput(1)) + (-1.5*Sqr(g2)*Sqr(Abs(MuInput))*Sqr(
      Yu(2,2)) + 1.5*Sqr(g2)*Sqr(AuInput(2,2))*Sqr(Yu(2,2)) - 3*Sqr(AuInput(2,2))*
      (Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2)))*Sqr(Yu(2,2)) + 3*Sqr(Abs(
      MuInput))*Sqr(Yu(2,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(
      msqInput(2),msqInput(2),msuInput(2)) + 3*Quad(Yd(0,0))*Sqr(Abs(MuInput))*Sqr
      (AdInput(0,0))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(0)) + 3*
      AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(
      msdInput(0),msdInput(0),msqInput(0),msqInput(1)) + 3*AdInput(0,0)*AdInput(0,
      2)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),
      msqInput(0),msqInput(2)) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*Sqr
      (Yd(0,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1),msqInput(0))
      + 3*Quad(Yd(0,1))*Sqr(Abs(MuInput))*Sqr(AdInput(0,1))*TCD0(msdInput(0),
      msdInput(0),msqInput(1),msqInput(1)) + 3*AdInput(0,1)*AdInput(0,2)*Sqr(Abs(
      MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(1)
      ,msqInput(2)) + 3*AdInput(0,0)*AdInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*
      Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),msqInput(0)) + 3*
      AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(
      msdInput(0),msdInput(0),msqInput(2),msqInput(1)) + 3*Quad(Yd(0,2))*Sqr(Abs(
      MuInput))*Sqr(AdInput(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),
      msqInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*Sqr(Yd(1,0
      ))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AdInput(0,1))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(0),msdInput(1),
      msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,2))*Sqr(Yd(0,2)
      )*Sqr(Yd(1,2))*TCD0(msdInput(0),msdInput(1),msqInput(2),msqInput(2)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(0),
      msdInput(2),msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,1))
      *Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(0),msdInput(2),msqInput(1),msqInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(0),msdInput(2),msqInput(2),msqInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))
      *Yd(0,0)*Yu(0,0)) + AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0))*TCD0(msdInput
      (0),msqInput(0),msqInput(0),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,0)
      *Yu(1,0)) + AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0))*TCD0(msdInput(0),
      msqInput(0),msqInput(0),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,0)*Yu(
      2,0)) + AdInput(0,0)*AuInput(2,0)*Yd(0,0)*Yu(2,0))*TCD0(msdInput(0),msqInput
      (0),msqInput(0),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,1)*Yu(0,1)) +
      AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1))*TCD0(msdInput(0),msqInput(1),
      msqInput(1),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)) +
      AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1))*TCD0(msdInput(0),msqInput(1),
      msqInput(1),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)) +
      AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu(2,1))*TCD0(msdInput(0),msqInput(1),
      msqInput(1),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,2)*Yu(0,2)) +
      AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2))*TCD0(msdInput(0),msqInput(2),
      msqInput(2),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) +
      AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2))*TCD0(msdInput(0),msqInput(2),
      msqInput(2),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) +
      AdInput(0,2)*AuInput(2,2)*Yd(0,2)*Yu(2,2))*TCD0(msdInput(0),msqInput(2),
      msqInput(2),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(1,0))*Sqr(Yd(0,0)
      )*Sqr(Yd(1,0))*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(0)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(1,1))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(1),
      msdInput(0),msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(1,2))
      *Sqr(Yd(0,2))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput
      (2)) + 3*Quad(Yd(1,0))*Sqr(Abs(MuInput))*Sqr(AdInput(1,0))*TCD0(msdInput(1),
      msdInput(1),msqInput(0),msqInput(0)) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(
      MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(0)
      ,msqInput(1)) + 3*AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*
      Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(0),msqInput(2)) + 3*
      AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(
      msdInput(1),msdInput(1),msqInput(1),msqInput(0)) + 3*Quad(Yd(1,1))*Sqr(Abs(
      MuInput))*Sqr(AdInput(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(1),
      msqInput(1)) + 3*AdInput(1,1)*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*
      Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(2)) + 3*
      AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(
      msdInput(1),msdInput(1),msqInput(2),msqInput(0)) + 3*AdInput(1,1)*AdInput(1,
      2)*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),
      msqInput(2),msqInput(1)) + 3*Quad(Yd(1,2))*Sqr(Abs(MuInput))*Sqr(AdInput(1,2
      ))*TCD0(msdInput(1),msdInput(1),msqInput(2),msqInput(2)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AdInput(1,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(1),msdInput(2),
      msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(1,1))*Sqr(Yd(1,1)
      )*Sqr(Yd(2,1))*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(1)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(1,2))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,0)*Yu(
      0,0)) + AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0))*TCD0(msdInput(1),msqInput
      (0),msqInput(0),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,0)*Yu(1,0)) +
      AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0))*TCD0(msdInput(1),msqInput(0),
      msqInput(0),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0)) +
      AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0))*TCD0(msdInput(1),msqInput(0),
      msqInput(0),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,1)*Yu(0,1)) +
      AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1))*TCD0(msdInput(1),msqInput(1),
      msqInput(1),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,1)*Yu(1,1)) +
      AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1))*TCD0(msdInput(1),msqInput(1),
      msqInput(1),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,1)*Yu(2,1)) +
      AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1))*TCD0(msdInput(1),msqInput(1),
      msqInput(1),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,2)*Yu(0,2)) +
      AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2))*TCD0(msdInput(1),msqInput(2),
      msqInput(2),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,2)*Yu(1,2)) +
      AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2))*TCD0(msdInput(1),msqInput(2),
      msqInput(2),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(1,2)*Yu(2,2)) +
      AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2))*TCD0(msdInput(1),msqInput(2),
      msqInput(2),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,0))*Sqr(Yd(0,0)
      )*Sqr(Yd(2,0))*TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(0)) + 3*Sqr
      (Abs(MuInput))*Sqr(AdInput(2,1))*Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(0),msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,2))
      *Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput
      (2)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(
      msdInput(2),msdInput(1),msqInput(0),msqInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(
      AdInput(2,1))*Sqr(Yd(1,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(1),
      msqInput(1),msqInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AdInput(2,2))*Sqr(Yd(1,2)
      )*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(2)) + 3*
      Quad(Yd(2,0))*Sqr(Abs(MuInput))*Sqr(AdInput(2,0))*TCD0(msdInput(2),msdInput(
      2),msqInput(0),msqInput(0)) + 3*AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*
      Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(
      1)) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,2)
      )*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(2)) + 3*AdInput(2,0)*
      AdInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(2),msqInput(1),msqInput(0)) + 3*Quad(Yd(2,1))*Sqr(Abs(MuInput))*Sqr
      (AdInput(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(1)) + 3*
      AdInput(2,1)*AdInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(2),msqInput(1),msqInput(2)) + 3*AdInput(2,0)*AdInput(2,
      2)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(2),msqInput(0)) + 3*AdInput(2,1)*AdInput(2,2)*Sqr(Abs(MuInput))*Sqr
      (Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(1))
      + 3*Quad(Yd(2,2))*Sqr(Abs(MuInput))*Sqr(AdInput(2,2))*TCD0(msdInput(2),
      msdInput(2),msqInput(2),msqInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,0)*Yu(
      0,0)) + AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0))*TCD0(msdInput(2),msqInput
      (0),msqInput(0),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,0)*Yu(1,0)) +
      AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0))*TCD0(msdInput(2),msqInput(0),
      msqInput(0),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,0)*Yu(2,0)) +
      AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0))*TCD0(msdInput(2),msqInput(0),
      msqInput(0),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,1)*Yu(0,1)) +
      AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1))*TCD0(msdInput(2),msqInput(1),
      msqInput(1),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,1)*Yu(1,1)) +
      AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1))*TCD0(msdInput(2),msqInput(1),
      msqInput(1),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,1)*Yu(2,1)) +
      AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1))*TCD0(msdInput(2),msqInput(1),
      msqInput(1),msuInput(2)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,2)*Yu(0,2)) +
      AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2))*TCD0(msdInput(2),msqInput(2),
      msqInput(2),msuInput(0)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,2)*Yu(1,2)) +
      AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2))*TCD0(msdInput(2),msqInput(2),
      msqInput(2),msuInput(1)) - 3*Sqr(-(Sqr(Abs(MuInput))*Yd(2,2)*Yu(2,2)) +
      AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2))*TCD0(msdInput(2),msqInput(2),
      msqInput(2),msuInput(2)) + Quad(Ye(0,0))*Sqr(Abs(MuInput))*Sqr(AeInput(0,0))
      *TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(0)) + AeInput(0,0)*
      AeInput(0,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(mseInput(0),
      mseInput(0),mslInput(0),mslInput(1)) + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(
      MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(0)
      ,mslInput(2)) + AeInput(0,0)*AeInput(0,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr
      (Ye(0,1))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(0)) + Quad(Ye(0,
      1))*Sqr(Abs(MuInput))*Sqr(AeInput(0,1))*TCD0(mseInput(0),mseInput(0),
      mslInput(1),mslInput(1)) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*Sqr(
      Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(2))
      + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0
      (mseInput(0),mseInput(0),mslInput(2),mslInput(0)) + AeInput(0,1)*AeInput(0,2
      )*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),
      mslInput(2),mslInput(1)) + Quad(Ye(0,2))*Sqr(Abs(MuInput))*Sqr(AeInput(0,2))
      *TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*
      Sqr(AeInput(0,0))*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(0),mseInput(1),
      mslInput(0),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(AeInput(0,1))*Sqr(Ye(0,1))*
      Sqr(Ye(1,1))*TCD0(mseInput(0),mseInput(1),mslInput(1),mslInput(1)) + Sqr(Abs
      (MuInput))*Sqr(AeInput(0,2))*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(0),
      mseInput(1),mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(0,0))*
      Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(
      0)) + Sqr(Abs(MuInput))*Sqr(AeInput(0,1))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(
      AeInput(0,2))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,0))*Sqr(Ye(0,0))*
      Sqr(Ye(1,0))*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(0)) + Sqr(Abs
      (MuInput))*Sqr(AeInput(1,1))*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0(mseInput(1),
      mseInput(0),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,2))*
      Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(0),mslInput(2),mslInput(
      2)) + Quad(Ye(1,0))*Sqr(Abs(MuInput))*Sqr(AeInput(1,0))*TCD0(mseInput(1),
      mseInput(1),mslInput(0),mslInput(0)) + AeInput(1,0)*AeInput(1,1)*Sqr(Abs(
      MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(0)
      ,mslInput(1)) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr
      (Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(2)) + AeInput(1,
      0)*AeInput(1,1)*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1)
      ,mseInput(1),mslInput(1),mslInput(0)) + Quad(Ye(1,1))*Sqr(Abs(MuInput))*Sqr(
      AeInput(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(1)) +
      AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(Ye(1,1))*Sqr(Ye(1,2))*TCD0(
      mseInput(1),mseInput(1),mslInput(1),mslInput(2)) + AeInput(1,0)*AeInput(1,2)
      *Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(2),mslInput(0)) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(
      Ye(1,1))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(2),mslInput(1))
      + Quad(Ye(1,2))*Sqr(Abs(MuInput))*Sqr(AeInput(1,2))*TCD0(mseInput(1),
      mseInput(1),mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,0))*
      Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(
      0)) + Sqr(Abs(MuInput))*Sqr(AeInput(1,1))*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(
      AeInput(1,2))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(mseInput(1),mseInput(2),
      mslInput(2),mslInput(2)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,0))*Sqr(Ye(0,0))*
      Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(0)) + Sqr(Abs
      (MuInput))*Sqr(AeInput(2,1))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(mseInput(2),
      mseInput(0),mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,2))*
      Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(
      2)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,0))*Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(
      mseInput(2),mseInput(1),mslInput(0),mslInput(0)) + Sqr(Abs(MuInput))*Sqr(
      AeInput(2,1))*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(1),
      mslInput(1),mslInput(1)) + Sqr(Abs(MuInput))*Sqr(AeInput(2,2))*Sqr(Ye(1,2))*
      Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(2)) + Quad(Ye
      (2,0))*Sqr(Abs(MuInput))*Sqr(AeInput(2,0))*TCD0(mseInput(2),mseInput(2),
      mslInput(0),mslInput(0)) + AeInput(2,0)*AeInput(2,1)*Sqr(Abs(MuInput))*Sqr(
      Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(0),mslInput(1))
      + AeInput(2,0)*AeInput(2,2)*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0
      (mseInput(2),mseInput(2),mslInput(0),mslInput(2)) + AeInput(2,0)*AeInput(2,1
      )*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),
      mslInput(1),mslInput(0)) + Quad(Ye(2,1))*Sqr(Abs(MuInput))*Sqr(AeInput(2,1))
      *TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(1)) + AeInput(2,1)*
      AeInput(2,2)*Sqr(Abs(MuInput))*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),
      mseInput(2),mslInput(1),mslInput(2)) + AeInput(2,0)*AeInput(2,2)*Sqr(Abs(
      MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2)
      ,mslInput(0)) + AeInput(2,1)*AeInput(2,2)*Sqr(Abs(MuInput))*Sqr(Ye(2,1))*Sqr
      (Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(1)) + Quad(Ye(2,
      2))*Sqr(Abs(MuInput))*Sqr(AeInput(2,2))*TCD0(mseInput(2),mseInput(2),
      mslInput(2),mslInput(2)) + 3*Quad(Yu(0,0))*Sqr(Abs(MuInput))*Sqr(AuInput(0,0
      ))*TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(0)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AuInput(1,0))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),
      msuInput(0),msuInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(2,0))*Sqr(Yu(0,0)
      )*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(2)) + 3*Sqr
      (Abs(MuInput))*Sqr(AuInput(0,0))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(1),msuInput(0)) + 3*Quad(Yu(1,0))*Sqr(Abs(MuInput))*Sqr
      (AuInput(1,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),msuInput(1)) + 3*Sqr
      (Abs(MuInput))*Sqr(AuInput(2,0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(1),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(0,0))
      *Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(2),msuInput
      (0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(
      msqInput(0),msqInput(0),msuInput(2),msuInput(1)) + 3*Quad(Yu(2,0))*Sqr(Abs(
      MuInput))*Sqr(AuInput(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(2),
      msuInput(2)) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*
      Sqr(Yu(0,1))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(0)) + 3*
      AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(
      msqInput(0),msqInput(1),msuInput(1),msuInput(1)) + 3*AuInput(2,0)*AuInput(2,
      1)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(0),msqInput(1),
      msuInput(2),msuInput(2)) + 3*AuInput(0,0)*AuInput(0,2)*Sqr(Abs(MuInput))*Sqr
      (Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(0),msqInput(2),msuInput(0),msuInput(0))
      + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*
      TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(1)) + 3*AuInput(2,0)*
      AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(msqInput(0),
      msqInput(2),msuInput(2),msuInput(2)) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(
      MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(1),msqInput(0),msuInput(0)
      ,msuInput(0)) + 3*AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*
      Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(1)) + 3*
      AuInput(2,0)*AuInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(
      msqInput(1),msqInput(0),msuInput(2),msuInput(2)) + 3*Quad(Yu(0,1))*Sqr(Abs(
      MuInput))*Sqr(AuInput(0,1))*TCD0(msqInput(1),msqInput(1),msuInput(0),
      msuInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,1))*Sqr(Yu(0,1))*Sqr(Yu(1,1
      ))*TCD0(msqInput(1),msqInput(1),msuInput(0),msuInput(1)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AuInput(2,1))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(0),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(0,1))*Sqr(Yu(0,1)
      )*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(0)) + 3*
      Quad(Yu(1,1))*Sqr(Abs(MuInput))*Sqr(AuInput(1,1))*TCD0(msqInput(1),msqInput(
      1),msuInput(1),msuInput(1)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(2,1))*Sqr(Yu(1
      ,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(2)) + 3*
      Sqr(Abs(MuInput))*Sqr(AuInput(0,1))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(
      1),msqInput(1),msuInput(2),msuInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,
      1))*Sqr(Yu(1,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(2),
      msuInput(1)) + 3*Quad(Yu(2,1))*Sqr(Abs(MuInput))*Sqr(AuInput(2,1))*TCD0(
      msqInput(1),msqInput(1),msuInput(2),msuInput(2)) + 3*AuInput(0,1)*AuInput(0,
      2)*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(1),msqInput(2),
      msuInput(0),msuInput(0)) + 3*AuInput(1,1)*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr
      (Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(1),msqInput(2),msuInput(1),msuInput(1))
      + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*Sqr(Yu(2,2))*
      TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(2)) + 3*AuInput(0,0)*
      AuInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(2),
      msqInput(0),msuInput(0),msuInput(0)) + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(
      MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(0),msuInput(1)
      ,msuInput(1)) + 3*AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*
      Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(2)) + 3*
      AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(
      msqInput(2),msqInput(1),msuInput(0),msuInput(0)) + 3*AuInput(1,1)*AuInput(1,
      2)*Sqr(Abs(MuInput))*Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(1),
      msuInput(1),msuInput(1)) + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr
      (Yu(2,1))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(2))
      + 3*Quad(Yu(0,2))*Sqr(Abs(MuInput))*Sqr(AuInput(0,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(0),msuInput(0)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(1,2))
      *Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput
      (1)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(2,2))*Sqr(Yu(0,2))*Sqr(Yu(2,2))*TCD0(
      msqInput(2),msqInput(2),msuInput(0),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(
      AuInput(0,2))*Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),
      msuInput(1),msuInput(0)) + 3*Quad(Yu(1,2))*Sqr(Abs(MuInput))*Sqr(AuInput(1,2
      ))*TCD0(msqInput(2),msqInput(2),msuInput(1),msuInput(1)) + 3*Sqr(Abs(MuInput
      ))*Sqr(AuInput(2,2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),
      msuInput(1),msuInput(2)) + 3*Sqr(Abs(MuInput))*Sqr(AuInput(0,2))*Sqr(Yu(0,2)
      )*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(0)) + 3*Sqr
      (Abs(MuInput))*Sqr(AuInput(1,2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(2),msuInput(1)) + 3*Quad(Yu(2,2))*Sqr(Abs(MuInput))*Sqr
      (AuInput(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(2)) - 0.36*
      Quad(g1)*Sqr(Abs(MuInput))*TCD2t(M1Input,M1Input,Abs(MuInput),Abs(MuInput))
      + 0.6*Sqr(g1)*Sqr(g2)*(1 + (-2*M1Input*M2Input + 2*Sqr(Abs(MuInput)))*TCD2t(
      M1Input,M2Input,Abs(MuInput),Abs(MuInput)) - 2*TCD4t(M1Input,M2Input,Abs(
      MuInput),Abs(MuInput),SCALE)) + Quad(g2)*((-2*Sqr(Abs(M2Input)) - 3*Sqr(Abs(
      MuInput)))*TCD2t(M2Input,M2Input,Abs(MuInput),Abs(MuInput)) + 2*TCD4t(
      M2Input,M2Input,Abs(MuInput),Abs(MuInput),SCALE)) + 3*AdInput(0,0)*AdInput(0
      ,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(1))*
      Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(
      MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0
      ,1)*Yd(1,0)*Yd(1,1) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(
      1,1) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*
      AdInput(0,0)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),
      msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*AdInput(0,0)*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),msqInput(2),
      msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*AdInput(1,0)*AdInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(2))*Yd(0
      ,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0
      )*Yd(1,2) + 3*AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(1),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*
      AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1),
      msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*AdInput(1,1)*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(1),
      msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*AdInput(1,1)*AdInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(1))*Yd(0
      ,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*
      TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0
      )*Yd(2,1) + 3*AdInput(0,0)*AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*
      AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),
      msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(2,0)*
      AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),
      msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(1,0)*AdInput(1,1)*
      Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(0),msqInput(1))*Yd(1
      ,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(1,0)*AdInput(1,1)*Sqr(Abs(MuInput))*
      TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0
      )*Yd(2,1) + 3*AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(1),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*
      AdInput(2,0)*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),
      msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*AdInput(0,0)*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(0),
      msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(0,0)*AdInput(0,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput(0))*Yd(0
      ,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0
      )*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),
      msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*
      AdInput(1,0)*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),
      msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(1,0)*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2),msqInput(2),
      msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(2))*Yd(1
      ,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*AdInput(2,0)*AdInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0
      )*Yd(2,2) + 3*AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*
      AdInput(0,1)*AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2),
      msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(2,1)*
      AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),
      msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(2,1)*AdInput(2,2)*
      Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput(1))*Yd(0
      ,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(1,1)*AdInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1
      )*Yd(2,2) + 3*AdInput(1,1)*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),
      msdInput(2),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*
      AdInput(2,1)*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),
      msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*AdInput(2,1)*
      AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(2),
      msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + AeInput(0,0)*AeInput(0,1)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,0)
      *Ye(0,1)*Ye(1,0)*Ye(1,1) + AeInput(0,0)*AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(
      1,1) + AeInput(1,0)*AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + AeInput(1,0)*
      AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + AeInput(0,0)*AeInput(0,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(2))*Ye(0,0)
      *Ye(0,2)*Ye(1,0)*Ye(1,2) + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(
      1,2) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + AeInput(1,0)*
      AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(2),
      mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + AeInput(0,1)*AeInput(0,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(1),mslInput(2))*Ye(0,1)
      *Ye(0,2)*Ye(1,1)*Ye(1,2) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(
      1,2) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + AeInput(1,1)*
      AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(2),
      mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + AeInput(0,0)*AeInput(0,1)*Sqr
      (Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(1))*Ye(0,0)
      *Ye(0,1)*Ye(2,0)*Ye(2,1) + AeInput(0,0)*AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(
      2,1) + AeInput(2,0)*AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput
      (0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + AeInput(2,0)*
      AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + AeInput(1,0)*AeInput(1,1)*Sqr
      (Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(1))*Ye(1,0)
      *Ye(1,1)*Ye(2,0)*Ye(2,1) + AeInput(1,0)*AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(
      2,1) + AeInput(2,0)*AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput
      (1),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + AeInput(2,0)*
      AeInput(2,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),
      mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + Sqr(Abs(MuInput))*TCC0(
      mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1
      ,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput
      (1),mslInput(0))*Ye(0,0)*Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0
      )*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput(0),mslInput(1))*Ye(
      1,0)*Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs
      (MuInput))*TCC0(mseInput(1),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*(Ye(0,0
      )*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(2),mslInput(0),mslInput(1))*Ye(2,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1
      ,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput
      (1),mslInput(0))*Ye(2,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0
      )*Ye(2,1)) + AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),
      mseInput(2),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) +
      AeInput(0,0)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(0),
      mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)
      *Ye(0,2)*Ye(2,0)*Ye(2,2) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(
      2,2) + AeInput(1,0)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (2),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(0),
      mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + AeInput(2,0)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)
      *Ye(1,2)*Ye(2,0)*Ye(2,2) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(
      2,2) + AeInput(0,1)*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput
      (2),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(1),
      mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)
      *Ye(0,2)*Ye(2,1)*Ye(2,2) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(
      2,2) + AeInput(1,1)*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (2),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*
      AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),
      mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + AeInput(2,1)*AeInput(2,2)*Sqr
      (Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(1))*Ye(1,1)
      *Ye(1,2)*Ye(2,1)*Ye(2,2) + Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput(0),
      mslInput(2))*Ye(0,0)*Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye
      (2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)
      *Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs(
      MuInput))*TCC0(mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*(Ye(0,0)
      *Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*(Ye(0,0)*Ye(0,2) + Ye(1
      ,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput
      (0),mslInput(2))*Ye(2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0
      )*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(2),mslInput(2),mslInput(0))*Ye(
      2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + Sqr(Abs
      (MuInput))*TCC0(mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*(Ye(0,1
      )*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*(Ye(0,1)*Ye(0,2) + Ye(1
      ,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput
      (1),mslInput(2))*Ye(1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1
      )*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(mseInput(1),mslInput(2),mslInput(1))*Ye(
      1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs
      (MuInput))*TCC0(mseInput(2),mslInput(1),mslInput(2))*Ye(2,1)*Ye(2,2)*(Ye(0,1
      )*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + Sqr(Abs(MuInput))*TCC0(
      mseInput(2),mslInput(2),mslInput(1))*Ye(2,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1
      ,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) - 3*TCD0(msdInput(0),msqInput(0),msqInput(1),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(0,1)*Yu(0,1)) - 3*TCD0(msdInput(0),msqInput(1),msqInput(0),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(0,1)*Yu(0,1)) - 3*TCD0(msdInput(1),msqInput(0),msqInput(1),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(0,1)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(0),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(0,1)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(1),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(0,1)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(0),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(0,1)) - 3*TCD0(msdInput(0),msqInput(0),msqInput(2),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(0),msqInput(2),msqInput(0),
      msuInput(0))*(AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(0,0)*Yu(0,0))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(0),msqInput(1),msqInput(2),
      msuInput(0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(0,1)*Yu(0,1))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(0),msqInput(2),msqInput(1),
      msuInput(0))*(AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(0,1)*Yu(0,1))*(AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(0,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(0),msqInput(2),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(0),
      msuInput(0))*(AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(0,0))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(2),
      msuInput(0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(0,1))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(1),
      msuInput(0))*(AdInput(1,1)*AuInput(0,1)*Yd(1,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(0,1))*(AdInput(1,2)*AuInput(0,2)*Yd(1,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(2),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(0),
      msuInput(0))*(AdInput(2,0)*AuInput(0,0)*Yd(2,0)*Yu(0,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(0,0))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(2),
      msuInput(0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(0,1))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(1),
      msuInput(0))*(AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(0,1))*(AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(0,2)) + TCC0(msdInput(0),msqInput(0),msuInput(0))*(-6*
      AdInput(0,0)*AuInput(0,0)*Yd(0,0)*Yu(0,0)*(Yd(0,0)*Yu(0,0) + Yd(0,1)*Yu(0,1)
      + Yd(0,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(0,0)*Yu(0,0)*(Yd(0,0)*Yu(0,0) +
      Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2))) + TCC0(msdInput(0),msqInput(1),msuInput(
      0))*(-6*AdInput(0,1)*AuInput(0,1)*Yd(0,1)*Yu(0,1)*(Yd(0,0)*Yu(0,0) + Yd(0,1)
      *Yu(0,1) + Yd(0,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(0,1)*Yu(0,1)*(Yd(0,0)*
      Yu(0,0) + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2))) + TCC0(msdInput(0),msqInput(2)
      ,msuInput(0))*(-6*AdInput(0,2)*AuInput(0,2)*Yd(0,2)*Yu(0,2)*(Yd(0,0)*Yu(0,0)
      + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(0,2)*Yu(0,2)*(
      Yd(0,0)*Yu(0,0) + Yd(0,1)*Yu(0,1) + Yd(0,2)*Yu(0,2))) + TCC0(msdInput(1),
      msqInput(0),msuInput(0))*(-6*AdInput(1,0)*AuInput(0,0)*Yd(1,0)*Yu(0,0)*(Yd(1
      ,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(1,
      0)*Yu(0,0)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2))) + TCC0(
      msdInput(1),msqInput(1),msuInput(0))*(-6*AdInput(1,1)*AuInput(0,1)*Yd(1,1)*
      Yu(0,1)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2)) + 6*Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(0,1)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,
      2))) + TCC0(msdInput(1),msqInput(2),msuInput(0))*(-6*AdInput(1,2)*AuInput(0,
      2)*Yd(1,2)*Yu(0,2)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1,2)*Yu(0,2)) + 6
      *Sqr(Abs(MuInput))*Yd(1,2)*Yu(0,2)*(Yd(1,0)*Yu(0,0) + Yd(1,1)*Yu(0,1) + Yd(1
      ,2)*Yu(0,2))) + TCC0(msdInput(2),msqInput(0),msuInput(0))*(-6*AdInput(2,0)*
      AuInput(0,0)*Yd(2,0)*Yu(0,0)*(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu
      (0,2)) + 6*Sqr(Abs(MuInput))*Yd(2,0)*Yu(0,0)*(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0
      ,1) + Yd(2,2)*Yu(0,2))) + TCC0(msdInput(2),msqInput(1),msuInput(0))*(-6*
      AdInput(2,1)*AuInput(0,1)*Yd(2,1)*Yu(0,1)*(Yd(2,0)*Yu(0,0) + Yd(2,1)*Yu(0,1)
      + Yd(2,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(2,1)*Yu(0,1)*(Yd(2,0)*Yu(0,0) +
      Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu(0,2))) + TCC0(msdInput(2),msqInput(2),msuInput(
      0))*(-6*AdInput(2,2)*AuInput(0,2)*Yd(2,2)*Yu(0,2)*(Yd(2,0)*Yu(0,0) + Yd(2,1)
      *Yu(0,1) + Yd(2,2)*Yu(0,2)) + 6*Sqr(Abs(MuInput))*Yd(2,2)*Yu(0,2)*(Yd(2,0)*
      Yu(0,0) + Yd(2,1)*Yu(0,1) + Yd(2,2)*Yu(0,2))) + 3*AuInput(1,0)*AuInput(1,1)*
      Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(1))*Yu(0
      ,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(1),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0
      )*Yu(1,1) + 3*AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),
      msqInput(0),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*
      AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),
      msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*TCD0(msdInput(0
      ),msqInput(0),msqInput(1),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu
      (1,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1
      )*Yu(1,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)) - 3*TCD0(msdInput(0),msqInput
      (1),msqInput(0),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0) -
      Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1
      ) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)) - 3*TCD0(msdInput(1),msqInput(0),
      msqInput(1),msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(
      Abs(MuInput))*Yd(1,0)*Yu(1,0))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) -
      Sqr(Abs(MuInput))*Yd(1,1)*Yu(1,1)) - 3*TCD0(msdInput(1),msqInput(1),msqInput
      (0),msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(Abs(
      MuInput))*Yd(1,0)*Yu(1,0))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) - Sqr(
      Abs(MuInput))*Yd(1,1)*Yu(1,1)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(1),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(1,1)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(0),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(
      MuInput))*Yd(2,1)*Yu(1,1)) + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(2),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0
      )*Yu(1,2) + 3*AuInput(0,0)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),
      msqInput(2),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*
      AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),
      msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*AuInput(0,0)*
      AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(1),
      msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*AuInput(1,1)*AuInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(0),msuInput(1))*Yu(0
      ,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(1),msqInput(2),msuInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1
      )*Yu(1,2) + 3*AuInput(1,1)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),
      msqInput(1),msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*
      AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),
      msuInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*TCD0(msdInput(0
      ),msqInput(0),msqInput(2),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu
      (1,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2
      )*Yu(1,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(0),msqInput
      (2),msqInput(0),msuInput(1))*(AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0) -
      Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2
      ) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(0),msqInput(1),
      msqInput(2),msuInput(1))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1) - Sqr(
      Abs(MuInput))*Yd(0,1)*Yu(1,1))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2) -
      Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(0),msqInput(2),msqInput
      (1),msuInput(1))*(AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1) - Sqr(Abs(
      MuInput))*Yd(0,1)*Yu(1,1))*(AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2) - Sqr(
      Abs(MuInput))*Yd(0,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(0),msqInput(2),
      msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(1,0))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(0),
      msuInput(1))*(AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(1,0)*Yu(1,0))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(2),
      msuInput(1))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(1,1))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(1),
      msuInput(1))*(AdInput(1,1)*AuInput(1,1)*Yd(1,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(1,1))*(AdInput(1,2)*AuInput(1,2)*Yd(1,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(2),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(0),
      msuInput(1))*(AdInput(2,0)*AuInput(1,0)*Yd(2,0)*Yu(1,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(1,0))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(2),
      msuInput(1))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(1,1))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(1),
      msuInput(1))*(AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(1,1))*(AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(1,2)) + TCC0(msdInput(0),msqInput(0),msuInput(1))*(-6*
      AdInput(0,0)*AuInput(1,0)*Yd(0,0)*Yu(1,0)*(Yd(0,0)*Yu(1,0) + Yd(0,1)*Yu(1,1)
      + Yd(0,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(0,0)*Yu(1,0)*(Yd(0,0)*Yu(1,0) +
      Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2))) + TCC0(msdInput(0),msqInput(1),msuInput(
      1))*(-6*AdInput(0,1)*AuInput(1,1)*Yd(0,1)*Yu(1,1)*(Yd(0,0)*Yu(1,0) + Yd(0,1)
      *Yu(1,1) + Yd(0,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(0,1)*Yu(1,1)*(Yd(0,0)*
      Yu(1,0) + Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2))) + TCC0(msdInput(0),msqInput(2)
      ,msuInput(1))*(-6*AdInput(0,2)*AuInput(1,2)*Yd(0,2)*Yu(1,2)*(Yd(0,0)*Yu(1,0)
      + Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(0,2)*Yu(1,2)*(
      Yd(0,0)*Yu(1,0) + Yd(0,1)*Yu(1,1) + Yd(0,2)*Yu(1,2))) + TCC0(msdInput(1),
      msqInput(0),msuInput(1))*(-6*AdInput(1,0)*AuInput(1,0)*Yd(1,0)*Yu(1,0)*(Yd(1
      ,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(1,
      0)*Yu(1,0)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2))) + TCC0(
      msdInput(1),msqInput(1),msuInput(1))*(-6*AdInput(1,1)*AuInput(1,1)*Yd(1,1)*
      Yu(1,1)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2)) + 6*Sqr(Abs(
      MuInput))*Yd(1,1)*Yu(1,1)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,
      2))) + TCC0(msdInput(1),msqInput(2),msuInput(1))*(-6*AdInput(1,2)*AuInput(1,
      2)*Yd(1,2)*Yu(1,2)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1,2)*Yu(1,2)) + 6
      *Sqr(Abs(MuInput))*Yd(1,2)*Yu(1,2)*(Yd(1,0)*Yu(1,0) + Yd(1,1)*Yu(1,1) + Yd(1
      ,2)*Yu(1,2))) + TCC0(msdInput(2),msqInput(0),msuInput(1))*(-6*AdInput(2,0)*
      AuInput(1,0)*Yd(2,0)*Yu(1,0)*(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu
      (1,2)) + 6*Sqr(Abs(MuInput))*Yd(2,0)*Yu(1,0)*(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1
      ,1) + Yd(2,2)*Yu(1,2))) + TCC0(msdInput(2),msqInput(1),msuInput(1))*(-6*
      AdInput(2,1)*AuInput(1,1)*Yd(2,1)*Yu(1,1)*(Yd(2,0)*Yu(1,0) + Yd(2,1)*Yu(1,1)
      + Yd(2,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(2,1)*Yu(1,1)*(Yd(2,0)*Yu(1,0) +
      Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu(1,2))) + TCC0(msdInput(2),msqInput(2),msuInput(
      1))*(-6*AdInput(2,2)*AuInput(1,2)*Yd(2,2)*Yu(1,2)*(Yd(2,0)*Yu(1,0) + Yd(2,1)
      *Yu(1,1) + Yd(2,2)*Yu(1,2)) + 6*Sqr(Abs(MuInput))*Yd(2,2)*Yu(1,2)*(Yd(2,0)*
      Yu(1,0) + Yd(2,1)*Yu(1,1) + Yd(2,2)*Yu(1,2))) + 3*AuInput(2,0)*AuInput(2,1)*
      Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(0),msuInput(2))*Yu(0
      ,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0
      )*Yu(2,1) + 3*AuInput(2,0)*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),
      msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*
      AuInput(0,0)*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0),
      msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(2,0)*
      AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(1),
      msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(1,0)*AuInput(1,1)*
      Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(1))*Yu(1
      ,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*AuInput(2,0)*AuInput(2,1)*Sqr(Abs(MuInput))*
      TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0
      )*Yu(2,1) + 3*AuInput(1,0)*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),
      msqInput(0),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*
      TCD0(msdInput(0),msqInput(0),msqInput(1),msuInput(2))*(AdInput(0,0)*AuInput(
      2,0)*Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,1)*
      AuInput(2,1)*Yd(0,1)*Yu(2,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)) - 3*TCD0(
      msdInput(0),msqInput(1),msqInput(0),msuInput(2))*(AdInput(0,0)*AuInput(2,0)*
      Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,1)*AuInput(2
      ,1)*Yd(0,1)*Yu(2,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)) - 3*TCD0(msdInput(1
      ),msqInput(0),msqInput(1),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu
      (2,0) - Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1
      )*Yu(2,1) - Sqr(Abs(MuInput))*Yd(1,1)*Yu(2,1)) - 3*TCD0(msdInput(1),msqInput
      (1),msqInput(0),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0) -
      Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1
      ) - Sqr(Abs(MuInput))*Yd(1,1)*Yu(2,1)) - 3*TCD0(msdInput(2),msqInput(0),
      msqInput(1),msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(
      Abs(MuInput))*Yd(2,0)*Yu(2,0))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) -
      Sqr(Abs(MuInput))*Yd(2,1)*Yu(2,1)) - 3*TCD0(msdInput(2),msqInput(1),msqInput
      (0),msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(Abs(
      MuInput))*Yd(2,0)*Yu(2,0))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) - Sqr(
      Abs(MuInput))*Yd(2,1)*Yu(2,1)) - 3*TCB0(msqInput(0),msqInput(1),SCALE)*(Yd(0
      ,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*(Yu(0,0)*Yu(0,1) + Yu(1,0)*
      Yu(1,1) + Yu(2,0)*Yu(2,1)) - 3*TCB0(msqInput(1),msqInput(0),SCALE)*(Yd(0,0)*
      Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1
      ,1) + Yu(2,0)*Yu(2,1)) + TCC0(msdInput(0),msqInput(0),msqInput(1))*(3*Sqr(
      Abs(MuInput))*Yd(0,0)*Yd(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*
      Yd(2,1)) - 3*AdInput(0,0)*AdInput(0,1)*Yd(0,0)*Yd(0,1)*(Yu(0,0)*Yu(0,1) + Yu
      (1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msdInput(0),msqInput(1),msqInput(0)
      )*(3*Sqr(Abs(MuInput))*Yd(0,0)*Yd(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) +
      Yd(2,0)*Yd(2,1)) - 3*AdInput(0,0)*AdInput(0,1)*Yd(0,0)*Yd(0,1)*(Yu(0,0)*Yu(0
      ,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msdInput(1),msqInput(0),
      msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)
      *Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*AdInput(1,0)*AdInput(1,1)*Yd(1,0)*Yd(1,1)*(
      Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msdInput(1),
      msqInput(1),msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,
      1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*AdInput(1,0)*AdInput(1,1)*Yd(1,0
      )*Yd(1,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(
      msdInput(2),msqInput(0),msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(2,0)*Yd(2,1)*(
      Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*AdInput(2,0)*
      AdInput(2,1)*Yd(2,0)*Yd(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu
      (2,1))) + TCC0(msdInput(2),msqInput(1),msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(
      2,0)*Yd(2,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) - 3*
      AdInput(2,0)*AdInput(2,1)*Yd(2,0)*Yd(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1)
      + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(0),msqInput(1),msuInput(0))*(-3*AuInput(
      0,0)*AuInput(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(0
      ,0)*Yu(0,1) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)
      *Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(1),msqInput(0),msuInput(0))*(-3
      *AuInput(0,0)*AuInput(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2
      ,1))*Yu(0,0)*Yu(0,1) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,1)*(Yu(0,0)*Yu(0,1)
      + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(0),msqInput(1),
      msuInput(1))*(-3*AuInput(1,0)*AuInput(1,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1
      ) + Yd(2,0)*Yd(2,1))*Yu(1,0)*Yu(1,1) + 3*Sqr(Abs(MuInput))*Yu(1,0)*Yu(1,1)*(
      Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(msqInput(1),
      msqInput(0),msuInput(1))*(-3*AuInput(1,0)*AuInput(1,1)*(Yd(0,0)*Yd(0,1) + Yd
      (1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(1,0)*Yu(1,1) + 3*Sqr(Abs(MuInput))*Yu(1,
      0)*Yu(1,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1))) + TCC0(
      msqInput(0),msqInput(1),msuInput(2))*(-3*AuInput(2,0)*AuInput(2,1)*(Yd(0,0)*
      Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(2,0)*Yu(2,1) + 3*Sqr(Abs(
      MuInput))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,
      1))) + TCC0(msqInput(1),msqInput(0),msuInput(2))*(-3*AuInput(2,0)*AuInput(2,
      1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1))*Yu(2,0)*Yu(2,1) + 3
      *Sqr(Abs(MuInput))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2
      ,0)*Yu(2,1))) + 3*AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(
      0),msqInput(2),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*
      AuInput(0,0)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),
      msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(2,0)*
      AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(0),
      msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(0,0)*AuInput(0,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(0))*Yu(0
      ,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0
      )*Yu(2,2) + 3*AuInput(1,0)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(0),
      msqInput(2),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*
      AuInput(2,0)*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),
      msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(1,0)*
      AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0),msuInput(2),
      msuInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*AuInput(2,1)*AuInput(2,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(0),msuInput(2))*Yu(0
      ,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1
      )*Yu(2,2) + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),
      msqInput(1),msuInput(0),msuInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*
      AuInput(0,1)*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1),
      msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(2,1)*
      AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(1),
      msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(1,1)*AuInput(1,2)*
      Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(1))*Yu(1
      ,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*AuInput(2,1)*AuInput(2,2)*Sqr(Abs(MuInput))*
      TCD0(msqInput(2),msqInput(1),msuInput(1),msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1
      )*Yu(2,2) + 3*AuInput(1,1)*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),
      msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*
      TCD0(msdInput(0),msqInput(0),msqInput(2),msuInput(2))*(AdInput(0,0)*AuInput(
      2,0)*Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,2)*
      AuInput(2,2)*Yd(0,2)*Yu(2,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(
      msdInput(0),msqInput(2),msqInput(0),msuInput(2))*(AdInput(0,0)*AuInput(2,0)*
      Yd(0,0)*Yu(2,0) - Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0))*(AdInput(0,2)*AuInput(2
      ,2)*Yd(0,2)*Yu(2,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(msdInput(0
      ),msqInput(1),msqInput(2),msuInput(2))*(AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu
      (2,1) - Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1))*(AdInput(0,2)*AuInput(2,2)*Yd(0,2
      )*Yu(2,2) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(msdInput(0),msqInput
      (2),msqInput(1),msuInput(2))*(AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu(2,1) -
      Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1))*(AdInput(0,2)*AuInput(2,2)*Yd(0,2)*Yu(2,2
      ) - Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(0),
      msqInput(2),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0) - Sqr(
      Abs(MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) -
      Sqr(Abs(MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput
      (0),msuInput(2))*(AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0) - Sqr(Abs(
      MuInput))*Yd(1,0)*Yu(2,0))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) - Sqr(
      Abs(MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(1),msqInput(2),
      msuInput(2))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(2,1))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(1),msqInput(2),msqInput(1),
      msuInput(2))*(AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(1,1)*Yu(2,1))*(AdInput(1,2)*AuInput(2,2)*Yd(1,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(0),msqInput(2),
      msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(2,0))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(0),
      msuInput(2))*(AdInput(2,0)*AuInput(2,0)*Yd(2,0)*Yu(2,0) - Sqr(Abs(MuInput))*
      Yd(2,0)*Yu(2,0))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(1),msqInput(2),
      msuInput(2))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(2,1))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCD0(msdInput(2),msqInput(2),msqInput(1),
      msuInput(2))*(AdInput(2,1)*AuInput(2,1)*Yd(2,1)*Yu(2,1) - Sqr(Abs(MuInput))*
      Yd(2,1)*Yu(2,1))*(AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2) - Sqr(Abs(
      MuInput))*Yd(2,2)*Yu(2,2)) - 3*TCB0(msqInput(0),msqInput(2),SCALE)*(Yd(0,0)*
      Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1
      ,2) + Yu(2,0)*Yu(2,2)) - 3*TCB0(msqInput(2),msqInput(0),SCALE)*(Yd(0,0)*Yd(0
      ,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2)
      + Yu(2,0)*Yu(2,2)) - 3*TCB0(msqInput(1),msqInput(2),SCALE)*(Yd(0,1)*Yd(0,2)
      + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu
      (2,1)*Yu(2,2)) - 3*TCB0(msqInput(2),msqInput(1),SCALE)*(Yd(0,1)*Yd(0,2) + Yd
      (1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1
      )*Yu(2,2)) + TCC0(msdInput(0),msqInput(0),msuInput(2))*(-6*AdInput(0,0)*
      AuInput(2,0)*Yd(0,0)*Yu(2,0)*(Yd(0,0)*Yu(2,0) + Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu
      (2,2)) + 6*Sqr(Abs(MuInput))*Yd(0,0)*Yu(2,0)*(Yd(0,0)*Yu(2,0) + Yd(0,1)*Yu(2
      ,1) + Yd(0,2)*Yu(2,2))) + TCC0(msdInput(0),msqInput(1),msuInput(2))*(-6*
      AdInput(0,1)*AuInput(2,1)*Yd(0,1)*Yu(2,1)*(Yd(0,0)*Yu(2,0) + Yd(0,1)*Yu(2,1)
      + Yd(0,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(0,1)*Yu(2,1)*(Yd(0,0)*Yu(2,0) +
      Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu(2,2))) + TCC0(msdInput(0),msqInput(2),msuInput(
      2))*(-6*AdInput(0,2)*AuInput(2,2)*Yd(0,2)*Yu(2,2)*(Yd(0,0)*Yu(2,0) + Yd(0,1)
      *Yu(2,1) + Yd(0,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(0,2)*Yu(2,2)*(Yd(0,0)*
      Yu(2,0) + Yd(0,1)*Yu(2,1) + Yd(0,2)*Yu(2,2))) + TCC0(msdInput(1),msqInput(0)
      ,msuInput(2))*(-6*AdInput(1,0)*AuInput(2,0)*Yd(1,0)*Yu(2,0)*(Yd(1,0)*Yu(2,0)
      + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(1,0)*Yu(2,0)*(
      Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2))) + TCC0(msdInput(1),
      msqInput(1),msuInput(2))*(-6*AdInput(1,1)*AuInput(2,1)*Yd(1,1)*Yu(2,1)*(Yd(1
      ,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(1,
      1)*Yu(2,1)*(Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2))) + TCC0(
      msdInput(1),msqInput(2),msuInput(2))*(-6*AdInput(1,2)*AuInput(2,2)*Yd(1,2)*
      Yu(2,2)*(Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,2)) + 6*Sqr(Abs(
      MuInput))*Yd(1,2)*Yu(2,2)*(Yd(1,0)*Yu(2,0) + Yd(1,1)*Yu(2,1) + Yd(1,2)*Yu(2,
      2))) + TCC0(msdInput(2),msqInput(0),msuInput(2))*(-6*AdInput(2,0)*AuInput(2,
      0)*Yd(2,0)*Yu(2,0)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1) + Yd(2,2)*Yu(2,2)) + 6
      *Sqr(Abs(MuInput))*Yd(2,0)*Yu(2,0)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1) + Yd(2
      ,2)*Yu(2,2))) + TCC0(msdInput(2),msqInput(1),msuInput(2))*(-6*AdInput(2,1)*
      AuInput(2,1)*Yd(2,1)*Yu(2,1)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1) + Yd(2,2)*Yu
      (2,2)) + 6*Sqr(Abs(MuInput))*Yd(2,1)*Yu(2,1)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2
      ,1) + Yd(2,2)*Yu(2,2))) + TCC0(msdInput(2),msqInput(2),msuInput(2))*(-6*
      AdInput(2,2)*AuInput(2,2)*Yd(2,2)*Yu(2,2)*(Yd(2,0)*Yu(2,0) + Yd(2,1)*Yu(2,1)
      + Yd(2,2)*Yu(2,2)) + 6*Sqr(Abs(MuInput))*Yd(2,2)*Yu(2,2)*(Yd(2,0)*Yu(2,0) +
      Yd(2,1)*Yu(2,1) + Yd(2,2)*Yu(2,2))) + TCC0(msdInput(0),msqInput(0),msqInput(
      2))*(3*Sqr(Abs(MuInput))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2)
      + Yd(2,0)*Yd(2,2)) - 3*AdInput(0,0)*AdInput(0,2)*Yd(0,0)*Yd(0,2)*(Yu(0,0)*Yu
      (0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(0),msqInput(2),
      msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)
      *Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*AdInput(0,0)*AdInput(0,2)*Yd(0,0)*Yd(0,2)*(
      Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(1),
      msqInput(0),msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,2)*(Yd(0,0)*Yd(0,
      2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*AdInput(1,0)*AdInput(1,2)*Yd(1,0
      )*Yd(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(
      msdInput(1),msqInput(2),msqInput(0))*(3*Sqr(Abs(MuInput))*Yd(1,0)*Yd(1,2)*(
      Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*AdInput(1,0)*
      AdInput(1,2)*Yd(1,0)*Yd(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu
      (2,2))) + TCC0(msdInput(2),msqInput(0),msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(
      2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) - 3*
      AdInput(2,0)*AdInput(2,2)*Yd(2,0)*Yd(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2)
      + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(2),msqInput(2),msqInput(0))*(3*Sqr(Abs(
      MuInput))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,
      2)) - 3*AdInput(2,0)*AdInput(2,2)*Yd(2,0)*Yd(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)
      *Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msqInput(0),msqInput(2),msuInput(0))*(-3
      *AuInput(0,0)*AuInput(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2
      ,2))*Yu(0,0)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,2)*(Yu(0,0)*Yu(0,2)
      + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msqInput(2),msqInput(0),
      msuInput(0))*(-3*AuInput(0,0)*AuInput(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2
      ) + Yd(2,0)*Yd(2,2))*Yu(0,0)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,0)*Yu(0,2)*(
      Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msqInput(0),
      msqInput(2),msuInput(1))*(-3*AuInput(1,0)*AuInput(1,2)*(Yd(0,0)*Yd(0,2) + Yd
      (1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(1,0)*Yu(1,2) + 3*Sqr(Abs(MuInput))*Yu(1,
      0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2))) + TCC0(
      msqInput(2),msqInput(0),msuInput(1))*(-3*AuInput(1,0)*AuInput(1,2)*(Yd(0,0)*
      Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(1,0)*Yu(1,2) + 3*Sqr(Abs(
      MuInput))*Yu(1,0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,
      2))) + TCC0(msqInput(0),msqInput(2),msuInput(2))*(-3*AuInput(2,0)*AuInput(2,
      2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(2,0)*Yu(2,2) + 3
      *Sqr(Abs(MuInput))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2
      ,0)*Yu(2,2))) + TCC0(msqInput(2),msqInput(0),msuInput(2))*(-3*AuInput(2,0)*
      AuInput(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2))*Yu(2,0)*
      Yu(2,2) + 3*Sqr(Abs(MuInput))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(
      1,2) + Yu(2,0)*Yu(2,2))) + TCC0(msdInput(0),msqInput(1),msqInput(2))*(3*Sqr(
      Abs(MuInput))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*
      Yd(2,2)) - 3*AdInput(0,1)*AdInput(0,2)*Yd(0,1)*Yd(0,2)*(Yu(0,1)*Yu(0,2) + Yu
      (1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msdInput(0),msqInput(2),msqInput(1)
      )*(3*Sqr(Abs(MuInput))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) +
      Yd(2,1)*Yd(2,2)) - 3*AdInput(0,1)*AdInput(0,2)*Yd(0,1)*Yd(0,2)*(Yu(0,1)*Yu(0
      ,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msdInput(1),msqInput(1),
      msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(1,1)*Yd(1,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)
      *Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*AdInput(1,1)*AdInput(1,2)*Yd(1,1)*Yd(1,2)*(
      Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msdInput(1),
      msqInput(2),msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(1,1)*Yd(1,2)*(Yd(0,1)*Yd(0,
      2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*AdInput(1,1)*AdInput(1,2)*Yd(1,1
      )*Yd(1,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(
      msdInput(2),msqInput(1),msqInput(2))*(3*Sqr(Abs(MuInput))*Yd(2,1)*Yd(2,2)*(
      Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*AdInput(2,1)*
      AdInput(2,2)*Yd(2,1)*Yd(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu
      (2,2))) + TCC0(msdInput(2),msqInput(2),msqInput(1))*(3*Sqr(Abs(MuInput))*Yd(
      2,1)*Yd(2,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) - 3*
      AdInput(2,1)*AdInput(2,2)*Yd(2,1)*Yd(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2)
      + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(1),msqInput(2),msuInput(0))*(-3*AuInput(
      0,1)*AuInput(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(0
      ,1)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)
      *Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(2),msqInput(1),msuInput(0))*(-3
      *AuInput(0,1)*AuInput(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2
      ,2))*Yu(0,1)*Yu(0,2) + 3*Sqr(Abs(MuInput))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2)
      + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(1),msqInput(2),
      msuInput(1))*(-3*AuInput(1,1)*AuInput(1,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2
      ) + Yd(2,1)*Yd(2,2))*Yu(1,1)*Yu(1,2) + 3*Sqr(Abs(MuInput))*Yu(1,1)*Yu(1,2)*(
      Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(msqInput(2),
      msqInput(1),msuInput(1))*(-3*AuInput(1,1)*AuInput(1,2)*(Yd(0,1)*Yd(0,2) + Yd
      (1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(1,1)*Yu(1,2) + 3*Sqr(Abs(MuInput))*Yu(1,
      1)*Yu(1,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2))) + TCC0(
      msqInput(1),msqInput(2),msuInput(2))*(-3*AuInput(2,1)*AuInput(2,2)*(Yd(0,1)*
      Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(2,1)*Yu(2,2) + 3*Sqr(Abs(
      MuInput))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,
      2))) + TCC0(msqInput(2),msqInput(1),msuInput(2))*(-3*AuInput(2,1)*AuInput(2,
      2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2))*Yu(2,1)*Yu(2,2) + 3
      *Sqr(Abs(MuInput))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2
      ,1)*Yu(2,2)))))));
   MODEL->set_Lambda5(Re(0.006332573977646111*UnitStep(-1 + LambdaLoopOrder)*(0.36
      *Quad(g1)*Sqr(M1Input)*Sqr(MuInput)*TCD0(M1Input,M1Input,Abs(MuInput),Abs(
      MuInput)) + 1.2*M1Input*M2Input*Sqr(g1)*Sqr(g2)*Sqr(MuInput)*TCD0(M1Input,
      M2Input,Abs(MuInput),Abs(MuInput)) + 3*Quad(g2)*Sqr(M2Input)*Sqr(MuInput)*
      TCD0(M2Input,M2Input,Abs(MuInput),Abs(MuInput)) - 3*Quad(Yd(0,0))*Sqr(
      MuInput)*Sqr(AdInput(0,0))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput
      (0)) - 3*AdInput(0,0)*AdInput(0,1)*Sqr(MuInput)*Sqr(Yd(0,0))*Sqr(Yd(0,1))*
      TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(1)) - 3*AdInput(0,0)*
      AdInput(0,2)*Sqr(MuInput)*Sqr(Yd(0,0))*Sqr(Yd(0,2))*TCD0(msdInput(0),
      msdInput(0),msqInput(0),msqInput(2)) - 3*AdInput(0,0)*AdInput(0,1)*Sqr(
      MuInput)*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1),
      msqInput(0)) - 3*Quad(Yd(0,1))*Sqr(MuInput)*Sqr(AdInput(0,1))*TCD0(msdInput(
      0),msdInput(0),msqInput(1),msqInput(1)) - 3*AdInput(0,1)*AdInput(0,2)*Sqr(
      MuInput)*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(1),
      msqInput(2)) - 3*AdInput(0,0)*AdInput(0,2)*Sqr(MuInput)*Sqr(Yd(0,0))*Sqr(Yd(
      0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),msqInput(0)) - 3*AdInput(0,1)
      *AdInput(0,2)*Sqr(MuInput)*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),
      msdInput(0),msqInput(2),msqInput(1)) - 3*Quad(Yd(0,2))*Sqr(MuInput)*Sqr(
      AdInput(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),msqInput(2)) - 3*
      AdInput(0,0)*AdInput(1,0)*Sqr(MuInput)*Sqr(Yd(0,0))*Sqr(Yd(1,0))*TCD0(
      msdInput(0),msdInput(1),msqInput(0),msqInput(0)) - 3*AdInput(0,1)*AdInput(1,
      1)*Sqr(MuInput)*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(0),msdInput(1),
      msqInput(1),msqInput(1)) - 3*AdInput(0,2)*AdInput(1,2)*Sqr(MuInput)*Sqr(Yd(0
      ,2))*Sqr(Yd(1,2))*TCD0(msdInput(0),msdInput(1),msqInput(2),msqInput(2)) - 3*
      AdInput(0,0)*AdInput(2,0)*Sqr(MuInput)*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(
      msdInput(0),msdInput(2),msqInput(0),msqInput(0)) - 3*AdInput(0,1)*AdInput(2,
      1)*Sqr(MuInput)*Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(0),msdInput(2),
      msqInput(1),msqInput(1)) - 3*AdInput(0,2)*AdInput(2,2)*Sqr(MuInput)*Sqr(Yd(0
      ,2))*Sqr(Yd(2,2))*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput(2)) - 3*
      AdInput(0,0)*AdInput(1,0)*Sqr(MuInput)*Sqr(Yd(0,0))*Sqr(Yd(1,0))*TCD0(
      msdInput(1),msdInput(0),msqInput(0),msqInput(0)) - 3*AdInput(0,1)*AdInput(1,
      1)*Sqr(MuInput)*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(0),
      msqInput(1),msqInput(1)) - 3*AdInput(0,2)*AdInput(1,2)*Sqr(MuInput)*Sqr(Yd(0
      ,2))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(2)) - 3*
      Quad(Yd(1,0))*Sqr(MuInput)*Sqr(AdInput(1,0))*TCD0(msdInput(1),msdInput(1),
      msqInput(0),msqInput(0)) - 3*AdInput(1,0)*AdInput(1,1)*Sqr(MuInput)*Sqr(Yd(1
      ,0))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(0),msqInput(1)) - 3*
      AdInput(1,0)*AdInput(1,2)*Sqr(MuInput)*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(
      msdInput(1),msdInput(1),msqInput(0),msqInput(2)) - 3*AdInput(1,0)*AdInput(1,
      1)*Sqr(MuInput)*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1),
      msqInput(1),msqInput(0)) - 3*Quad(Yd(1,1))*Sqr(MuInput)*Sqr(AdInput(1,1))*
      TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(1)) - 3*AdInput(1,1)*
      AdInput(1,2)*Sqr(MuInput)*Sqr(Yd(1,1))*Sqr(Yd(1,2))*TCD0(msdInput(1),
      msdInput(1),msqInput(1),msqInput(2)) - 3*AdInput(1,0)*AdInput(1,2)*Sqr(
      MuInput)*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2),
      msqInput(0)) - 3*AdInput(1,1)*AdInput(1,2)*Sqr(MuInput)*Sqr(Yd(1,1))*Sqr(Yd(
      1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2),msqInput(1)) - 3*Quad(Yd(1,2)
      )*Sqr(MuInput)*Sqr(AdInput(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2),
      msqInput(2)) - 3*AdInput(1,0)*AdInput(2,0)*Sqr(MuInput)*Sqr(Yd(1,0))*Sqr(Yd(
      2,0))*TCD0(msdInput(1),msdInput(2),msqInput(0),msqInput(0)) - 3*AdInput(1,1)
      *AdInput(2,1)*Sqr(MuInput)*Sqr(Yd(1,1))*Sqr(Yd(2,1))*TCD0(msdInput(1),
      msdInput(2),msqInput(1),msqInput(1)) - 3*AdInput(1,2)*AdInput(2,2)*Sqr(
      MuInput)*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(msdInput(1),msdInput(2),msqInput(2),
      msqInput(2)) - 3*AdInput(0,0)*AdInput(2,0)*Sqr(MuInput)*Sqr(Yd(0,0))*Sqr(Yd(
      2,0))*TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(0)) - 3*AdInput(0,1)
      *AdInput(2,1)*Sqr(MuInput)*Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(0),msqInput(1),msqInput(1)) - 3*AdInput(0,2)*AdInput(2,2)*Sqr(
      MuInput)*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(0),msqInput(2),
      msqInput(2)) - 3*AdInput(1,0)*AdInput(2,0)*Sqr(MuInput)*Sqr(Yd(1,0))*Sqr(Yd(
      2,0))*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(0)) - 3*AdInput(1,1)
      *AdInput(2,1)*Sqr(MuInput)*Sqr(Yd(1,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(1),msqInput(1),msqInput(1)) - 3*AdInput(1,2)*AdInput(2,2)*Sqr(
      MuInput)*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(1),msqInput(2),
      msqInput(2)) - 3*Quad(Yd(2,0))*Sqr(MuInput)*Sqr(AdInput(2,0))*TCD0(msdInput(
      2),msdInput(2),msqInput(0),msqInput(0)) - 3*AdInput(2,0)*AdInput(2,1)*Sqr(
      MuInput)*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(0),
      msqInput(1)) - 3*AdInput(2,0)*AdInput(2,2)*Sqr(MuInput)*Sqr(Yd(2,0))*Sqr(Yd(
      2,2))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(2)) - 3*AdInput(2,0)
      *AdInput(2,1)*Sqr(MuInput)*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),
      msdInput(2),msqInput(1),msqInput(0)) - 3*Quad(Yd(2,1))*Sqr(MuInput)*Sqr(
      AdInput(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(1)) - 3*
      AdInput(2,1)*AdInput(2,2)*Sqr(MuInput)*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(2),msqInput(1),msqInput(2)) - 3*AdInput(2,0)*AdInput(2,
      2)*Sqr(MuInput)*Sqr(Yd(2,0))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(2),msqInput(0)) - 3*AdInput(2,1)*AdInput(2,2)*Sqr(MuInput)*Sqr(Yd(2
      ,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(1)) - 3*
      Quad(Yd(2,2))*Sqr(MuInput)*Sqr(AdInput(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(2),msqInput(2)) - Quad(Ye(0,0))*Sqr(MuInput)*Sqr(AeInput(0,0))*TCD0
      (mseInput(0),mseInput(0),mslInput(0),mslInput(0)) - AeInput(0,0)*AeInput(0,1
      )*Sqr(MuInput)*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(mseInput(0),mseInput(0),
      mslInput(0),mslInput(1)) - AeInput(0,0)*AeInput(0,2)*Sqr(MuInput)*Sqr(Ye(0,0
      ))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(2)) -
      AeInput(0,0)*AeInput(0,1)*Sqr(MuInput)*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(
      mseInput(0),mseInput(0),mslInput(1),mslInput(0)) - Quad(Ye(0,1))*Sqr(MuInput
      )*Sqr(AeInput(0,1))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(1)) -
      AeInput(0,1)*AeInput(0,2)*Sqr(MuInput)*Sqr(Ye(0,1))*Sqr(Ye(0,2))*TCD0(
      mseInput(0),mseInput(0),mslInput(1),mslInput(2)) - AeInput(0,0)*AeInput(0,2)
      *Sqr(MuInput)*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),
      mslInput(2),mslInput(0)) - AeInput(0,1)*AeInput(0,2)*Sqr(MuInput)*Sqr(Ye(0,1
      ))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(1)) - Quad
      (Ye(0,2))*Sqr(MuInput)*Sqr(AeInput(0,2))*TCD0(mseInput(0),mseInput(0),
      mslInput(2),mslInput(2)) - AeInput(0,0)*AeInput(1,0)*Sqr(MuInput)*Sqr(Ye(0,0
      ))*Sqr(Ye(1,0))*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(0)) -
      AeInput(0,1)*AeInput(1,1)*Sqr(MuInput)*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0(
      mseInput(0),mseInput(1),mslInput(1),mslInput(1)) - AeInput(0,2)*AeInput(1,2)
      *Sqr(MuInput)*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(0),mseInput(1),
      mslInput(2),mslInput(2)) - AeInput(0,0)*AeInput(2,0)*Sqr(MuInput)*Sqr(Ye(0,0
      ))*Sqr(Ye(2,0))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(0)) -
      AeInput(0,1)*AeInput(2,1)*Sqr(MuInput)*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(1)) - AeInput(0,2)*AeInput(2,2)
      *Sqr(MuInput)*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(2)) - AeInput(0,0)*AeInput(1,0)*Sqr(MuInput)*Sqr(Ye(0,0
      ))*Sqr(Ye(1,0))*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(0)) -
      AeInput(0,1)*AeInput(1,1)*Sqr(MuInput)*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0(
      mseInput(1),mseInput(0),mslInput(1),mslInput(1)) - AeInput(0,2)*AeInput(1,2)
      *Sqr(MuInput)*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(0),
      mslInput(2),mslInput(2)) - Quad(Ye(1,0))*Sqr(MuInput)*Sqr(AeInput(1,0))*TCD0
      (mseInput(1),mseInput(1),mslInput(0),mslInput(0)) - AeInput(1,0)*AeInput(1,1
      )*Sqr(MuInput)*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1),mseInput(1),
      mslInput(0),mslInput(1)) - AeInput(1,0)*AeInput(1,2)*Sqr(MuInput)*Sqr(Ye(1,0
      ))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(2)) -
      AeInput(1,0)*AeInput(1,1)*Sqr(MuInput)*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(
      mseInput(1),mseInput(1),mslInput(1),mslInput(0)) - Quad(Ye(1,1))*Sqr(MuInput
      )*Sqr(AeInput(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(1)) -
      AeInput(1,1)*AeInput(1,2)*Sqr(MuInput)*Sqr(Ye(1,1))*Sqr(Ye(1,2))*TCD0(
      mseInput(1),mseInput(1),mslInput(1),mslInput(2)) - AeInput(1,0)*AeInput(1,2)
      *Sqr(MuInput)*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(2),mslInput(0)) - AeInput(1,1)*AeInput(1,2)*Sqr(MuInput)*Sqr(Ye(1,1
      ))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(2),mslInput(1)) - Quad
      (Ye(1,2))*Sqr(MuInput)*Sqr(AeInput(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(2),mslInput(2)) - AeInput(1,0)*AeInput(2,0)*Sqr(MuInput)*Sqr(Ye(1,0
      ))*Sqr(Ye(2,0))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(0)) -
      AeInput(1,1)*AeInput(2,1)*Sqr(MuInput)*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(1),mseInput(2),mslInput(1),mslInput(1)) - AeInput(1,2)*AeInput(2,2)
      *Sqr(MuInput)*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(mseInput(1),mseInput(2),
      mslInput(2),mslInput(2)) - AeInput(0,0)*AeInput(2,0)*Sqr(MuInput)*Sqr(Ye(0,0
      ))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(0)) -
      AeInput(0,1)*AeInput(2,1)*Sqr(MuInput)*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(2),mseInput(0),mslInput(1),mslInput(1)) - AeInput(0,2)*AeInput(2,2)
      *Sqr(MuInput)*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(0),
      mslInput(2),mslInput(2)) - AeInput(1,0)*AeInput(2,0)*Sqr(MuInput)*Sqr(Ye(1,0
      ))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(1),mslInput(0),mslInput(0)) -
      AeInput(1,1)*AeInput(2,1)*Sqr(MuInput)*Sqr(Ye(1,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(2),mseInput(1),mslInput(1),mslInput(1)) - AeInput(1,2)*AeInput(2,2)
      *Sqr(MuInput)*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(1),
      mslInput(2),mslInput(2)) - Quad(Ye(2,0))*Sqr(MuInput)*Sqr(AeInput(2,0))*TCD0
      (mseInput(2),mseInput(2),mslInput(0),mslInput(0)) - AeInput(2,0)*AeInput(2,1
      )*Sqr(MuInput)*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),
      mslInput(0),mslInput(1)) - AeInput(2,0)*AeInput(2,2)*Sqr(MuInput)*Sqr(Ye(2,0
      ))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(0),mslInput(2)) -
      AeInput(2,0)*AeInput(2,1)*Sqr(MuInput)*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(
      mseInput(2),mseInput(2),mslInput(1),mslInput(0)) - Quad(Ye(2,1))*Sqr(MuInput
      )*Sqr(AeInput(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(1)) -
      AeInput(2,1)*AeInput(2,2)*Sqr(MuInput)*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(
      mseInput(2),mseInput(2),mslInput(1),mslInput(2)) - AeInput(2,0)*AeInput(2,2)
      *Sqr(MuInput)*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),
      mslInput(2),mslInput(0)) - AeInput(2,1)*AeInput(2,2)*Sqr(MuInput)*Sqr(Ye(2,1
      ))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(1)) - Quad
      (Ye(2,2))*Sqr(MuInput)*Sqr(AeInput(2,2))*TCD0(mseInput(2),mseInput(2),
      mslInput(2),mslInput(2)) - 3*Quad(Yu(0,0))*Sqr(MuInput)*Sqr(AuInput(0,0))*
      TCD0(msqInput(0),msqInput(0),msuInput(0),msuInput(0)) - 3*AuInput(0,0)*
      AuInput(1,0)*Sqr(MuInput)*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(0),msuInput(1)) - 3*AuInput(0,0)*AuInput(2,0)*Sqr(
      MuInput)*Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(0),
      msuInput(2)) - 3*AuInput(0,0)*AuInput(1,0)*Sqr(MuInput)*Sqr(Yu(0,0))*Sqr(Yu(
      1,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),msuInput(0)) - 3*Quad(Yu(1,0)
      )*Sqr(MuInput)*Sqr(AuInput(1,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),
      msuInput(1)) - 3*AuInput(1,0)*AuInput(2,0)*Sqr(MuInput)*Sqr(Yu(1,0))*Sqr(Yu(
      2,0))*TCD0(msqInput(0),msqInput(0),msuInput(1),msuInput(2)) - 3*AuInput(0,0)
      *AuInput(2,0)*Sqr(MuInput)*Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),
      msqInput(0),msuInput(2),msuInput(0)) - 3*AuInput(1,0)*AuInput(2,0)*Sqr(
      MuInput)*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msuInput(2),
      msuInput(1)) - 3*Quad(Yu(2,0))*Sqr(MuInput)*Sqr(AuInput(2,0))*TCD0(msqInput(
      0),msqInput(0),msuInput(2),msuInput(2)) - 3*AuInput(0,0)*AuInput(0,1)*Sqr(
      MuInput)*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(0),msqInput(1),msuInput(0),
      msuInput(0)) - 3*AuInput(1,0)*AuInput(1,1)*Sqr(MuInput)*Sqr(Yu(1,0))*Sqr(Yu(
      1,1))*TCD0(msqInput(0),msqInput(1),msuInput(1),msuInput(1)) - 3*AuInput(2,0)
      *AuInput(2,1)*Sqr(MuInput)*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(0),
      msqInput(1),msuInput(2),msuInput(2)) - 3*AuInput(0,0)*AuInput(0,2)*Sqr(
      MuInput)*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(0),msqInput(2),msuInput(0),
      msuInput(0)) - 3*AuInput(1,0)*AuInput(1,2)*Sqr(MuInput)*Sqr(Yu(1,0))*Sqr(Yu(
      1,2))*TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(1)) - 3*AuInput(2,0)
      *AuInput(2,2)*Sqr(MuInput)*Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(msqInput(0),
      msqInput(2),msuInput(2),msuInput(2)) - 3*AuInput(0,0)*AuInput(0,1)*Sqr(
      MuInput)*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(1),msqInput(0),msuInput(0),
      msuInput(0)) - 3*AuInput(1,0)*AuInput(1,1)*Sqr(MuInput)*Sqr(Yu(1,0))*Sqr(Yu(
      1,1))*TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(1)) - 3*AuInput(2,0)
      *AuInput(2,1)*Sqr(MuInput)*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(1),
      msqInput(0),msuInput(2),msuInput(2)) - 3*Quad(Yu(0,1))*Sqr(MuInput)*Sqr(
      AuInput(0,1))*TCD0(msqInput(1),msqInput(1),msuInput(0),msuInput(0)) - 3*
      AuInput(0,1)*AuInput(1,1)*Sqr(MuInput)*Sqr(Yu(0,1))*Sqr(Yu(1,1))*TCD0(
      msqInput(1),msqInput(1),msuInput(0),msuInput(1)) - 3*AuInput(0,1)*AuInput(2,
      1)*Sqr(MuInput)*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(0),msuInput(2)) - 3*AuInput(0,1)*AuInput(1,1)*Sqr(MuInput)*Sqr(Yu(0
      ,1))*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(0)) - 3*
      Quad(Yu(1,1))*Sqr(MuInput)*Sqr(AuInput(1,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(1),msuInput(1)) - 3*AuInput(1,1)*AuInput(2,1)*Sqr(MuInput)*Sqr(Yu(1
      ,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msuInput(1),msuInput(2)) - 3*
      AuInput(0,1)*AuInput(2,1)*Sqr(MuInput)*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(
      msqInput(1),msqInput(1),msuInput(2),msuInput(0)) - 3*AuInput(1,1)*AuInput(2,
      1)*Sqr(MuInput)*Sqr(Yu(1,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),
      msuInput(2),msuInput(1)) - 3*Quad(Yu(2,1))*Sqr(MuInput)*Sqr(AuInput(2,1))*
      TCD0(msqInput(1),msqInput(1),msuInput(2),msuInput(2)) - 3*AuInput(0,1)*
      AuInput(0,2)*Sqr(MuInput)*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(1),
      msqInput(2),msuInput(0),msuInput(0)) - 3*AuInput(1,1)*AuInput(1,2)*Sqr(
      MuInput)*Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(1),msqInput(2),msuInput(1),
      msuInput(1)) - 3*AuInput(2,1)*AuInput(2,2)*Sqr(MuInput)*Sqr(Yu(2,1))*Sqr(Yu(
      2,2))*TCD0(msqInput(1),msqInput(2),msuInput(2),msuInput(2)) - 3*AuInput(0,0)
      *AuInput(0,2)*Sqr(MuInput)*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(2),
      msqInput(0),msuInput(0),msuInput(0)) - 3*AuInput(1,0)*AuInput(1,2)*Sqr(
      MuInput)*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(0),msuInput(1),
      msuInput(1)) - 3*AuInput(2,0)*AuInput(2,2)*Sqr(MuInput)*Sqr(Yu(2,0))*Sqr(Yu(
      2,2))*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(2)) - 3*AuInput(0,1)
      *AuInput(0,2)*Sqr(MuInput)*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(2),
      msqInput(1),msuInput(0),msuInput(0)) - 3*AuInput(1,1)*AuInput(1,2)*Sqr(
      MuInput)*Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(1),msuInput(1),
      msuInput(1)) - 3*AuInput(2,1)*AuInput(2,2)*Sqr(MuInput)*Sqr(Yu(2,1))*Sqr(Yu(
      2,2))*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(2)) - 3*Quad(Yu(0,2)
      )*Sqr(MuInput)*Sqr(AuInput(0,2))*TCD0(msqInput(2),msqInput(2),msuInput(0),
      msuInput(0)) - 3*AuInput(0,2)*AuInput(1,2)*Sqr(MuInput)*Sqr(Yu(0,2))*Sqr(Yu(
      1,2))*TCD0(msqInput(2),msqInput(2),msuInput(0),msuInput(1)) - 3*AuInput(0,2)
      *AuInput(2,2)*Sqr(MuInput)*Sqr(Yu(0,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(0),msuInput(2)) - 3*AuInput(0,2)*AuInput(1,2)*Sqr(
      MuInput)*Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),msuInput(1),
      msuInput(0)) - 3*Quad(Yu(1,2))*Sqr(MuInput)*Sqr(AuInput(1,2))*TCD0(msqInput(
      2),msqInput(2),msuInput(1),msuInput(1)) - 3*AuInput(1,2)*AuInput(2,2)*Sqr(
      MuInput)*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(1),
      msuInput(2)) - 3*AuInput(0,2)*AuInput(2,2)*Sqr(MuInput)*Sqr(Yu(0,2))*Sqr(Yu(
      2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(0)) - 3*AuInput(1,2)
      *AuInput(2,2)*Sqr(MuInput)*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msuInput(2),msuInput(1)) - 3*Quad(Yu(2,2))*Sqr(MuInput)*Sqr(
      AuInput(2,2))*TCD0(msqInput(2),msqInput(2),msuInput(2),msuInput(2)) - 3*
      AdInput(0,1)*AdInput(1,0)*Sqr(MuInput)*TCD0(msdInput(0),msdInput(1),msqInput
      (0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*AdInput(0,0)*AdInput(1,
      1)*Sqr(MuInput)*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(0))*Yd(0,0
      )*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*AdInput(0,0)*AdInput(1,1)*Sqr(MuInput)*TCD0(
      msdInput(1),msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(
      1,1) - 3*AdInput(0,1)*AdInput(1,0)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(0)
      ,msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) - 3*AdInput(0,2)*
      AdInput(1,0)*Sqr(MuInput)*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(
      2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*AdInput(0,0)*AdInput(1,2)*Sqr(
      MuInput)*TCD0(msdInput(0),msdInput(1),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,
      2)*Yd(1,0)*Yd(1,2) - 3*AdInput(0,0)*AdInput(1,2)*Sqr(MuInput)*TCD0(msdInput(
      1),msdInput(0),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*
      AdInput(0,2)*AdInput(1,0)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(0),msqInput
      (2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) - 3*AdInput(0,2)*AdInput(1,
      1)*Sqr(MuInput)*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(2))*Yd(0,1
      )*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*AdInput(0,1)*AdInput(1,2)*Sqr(MuInput)*TCD0(
      msdInput(0),msdInput(1),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(
      1,2) - 3*AdInput(0,1)*AdInput(1,2)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(0)
      ,msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*AdInput(0,2)*
      AdInput(1,1)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(
      1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) - 3*AdInput(0,1)*AdInput(2,0)*Sqr(
      MuInput)*TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,
      1)*Yd(2,0)*Yd(2,1) - 3*AdInput(0,0)*AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(
      0),msdInput(2),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*
      AdInput(0,0)*AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(2),msdInput(0),msqInput
      (0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*AdInput(0,1)*AdInput(2,
      0)*Sqr(MuInput)*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(0))*Yd(0,0
      )*Yd(0,1)*Yd(2,0)*Yd(2,1) - 3*AdInput(1,1)*AdInput(2,0)*Sqr(MuInput)*TCD0(
      msdInput(1),msdInput(2),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(
      2,1) - 3*AdInput(1,0)*AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(2)
      ,msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*AdInput(1,0)*
      AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(
      1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) - 3*AdInput(1,1)*AdInput(2,0)*Sqr(
      MuInput)*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,
      1)*Yd(2,0)*Yd(2,1) - 3*AdInput(0,2)*AdInput(2,0)*Sqr(MuInput)*TCD0(msdInput(
      0),msdInput(2),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(0,0)*AdInput(2,2)*Sqr(MuInput)*TCD0(msdInput(0),msdInput(2),msqInput
      (2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*AdInput(0,0)*AdInput(2,
      2)*Sqr(MuInput)*TCD0(msdInput(2),msdInput(0),msqInput(0),msqInput(2))*Yd(0,0
      )*Yd(0,2)*Yd(2,0)*Yd(2,2) - 3*AdInput(0,2)*AdInput(2,0)*Sqr(MuInput)*TCD0(
      msdInput(2),msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(
      2,2) - 3*AdInput(1,2)*AdInput(2,0)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(2)
      ,msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*AdInput(1,0)*
      AdInput(2,2)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput(
      0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*AdInput(1,0)*AdInput(2,2)*Sqr(
      MuInput)*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,
      2)*Yd(2,0)*Yd(2,2) - 3*AdInput(1,2)*AdInput(2,0)*Sqr(MuInput)*TCD0(msdInput(
      2),msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) - 3*
      AdInput(0,2)*AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(0),msdInput(2),msqInput
      (1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*AdInput(0,1)*AdInput(2,
      2)*Sqr(MuInput)*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput(1))*Yd(0,1
      )*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*AdInput(0,1)*AdInput(2,2)*Sqr(MuInput)*TCD0(
      msdInput(2),msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(
      2,2) - 3*AdInput(0,2)*AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(2),msdInput(0)
      ,msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) - 3*AdInput(1,2)*
      AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(
      2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*AdInput(1,1)*AdInput(2,2)*Sqr(
      MuInput)*TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,
      2)*Yd(2,1)*Yd(2,2) - 3*AdInput(1,1)*AdInput(2,2)*Sqr(MuInput)*TCD0(msdInput(
      2),msdInput(1),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - 3*
      AdInput(1,2)*AdInput(2,1)*Sqr(MuInput)*TCD0(msdInput(2),msdInput(1),msqInput
      (2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) - AeInput(0,1)*AeInput(1,0)
      *Sqr(MuInput)*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,0)*
      Ye(0,1)*Ye(1,0)*Ye(1,1) - AeInput(0,0)*AeInput(1,1)*Sqr(MuInput)*TCD0(
      mseInput(0),mseInput(1),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(
      1,1) - AeInput(0,0)*AeInput(1,1)*Sqr(MuInput)*TCD0(mseInput(1),mseInput(0),
      mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) - AeInput(0,1)*
      AeInput(1,0)*Sqr(MuInput)*TCD0(mseInput(1),mseInput(0),mslInput(1),mslInput(
      0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) - AeInput(0,2)*AeInput(1,0)*Sqr(MuInput)
      *TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,
      0)*Ye(1,2) - AeInput(0,0)*AeInput(1,2)*Sqr(MuInput)*TCD0(mseInput(0),
      mseInput(1),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) -
      AeInput(0,0)*AeInput(1,2)*Sqr(MuInput)*TCD0(mseInput(1),mseInput(0),mslInput
      (0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) - AeInput(0,2)*AeInput(1,0)
      *Sqr(MuInput)*TCD0(mseInput(1),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*
      Ye(0,2)*Ye(1,0)*Ye(1,2) - AeInput(0,2)*AeInput(1,1)*Sqr(MuInput)*TCD0(
      mseInput(0),mseInput(1),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(
      1,2) - AeInput(0,1)*AeInput(1,2)*Sqr(MuInput)*TCD0(mseInput(0),mseInput(1),
      mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) - AeInput(0,1)*
      AeInput(1,2)*Sqr(MuInput)*TCD0(mseInput(1),mseInput(0),mslInput(1),mslInput(
      2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) - AeInput(0,2)*AeInput(1,1)*Sqr(MuInput)
      *TCD0(mseInput(1),mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,
      1)*Ye(1,2) - AeInput(0,1)*AeInput(2,0)*Sqr(MuInput)*TCD0(mseInput(0),
      mseInput(2),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) -
      AeInput(0,0)*AeInput(2,1)*Sqr(MuInput)*TCD0(mseInput(0),mseInput(2),mslInput
      (1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) - AeInput(0,0)*AeInput(2,1)
      *Sqr(MuInput)*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*
      Ye(0,1)*Ye(2,0)*Ye(2,1) - AeInput(0,1)*AeInput(2,0)*Sqr(MuInput)*TCD0(
      mseInput(2),mseInput(0),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(
      2,1) - AeInput(1,1)*AeInput(2,0)*Sqr(MuInput)*TCD0(mseInput(1),mseInput(2),
      mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) - AeInput(1,0)*
      AeInput(2,1)*Sqr(MuInput)*TCD0(mseInput(1),mseInput(2),mslInput(1),mslInput(
      0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) - AeInput(1,0)*AeInput(2,1)*Sqr(MuInput)
      *TCD0(mseInput(2),mseInput(1),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,
      0)*Ye(2,1) - AeInput(1,1)*AeInput(2,0)*Sqr(MuInput)*TCD0(mseInput(2),
      mseInput(1),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) -
      AeInput(0,2)*AeInput(2,0)*Sqr(MuInput)*TCD0(mseInput(0),mseInput(2),mslInput
      (0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) - AeInput(0,0)*AeInput(2,2)
      *Sqr(MuInput)*TCD0(mseInput(0),mseInput(2),mslInput(2),mslInput(0))*Ye(0,0)*
      Ye(0,2)*Ye(2,0)*Ye(2,2) - AeInput(0,0)*AeInput(2,2)*Sqr(MuInput)*TCD0(
      mseInput(2),mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(
      2,2) - AeInput(0,2)*AeInput(2,0)*Sqr(MuInput)*TCD0(mseInput(2),mseInput(0),
      mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) - AeInput(1,2)*
      AeInput(2,0)*Sqr(MuInput)*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(
      2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) - AeInput(1,0)*AeInput(2,2)*Sqr(MuInput)
      *TCD0(mseInput(1),mseInput(2),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,
      0)*Ye(2,2) - AeInput(1,0)*AeInput(2,2)*Sqr(MuInput)*TCD0(mseInput(2),
      mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) -
      AeInput(1,2)*AeInput(2,0)*Sqr(MuInput)*TCD0(mseInput(2),mseInput(1),mslInput
      (2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) - AeInput(0,2)*AeInput(2,1)
      *Sqr(MuInput)*TCD0(mseInput(0),mseInput(2),mslInput(1),mslInput(2))*Ye(0,1)*
      Ye(0,2)*Ye(2,1)*Ye(2,2) - AeInput(0,1)*AeInput(2,2)*Sqr(MuInput)*TCD0(
      mseInput(0),mseInput(2),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(
      2,2) - AeInput(0,1)*AeInput(2,2)*Sqr(MuInput)*TCD0(mseInput(2),mseInput(0),
      mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) - AeInput(0,2)*
      AeInput(2,1)*Sqr(MuInput)*TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(
      1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) - AeInput(1,2)*AeInput(2,1)*Sqr(MuInput)
      *TCD0(mseInput(1),mseInput(2),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,
      1)*Ye(2,2) - AeInput(1,1)*AeInput(2,2)*Sqr(MuInput)*TCD0(mseInput(1),
      mseInput(2),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) -
      AeInput(1,1)*AeInput(2,2)*Sqr(MuInput)*TCD0(mseInput(2),mseInput(1),mslInput
      (1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) - AeInput(1,2)*AeInput(2,1)
      *Sqr(MuInput)*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(1))*Ye(1,1)*
      Ye(1,2)*Ye(2,1)*Ye(2,2) - 3*AuInput(0,0)*AuInput(1,1)*Sqr(MuInput)*TCD0(
      msqInput(0),msqInput(1),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(
      1,1) - 3*AuInput(0,1)*AuInput(1,0)*Sqr(MuInput)*TCD0(msqInput(0),msqInput(1)
      ,msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*AuInput(0,1)*
      AuInput(1,0)*Sqr(MuInput)*TCD0(msqInput(1),msqInput(0),msuInput(0),msuInput(
      1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) - 3*AuInput(0,0)*AuInput(1,1)*Sqr(
      MuInput)*TCD0(msqInput(1),msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,
      1)*Yu(1,0)*Yu(1,1) - 3*AuInput(0,0)*AuInput(1,2)*Sqr(MuInput)*TCD0(msqInput(
      0),msqInput(2),msuInput(0),msuInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*
      AuInput(0,2)*AuInput(1,0)*Sqr(MuInput)*TCD0(msqInput(0),msqInput(2),msuInput
      (1),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*AuInput(0,2)*AuInput(1,
      0)*Sqr(MuInput)*TCD0(msqInput(2),msqInput(0),msuInput(0),msuInput(1))*Yu(0,0
      )*Yu(0,2)*Yu(1,0)*Yu(1,2) - 3*AuInput(0,0)*AuInput(1,2)*Sqr(MuInput)*TCD0(
      msqInput(2),msqInput(0),msuInput(1),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(
      1,2) - 3*AuInput(0,1)*AuInput(1,2)*Sqr(MuInput)*TCD0(msqInput(1),msqInput(2)
      ,msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*AuInput(0,2)*
      AuInput(1,1)*Sqr(MuInput)*TCD0(msqInput(1),msqInput(2),msuInput(1),msuInput(
      0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*AuInput(0,2)*AuInput(1,1)*Sqr(
      MuInput)*TCD0(msqInput(2),msqInput(1),msuInput(0),msuInput(1))*Yu(0,1)*Yu(0,
      2)*Yu(1,1)*Yu(1,2) - 3*AuInput(0,1)*AuInput(1,2)*Sqr(MuInput)*TCD0(msqInput(
      2),msqInput(1),msuInput(1),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) - 3*
      AuInput(0,0)*AuInput(2,1)*Sqr(MuInput)*TCD0(msqInput(0),msqInput(1),msuInput
      (0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) - 3*AuInput(0,1)*AuInput(2,
      0)*Sqr(MuInput)*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(0))*Yu(0,0
      )*Yu(0,1)*Yu(2,0)*Yu(2,1) - 3*AuInput(0,1)*AuInput(2,0)*Sqr(MuInput)*TCD0(
      msqInput(1),msqInput(0),msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(
      2,1) - 3*AuInput(0,0)*AuInput(2,1)*Sqr(MuInput)*TCD0(msqInput(1),msqInput(0)
      ,msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) - 3*AuInput(1,0)*
      AuInput(2,1)*Sqr(MuInput)*TCD0(msqInput(0),msqInput(1),msuInput(1),msuInput(
      2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*AuInput(1,1)*AuInput(2,0)*Sqr(
      MuInput)*TCD0(msqInput(0),msqInput(1),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,
      1)*Yu(2,0)*Yu(2,1) - 3*AuInput(1,1)*AuInput(2,0)*Sqr(MuInput)*TCD0(msqInput(
      1),msqInput(0),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*
      AuInput(1,0)*AuInput(2,1)*Sqr(MuInput)*TCD0(msqInput(1),msqInput(0),msuInput
      (2),msuInput(1))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) - 3*AuInput(0,0)*AuInput(2,
      2)*Sqr(MuInput)*TCD0(msqInput(0),msqInput(2),msuInput(0),msuInput(2))*Yu(0,0
      )*Yu(0,2)*Yu(2,0)*Yu(2,2) - 3*AuInput(0,2)*AuInput(2,0)*Sqr(MuInput)*TCD0(
      msqInput(0),msqInput(2),msuInput(2),msuInput(0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(
      2,2) - 3*AuInput(0,2)*AuInput(2,0)*Sqr(MuInput)*TCD0(msqInput(2),msqInput(0)
      ,msuInput(0),msuInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) - 3*AuInput(0,0)*
      AuInput(2,2)*Sqr(MuInput)*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(
      0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) - 3*AuInput(1,0)*AuInput(2,2)*Sqr(
      MuInput)*TCD0(msqInput(0),msqInput(2),msuInput(1),msuInput(2))*Yu(1,0)*Yu(1,
      2)*Yu(2,0)*Yu(2,2) - 3*AuInput(1,2)*AuInput(2,0)*Sqr(MuInput)*TCD0(msqInput(
      0),msqInput(2),msuInput(2),msuInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*
      AuInput(1,2)*AuInput(2,0)*Sqr(MuInput)*TCD0(msqInput(2),msqInput(0),msuInput
      (1),msuInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*AuInput(1,0)*AuInput(2,
      2)*Sqr(MuInput)*TCD0(msqInput(2),msqInput(0),msuInput(2),msuInput(1))*Yu(1,0
      )*Yu(1,2)*Yu(2,0)*Yu(2,2) - 3*AuInput(0,1)*AuInput(2,2)*Sqr(MuInput)*TCD0(
      msqInput(1),msqInput(2),msuInput(0),msuInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(
      2,2) - 3*AuInput(0,2)*AuInput(2,1)*Sqr(MuInput)*TCD0(msqInput(1),msqInput(2)
      ,msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*AuInput(0,2)*
      AuInput(2,1)*Sqr(MuInput)*TCD0(msqInput(2),msqInput(1),msuInput(0),msuInput(
      2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) - 3*AuInput(0,1)*AuInput(2,2)*Sqr(
      MuInput)*TCD0(msqInput(2),msqInput(1),msuInput(2),msuInput(0))*Yu(0,1)*Yu(0,
      2)*Yu(2,1)*Yu(2,2) - 3*AuInput(1,1)*AuInput(2,2)*Sqr(MuInput)*TCD0(msqInput(
      1),msqInput(2),msuInput(1),msuInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*
      AuInput(1,2)*AuInput(2,1)*Sqr(MuInput)*TCD0(msqInput(1),msqInput(2),msuInput
      (2),msuInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*AuInput(1,2)*AuInput(2,
      1)*Sqr(MuInput)*TCD0(msqInput(2),msqInput(1),msuInput(1),msuInput(2))*Yu(1,1
      )*Yu(1,2)*Yu(2,1)*Yu(2,2) - 3*AuInput(1,1)*AuInput(2,2)*Sqr(MuInput)*TCD0(
      msqInput(2),msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(
      2,2))));
   MODEL->set_Lambda6(Re((0.000053468657576480914*(-1 + AuInput(2,2)/MSUSY)*Cube(
      MuInput)*Quad(Yu(2,2))*Sqr(g3)*UnitStep(-2 + LambdaLoopOrder))/Cube(MSUSY) +
      0.006332573977646111*UnitStep(-1 + LambdaLoopOrder)*((-0.3*MuInput*AdInput(0
      ,0)*Sqr(g1)*Sqr(Yd(0,0)) + 3*MuInput*AdInput(0,0)*Sqr(Yd(0,0))*(Sqr(Yd(0,0))
      + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))*TCC0(msdInput(0),msdInput(0),msqInput(0)) +
      (-0.3*MuInput*AdInput(0,1)*Sqr(g1)*Sqr(Yd(0,1)) + 3*MuInput*AdInput(0,1)*Sqr
      (Yd(0,1))*(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,2))))*TCC0(msdInput(0),
      msdInput(0),msqInput(1)) + (-0.3*MuInput*AdInput(0,2)*Sqr(g1)*Sqr(Yd(0,2)) +
      3*MuInput*AdInput(0,2)*Sqr(Yd(0,2))*(Sqr(Yd(0,0)) + Sqr(Yd(0,1)) + Sqr(Yd(0,
      2))))*TCC0(msdInput(0),msdInput(0),msqInput(2)) + (0.25*MuInput*AdInput(0,0)
      *(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yd(0,0)) + 3*MuInput*AdInput(0,0)*Sqr(Yd(0,0
      ))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0))))*TCC0(msdInput(0),msqInput(0
      ),msqInput(0)) + (0.25*MuInput*AdInput(0,1)*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(
      Yd(0,1)) + 3*MuInput*AdInput(0,1)*Sqr(Yd(0,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1))
      + Sqr(Yd(2,1))))*TCC0(msdInput(0),msqInput(1),msqInput(1)) + (0.25*MuInput*
      AdInput(0,2)*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yd(0,2)) + 3*MuInput*AdInput(0,2
      )*Sqr(Yd(0,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))))*TCC0(msdInput(0
      ),msqInput(2),msqInput(2)) + (-0.3*MuInput*AdInput(1,0)*Sqr(g1)*Sqr(Yd(1,0))
      + 3*MuInput*AdInput(1,0)*Sqr(Yd(1,0))*(Sqr(Yd(1,0)) + Sqr(Yd(1,1)) + Sqr(Yd(
      1,2))))*TCC0(msdInput(1),msdInput(1),msqInput(0)) + (-0.3*MuInput*AdInput(1,
      1)*Sqr(g1)*Sqr(Yd(1,1)) + 3*MuInput*AdInput(1,1)*Sqr(Yd(1,1))*(Sqr(Yd(1,0))
      + Sqr(Yd(1,1)) + Sqr(Yd(1,2))))*TCC0(msdInput(1),msdInput(1),msqInput(1)) +
      (-0.3*MuInput*AdInput(1,2)*Sqr(g1)*Sqr(Yd(1,2)) + 3*MuInput*AdInput(1,2)*Sqr
      (Yd(1,2))*(Sqr(Yd(1,0)) + Sqr(Yd(1,1)) + Sqr(Yd(1,2))))*TCC0(msdInput(1),
      msdInput(1),msqInput(2)) + (0.25*MuInput*AdInput(1,0)*(-0.6*Sqr(g1) - 3*Sqr(
      g2))*Sqr(Yd(1,0)) + 3*MuInput*AdInput(1,0)*Sqr(Yd(1,0))*(Sqr(Yd(0,0)) + Sqr(
      Yd(1,0)) + Sqr(Yd(2,0))))*TCC0(msdInput(1),msqInput(0),msqInput(0)) + (0.25*
      MuInput*AdInput(1,1)*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yd(1,1)) + 3*MuInput*
      AdInput(1,1)*Sqr(Yd(1,1))*(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*TCC0
      (msdInput(1),msqInput(1),msqInput(1)) + (0.25*MuInput*AdInput(1,2)*(-0.6*Sqr
      (g1) - 3*Sqr(g2))*Sqr(Yd(1,2)) + 3*MuInput*AdInput(1,2)*Sqr(Yd(1,2))*(Sqr(Yd
      (0,2)) + Sqr(Yd(1,2)) + Sqr(Yd(2,2))))*TCC0(msdInput(1),msqInput(2),msqInput
      (2)) + (-0.3*MuInput*AdInput(2,0)*Sqr(g1)*Sqr(Yd(2,0)) + 3*MuInput*AdInput(2
      ,0)*Sqr(Yd(2,0))*(Sqr(Yd(2,0)) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))))*TCC0(msdInput
      (2),msdInput(2),msqInput(0)) + (-0.3*MuInput*AdInput(2,1)*Sqr(g1)*Sqr(Yd(2,1
      )) + 3*MuInput*AdInput(2,1)*Sqr(Yd(2,1))*(Sqr(Yd(2,0)) + Sqr(Yd(2,1)) + Sqr(
      Yd(2,2))))*TCC0(msdInput(2),msdInput(2),msqInput(1)) + (-0.3*MuInput*AdInput
      (2,2)*Sqr(g1)*Sqr(Yd(2,2)) + 3*MuInput*AdInput(2,2)*Sqr(Yd(2,2))*(Sqr(Yd(2,0
      )) + Sqr(Yd(2,1)) + Sqr(Yd(2,2))))*TCC0(msdInput(2),msdInput(2),msqInput(2))
      + (0.25*MuInput*AdInput(2,0)*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yd(2,0)) + 3*
      MuInput*AdInput(2,0)*Sqr(Yd(2,0))*(Sqr(Yd(0,0)) + Sqr(Yd(1,0)) + Sqr(Yd(2,0)
      )))*TCC0(msdInput(2),msqInput(0),msqInput(0)) + (0.25*MuInput*AdInput(2,1)*(
      -0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yd(2,1)) + 3*MuInput*AdInput(2,1)*Sqr(Yd(2,1))
      *(Sqr(Yd(0,1)) + Sqr(Yd(1,1)) + Sqr(Yd(2,1))))*TCC0(msdInput(2),msqInput(1),
      msqInput(1)) + (0.25*MuInput*AdInput(2,2)*(-0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yd(
      2,2)) + 3*MuInput*AdInput(2,2)*Sqr(Yd(2,2))*(Sqr(Yd(0,2)) + Sqr(Yd(1,2)) +
      Sqr(Yd(2,2))))*TCC0(msdInput(2),msqInput(2),msqInput(2)) + (-0.3*MuInput*
      AeInput(0,0)*Sqr(g1)*Sqr(Ye(0,0)) + MuInput*AeInput(0,0)*Sqr(Ye(0,0))*(Sqr(
      Ye(0,0)) + Sqr(Ye(0,1)) + Sqr(Ye(0,2))))*TCC0(mseInput(0),mseInput(0),
      mslInput(0)) + (-0.3*MuInput*AeInput(0,1)*Sqr(g1)*Sqr(Ye(0,1)) + MuInput*
      AeInput(0,1)*Sqr(Ye(0,1))*(Sqr(Ye(0,0)) + Sqr(Ye(0,1)) + Sqr(Ye(0,2))))*TCC0
      (mseInput(0),mseInput(0),mslInput(1)) + (-0.3*MuInput*AeInput(0,2)*Sqr(g1)*
      Sqr(Ye(0,2)) + MuInput*AeInput(0,2)*Sqr(Ye(0,2))*(Sqr(Ye(0,0)) + Sqr(Ye(0,1)
      ) + Sqr(Ye(0,2))))*TCC0(mseInput(0),mseInput(0),mslInput(2)) + (0.25*MuInput
      *AeInput(0,0)*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Ye(0,0)) + MuInput*AeInput(0,0)*
      Sqr(Ye(0,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(0),
      mslInput(0),mslInput(0)) + (0.25*MuInput*AeInput(0,1)*(0.6*Sqr(g1) - Sqr(g2)
      )*Sqr(Ye(0,1)) + MuInput*AeInput(0,1)*Sqr(Ye(0,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,
      1)) + Sqr(Ye(2,1))))*TCC0(mseInput(0),mslInput(1),mslInput(1)) + (0.25*
      MuInput*AeInput(0,2)*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Ye(0,2)) + MuInput*AeInput(
      0,2)*Sqr(Ye(0,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(
      mseInput(0),mslInput(2),mslInput(2)) + (-0.3*MuInput*AeInput(1,0)*Sqr(g1)*
      Sqr(Ye(1,0)) + MuInput*AeInput(1,0)*Sqr(Ye(1,0))*(Sqr(Ye(1,0)) + Sqr(Ye(1,1)
      ) + Sqr(Ye(1,2))))*TCC0(mseInput(1),mseInput(1),mslInput(0)) + (-0.3*MuInput
      *AeInput(1,1)*Sqr(g1)*Sqr(Ye(1,1)) + MuInput*AeInput(1,1)*Sqr(Ye(1,1))*(Sqr(
      Ye(1,0)) + Sqr(Ye(1,1)) + Sqr(Ye(1,2))))*TCC0(mseInput(1),mseInput(1),
      mslInput(1)) + (-0.3*MuInput*AeInput(1,2)*Sqr(g1)*Sqr(Ye(1,2)) + MuInput*
      AeInput(1,2)*Sqr(Ye(1,2))*(Sqr(Ye(1,0)) + Sqr(Ye(1,1)) + Sqr(Ye(1,2))))*TCC0
      (mseInput(1),mseInput(1),mslInput(2)) + (0.25*MuInput*AeInput(1,0)*(0.6*Sqr(
      g1) - Sqr(g2))*Sqr(Ye(1,0)) + MuInput*AeInput(1,0)*Sqr(Ye(1,0))*(Sqr(Ye(0,0)
      ) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(mseInput(1),mslInput(0),mslInput(0))
      + (0.25*MuInput*AeInput(1,1)*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Ye(1,1)) + MuInput*
      AeInput(1,1)*Sqr(Ye(1,1))*(Sqr(Ye(0,1)) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0
      (mseInput(1),mslInput(1),mslInput(1)) + (0.25*MuInput*AeInput(1,2)*(0.6*Sqr(
      g1) - Sqr(g2))*Sqr(Ye(1,2)) + MuInput*AeInput(1,2)*Sqr(Ye(1,2))*(Sqr(Ye(0,2)
      ) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0(mseInput(1),mslInput(2),mslInput(2))
      + (-0.3*MuInput*AeInput(2,0)*Sqr(g1)*Sqr(Ye(2,0)) + MuInput*AeInput(2,0)*Sqr
      (Ye(2,0))*(Sqr(Ye(2,0)) + Sqr(Ye(2,1)) + Sqr(Ye(2,2))))*TCC0(mseInput(2),
      mseInput(2),mslInput(0)) + (-0.3*MuInput*AeInput(2,1)*Sqr(g1)*Sqr(Ye(2,1)) +
      MuInput*AeInput(2,1)*Sqr(Ye(2,1))*(Sqr(Ye(2,0)) + Sqr(Ye(2,1)) + Sqr(Ye(2,2)
      )))*TCC0(mseInput(2),mseInput(2),mslInput(1)) + (-0.3*MuInput*AeInput(2,2)*
      Sqr(g1)*Sqr(Ye(2,2)) + MuInput*AeInput(2,2)*Sqr(Ye(2,2))*(Sqr(Ye(2,0)) + Sqr
      (Ye(2,1)) + Sqr(Ye(2,2))))*TCC0(mseInput(2),mseInput(2),mslInput(2)) + (0.25
      *MuInput*AeInput(2,0)*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Ye(2,0)) + MuInput*AeInput
      (2,0)*Sqr(Ye(2,0))*(Sqr(Ye(0,0)) + Sqr(Ye(1,0)) + Sqr(Ye(2,0))))*TCC0(
      mseInput(2),mslInput(0),mslInput(0)) + (0.25*MuInput*AeInput(2,1)*(0.6*Sqr(
      g1) - Sqr(g2))*Sqr(Ye(2,1)) + MuInput*AeInput(2,1)*Sqr(Ye(2,1))*(Sqr(Ye(0,1)
      ) + Sqr(Ye(1,1)) + Sqr(Ye(2,1))))*TCC0(mseInput(2),mslInput(1),mslInput(1))
      + (0.25*MuInput*AeInput(2,2)*(0.6*Sqr(g1) - Sqr(g2))*Sqr(Ye(2,2)) + MuInput*
      AeInput(2,2)*Sqr(Ye(2,2))*(Sqr(Ye(0,2)) + Sqr(Ye(1,2)) + Sqr(Ye(2,2))))*TCC0
      (mseInput(2),mslInput(2),mslInput(2)) + 0.25*MuInput*AuInput(0,0)*(-0.6*Sqr(
      g1) + 3*Sqr(g2))*Sqr(Yu(0,0))*TCC0(msqInput(0),msqInput(0),msuInput(0)) +
      0.25*MuInput*AuInput(1,0)*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yu(1,0))*TCC0(
      msqInput(0),msqInput(0),msuInput(1)) + 0.25*MuInput*AuInput(2,0)*(-0.6*Sqr(
      g1) + 3*Sqr(g2))*Sqr(Yu(2,0))*TCC0(msqInput(0),msqInput(0),msuInput(2)) +
      0.6*MuInput*AuInput(0,0)*Sqr(g1)*Sqr(Yu(0,0))*TCC0(msqInput(0),msuInput(0),
      msuInput(0)) + 0.6*MuInput*AuInput(1,0)*Sqr(g1)*Sqr(Yu(1,0))*TCC0(msqInput(0
      ),msuInput(1),msuInput(1)) + 0.6*MuInput*AuInput(2,0)*Sqr(g1)*Sqr(Yu(2,0))*
      TCC0(msqInput(0),msuInput(2),msuInput(2)) + 0.25*MuInput*AuInput(0,1)*(-0.6*
      Sqr(g1) + 3*Sqr(g2))*Sqr(Yu(0,1))*TCC0(msqInput(1),msqInput(1),msuInput(0))
      + 0.25*MuInput*AuInput(1,1)*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yu(1,1))*TCC0(
      msqInput(1),msqInput(1),msuInput(1)) + 0.25*MuInput*AuInput(2,1)*(-0.6*Sqr(
      g1) + 3*Sqr(g2))*Sqr(Yu(2,1))*TCC0(msqInput(1),msqInput(1),msuInput(2)) +
      0.6*MuInput*AuInput(0,1)*Sqr(g1)*Sqr(Yu(0,1))*TCC0(msqInput(1),msuInput(0),
      msuInput(0)) + 0.6*MuInput*AuInput(1,1)*Sqr(g1)*Sqr(Yu(1,1))*TCC0(msqInput(1
      ),msuInput(1),msuInput(1)) + 0.6*MuInput*AuInput(2,1)*Sqr(g1)*Sqr(Yu(2,1))*
      TCC0(msqInput(1),msuInput(2),msuInput(2)) + 0.25*MuInput*AuInput(0,2)*(-0.6*
      Sqr(g1) + 3*Sqr(g2))*Sqr(Yu(0,2))*TCC0(msqInput(2),msqInput(2),msuInput(0))
      + 0.25*MuInput*AuInput(1,2)*(-0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yu(1,2))*TCC0(
      msqInput(2),msqInput(2),msuInput(1)) + 0.25*MuInput*AuInput(2,2)*(-0.6*Sqr(
      g1) + 3*Sqr(g2))*Sqr(Yu(2,2))*TCC0(msqInput(2),msqInput(2),msuInput(2)) +
      0.6*MuInput*AuInput(0,2)*Sqr(g1)*Sqr(Yu(0,2))*TCC0(msqInput(2),msuInput(0),
      msuInput(0)) + 0.6*MuInput*AuInput(1,2)*Sqr(g1)*Sqr(Yu(1,2))*TCC0(msqInput(2
      ),msuInput(1),msuInput(1)) + 0.6*MuInput*AuInput(2,2)*Sqr(g1)*Sqr(Yu(2,2))*
      TCC0(msqInput(2),msuInput(2),msuInput(2)) + 3*MuInput*Cube(AdInput(0,0))*
      Quad(Yd(0,0))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(0)) + 3*
      MuInput*AdInput(0,0)*Sqr(AdInput(0,1))*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(
      msdInput(0),msdInput(0),msqInput(0),msqInput(1)) + 3*MuInput*AdInput(0,0)*
      Sqr(AdInput(0,2))*Sqr(Yd(0,0))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),
      msqInput(0),msqInput(2)) + 3*MuInput*AdInput(0,1)*Sqr(AdInput(0,0))*Sqr(Yd(0
      ,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1),msqInput(0)) + 3*
      MuInput*Cube(AdInput(0,1))*Quad(Yd(0,1))*TCD0(msdInput(0),msdInput(0),
      msqInput(1),msqInput(1)) + 3*MuInput*AdInput(0,1)*Sqr(AdInput(0,2))*Sqr(Yd(0
      ,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(1),msqInput(2)) + 3*
      MuInput*AdInput(0,2)*Sqr(AdInput(0,0))*Sqr(Yd(0,0))*Sqr(Yd(0,2))*TCD0(
      msdInput(0),msdInput(0),msqInput(2),msqInput(0)) + 3*MuInput*AdInput(0,2)*
      Sqr(AdInput(0,1))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),
      msqInput(2),msqInput(1)) + 3*MuInput*Cube(AdInput(0,2))*Quad(Yd(0,2))*TCD0(
      msdInput(0),msdInput(0),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(1,0)*
      Sqr(AdInput(0,0))*Sqr(Yd(0,0))*Sqr(Yd(1,0))*TCD0(msdInput(0),msdInput(1),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(1,1)*Sqr(AdInput(0,1))*Sqr(Yd(0
      ,1))*Sqr(Yd(1,1))*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(1,2)*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*Sqr(Yd(1,2))*TCD0(
      msdInput(0),msdInput(1),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(2,0)*
      Sqr(AdInput(0,0))*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(0),msdInput(2),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(2,1)*Sqr(AdInput(0,1))*Sqr(Yd(0
      ,1))*Sqr(Yd(2,1))*TCD0(msdInput(0),msdInput(2),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(2,2)*Sqr(AdInput(0,2))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(0),msdInput(2),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(0,0)*
      Sqr(AdInput(1,0))*Sqr(Yd(0,0))*Sqr(Yd(1,0))*TCD0(msdInput(1),msdInput(0),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(0,1)*Sqr(AdInput(1,1))*Sqr(Yd(0
      ,1))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(0),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(0,2)*Sqr(AdInput(1,2))*Sqr(Yd(0,2))*Sqr(Yd(1,2))*TCD0(
      msdInput(1),msdInput(0),msqInput(2),msqInput(2)) + 3*MuInput*Cube(AdInput(1,
      0))*Quad(Yd(1,0))*TCD0(msdInput(1),msdInput(1),msqInput(0),msqInput(0)) + 3*
      MuInput*AdInput(1,0)*Sqr(AdInput(1,1))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(
      msdInput(1),msdInput(1),msqInput(0),msqInput(1)) + 3*MuInput*AdInput(1,0)*
      Sqr(AdInput(1,2))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),
      msqInput(0),msqInput(2)) + 3*MuInput*AdInput(1,1)*Sqr(AdInput(1,0))*Sqr(Yd(1
      ,0))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(0)) + 3*
      MuInput*Cube(AdInput(1,1))*Quad(Yd(1,1))*TCD0(msdInput(1),msdInput(1),
      msqInput(1),msqInput(1)) + 3*MuInput*AdInput(1,1)*Sqr(AdInput(1,2))*Sqr(Yd(1
      ,1))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(2)) + 3*
      MuInput*AdInput(1,2)*Sqr(AdInput(1,0))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(
      msdInput(1),msdInput(1),msqInput(2),msqInput(0)) + 3*MuInput*AdInput(1,2)*
      Sqr(AdInput(1,1))*Sqr(Yd(1,1))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),
      msqInput(2),msqInput(1)) + 3*MuInput*Cube(AdInput(1,2))*Quad(Yd(1,2))*TCD0(
      msdInput(1),msdInput(1),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(2,0)*
      Sqr(AdInput(1,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(1),msdInput(2),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(2,1)*Sqr(AdInput(1,1))*Sqr(Yd(1
      ,1))*Sqr(Yd(2,1))*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(2,2)*Sqr(AdInput(1,2))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(1),msdInput(2),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(0,0)*
      Sqr(AdInput(2,0))*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(2),msdInput(0),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(0,1)*Sqr(AdInput(2,1))*Sqr(Yd(0
      ,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(0,2)*Sqr(AdInput(2,2))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(0),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(1,0)*
      Sqr(AdInput(2,0))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(2),msdInput(1),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(1,1)*Sqr(AdInput(2,1))*Sqr(Yd(1
      ,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(1,2)*Sqr(AdInput(2,2))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(1),msqInput(2),msqInput(2)) + 3*MuInput*Cube(AdInput(2,
      0))*Quad(Yd(2,0))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(0)) + 3*
      MuInput*AdInput(2,0)*Sqr(AdInput(2,1))*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(
      msdInput(2),msdInput(2),msqInput(0),msqInput(1)) + 3*MuInput*AdInput(2,0)*
      Sqr(AdInput(2,2))*Sqr(Yd(2,0))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(0),msqInput(2)) + 3*MuInput*AdInput(2,1)*Sqr(AdInput(2,0))*Sqr(Yd(2
      ,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(0)) + 3*
      MuInput*Cube(AdInput(2,1))*Quad(Yd(2,1))*TCD0(msdInput(2),msdInput(2),
      msqInput(1),msqInput(1)) + 3*MuInput*AdInput(2,1)*Sqr(AdInput(2,2))*Sqr(Yd(2
      ,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(1),msqInput(2)) + 3*
      MuInput*AdInput(2,2)*Sqr(AdInput(2,0))*Sqr(Yd(2,0))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(2),msqInput(2),msqInput(0)) + 3*MuInput*AdInput(2,2)*
      Sqr(AdInput(2,1))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(2),msqInput(1)) + 3*MuInput*Cube(AdInput(2,2))*Quad(Yd(2,2))*TCD0(
      msdInput(2),msdInput(2),msqInput(2),msqInput(2)) + MuInput*Cube(AeInput(0,0)
      )*Quad(Ye(0,0))*TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(0)) +
      MuInput*AeInput(0,0)*Sqr(AeInput(0,1))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(
      mseInput(0),mseInput(0),mslInput(0),mslInput(1)) + MuInput*AeInput(0,0)*Sqr(
      AeInput(0,2))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),
      mslInput(0),mslInput(2)) + MuInput*AeInput(0,1)*Sqr(AeInput(0,0))*Sqr(Ye(0,0
      ))*Sqr(Ye(0,1))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(0)) +
      MuInput*Cube(AeInput(0,1))*Quad(Ye(0,1))*TCD0(mseInput(0),mseInput(0),
      mslInput(1),mslInput(1)) + MuInput*AeInput(0,1)*Sqr(AeInput(0,2))*Sqr(Ye(0,1
      ))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(2)) +
      MuInput*AeInput(0,2)*Sqr(AeInput(0,0))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(
      mseInput(0),mseInput(0),mslInput(2),mslInput(0)) + MuInput*AeInput(0,2)*Sqr(
      AeInput(0,1))*Sqr(Ye(0,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),
      mslInput(2),mslInput(1)) + MuInput*Cube(AeInput(0,2))*Quad(Ye(0,2))*TCD0(
      mseInput(0),mseInput(0),mslInput(2),mslInput(2)) + MuInput*AeInput(1,0)*Sqr(
      AeInput(0,0))*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(0),mseInput(1),
      mslInput(0),mslInput(0)) + MuInput*AeInput(1,1)*Sqr(AeInput(0,1))*Sqr(Ye(0,1
      ))*Sqr(Ye(1,1))*TCD0(mseInput(0),mseInput(1),mslInput(1),mslInput(1)) +
      MuInput*AeInput(1,2)*Sqr(AeInput(0,2))*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(
      mseInput(0),mseInput(1),mslInput(2),mslInput(2)) + MuInput*AeInput(2,0)*Sqr(
      AeInput(0,0))*Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(0),mseInput(2),
      mslInput(0),mslInput(0)) + MuInput*AeInput(2,1)*Sqr(AeInput(0,1))*Sqr(Ye(0,1
      ))*Sqr(Ye(2,1))*TCD0(mseInput(0),mseInput(2),mslInput(1),mslInput(1)) +
      MuInput*AeInput(2,2)*Sqr(AeInput(0,2))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(0),mseInput(2),mslInput(2),mslInput(2)) + MuInput*AeInput(0,0)*Sqr(
      AeInput(1,0))*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(1),mseInput(0),
      mslInput(0),mslInput(0)) + MuInput*AeInput(0,1)*Sqr(AeInput(1,1))*Sqr(Ye(0,1
      ))*Sqr(Ye(1,1))*TCD0(mseInput(1),mseInput(0),mslInput(1),mslInput(1)) +
      MuInput*AeInput(0,2)*Sqr(AeInput(1,2))*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(
      mseInput(1),mseInput(0),mslInput(2),mslInput(2)) + MuInput*Cube(AeInput(1,0)
      )*Quad(Ye(1,0))*TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(0)) +
      MuInput*AeInput(1,0)*Sqr(AeInput(1,1))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(
      mseInput(1),mseInput(1),mslInput(0),mslInput(1)) + MuInput*AeInput(1,0)*Sqr(
      AeInput(1,2))*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(0),mslInput(2)) + MuInput*AeInput(1,1)*Sqr(AeInput(1,0))*Sqr(Ye(1,0
      ))*Sqr(Ye(1,1))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(0)) +
      MuInput*Cube(AeInput(1,1))*Quad(Ye(1,1))*TCD0(mseInput(1),mseInput(1),
      mslInput(1),mslInput(1)) + MuInput*AeInput(1,1)*Sqr(AeInput(1,2))*Sqr(Ye(1,1
      ))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(2)) +
      MuInput*AeInput(1,2)*Sqr(AeInput(1,0))*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(
      mseInput(1),mseInput(1),mslInput(2),mslInput(0)) + MuInput*AeInput(1,2)*Sqr(
      AeInput(1,1))*Sqr(Ye(1,1))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),
      mslInput(2),mslInput(1)) + MuInput*Cube(AeInput(1,2))*Quad(Ye(1,2))*TCD0(
      mseInput(1),mseInput(1),mslInput(2),mslInput(2)) + MuInput*AeInput(2,0)*Sqr(
      AeInput(1,0))*Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(1),mseInput(2),
      mslInput(0),mslInput(0)) + MuInput*AeInput(2,1)*Sqr(AeInput(1,1))*Sqr(Ye(1,1
      ))*Sqr(Ye(2,1))*TCD0(mseInput(1),mseInput(2),mslInput(1),mslInput(1)) +
      MuInput*AeInput(2,2)*Sqr(AeInput(1,2))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(1),mseInput(2),mslInput(2),mslInput(2)) + MuInput*AeInput(0,0)*Sqr(
      AeInput(2,0))*Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(0),
      mslInput(0),mslInput(0)) + MuInput*AeInput(0,1)*Sqr(AeInput(2,1))*Sqr(Ye(0,1
      ))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(0),mslInput(1),mslInput(1)) +
      MuInput*AeInput(0,2)*Sqr(AeInput(2,2))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(2),mseInput(0),mslInput(2),mslInput(2)) + MuInput*AeInput(1,0)*Sqr(
      AeInput(2,0))*Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(1),
      mslInput(0),mslInput(0)) + MuInput*AeInput(1,1)*Sqr(AeInput(2,1))*Sqr(Ye(1,1
      ))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(1),mslInput(1),mslInput(1)) +
      MuInput*AeInput(1,2)*Sqr(AeInput(2,2))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(2),mseInput(1),mslInput(2),mslInput(2)) + MuInput*Cube(AeInput(2,0)
      )*Quad(Ye(2,0))*TCD0(mseInput(2),mseInput(2),mslInput(0),mslInput(0)) +
      MuInput*AeInput(2,0)*Sqr(AeInput(2,1))*Sqr(Ye(2,0))*Sqr(Ye(2,1))*TCD0(
      mseInput(2),mseInput(2),mslInput(0),mslInput(1)) + MuInput*AeInput(2,0)*Sqr(
      AeInput(2,2))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),
      mslInput(0),mslInput(2)) + MuInput*AeInput(2,1)*Sqr(AeInput(2,0))*Sqr(Ye(2,0
      ))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(0)) +
      MuInput*Cube(AeInput(2,1))*Quad(Ye(2,1))*TCD0(mseInput(2),mseInput(2),
      mslInput(1),mslInput(1)) + MuInput*AeInput(2,1)*Sqr(AeInput(2,2))*Sqr(Ye(2,1
      ))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(2)) +
      MuInput*AeInput(2,2)*Sqr(AeInput(2,0))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0(
      mseInput(2),mseInput(2),mslInput(2),mslInput(0)) + MuInput*AeInput(2,2)*Sqr(
      AeInput(2,1))*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),
      mslInput(2),mslInput(1)) + MuInput*Cube(AeInput(2,2))*Quad(Ye(2,2))*TCD0(
      mseInput(2),mseInput(2),mslInput(2),mslInput(2)) + 3*MuInput*AuInput(0,0)*
      Quad(Yu(0,0))*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(0),msqInput(0),
      msqInput(0)) + 3*MuInput*AuInput(1,0)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(
      1,0))*TCD0(msqInput(0),msqInput(0),msqInput(0),msqInput(1)) + 3*MuInput*
      AuInput(2,0)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),
      msqInput(0),msqInput(0),msqInput(2)) + 3*MuInput*AuInput(0,0)*Sqr(Abs(
      MuInput))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),msqInput(1)
      ,msqInput(0)) + 3*MuInput*AuInput(1,0)*Quad(Yu(1,0))*Sqr(Abs(MuInput))*TCD0(
      msqInput(0),msqInput(0),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(2,0)*
      Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),
      msqInput(1),msqInput(2)) + 3*MuInput*AuInput(0,0)*Sqr(Abs(MuInput))*Sqr(Yu(0
      ,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msqInput(2),msqInput(0)) + 3*
      MuInput*AuInput(1,0)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(
      msqInput(0),msqInput(0),msqInput(2),msqInput(1)) + 3*MuInput*AuInput(2,0)*
      Quad(Yu(2,0))*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(0),msqInput(2),
      msqInput(2)) + 3*MuInput*AuInput(0,0)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(
      0,1))*TCD0(msqInput(0),msqInput(1),msqInput(0),msqInput(0)) + 3*MuInput*
      AuInput(1,0)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(msqInput(0),
      msqInput(1),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(2,0)*Sqr(Abs(
      MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(0),msqInput(1),msqInput(2)
      ,msqInput(2)) + 3*MuInput*AuInput(0,0)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu
      (0,2))*TCD0(msqInput(0),msqInput(2),msqInput(0),msqInput(0)) + 3*MuInput*
      AuInput(1,0)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(0),
      msqInput(2),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(2,0)*Sqr(Abs(
      MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(msqInput(0),msqInput(2),msqInput(2)
      ,msqInput(2)) + 3*MuInput*AuInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu
      (0,1))*TCD0(msqInput(1),msqInput(0),msqInput(0),msqInput(0)) + 3*MuInput*
      AuInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(msqInput(1),
      msqInput(0),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(2,1)*Sqr(Abs(
      MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(0),msqInput(2)
      ,msqInput(2)) + 3*MuInput*AuInput(0,1)*Quad(Yu(0,1))*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(1),msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,1)*
      Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(1),
      msqInput(0),msqInput(1)) + 3*MuInput*AuInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yu(0
      ,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msqInput(0),msqInput(2)) + 3*
      MuInput*AuInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(1,1))*TCD0(
      msqInput(1),msqInput(1),msqInput(1),msqInput(0)) + 3*MuInput*AuInput(1,1)*
      Quad(Yu(1,1))*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(1),msqInput(1),
      msqInput(1)) + 3*MuInput*AuInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yu(1,1))*Sqr(Yu(
      2,1))*TCD0(msqInput(1),msqInput(1),msqInput(1),msqInput(2)) + 3*MuInput*
      AuInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),
      msqInput(1),msqInput(2),msqInput(0)) + 3*MuInput*AuInput(1,1)*Sqr(Abs(
      MuInput))*Sqr(Yu(1,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msqInput(2)
      ,msqInput(1)) + 3*MuInput*AuInput(2,1)*Quad(Yu(2,1))*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(1),msqInput(2),msqInput(2)) + 3*MuInput*AuInput(0,1)*
      Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(1),msqInput(2),
      msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yu(1
      ,1))*Sqr(Yu(1,2))*TCD0(msqInput(1),msqInput(2),msqInput(1),msqInput(1)) + 3*
      MuInput*AuInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*Sqr(Yu(2,2))*TCD0(
      msqInput(1),msqInput(2),msqInput(2),msqInput(2)) + 3*MuInput*AuInput(0,2)*
      Sqr(Abs(MuInput))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(2),msqInput(0),
      msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yu(1
      ,0))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(0),msqInput(1),msqInput(1)) + 3*
      MuInput*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,0))*Sqr(Yu(2,2))*TCD0(
      msqInput(2),msqInput(0),msqInput(2),msqInput(2)) + 3*MuInput*AuInput(0,2)*
      Sqr(Abs(MuInput))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(2),msqInput(1),
      msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yu(1
      ,1))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(1),msqInput(1),msqInput(1)) + 3*
      MuInput*AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(2,1))*Sqr(Yu(2,2))*TCD0(
      msqInput(2),msqInput(1),msqInput(2),msqInput(2)) + 3*MuInput*AuInput(0,2)*
      Quad(Yu(0,2))*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(2),msqInput(0),
      msqInput(0)) + 3*MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*Sqr(Yu(
      1,2))*TCD0(msqInput(2),msqInput(2),msqInput(0),msqInput(1)) + 3*MuInput*
      AuInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yu(0,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msqInput(0),msqInput(2)) + 3*MuInput*AuInput(0,2)*Sqr(Abs(
      MuInput))*Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),msqInput(1)
      ,msqInput(0)) + 3*MuInput*AuInput(1,2)*Quad(Yu(1,2))*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(2),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(2,2)*
      Sqr(Abs(MuInput))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),
      msqInput(1),msqInput(2)) + 3*MuInput*AuInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yu(0
      ,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msqInput(2),msqInput(0)) + 3*
      MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(
      msqInput(2),msqInput(2),msqInput(2),msqInput(1)) + 3*MuInput*AuInput(2,2)*
      Quad(Yu(2,2))*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(2),msqInput(2),
      msqInput(2)) + 0.36*M1Input*MuInput*Quad(g1)*TCD2t(M1Input,M1Input,Abs(
      MuInput),Abs(MuInput)) + 0.6*(M1Input + M2Input)*MuInput*Sqr(g1)*Sqr(g2)*
      TCD2t(M1Input,M2Input,Abs(MuInput),Abs(MuInput)) + 3*M2Input*MuInput*Quad(g2
      )*TCD2t(M2Input,M2Input,Abs(MuInput),Abs(MuInput)) + 3*MuInput*AdInput(0,0)*
      AdInput(0,1)*AdInput(1,1)*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(
      1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*MuInput*AdInput(0,0)*AdInput(0,1)*
      AdInput(1,0)*TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(0))*Yd(0,0)*
      Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*MuInput*AdInput(0,1)*AdInput(1,0)*AdInput(1,1)*
      TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0
      )*Yd(1,1) + 3*MuInput*AdInput(0,0)*AdInput(1,0)*AdInput(1,1)*TCD0(msdInput(1
      ),msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*
      MuInput*AdInput(0,0)*AdInput(0,2)*AdInput(1,2)*TCD0(msdInput(0),msdInput(1),
      msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*MuInput*AdInput
      (0,0)*AdInput(0,2)*AdInput(1,0)*TCD0(msdInput(0),msdInput(1),msqInput(2),
      msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*MuInput*AdInput(0,2)*
      AdInput(1,0)*AdInput(1,2)*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(
      2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*MuInput*AdInput(0,0)*AdInput(1,0)*
      AdInput(1,2)*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*
      Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*MuInput*AdInput(0,1)*AdInput(0,2)*AdInput(1,2)*
      TCD0(msdInput(0),msdInput(1),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1
      )*Yd(1,2) + 3*MuInput*AdInput(0,1)*AdInput(0,2)*AdInput(1,1)*TCD0(msdInput(0
      ),msdInput(1),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*
      MuInput*AdInput(0,2)*AdInput(1,1)*AdInput(1,2)*TCD0(msdInput(1),msdInput(0),
      msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*MuInput*AdInput
      (0,1)*AdInput(1,1)*AdInput(1,2)*TCD0(msdInput(1),msdInput(0),msqInput(2),
      msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*MuInput*AdInput(0,0)*TCC0(
      msdInput(0),msdInput(1),msqInput(0))*Yd(0,0)*Yd(1,0)*(Yd(0,0)*Yd(1,0) + Yd(0
      ,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) + 3*MuInput*AdInput(1,0)*TCC0(msdInput(1),
      msdInput(0),msqInput(0))*Yd(0,0)*Yd(1,0)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1)
      + Yd(0,2)*Yd(1,2)) + 3*MuInput*AdInput(0,1)*TCC0(msdInput(0),msdInput(1),
      msqInput(1))*Yd(0,1)*Yd(1,1)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd
      (1,2)) + 3*MuInput*AdInput(1,1)*TCC0(msdInput(1),msdInput(0),msqInput(1))*Yd
      (0,1)*Yd(1,1)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) + 3*
      MuInput*AdInput(0,2)*TCC0(msdInput(0),msdInput(1),msqInput(2))*Yd(0,2)*Yd(1,
      2)*(Yd(0,0)*Yd(1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) + 3*MuInput*AdInput
      (1,2)*TCC0(msdInput(1),msdInput(0),msqInput(2))*Yd(0,2)*Yd(1,2)*(Yd(0,0)*Yd(
      1,0) + Yd(0,1)*Yd(1,1) + Yd(0,2)*Yd(1,2)) + 3*MuInput*AdInput(0,0)*AdInput(0
      ,1)*AdInput(2,1)*TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(1))*Yd(0,
      0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(0,0)*AdInput(0,1)*AdInput(2,0
      )*TCD0(msdInput(0),msdInput(2),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2
      ,0)*Yd(2,1) + 3*MuInput*AdInput(0,1)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput
      (2),msdInput(0),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3
      *MuInput*AdInput(0,0)*AdInput(2,0)*AdInput(2,1)*TCD0(msdInput(2),msdInput(0)
      ,msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*
      AdInput(1,0)*AdInput(1,1)*AdInput(2,1)*TCD0(msdInput(1),msdInput(2),msqInput
      (0),msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(1,0)*
      AdInput(1,1)*AdInput(2,0)*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(
      0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(1,1)*AdInput(2,0)*
      AdInput(2,1)*TCD0(msdInput(2),msdInput(1),msqInput(0),msqInput(1))*Yd(1,0)*
      Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(1,0)*AdInput(2,0)*AdInput(2,1)*
      TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0
      )*Yd(2,1) + 3*MuInput*AdInput(0,0)*TCC0(msdInput(0),msqInput(0),msqInput(1))
      *Yd(0,0)*Yd(0,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) + 3*
      MuInput*AdInput(0,1)*TCC0(msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,
      1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) + 3*MuInput*AdInput
      (1,0)*TCC0(msdInput(1),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(
      0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) + 3*MuInput*AdInput(1,1)*TCC0(
      msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*(Yd(0,0)*Yd(0,1) + Yd(1
      ,0)*Yd(1,1) + Yd(2,0)*Yd(2,1)) + 3*MuInput*AdInput(2,0)*TCC0(msdInput(2),
      msqInput(0),msqInput(1))*Yd(2,0)*Yd(2,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1)
      + Yd(2,0)*Yd(2,1)) + 3*MuInput*AdInput(2,1)*TCC0(msdInput(2),msqInput(1),
      msqInput(0))*Yd(2,0)*Yd(2,1)*(Yd(0,0)*Yd(0,1) + Yd(1,0)*Yd(1,1) + Yd(2,0)*Yd
      (2,1)) + 3*MuInput*AdInput(0,0)*AdInput(0,2)*AdInput(2,2)*TCD0(msdInput(0),
      msdInput(2),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*
      MuInput*AdInput(0,0)*AdInput(0,2)*AdInput(2,0)*TCD0(msdInput(0),msdInput(2),
      msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput
      (0,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2),msdInput(0),msqInput(0),
      msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput(0,0)*
      AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput(
      0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput(1,0)*AdInput(1,2)*
      AdInput(2,2)*TCD0(msdInput(1),msdInput(2),msqInput(0),msqInput(2))*Yd(1,0)*
      Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput(1,0)*AdInput(1,2)*AdInput(2,0)*
      TCD0(msdInput(1),msdInput(2),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0
      )*Yd(2,2) + 3*MuInput*AdInput(1,2)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2
      ),msdInput(1),msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*
      MuInput*AdInput(1,0)*AdInput(2,0)*AdInput(2,2)*TCD0(msdInput(2),msdInput(1),
      msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput
      (0,1)*AdInput(0,2)*AdInput(2,2)*TCD0(msdInput(0),msdInput(2),msqInput(1),
      msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(0,1)*
      AdInput(0,2)*AdInput(2,1)*TCD0(msdInput(0),msdInput(2),msqInput(2),msqInput(
      1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(0,2)*AdInput(2,1)*
      AdInput(2,2)*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(2))*Yd(0,1)*
      Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(0,1)*AdInput(2,1)*AdInput(2,2)*
      TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1
      )*Yd(2,2) + 3*MuInput*AdInput(1,1)*AdInput(1,2)*AdInput(2,2)*TCD0(msdInput(1
      ),msdInput(2),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*
      MuInput*AdInput(1,1)*AdInput(1,2)*AdInput(2,1)*TCD0(msdInput(1),msdInput(2),
      msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput
      (1,2)*AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(2),msdInput(1),msqInput(1),
      msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(1,1)*
      AdInput(2,1)*AdInput(2,2)*TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(
      1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(0,0)*TCC0(msdInput(0
      ),msdInput(2),msqInput(0))*Yd(0,0)*Yd(2,0)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1
      ) + Yd(0,2)*Yd(2,2)) + 3*MuInput*AdInput(2,0)*TCC0(msdInput(2),msdInput(0),
      msqInput(0))*Yd(0,0)*Yd(2,0)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd
      (2,2)) + 3*MuInput*AdInput(0,1)*TCC0(msdInput(0),msdInput(2),msqInput(1))*Yd
      (0,1)*Yd(2,1)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) + 3*
      MuInput*AdInput(2,1)*TCC0(msdInput(2),msdInput(0),msqInput(1))*Yd(0,1)*Yd(2,
      1)*(Yd(0,0)*Yd(2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) + 3*MuInput*AdInput
      (0,2)*TCC0(msdInput(0),msdInput(2),msqInput(2))*Yd(0,2)*Yd(2,2)*(Yd(0,0)*Yd(
      2,0) + Yd(0,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) + 3*MuInput*AdInput(2,2)*TCC0(
      msdInput(2),msdInput(0),msqInput(2))*Yd(0,2)*Yd(2,2)*(Yd(0,0)*Yd(2,0) + Yd(0
      ,1)*Yd(2,1) + Yd(0,2)*Yd(2,2)) + 3*MuInput*AdInput(1,0)*TCC0(msdInput(1),
      msdInput(2),msqInput(0))*Yd(1,0)*Yd(2,0)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1)
      + Yd(1,2)*Yd(2,2)) + 3*MuInput*AdInput(2,0)*TCC0(msdInput(2),msdInput(1),
      msqInput(0))*Yd(1,0)*Yd(2,0)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd
      (2,2)) + 3*MuInput*AdInput(1,1)*TCC0(msdInput(1),msdInput(2),msqInput(1))*Yd
      (1,1)*Yd(2,1)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) + 3*
      MuInput*AdInput(2,1)*TCC0(msdInput(2),msdInput(1),msqInput(1))*Yd(1,1)*Yd(2,
      1)*(Yd(1,0)*Yd(2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) + 3*MuInput*AdInput
      (1,2)*TCC0(msdInput(1),msdInput(2),msqInput(2))*Yd(1,2)*Yd(2,2)*(Yd(1,0)*Yd(
      2,0) + Yd(1,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) + 3*MuInput*AdInput(2,2)*TCC0(
      msdInput(2),msdInput(1),msqInput(2))*Yd(1,2)*Yd(2,2)*(Yd(1,0)*Yd(2,0) + Yd(1
      ,1)*Yd(2,1) + Yd(1,2)*Yd(2,2)) + 3*MuInput*AdInput(0,0)*TCC0(msdInput(0),
      msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2)
      + Yd(2,0)*Yd(2,2)) + 3*MuInput*AdInput(0,2)*TCC0(msdInput(0),msqInput(2),
      msqInput(0))*Yd(0,0)*Yd(0,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd
      (2,2)) + 3*MuInput*AdInput(1,0)*TCC0(msdInput(1),msqInput(0),msqInput(2))*Yd
      (1,0)*Yd(1,2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) + 3*
      MuInput*AdInput(1,2)*TCC0(msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,
      2)*(Yd(0,0)*Yd(0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) + 3*MuInput*AdInput
      (2,0)*TCC0(msdInput(2),msqInput(0),msqInput(2))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(
      0,2) + Yd(1,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) + 3*MuInput*AdInput(2,2)*TCC0(
      msdInput(2),msqInput(2),msqInput(0))*Yd(2,0)*Yd(2,2)*(Yd(0,0)*Yd(0,2) + Yd(1
      ,0)*Yd(1,2) + Yd(2,0)*Yd(2,2)) + 3*MuInput*AdInput(0,1)*TCC0(msdInput(0),
      msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2)
      + Yd(2,1)*Yd(2,2)) + 3*MuInput*AdInput(0,2)*TCC0(msdInput(0),msqInput(2),
      msqInput(1))*Yd(0,1)*Yd(0,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd
      (2,2)) + 3*MuInput*AdInput(1,1)*TCC0(msdInput(1),msqInput(1),msqInput(2))*Yd
      (1,1)*Yd(1,2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) + 3*
      MuInput*AdInput(1,2)*TCC0(msdInput(1),msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,
      2)*(Yd(0,1)*Yd(0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) + 3*MuInput*AdInput
      (2,1)*TCC0(msdInput(2),msqInput(1),msqInput(2))*Yd(2,1)*Yd(2,2)*(Yd(0,1)*Yd(
      0,2) + Yd(1,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) + 3*MuInput*AdInput(2,2)*TCC0(
      msdInput(2),msqInput(2),msqInput(1))*Yd(2,1)*Yd(2,2)*(Yd(0,1)*Yd(0,2) + Yd(1
      ,1)*Yd(1,2) + Yd(2,1)*Yd(2,2)) + MuInput*AeInput(0,0)*AeInput(0,1)*AeInput(1
      ,1)*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye
      (1,0)*Ye(1,1) + MuInput*AeInput(0,0)*AeInput(0,1)*AeInput(1,0)*TCD0(mseInput
      (0),mseInput(1),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) +
      MuInput*AeInput(0,1)*AeInput(1,0)*AeInput(1,1)*TCD0(mseInput(1),mseInput(0),
      mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + MuInput*AeInput(0
      ,0)*AeInput(1,0)*AeInput(1,1)*TCD0(mseInput(1),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + MuInput*AeInput(0,0)*AeInput(
      0,2)*AeInput(1,2)*TCD0(mseInput(0),mseInput(1),mslInput(0),mslInput(2))*Ye(0
      ,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + MuInput*AeInput(0,0)*AeInput(0,2)*AeInput(1,0)
      *TCD0(mseInput(0),mseInput(1),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,
      0)*Ye(1,2) + MuInput*AeInput(0,2)*AeInput(1,0)*AeInput(1,2)*TCD0(mseInput(1)
      ,mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) +
      MuInput*AeInput(0,0)*AeInput(1,0)*AeInput(1,2)*TCD0(mseInput(1),mseInput(0),
      mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + MuInput*AeInput(0
      ,1)*AeInput(0,2)*AeInput(1,2)*TCD0(mseInput(0),mseInput(1),mslInput(1),
      mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + MuInput*AeInput(0,1)*AeInput(
      0,2)*AeInput(1,1)*TCD0(mseInput(0),mseInput(1),mslInput(2),mslInput(1))*Ye(0
      ,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + MuInput*AeInput(0,2)*AeInput(1,1)*AeInput(1,2)
      *TCD0(mseInput(1),mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,
      1)*Ye(1,2) + MuInput*AeInput(0,1)*AeInput(1,1)*AeInput(1,2)*TCD0(mseInput(1)
      ,mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) +
      MuInput*AeInput(0,0)*TCC0(mseInput(0),mseInput(1),mslInput(0))*Ye(0,0)*Ye(1,
      0)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye(1,2)) + MuInput*AeInput(1
      ,0)*TCC0(mseInput(1),mseInput(0),mslInput(0))*Ye(0,0)*Ye(1,0)*(Ye(0,0)*Ye(1,
      0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye(1,2)) + MuInput*AeInput(0,1)*TCC0(mseInput
      (0),mseInput(1),mslInput(1))*Ye(0,1)*Ye(1,1)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1
      ,1) + Ye(0,2)*Ye(1,2)) + MuInput*AeInput(1,1)*TCC0(mseInput(1),mseInput(0),
      mslInput(1))*Ye(0,1)*Ye(1,1)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye
      (1,2)) + MuInput*AeInput(0,2)*TCC0(mseInput(0),mseInput(1),mslInput(2))*Ye(0
      ,2)*Ye(1,2)*(Ye(0,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye(1,2)) + MuInput*
      AeInput(1,2)*TCC0(mseInput(1),mseInput(0),mslInput(2))*Ye(0,2)*Ye(1,2)*(Ye(0
      ,0)*Ye(1,0) + Ye(0,1)*Ye(1,1) + Ye(0,2)*Ye(1,2)) + MuInput*AeInput(0,0)*
      AeInput(0,1)*AeInput(2,1)*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(
      1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(0,0)*AeInput(0,1)*
      AeInput(2,0)*TCD0(mseInput(0),mseInput(2),mslInput(1),mslInput(0))*Ye(0,0)*
      Ye(0,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(0,1)*AeInput(2,0)*AeInput(2,1)*
      TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(2,0
      )*Ye(2,1) + MuInput*AeInput(0,0)*AeInput(2,0)*AeInput(2,1)*TCD0(mseInput(2),
      mseInput(0),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) +
      MuInput*AeInput(1,0)*AeInput(1,1)*AeInput(2,1)*TCD0(mseInput(1),mseInput(2),
      mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(1
      ,0)*AeInput(1,1)*AeInput(2,0)*TCD0(mseInput(1),mseInput(2),mslInput(1),
      mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(1,1)*AeInput(
      2,0)*AeInput(2,1)*TCD0(mseInput(2),mseInput(1),mslInput(0),mslInput(1))*Ye(1
      ,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(1,0)*AeInput(2,0)*AeInput(2,1)
      *TCD0(mseInput(2),mseInput(1),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,
      0)*Ye(2,1) + MuInput*AeInput(0,0)*TCC0(mseInput(0),mslInput(0),mslInput(1))*
      Ye(0,0)*Ye(0,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) +
      MuInput*AeInput(0,1)*TCC0(mseInput(0),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,
      1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + MuInput*AeInput(1
      ,0)*TCC0(mseInput(1),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1,1)*(Ye(0,0)*Ye(0,
      1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + MuInput*AeInput(1,1)*TCC0(mseInput
      (1),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1
      ,1) + Ye(2,0)*Ye(2,1)) + MuInput*AeInput(2,0)*TCC0(mseInput(2),mslInput(0),
      mslInput(1))*Ye(2,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye
      (2,1)) + MuInput*AeInput(2,1)*TCC0(mseInput(2),mslInput(1),mslInput(0))*Ye(2
      ,0)*Ye(2,1)*(Ye(0,0)*Ye(0,1) + Ye(1,0)*Ye(1,1) + Ye(2,0)*Ye(2,1)) + MuInput*
      AeInput(0,0)*AeInput(0,2)*AeInput(2,2)*TCD0(mseInput(0),mseInput(2),mslInput
      (0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(0,0)*
      AeInput(0,2)*AeInput(2,0)*TCD0(mseInput(0),mseInput(2),mslInput(2),mslInput(
      0))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(0,2)*AeInput(2,0)*
      AeInput(2,2)*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*
      Ye(0,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(0,0)*AeInput(2,0)*AeInput(2,2)*
      TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2,0
      )*Ye(2,2) + MuInput*AeInput(1,0)*AeInput(1,2)*AeInput(2,2)*TCD0(mseInput(1),
      mseInput(2),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) +
      MuInput*AeInput(1,0)*AeInput(1,2)*AeInput(2,0)*TCD0(mseInput(1),mseInput(2),
      mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(1
      ,2)*AeInput(2,0)*AeInput(2,2)*TCD0(mseInput(2),mseInput(1),mslInput(0),
      mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(1,0)*AeInput(
      2,0)*AeInput(2,2)*TCD0(mseInput(2),mseInput(1),mslInput(2),mslInput(0))*Ye(1
      ,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(0,1)*AeInput(0,2)*AeInput(2,2)
      *TCD0(mseInput(0),mseInput(2),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,
      1)*Ye(2,2) + MuInput*AeInput(0,1)*AeInput(0,2)*AeInput(2,1)*TCD0(mseInput(0)
      ,mseInput(2),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) +
      MuInput*AeInput(0,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2),mseInput(0),
      mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(0
      ,1)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2),mseInput(0),mslInput(2),
      mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(1,1)*AeInput(
      1,2)*AeInput(2,2)*TCD0(mseInput(1),mseInput(2),mslInput(1),mslInput(2))*Ye(1
      ,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(1,1)*AeInput(1,2)*AeInput(2,1)
      *TCD0(mseInput(1),mseInput(2),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,
      1)*Ye(2,2) + MuInput*AeInput(1,2)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2)
      ,mseInput(1),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) +
      MuInput*AeInput(1,1)*AeInput(2,1)*AeInput(2,2)*TCD0(mseInput(2),mseInput(1),
      mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(0
      ,0)*TCC0(mseInput(0),mseInput(2),mslInput(0))*Ye(0,0)*Ye(2,0)*(Ye(0,0)*Ye(2,
      0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) + MuInput*AeInput(2,0)*TCC0(mseInput
      (2),mseInput(0),mslInput(0))*Ye(0,0)*Ye(2,0)*(Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2
      ,1) + Ye(0,2)*Ye(2,2)) + MuInput*AeInput(0,1)*TCC0(mseInput(0),mseInput(2),
      mslInput(1))*Ye(0,1)*Ye(2,1)*(Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye
      (2,2)) + MuInput*AeInput(2,1)*TCC0(mseInput(2),mseInput(0),mslInput(1))*Ye(0
      ,1)*Ye(2,1)*(Ye(0,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) + MuInput*
      AeInput(0,2)*TCC0(mseInput(0),mseInput(2),mslInput(2))*Ye(0,2)*Ye(2,2)*(Ye(0
      ,0)*Ye(2,0) + Ye(0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) + MuInput*AeInput(2,2)*TCC0
      (mseInput(2),mseInput(0),mslInput(2))*Ye(0,2)*Ye(2,2)*(Ye(0,0)*Ye(2,0) + Ye(
      0,1)*Ye(2,1) + Ye(0,2)*Ye(2,2)) + MuInput*AeInput(1,0)*TCC0(mseInput(1),
      mseInput(2),mslInput(0))*Ye(1,0)*Ye(2,0)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1)
      + Ye(1,2)*Ye(2,2)) + MuInput*AeInput(2,0)*TCC0(mseInput(2),mseInput(1),
      mslInput(0))*Ye(1,0)*Ye(2,0)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2)*Ye
      (2,2)) + MuInput*AeInput(1,1)*TCC0(mseInput(1),mseInput(2),mslInput(1))*Ye(1
      ,1)*Ye(2,1)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2)*Ye(2,2)) + MuInput*
      AeInput(2,1)*TCC0(mseInput(2),mseInput(1),mslInput(1))*Ye(1,1)*Ye(2,1)*(Ye(1
      ,0)*Ye(2,0) + Ye(1,1)*Ye(2,1) + Ye(1,2)*Ye(2,2)) + MuInput*AeInput(1,2)*TCC0
      (mseInput(1),mseInput(2),mslInput(2))*Ye(1,2)*Ye(2,2)*(Ye(1,0)*Ye(2,0) + Ye(
      1,1)*Ye(2,1) + Ye(1,2)*Ye(2,2)) + MuInput*AeInput(2,2)*TCC0(mseInput(2),
      mseInput(1),mslInput(2))*Ye(1,2)*Ye(2,2)*(Ye(1,0)*Ye(2,0) + Ye(1,1)*Ye(2,1)
      + Ye(1,2)*Ye(2,2)) + MuInput*AeInput(0,0)*TCC0(mseInput(0),mslInput(0),
      mslInput(2))*Ye(0,0)*Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye
      (2,2)) + MuInput*AeInput(0,2)*TCC0(mseInput(0),mslInput(2),mslInput(0))*Ye(0
      ,0)*Ye(0,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + MuInput*
      AeInput(1,0)*TCC0(mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*(Ye(0
      ,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + MuInput*AeInput(1,2)*TCC0
      (mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*(Ye(0,0)*Ye(0,2) + Ye(
      1,0)*Ye(1,2) + Ye(2,0)*Ye(2,2)) + MuInput*AeInput(2,0)*TCC0(mseInput(2),
      mslInput(0),mslInput(2))*Ye(2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2)
      + Ye(2,0)*Ye(2,2)) + MuInput*AeInput(2,2)*TCC0(mseInput(2),mslInput(2),
      mslInput(0))*Ye(2,0)*Ye(2,2)*(Ye(0,0)*Ye(0,2) + Ye(1,0)*Ye(1,2) + Ye(2,0)*Ye
      (2,2)) + MuInput*AeInput(0,1)*TCC0(mseInput(0),mslInput(1),mslInput(2))*Ye(0
      ,1)*Ye(0,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + MuInput*
      AeInput(0,2)*TCC0(mseInput(0),mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*(Ye(0
      ,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + MuInput*AeInput(1,1)*TCC0
      (mseInput(1),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(
      1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + MuInput*AeInput(1,2)*TCC0(mseInput(1),
      mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2)
      + Ye(2,1)*Ye(2,2)) + MuInput*AeInput(2,1)*TCC0(mseInput(2),mslInput(1),
      mslInput(2))*Ye(2,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye
      (2,2)) + MuInput*AeInput(2,2)*TCC0(mseInput(2),mslInput(2),mslInput(1))*Ye(2
      ,1)*Ye(2,2)*(Ye(0,1)*Ye(0,2) + Ye(1,1)*Ye(1,2) + Ye(2,1)*Ye(2,2)) + 3*
      MuInput*AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msqInput
      (0),msqInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*MuInput*AuInput(0,0)*
      Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msqInput(1),msqInput(0))*Yu(0
      ,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*MuInput*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(0),msqInput(0),msqInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(
      1,1) + 3*MuInput*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0)
      ,msqInput(1),msqInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*MuInput*
      AuInput(1,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msqInput(0),
      msqInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*MuInput*AuInput(0,0)*Sqr(
      Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msqInput(1),msqInput(0))*Yu(0,0)*
      Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(0),msqInput(0),msqInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(
      1,2) + 3*MuInput*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0)
      ,msqInput(1),msqInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*MuInput*
      AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msqInput(0),
      msqInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*MuInput*AuInput(0,1)*Sqr(
      Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msqInput(1),msqInput(0))*Yu(0,1)*
      Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(1),msqInput(0),msqInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(
      1,2) + 3*MuInput*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1)
      ,msqInput(1),msqInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*MuInput*
      AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msqInput(0),
      msqInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(0,0)*Sqr(
      Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msqInput(2),msqInput(0))*Yu(0,0)*
      Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(0),msqInput(0),msqInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(
      2,1) + 3*MuInput*AuInput(0,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0)
      ,msqInput(2),msqInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*
      AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msqInput(1),
      msqInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(1,0)*Sqr(
      Abs(MuInput))*TCD0(msqInput(0),msqInput(1),msqInput(2),msqInput(1))*Yu(1,0)*
      Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(
      msqInput(1),msqInput(0),msqInput(1),msqInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(
      2,1) + 3*MuInput*AuInput(1,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(0)
      ,msqInput(2),msqInput(1))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*
      AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msqInput(0),
      msqInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput(0,0)*Sqr(
      Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msqInput(2),msqInput(0))*Yu(0,0)*
      Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(0),msqInput(0),msqInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(
      2,2) + 3*MuInput*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0)
      ,msqInput(2),msqInput(0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*
      AuInput(2,0)*Sqr(Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msqInput(1),
      msqInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput(1,0)*Sqr(
      Abs(MuInput))*TCD0(msqInput(0),msqInput(2),msqInput(2),msqInput(1))*Yu(1,0)*
      Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(0),msqInput(1),msqInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(
      2,2) + 3*MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(0)
      ,msqInput(2),msqInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*
      AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msqInput(0),
      msqInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(0,1)*Sqr(
      Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msqInput(2),msqInput(0))*Yu(0,1)*
      Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(1),msqInput(0),msqInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(
      2,2) + 3*MuInput*AuInput(0,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1)
      ,msqInput(2),msqInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*
      AuInput(2,1)*Sqr(Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msqInput(1),
      msqInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(1,1)*Sqr(
      Abs(MuInput))*TCD0(msqInput(1),msqInput(2),msqInput(2),msqInput(1))*Yu(1,1)*
      Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msqInput(2),msqInput(1),msqInput(1),msqInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(
      2,2) + 3*MuInput*AuInput(1,2)*Sqr(Abs(MuInput))*TCD0(msqInput(2),msqInput(1)
      ,msqInput(2),msqInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2))));
   MODEL->set_Lambda7(Re((0.00016040597272944275*MuInput*Quad(Yu(2,2))*Sqr(g3)*(2
      + (0.3333333333333333*Cube(AuInput(2,2)))/Cube(MSUSY) - Sqr(AuInput(2,2))/
      Sqr(MSUSY))*UnitStep(-2 + LambdaLoopOrder))/MSUSY + 0.006332573977646111*
      UnitStep(-1 + LambdaLoopOrder)*(0.3*MuInput*AdInput(0,0)*Sqr(g1)*Sqr(Yd(0,0)
      )*TCC0(msdInput(0),msdInput(0),msqInput(0)) + 0.3*MuInput*AdInput(0,1)*Sqr(
      g1)*Sqr(Yd(0,1))*TCC0(msdInput(0),msdInput(0),msqInput(1)) + 0.3*MuInput*
      AdInput(0,2)*Sqr(g1)*Sqr(Yd(0,2))*TCC0(msdInput(0),msdInput(0),msqInput(2))
      + 0.25*MuInput*AdInput(0,0)*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yd(0,0))*TCC0(
      msdInput(0),msqInput(0),msqInput(0)) + 0.25*MuInput*AdInput(0,1)*(0.6*Sqr(g1
      ) + 3*Sqr(g2))*Sqr(Yd(0,1))*TCC0(msdInput(0),msqInput(1),msqInput(1)) + 0.25
      *MuInput*AdInput(0,2)*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yd(0,2))*TCC0(msdInput(0
      ),msqInput(2),msqInput(2)) + 0.3*MuInput*AdInput(1,0)*Sqr(g1)*Sqr(Yd(1,0))*
      TCC0(msdInput(1),msdInput(1),msqInput(0)) + 0.3*MuInput*AdInput(1,1)*Sqr(g1)
      *Sqr(Yd(1,1))*TCC0(msdInput(1),msdInput(1),msqInput(1)) + 0.3*MuInput*
      AdInput(1,2)*Sqr(g1)*Sqr(Yd(1,2))*TCC0(msdInput(1),msdInput(1),msqInput(2))
      + 0.25*MuInput*AdInput(1,0)*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yd(1,0))*TCC0(
      msdInput(1),msqInput(0),msqInput(0)) + 0.25*MuInput*AdInput(1,1)*(0.6*Sqr(g1
      ) + 3*Sqr(g2))*Sqr(Yd(1,1))*TCC0(msdInput(1),msqInput(1),msqInput(1)) + 0.25
      *MuInput*AdInput(1,2)*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yd(1,2))*TCC0(msdInput(1
      ),msqInput(2),msqInput(2)) + 0.3*MuInput*AdInput(2,0)*Sqr(g1)*Sqr(Yd(2,0))*
      TCC0(msdInput(2),msdInput(2),msqInput(0)) + 0.3*MuInput*AdInput(2,1)*Sqr(g1)
      *Sqr(Yd(2,1))*TCC0(msdInput(2),msdInput(2),msqInput(1)) + 0.3*MuInput*
      AdInput(2,2)*Sqr(g1)*Sqr(Yd(2,2))*TCC0(msdInput(2),msdInput(2),msqInput(2))
      + 0.25*MuInput*AdInput(2,0)*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yd(2,0))*TCC0(
      msdInput(2),msqInput(0),msqInput(0)) + 0.25*MuInput*AdInput(2,1)*(0.6*Sqr(g1
      ) + 3*Sqr(g2))*Sqr(Yd(2,1))*TCC0(msdInput(2),msqInput(1),msqInput(1)) + 0.25
      *MuInput*AdInput(2,2)*(0.6*Sqr(g1) + 3*Sqr(g2))*Sqr(Yd(2,2))*TCC0(msdInput(2
      ),msqInput(2),msqInput(2)) + 0.3*MuInput*AeInput(0,0)*Sqr(g1)*Sqr(Ye(0,0))*
      TCC0(mseInput(0),mseInput(0),mslInput(0)) + 0.3*MuInput*AeInput(0,1)*Sqr(g1)
      *Sqr(Ye(0,1))*TCC0(mseInput(0),mseInput(0),mslInput(1)) + 0.3*MuInput*
      AeInput(0,2)*Sqr(g1)*Sqr(Ye(0,2))*TCC0(mseInput(0),mseInput(0),mslInput(2))
      + 0.25*MuInput*AeInput(0,0)*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Ye(0,0))*TCC0(
      mseInput(0),mslInput(0),mslInput(0)) + 0.25*MuInput*AeInput(0,1)*(-0.6*Sqr(
      g1) + Sqr(g2))*Sqr(Ye(0,1))*TCC0(mseInput(0),mslInput(1),mslInput(1)) + 0.25
      *MuInput*AeInput(0,2)*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Ye(0,2))*TCC0(mseInput(0)
      ,mslInput(2),mslInput(2)) + 0.3*MuInput*AeInput(1,0)*Sqr(g1)*Sqr(Ye(1,0))*
      TCC0(mseInput(1),mseInput(1),mslInput(0)) + 0.3*MuInput*AeInput(1,1)*Sqr(g1)
      *Sqr(Ye(1,1))*TCC0(mseInput(1),mseInput(1),mslInput(1)) + 0.3*MuInput*
      AeInput(1,2)*Sqr(g1)*Sqr(Ye(1,2))*TCC0(mseInput(1),mseInput(1),mslInput(2))
      + 0.25*MuInput*AeInput(1,0)*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Ye(1,0))*TCC0(
      mseInput(1),mslInput(0),mslInput(0)) + 0.25*MuInput*AeInput(1,1)*(-0.6*Sqr(
      g1) + Sqr(g2))*Sqr(Ye(1,1))*TCC0(mseInput(1),mslInput(1),mslInput(1)) + 0.25
      *MuInput*AeInput(1,2)*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Ye(1,2))*TCC0(mseInput(1)
      ,mslInput(2),mslInput(2)) + 0.3*MuInput*AeInput(2,0)*Sqr(g1)*Sqr(Ye(2,0))*
      TCC0(mseInput(2),mseInput(2),mslInput(0)) + 0.3*MuInput*AeInput(2,1)*Sqr(g1)
      *Sqr(Ye(2,1))*TCC0(mseInput(2),mseInput(2),mslInput(1)) + 0.3*MuInput*
      AeInput(2,2)*Sqr(g1)*Sqr(Ye(2,2))*TCC0(mseInput(2),mseInput(2),mslInput(2))
      + 0.25*MuInput*AeInput(2,0)*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Ye(2,0))*TCC0(
      mseInput(2),mslInput(0),mslInput(0)) + 0.25*MuInput*AeInput(2,1)*(-0.6*Sqr(
      g1) + Sqr(g2))*Sqr(Ye(2,1))*TCC0(mseInput(2),mslInput(1),mslInput(1)) + 0.25
      *MuInput*AeInput(2,2)*(-0.6*Sqr(g1) + Sqr(g2))*Sqr(Ye(2,2))*TCC0(mseInput(2)
      ,mslInput(2),mslInput(2)) + (0.25*MuInput*AuInput(0,0)*(0.6*Sqr(g1) - 3*Sqr(
      g2))*Sqr(Yu(0,0)) + 3*MuInput*AuInput(0,0)*Sqr(Yu(0,0))*(Sqr(Yu(0,0)) + Sqr(
      Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(0)) + (0.25*
      MuInput*AuInput(1,0)*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yu(1,0)) + 3*MuInput*
      AuInput(1,0)*Sqr(Yu(1,0))*(Sqr(Yu(0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCC0
      (msqInput(0),msqInput(0),msuInput(1)) + (0.25*MuInput*AuInput(2,0)*(0.6*Sqr(
      g1) - 3*Sqr(g2))*Sqr(Yu(2,0)) + 3*MuInput*AuInput(2,0)*Sqr(Yu(2,0))*(Sqr(Yu(
      0,0)) + Sqr(Yu(1,0)) + Sqr(Yu(2,0))))*TCC0(msqInput(0),msqInput(0),msuInput(
      2)) + (-0.6*MuInput*AuInput(0,0)*Sqr(g1)*Sqr(Yu(0,0)) + 3*MuInput*AuInput(0,
      0)*Sqr(Yu(0,0))*(Sqr(Yu(0,0)) + Sqr(Yu(0,1)) + Sqr(Yu(0,2))))*TCC0(msqInput(
      0),msuInput(0),msuInput(0)) + (-0.6*MuInput*AuInput(1,0)*Sqr(g1)*Sqr(Yu(1,0)
      ) + 3*MuInput*AuInput(1,0)*Sqr(Yu(1,0))*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(
      Yu(1,2))))*TCC0(msqInput(0),msuInput(1),msuInput(1)) + (-0.6*MuInput*AuInput
      (2,0)*Sqr(g1)*Sqr(Yu(2,0)) + 3*MuInput*AuInput(2,0)*Sqr(Yu(2,0))*(Sqr(Yu(2,0
      )) + Sqr(Yu(2,1)) + Sqr(Yu(2,2))))*TCC0(msqInput(0),msuInput(2),msuInput(2))
      + (0.25*MuInput*AuInput(0,1)*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yu(0,1)) + 3*
      MuInput*AuInput(0,1)*Sqr(Yu(0,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1)
      )))*TCC0(msqInput(1),msqInput(1),msuInput(0)) + (0.25*MuInput*AuInput(1,1)*(
      0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yu(1,1)) + 3*MuInput*AuInput(1,1)*Sqr(Yu(1,1))*
      (Sqr(Yu(0,1)) + Sqr(Yu(1,1)) + Sqr(Yu(2,1))))*TCC0(msqInput(1),msqInput(1),
      msuInput(1)) + (0.25*MuInput*AuInput(2,1)*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yu(2
      ,1)) + 3*MuInput*AuInput(2,1)*Sqr(Yu(2,1))*(Sqr(Yu(0,1)) + Sqr(Yu(1,1)) +
      Sqr(Yu(2,1))))*TCC0(msqInput(1),msqInput(1),msuInput(2)) + (-0.6*MuInput*
      AuInput(0,1)*Sqr(g1)*Sqr(Yu(0,1)) + 3*MuInput*AuInput(0,1)*Sqr(Yu(0,1))*(Sqr
      (Yu(0,0)) + Sqr(Yu(0,1)) + Sqr(Yu(0,2))))*TCC0(msqInput(1),msuInput(0),
      msuInput(0)) + (-0.6*MuInput*AuInput(1,1)*Sqr(g1)*Sqr(Yu(1,1)) + 3*MuInput*
      AuInput(1,1)*Sqr(Yu(1,1))*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))*TCC0
      (msqInput(1),msuInput(1),msuInput(1)) + (-0.6*MuInput*AuInput(2,1)*Sqr(g1)*
      Sqr(Yu(2,1)) + 3*MuInput*AuInput(2,1)*Sqr(Yu(2,1))*(Sqr(Yu(2,0)) + Sqr(Yu(2,
      1)) + Sqr(Yu(2,2))))*TCC0(msqInput(1),msuInput(2),msuInput(2)) + (0.25*
      MuInput*AuInput(0,2)*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yu(0,2)) + 3*MuInput*
      AuInput(0,2)*Sqr(Yu(0,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0
      (msqInput(2),msqInput(2),msuInput(0)) + (0.25*MuInput*AuInput(1,2)*(0.6*Sqr(
      g1) - 3*Sqr(g2))*Sqr(Yu(1,2)) + 3*MuInput*AuInput(1,2)*Sqr(Yu(1,2))*(Sqr(Yu(
      0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2))))*TCC0(msqInput(2),msqInput(2),msuInput(
      1)) + (0.25*MuInput*AuInput(2,2)*(0.6*Sqr(g1) - 3*Sqr(g2))*Sqr(Yu(2,2)) + 3*
      MuInput*AuInput(2,2)*Sqr(Yu(2,2))*(Sqr(Yu(0,2)) + Sqr(Yu(1,2)) + Sqr(Yu(2,2)
      )))*TCC0(msqInput(2),msqInput(2),msuInput(2)) + (-0.6*MuInput*AuInput(0,2)*
      Sqr(g1)*Sqr(Yu(0,2)) + 3*MuInput*AuInput(0,2)*Sqr(Yu(0,2))*(Sqr(Yu(0,0)) +
      Sqr(Yu(0,1)) + Sqr(Yu(0,2))))*TCC0(msqInput(2),msuInput(0),msuInput(0)) + (-
      0.6*MuInput*AuInput(1,2)*Sqr(g1)*Sqr(Yu(1,2)) + 3*MuInput*AuInput(1,2)*Sqr(
      Yu(1,2))*(Sqr(Yu(1,0)) + Sqr(Yu(1,1)) + Sqr(Yu(1,2))))*TCC0(msqInput(2),
      msuInput(1),msuInput(1)) + (-0.6*MuInput*AuInput(2,2)*Sqr(g1)*Sqr(Yu(2,2)) +
      3*MuInput*AuInput(2,2)*Sqr(Yu(2,2))*(Sqr(Yu(2,0)) + Sqr(Yu(2,1)) + Sqr(Yu(2,
      2))))*TCC0(msqInput(2),msuInput(2),msuInput(2)) + 3*MuInput*AdInput(0,0)*
      Quad(Yd(0,0))*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(0),msqInput(0),
      msqInput(0)) + 3*MuInput*AdInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(
      0,1))*TCD0(msdInput(0),msdInput(0),msqInput(0),msqInput(1)) + 3*MuInput*
      AdInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,2))*TCD0(msdInput(0),
      msdInput(0),msqInput(0),msqInput(2)) + 3*MuInput*AdInput(0,0)*Sqr(Abs(
      MuInput))*Sqr(Yd(0,0))*Sqr(Yd(0,1))*TCD0(msdInput(0),msdInput(0),msqInput(1)
      ,msqInput(0)) + 3*MuInput*AdInput(0,1)*Quad(Yd(0,1))*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(0),msqInput(1),msqInput(1)) + 3*MuInput*AdInput(0,2)*
      Sqr(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),
      msqInput(1),msqInput(2)) + 3*MuInput*AdInput(0,0)*Sqr(Abs(MuInput))*Sqr(Yd(0
      ,0))*Sqr(Yd(0,2))*TCD0(msdInput(0),msdInput(0),msqInput(2),msqInput(0)) + 3*
      MuInput*AdInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(0,2))*TCD0(
      msdInput(0),msdInput(0),msqInput(2),msqInput(1)) + 3*MuInput*AdInput(0,2)*
      Quad(Yd(0,2))*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(0),msqInput(2),
      msqInput(2)) + 3*MuInput*AdInput(1,0)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(
      1,0))*TCD0(msdInput(0),msdInput(1),msqInput(0),msqInput(0)) + 3*MuInput*
      AdInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(0),
      msdInput(1),msqInput(1),msqInput(1)) + 3*MuInput*AdInput(1,2)*Sqr(Abs(
      MuInput))*Sqr(Yd(0,2))*Sqr(Yd(1,2))*TCD0(msdInput(0),msdInput(1),msqInput(2)
      ,msqInput(2)) + 3*MuInput*AdInput(2,0)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd
      (2,0))*TCD0(msdInput(0),msdInput(2),msqInput(0),msqInput(0)) + 3*MuInput*
      AdInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(2,1))*TCD0(msdInput(0),
      msdInput(2),msqInput(1),msqInput(1)) + 3*MuInput*AdInput(2,2)*Sqr(Abs(
      MuInput))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(msdInput(0),msdInput(2),msqInput(2)
      ,msqInput(2)) + 3*MuInput*AdInput(0,0)*Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd
      (1,0))*TCD0(msdInput(1),msdInput(0),msqInput(0),msqInput(0)) + 3*MuInput*
      AdInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yd(0,1))*Sqr(Yd(1,1))*TCD0(msdInput(1),
      msdInput(0),msqInput(1),msqInput(1)) + 3*MuInput*AdInput(0,2)*Sqr(Abs(
      MuInput))*Sqr(Yd(0,2))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(0),msqInput(2)
      ,msqInput(2)) + 3*MuInput*AdInput(1,0)*Quad(Yd(1,0))*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(1),msqInput(0),msqInput(0)) + 3*MuInput*AdInput(1,1)*
      Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(msdInput(1),msdInput(1),
      msqInput(0),msqInput(1)) + 3*MuInput*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1
      ,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(0),msqInput(2)) + 3*
      MuInput*AdInput(1,0)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,1))*TCD0(
      msdInput(1),msdInput(1),msqInput(1),msqInput(0)) + 3*MuInput*AdInput(1,1)*
      Quad(Yd(1,1))*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(1),msqInput(1),
      msqInput(1)) + 3*MuInput*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,1))*Sqr(Yd(
      1,2))*TCD0(msdInput(1),msdInput(1),msqInput(1),msqInput(2)) + 3*MuInput*
      AdInput(1,0)*Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(1,2))*TCD0(msdInput(1),
      msdInput(1),msqInput(2),msqInput(0)) + 3*MuInput*AdInput(1,1)*Sqr(Abs(
      MuInput))*Sqr(Yd(1,1))*Sqr(Yd(1,2))*TCD0(msdInput(1),msdInput(1),msqInput(2)
      ,msqInput(1)) + 3*MuInput*AdInput(1,2)*Quad(Yd(1,2))*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(1),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(2,0)*
      Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(1),msdInput(2),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yd(1
      ,1))*Sqr(Yd(2,1))*TCD0(msdInput(1),msdInput(2),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(1),msdInput(2),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(0,0)*
      Sqr(Abs(MuInput))*Sqr(Yd(0,0))*Sqr(Yd(2,0))*TCD0(msdInput(2),msdInput(0),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(0,1)*Sqr(Abs(MuInput))*Sqr(Yd(0
      ,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(0,2)*Sqr(Abs(MuInput))*Sqr(Yd(0,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(0),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(1,0)*
      Sqr(Abs(MuInput))*Sqr(Yd(1,0))*Sqr(Yd(2,0))*TCD0(msdInput(2),msdInput(1),
      msqInput(0),msqInput(0)) + 3*MuInput*AdInput(1,1)*Sqr(Abs(MuInput))*Sqr(Yd(1
      ,1))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(1)) + 3*
      MuInput*AdInput(1,2)*Sqr(Abs(MuInput))*Sqr(Yd(1,2))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(1),msqInput(2),msqInput(2)) + 3*MuInput*AdInput(2,0)*
      Quad(Yd(2,0))*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(2),msqInput(0),
      msqInput(0)) + 3*MuInput*AdInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(
      2,1))*TCD0(msdInput(2),msdInput(2),msqInput(0),msqInput(1)) + 3*MuInput*
      AdInput(2,2)*Sqr(Abs(MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,2))*TCD0(msdInput(2),
      msdInput(2),msqInput(0),msqInput(2)) + 3*MuInput*AdInput(2,0)*Sqr(Abs(
      MuInput))*Sqr(Yd(2,0))*Sqr(Yd(2,1))*TCD0(msdInput(2),msdInput(2),msqInput(1)
      ,msqInput(0)) + 3*MuInput*AdInput(2,1)*Quad(Yd(2,1))*Sqr(Abs(MuInput))*TCD0(
      msdInput(2),msdInput(2),msqInput(1),msqInput(1)) + 3*MuInput*AdInput(2,2)*
      Sqr(Abs(MuInput))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),
      msqInput(1),msqInput(2)) + 3*MuInput*AdInput(2,0)*Sqr(Abs(MuInput))*Sqr(Yd(2
      ,0))*Sqr(Yd(2,2))*TCD0(msdInput(2),msdInput(2),msqInput(2),msqInput(0)) + 3*
      MuInput*AdInput(2,1)*Sqr(Abs(MuInput))*Sqr(Yd(2,1))*Sqr(Yd(2,2))*TCD0(
      msdInput(2),msdInput(2),msqInput(2),msqInput(1)) + 3*MuInput*AdInput(2,2)*
      Quad(Yd(2,2))*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(2),msqInput(2),
      msqInput(2)) + MuInput*AeInput(0,0)*Quad(Ye(0,0))*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(0),mslInput(0),mslInput(0)) + MuInput*AeInput(0,1)*Sqr(
      Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(mseInput(0),mseInput(0),
      mslInput(0),mslInput(1)) + MuInput*AeInput(0,2)*Sqr(Abs(MuInput))*Sqr(Ye(0,0
      ))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(0),mslInput(2)) +
      MuInput*AeInput(0,0)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,1))*TCD0(
      mseInput(0),mseInput(0),mslInput(1),mslInput(0)) + MuInput*AeInput(0,1)*Quad
      (Ye(0,1))*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(0),mslInput(1),
      mslInput(1)) + MuInput*AeInput(0,2)*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(0,
      2))*TCD0(mseInput(0),mseInput(0),mslInput(1),mslInput(2)) + MuInput*AeInput(
      0,0)*Sqr(Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0
      ),mslInput(2),mslInput(0)) + MuInput*AeInput(0,1)*Sqr(Abs(MuInput))*Sqr(Ye(0
      ,1))*Sqr(Ye(0,2))*TCD0(mseInput(0),mseInput(0),mslInput(2),mslInput(1)) +
      MuInput*AeInput(0,2)*Quad(Ye(0,2))*Sqr(Abs(MuInput))*TCD0(mseInput(0),
      mseInput(0),mslInput(2),mslInput(2)) + MuInput*AeInput(1,0)*Sqr(Abs(MuInput)
      )*Sqr(Ye(0,0))*Sqr(Ye(1,0))*TCD0(mseInput(0),mseInput(1),mslInput(0),
      mslInput(0)) + MuInput*AeInput(1,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(1,
      1))*TCD0(mseInput(0),mseInput(1),mslInput(1),mslInput(1)) + MuInput*AeInput(
      1,2)*Sqr(Abs(MuInput))*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(0),mseInput(1
      ),mslInput(2),mslInput(2)) + MuInput*AeInput(2,0)*Sqr(Abs(MuInput))*Sqr(Ye(0
      ,0))*Sqr(Ye(2,0))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(0)) +
      MuInput*AeInput(2,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(2,1))*TCD0(
      mseInput(0),mseInput(2),mslInput(1),mslInput(1)) + MuInput*AeInput(2,2)*Sqr(
      Abs(MuInput))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(2)) + MuInput*AeInput(0,0)*Sqr(Abs(MuInput))*Sqr(Ye(0,0
      ))*Sqr(Ye(1,0))*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(0)) +
      MuInput*AeInput(0,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,1))*Sqr(Ye(1,1))*TCD0(
      mseInput(1),mseInput(0),mslInput(1),mslInput(1)) + MuInput*AeInput(0,2)*Sqr(
      Abs(MuInput))*Sqr(Ye(0,2))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(0),
      mslInput(2),mslInput(2)) + MuInput*AeInput(1,0)*Quad(Ye(1,0))*Sqr(Abs(
      MuInput))*TCD0(mseInput(1),mseInput(1),mslInput(0),mslInput(0)) + MuInput*
      AeInput(1,1)*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,1))*TCD0(mseInput(1),
      mseInput(1),mslInput(0),mslInput(1)) + MuInput*AeInput(1,2)*Sqr(Abs(MuInput)
      )*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(0),
      mslInput(2)) + MuInput*AeInput(1,0)*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,
      1))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(0)) + MuInput*AeInput(
      1,1)*Quad(Ye(1,1))*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(1),mslInput(1
      ),mslInput(1)) + MuInput*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(Ye(1,1))*Sqr(Ye(
      1,2))*TCD0(mseInput(1),mseInput(1),mslInput(1),mslInput(2)) + MuInput*
      AeInput(1,0)*Sqr(Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(1,2))*TCD0(mseInput(1),
      mseInput(1),mslInput(2),mslInput(0)) + MuInput*AeInput(1,1)*Sqr(Abs(MuInput)
      )*Sqr(Ye(1,1))*Sqr(Ye(1,2))*TCD0(mseInput(1),mseInput(1),mslInput(2),
      mslInput(1)) + MuInput*AeInput(1,2)*Quad(Ye(1,2))*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(1),mslInput(2),mslInput(2)) + MuInput*AeInput(2,0)*Sqr(
      Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(1),mseInput(2),
      mslInput(0),mslInput(0)) + MuInput*AeInput(2,1)*Sqr(Abs(MuInput))*Sqr(Ye(1,1
      ))*Sqr(Ye(2,1))*TCD0(mseInput(1),mseInput(2),mslInput(1),mslInput(1)) +
      MuInput*AeInput(2,2)*Sqr(Abs(MuInput))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(1),mseInput(2),mslInput(2),mslInput(2)) + MuInput*AeInput(0,0)*Sqr(
      Abs(MuInput))*Sqr(Ye(0,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(0),
      mslInput(0),mslInput(0)) + MuInput*AeInput(0,1)*Sqr(Abs(MuInput))*Sqr(Ye(0,1
      ))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(0),mslInput(1),mslInput(1)) +
      MuInput*AeInput(0,2)*Sqr(Abs(MuInput))*Sqr(Ye(0,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(2),mseInput(0),mslInput(2),mslInput(2)) + MuInput*AeInput(1,0)*Sqr(
      Abs(MuInput))*Sqr(Ye(1,0))*Sqr(Ye(2,0))*TCD0(mseInput(2),mseInput(1),
      mslInput(0),mslInput(0)) + MuInput*AeInput(1,1)*Sqr(Abs(MuInput))*Sqr(Ye(1,1
      ))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(1),mslInput(1),mslInput(1)) +
      MuInput*AeInput(1,2)*Sqr(Abs(MuInput))*Sqr(Ye(1,2))*Sqr(Ye(2,2))*TCD0(
      mseInput(2),mseInput(1),mslInput(2),mslInput(2)) + MuInput*AeInput(2,0)*Quad
      (Ye(2,0))*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(2),mslInput(0),
      mslInput(0)) + MuInput*AeInput(2,1)*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,
      1))*TCD0(mseInput(2),mseInput(2),mslInput(0),mslInput(1)) + MuInput*AeInput(
      2,2)*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2
      ),mslInput(0),mslInput(2)) + MuInput*AeInput(2,0)*Sqr(Abs(MuInput))*Sqr(Ye(2
      ,0))*Sqr(Ye(2,1))*TCD0(mseInput(2),mseInput(2),mslInput(1),mslInput(0)) +
      MuInput*AeInput(2,1)*Quad(Ye(2,1))*Sqr(Abs(MuInput))*TCD0(mseInput(2),
      mseInput(2),mslInput(1),mslInput(1)) + MuInput*AeInput(2,2)*Sqr(Abs(MuInput)
      )*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2),mslInput(1),
      mslInput(2)) + MuInput*AeInput(2,0)*Sqr(Abs(MuInput))*Sqr(Ye(2,0))*Sqr(Ye(2,
      2))*TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(0)) + MuInput*AeInput(
      2,1)*Sqr(Abs(MuInput))*Sqr(Ye(2,1))*Sqr(Ye(2,2))*TCD0(mseInput(2),mseInput(2
      ),mslInput(2),mslInput(1)) + MuInput*AeInput(2,2)*Quad(Ye(2,2))*Sqr(Abs(
      MuInput))*TCD0(mseInput(2),mseInput(2),mslInput(2),mslInput(2)) + 3*MuInput*
      Cube(AuInput(0,0))*Quad(Yu(0,0))*TCD0(msqInput(0),msqInput(0),msqInput(0),
      msqInput(0)) + 3*MuInput*AuInput(0,0)*Sqr(AuInput(1,0))*Sqr(Yu(0,0))*Sqr(Yu(
      1,0))*TCD0(msqInput(0),msqInput(0),msqInput(0),msqInput(1)) + 3*MuInput*
      AuInput(0,0)*Sqr(AuInput(2,0))*Sqr(Yu(0,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),
      msqInput(0),msqInput(0),msqInput(2)) + 3*MuInput*AuInput(1,0)*Sqr(AuInput(0,
      0))*Sqr(Yu(0,0))*Sqr(Yu(1,0))*TCD0(msqInput(0),msqInput(0),msqInput(1),
      msqInput(0)) + 3*MuInput*Cube(AuInput(1,0))*Quad(Yu(1,0))*TCD0(msqInput(0),
      msqInput(0),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(1,0)*Sqr(AuInput(2,
      0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),msqInput(0),msqInput(1),
      msqInput(2)) + 3*MuInput*AuInput(2,0)*Sqr(AuInput(0,0))*Sqr(Yu(0,0))*Sqr(Yu(
      2,0))*TCD0(msqInput(0),msqInput(0),msqInput(2),msqInput(0)) + 3*MuInput*
      AuInput(2,0)*Sqr(AuInput(1,0))*Sqr(Yu(1,0))*Sqr(Yu(2,0))*TCD0(msqInput(0),
      msqInput(0),msqInput(2),msqInput(1)) + 3*MuInput*Cube(AuInput(2,0))*Quad(Yu(
      2,0))*TCD0(msqInput(0),msqInput(0),msqInput(2),msqInput(2)) + 3*MuInput*
      AuInput(0,0)*Sqr(AuInput(0,1))*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(0),
      msqInput(1),msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,0)*Sqr(AuInput(1,
      1))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(msqInput(0),msqInput(1),msqInput(1),
      msqInput(1)) + 3*MuInput*AuInput(2,0)*Sqr(AuInput(2,1))*Sqr(Yu(2,0))*Sqr(Yu(
      2,1))*TCD0(msqInput(0),msqInput(1),msqInput(2),msqInput(2)) + 3*MuInput*
      AuInput(0,0)*Sqr(AuInput(0,2))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(0),
      msqInput(2),msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,0)*Sqr(AuInput(1,
      2))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(0),msqInput(2),msqInput(1),
      msqInput(1)) + 3*MuInput*AuInput(2,0)*Sqr(AuInput(2,2))*Sqr(Yu(2,0))*Sqr(Yu(
      2,2))*TCD0(msqInput(0),msqInput(2),msqInput(2),msqInput(2)) + 3*MuInput*
      AuInput(0,1)*Sqr(AuInput(0,0))*Sqr(Yu(0,0))*Sqr(Yu(0,1))*TCD0(msqInput(1),
      msqInput(0),msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,1)*Sqr(AuInput(1,
      0))*Sqr(Yu(1,0))*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(0),msqInput(1),
      msqInput(1)) + 3*MuInput*AuInput(2,1)*Sqr(AuInput(2,0))*Sqr(Yu(2,0))*Sqr(Yu(
      2,1))*TCD0(msqInput(1),msqInput(0),msqInput(2),msqInput(2)) + 3*MuInput*Cube
      (AuInput(0,1))*Quad(Yu(0,1))*TCD0(msqInput(1),msqInput(1),msqInput(0),
      msqInput(0)) + 3*MuInput*AuInput(0,1)*Sqr(AuInput(1,1))*Sqr(Yu(0,1))*Sqr(Yu(
      1,1))*TCD0(msqInput(1),msqInput(1),msqInput(0),msqInput(1)) + 3*MuInput*
      AuInput(0,1)*Sqr(AuInput(2,1))*Sqr(Yu(0,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),
      msqInput(1),msqInput(0),msqInput(2)) + 3*MuInput*AuInput(1,1)*Sqr(AuInput(0,
      1))*Sqr(Yu(0,1))*Sqr(Yu(1,1))*TCD0(msqInput(1),msqInput(1),msqInput(1),
      msqInput(0)) + 3*MuInput*Cube(AuInput(1,1))*Quad(Yu(1,1))*TCD0(msqInput(1),
      msqInput(1),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(1,1)*Sqr(AuInput(2,
      1))*Sqr(Yu(1,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),msqInput(1),msqInput(1),
      msqInput(2)) + 3*MuInput*AuInput(2,1)*Sqr(AuInput(0,1))*Sqr(Yu(0,1))*Sqr(Yu(
      2,1))*TCD0(msqInput(1),msqInput(1),msqInput(2),msqInput(0)) + 3*MuInput*
      AuInput(2,1)*Sqr(AuInput(1,1))*Sqr(Yu(1,1))*Sqr(Yu(2,1))*TCD0(msqInput(1),
      msqInput(1),msqInput(2),msqInput(1)) + 3*MuInput*Cube(AuInput(2,1))*Quad(Yu(
      2,1))*TCD0(msqInput(1),msqInput(1),msqInput(2),msqInput(2)) + 3*MuInput*
      AuInput(0,1)*Sqr(AuInput(0,2))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(1),
      msqInput(2),msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,1)*Sqr(AuInput(1,
      2))*Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(1),msqInput(2),msqInput(1),
      msqInput(1)) + 3*MuInput*AuInput(2,1)*Sqr(AuInput(2,2))*Sqr(Yu(2,1))*Sqr(Yu(
      2,2))*TCD0(msqInput(1),msqInput(2),msqInput(2),msqInput(2)) + 3*MuInput*
      AuInput(0,2)*Sqr(AuInput(0,0))*Sqr(Yu(0,0))*Sqr(Yu(0,2))*TCD0(msqInput(2),
      msqInput(0),msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,2)*Sqr(AuInput(1,
      0))*Sqr(Yu(1,0))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(0),msqInput(1),
      msqInput(1)) + 3*MuInput*AuInput(2,2)*Sqr(AuInput(2,0))*Sqr(Yu(2,0))*Sqr(Yu(
      2,2))*TCD0(msqInput(2),msqInput(0),msqInput(2),msqInput(2)) + 3*MuInput*
      AuInput(0,2)*Sqr(AuInput(0,1))*Sqr(Yu(0,1))*Sqr(Yu(0,2))*TCD0(msqInput(2),
      msqInput(1),msqInput(0),msqInput(0)) + 3*MuInput*AuInput(1,2)*Sqr(AuInput(1,
      1))*Sqr(Yu(1,1))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(1),msqInput(1),
      msqInput(1)) + 3*MuInput*AuInput(2,2)*Sqr(AuInput(2,1))*Sqr(Yu(2,1))*Sqr(Yu(
      2,2))*TCD0(msqInput(2),msqInput(1),msqInput(2),msqInput(2)) + 3*MuInput*Cube
      (AuInput(0,2))*Quad(Yu(0,2))*TCD0(msqInput(2),msqInput(2),msqInput(0),
      msqInput(0)) + 3*MuInput*AuInput(0,2)*Sqr(AuInput(1,2))*Sqr(Yu(0,2))*Sqr(Yu(
      1,2))*TCD0(msqInput(2),msqInput(2),msqInput(0),msqInput(1)) + 3*MuInput*
      AuInput(0,2)*Sqr(AuInput(2,2))*Sqr(Yu(0,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msqInput(0),msqInput(2)) + 3*MuInput*AuInput(1,2)*Sqr(AuInput(0,
      2))*Sqr(Yu(0,2))*Sqr(Yu(1,2))*TCD0(msqInput(2),msqInput(2),msqInput(1),
      msqInput(0)) + 3*MuInput*Cube(AuInput(1,2))*Quad(Yu(1,2))*TCD0(msqInput(2),
      msqInput(2),msqInput(1),msqInput(1)) + 3*MuInput*AuInput(1,2)*Sqr(AuInput(2,
      2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),msqInput(2),msqInput(1),
      msqInput(2)) + 3*MuInput*AuInput(2,2)*Sqr(AuInput(0,2))*Sqr(Yu(0,2))*Sqr(Yu(
      2,2))*TCD0(msqInput(2),msqInput(2),msqInput(2),msqInput(0)) + 3*MuInput*
      AuInput(2,2)*Sqr(AuInput(1,2))*Sqr(Yu(1,2))*Sqr(Yu(2,2))*TCD0(msqInput(2),
      msqInput(2),msqInput(2),msqInput(1)) + 3*MuInput*Cube(AuInput(2,2))*Quad(Yu(
      2,2))*TCD0(msqInput(2),msqInput(2),msqInput(2),msqInput(2)) + 0.36*M1Input*
      MuInput*Quad(g1)*TCD2t(M1Input,M1Input,Abs(MuInput),Abs(MuInput)) + 0.6*(
      M1Input + M2Input)*MuInput*Sqr(g1)*Sqr(g2)*TCD2t(M1Input,M2Input,Abs(MuInput
      ),Abs(MuInput)) + 3*M2Input*MuInput*Quad(g2)*TCD2t(M2Input,M2Input,Abs(
      MuInput),Abs(MuInput)) + 3*MuInput*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(1),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(
      1,1) + 3*MuInput*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1)
      ,msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*MuInput*
      AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),
      msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*MuInput*AdInput(0,0)*Sqr(
      Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*
      Yd(0,1)*Yd(1,0)*Yd(1,1) + 3*MuInput*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(1),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(
      1,2) + 3*MuInput*AdInput(1,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1)
      ,msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*MuInput*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(0),
      msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*MuInput*AdInput(0,0)*Sqr(
      Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*
      Yd(0,2)*Yd(1,0)*Yd(1,2) + 3*MuInput*AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(1),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(
      1,2) + 3*MuInput*AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(1)
      ,msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*MuInput*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(1),
      msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*MuInput*AdInput(0,1)*Sqr(
      Abs(MuInput))*TCD0(msdInput(1),msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*
      Yd(0,2)*Yd(1,1)*Yd(1,2) + 3*MuInput*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(2),msqInput(0),msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(
      2,1) + 3*MuInput*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2)
      ,msqInput(1),msqInput(0))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*
      AdInput(0,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(0),
      msqInput(1))*Yd(0,0)*Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(0,0)*Sqr(
      Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),msqInput(0))*Yd(0,0)*
      Yd(0,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(2),msqInput(0),msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(
      2,1) + 3*MuInput*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2)
      ,msqInput(1),msqInput(0))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*
      AdInput(1,1)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),
      msqInput(1))*Yd(1,0)*Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(1,0)*Sqr(
      Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(1),msqInput(0))*Yd(1,0)*
      Yd(1,1)*Yd(2,0)*Yd(2,1) + 3*MuInput*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(2),msqInput(0),msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(
      2,2) + 3*MuInput*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2)
      ,msqInput(2),msqInput(0))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(0),
      msqInput(2))*Yd(0,0)*Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput(0,0)*Sqr(
      Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput(0))*Yd(0,0)*
      Yd(0,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(2),msqInput(0),msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(
      2,2) + 3*MuInput*AdInput(2,0)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2)
      ,msqInput(2),msqInput(0))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(0),
      msqInput(2))*Yd(1,0)*Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput(1,0)*Sqr(
      Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(0))*Yd(1,0)*
      Yd(1,2)*Yd(2,0)*Yd(2,2) + 3*MuInput*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(0),msdInput(2),msqInput(1),msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(
      2,2) + 3*MuInput*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(0),msdInput(2)
      ,msqInput(2),msqInput(1))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*
      AdInput(0,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(1),
      msqInput(2))*Yd(0,1)*Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(0,1)*Sqr(
      Abs(MuInput))*TCD0(msdInput(2),msdInput(0),msqInput(2),msqInput(1))*Yd(0,1)*
      Yd(0,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      msdInput(1),msdInput(2),msqInput(1),msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(
      2,2) + 3*MuInput*AdInput(2,1)*Sqr(Abs(MuInput))*TCD0(msdInput(1),msdInput(2)
      ,msqInput(2),msqInput(1))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*
      AdInput(1,2)*Sqr(Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(1),
      msqInput(2))*Yd(1,1)*Yd(1,2)*Yd(2,1)*Yd(2,2) + 3*MuInput*AdInput(1,1)*Sqr(
      Abs(MuInput))*TCD0(msdInput(2),msdInput(1),msqInput(2),msqInput(1))*Yd(1,1)*
      Yd(1,2)*Yd(2,1)*Yd(2,2) + MuInput*AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(1),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(
      1,1) + MuInput*AeInput(1,0)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),
      mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + MuInput*AeInput(0
      ,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(1))*
      Ye(0,0)*Ye(0,1)*Ye(1,0)*Ye(1,1) + MuInput*AeInput(0,0)*Sqr(Abs(MuInput))*
      TCD0(mseInput(1),mseInput(0),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(1,0
      )*Ye(1,1) + MuInput*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput
      (1),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + MuInput*
      AeInput(1,0)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(2),
      mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) + MuInput*AeInput(0,2)*Sqr(Abs(
      MuInput))*TCD0(mseInput(1),mseInput(0),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0
      ,2)*Ye(1,0)*Ye(1,2) + MuInput*AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput(1
      ),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(1,0)*Ye(1,2) +
      MuInput*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput
      (1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + MuInput*AeInput(1,1)*Sqr(
      Abs(MuInput))*TCD0(mseInput(0),mseInput(1),mslInput(2),mslInput(1))*Ye(0,1)*
      Ye(0,2)*Ye(1,1)*Ye(1,2) + MuInput*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(1),mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(
      1,2) + MuInput*AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(0),
      mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(1,1)*Ye(1,2) + MuInput*AeInput(2
      ,1)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(0),mslInput(1))*
      Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(2,0)*Sqr(Abs(MuInput))*
      TCD0(mseInput(0),mseInput(2),mslInput(1),mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0
      )*Ye(2,1) + MuInput*AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput
      (0),mslInput(0),mslInput(1))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + MuInput*
      AeInput(0,0)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(1),
      mslInput(0))*Ye(0,0)*Ye(0,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(2,1)*Sqr(Abs(
      MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(0),mslInput(1))*Ye(1,0)*Ye(1
      ,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput(1
      ),mseInput(2),mslInput(1),mslInput(0))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) +
      MuInput*AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput
      (0),mslInput(1))*Ye(1,0)*Ye(1,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(1,0)*Sqr(
      Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(1),mslInput(0))*Ye(1,0)*
      Ye(1,1)*Ye(2,0)*Ye(2,1) + MuInput*AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(0),mseInput(2),mslInput(0),mslInput(2))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(
      2,2) + MuInput*AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),
      mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(0
      ,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),mslInput(0),mslInput(2))*
      Ye(0,0)*Ye(0,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(0,0)*Sqr(Abs(MuInput))*
      TCD0(mseInput(2),mseInput(0),mslInput(2),mslInput(0))*Ye(0,0)*Ye(0,2)*Ye(2,0
      )*Ye(2,2) + MuInput*AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput
      (2),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + MuInput*
      AeInput(2,0)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(2),
      mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(1,2)*Sqr(Abs(
      MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(0),mslInput(2))*Ye(1,0)*Ye(1
      ,2)*Ye(2,0)*Ye(2,2) + MuInput*AeInput(1,0)*Sqr(Abs(MuInput))*TCD0(mseInput(2
      ),mseInput(1),mslInput(2),mslInput(0))*Ye(1,0)*Ye(1,2)*Ye(2,0)*Ye(2,2) +
      MuInput*AeInput(2,2)*Sqr(Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput
      (1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(2,1)*Sqr(
      Abs(MuInput))*TCD0(mseInput(0),mseInput(2),mslInput(2),mslInput(1))*Ye(0,1)*
      Ye(0,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(0,2)*Sqr(Abs(MuInput))*TCD0(
      mseInput(2),mseInput(0),mslInput(1),mslInput(2))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(
      2,2) + MuInput*AeInput(0,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(0),
      mslInput(2),mslInput(1))*Ye(0,1)*Ye(0,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(2
      ,2)*Sqr(Abs(MuInput))*TCD0(mseInput(1),mseInput(2),mslInput(1),mslInput(2))*
      Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + MuInput*AeInput(2,1)*Sqr(Abs(MuInput))*
      TCD0(mseInput(1),mseInput(2),mslInput(2),mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1
      )*Ye(2,2) + MuInput*AeInput(1,2)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput
      (1),mslInput(1),mslInput(2))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + MuInput*
      AeInput(1,1)*Sqr(Abs(MuInput))*TCD0(mseInput(2),mseInput(1),mslInput(2),
      mslInput(1))*Ye(1,1)*Ye(1,2)*Ye(2,1)*Ye(2,2) + 3*MuInput*AuInput(0,1)*
      AuInput(1,0)*AuInput(1,1)*TCD0(msqInput(0),msqInput(1),msqInput(0),msqInput(
      1))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*MuInput*AuInput(0,0)*AuInput(0,1)*
      AuInput(1,1)*TCD0(msqInput(0),msqInput(1),msqInput(1),msqInput(0))*Yu(0,0)*
      Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*MuInput*AuInput(0,0)*AuInput(1,0)*AuInput(1,1)*
      TCD0(msqInput(1),msqInput(0),msqInput(0),msqInput(1))*Yu(0,0)*Yu(0,1)*Yu(1,0
      )*Yu(1,1) + 3*MuInput*AuInput(0,0)*AuInput(0,1)*AuInput(1,0)*TCD0(msqInput(1
      ),msqInput(0),msqInput(1),msqInput(0))*Yu(0,0)*Yu(0,1)*Yu(1,0)*Yu(1,1) + 3*
      MuInput*AuInput(0,2)*AuInput(1,0)*AuInput(1,2)*TCD0(msqInput(0),msqInput(2),
      msqInput(0),msqInput(1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*MuInput*AuInput
      (0,0)*AuInput(0,2)*AuInput(1,2)*TCD0(msqInput(0),msqInput(2),msqInput(1),
      msqInput(0))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*MuInput*AuInput(0,0)*
      AuInput(1,0)*AuInput(1,2)*TCD0(msqInput(2),msqInput(0),msqInput(0),msqInput(
      1))*Yu(0,0)*Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*MuInput*AuInput(0,0)*AuInput(0,2)*
      AuInput(1,0)*TCD0(msqInput(2),msqInput(0),msqInput(1),msqInput(0))*Yu(0,0)*
      Yu(0,2)*Yu(1,0)*Yu(1,2) + 3*MuInput*AuInput(0,2)*AuInput(1,1)*AuInput(1,2)*
      TCD0(msqInput(1),msqInput(2),msqInput(0),msqInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1
      )*Yu(1,2) + 3*MuInput*AuInput(0,1)*AuInput(0,2)*AuInput(1,2)*TCD0(msqInput(1
      ),msqInput(2),msqInput(1),msqInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*
      MuInput*AuInput(0,1)*AuInput(1,1)*AuInput(1,2)*TCD0(msqInput(2),msqInput(1),
      msqInput(0),msqInput(1))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*MuInput*AuInput
      (0,1)*AuInput(0,2)*AuInput(1,1)*TCD0(msqInput(2),msqInput(1),msqInput(1),
      msqInput(0))*Yu(0,1)*Yu(0,2)*Yu(1,1)*Yu(1,2) + 3*MuInput*AuInput(0,0)*TCC0(
      msqInput(0),msuInput(0),msuInput(1))*Yu(0,0)*Yu(1,0)*(Yu(0,0)*Yu(1,0) + Yu(0
      ,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) + 3*MuInput*AuInput(1,0)*TCC0(msqInput(0),
      msuInput(1),msuInput(0))*Yu(0,0)*Yu(1,0)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1)
      + Yu(0,2)*Yu(1,2)) + 3*MuInput*AuInput(0,1)*TCC0(msqInput(1),msuInput(0),
      msuInput(1))*Yu(0,1)*Yu(1,1)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu
      (1,2)) + 3*MuInput*AuInput(1,1)*TCC0(msqInput(1),msuInput(1),msuInput(0))*Yu
      (0,1)*Yu(1,1)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) + 3*
      MuInput*AuInput(0,2)*TCC0(msqInput(2),msuInput(0),msuInput(1))*Yu(0,2)*Yu(1,
      2)*(Yu(0,0)*Yu(1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) + 3*MuInput*AuInput
      (1,2)*TCC0(msqInput(2),msuInput(1),msuInput(0))*Yu(0,2)*Yu(1,2)*(Yu(0,0)*Yu(
      1,0) + Yu(0,1)*Yu(1,1) + Yu(0,2)*Yu(1,2)) + 3*MuInput*AuInput(0,1)*AuInput(2
      ,0)*AuInput(2,1)*TCD0(msqInput(0),msqInput(1),msqInput(0),msqInput(2))*Yu(0,
      0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(0,0)*AuInput(0,1)*AuInput(2,1
      )*TCD0(msqInput(0),msqInput(1),msqInput(2),msqInput(0))*Yu(0,0)*Yu(0,1)*Yu(2
      ,0)*Yu(2,1) + 3*MuInput*AuInput(0,0)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput
      (1),msqInput(0),msqInput(0),msqInput(2))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3
      *MuInput*AuInput(0,0)*AuInput(0,1)*AuInput(2,0)*TCD0(msqInput(1),msqInput(0)
      ,msqInput(2),msqInput(0))*Yu(0,0)*Yu(0,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*
      AuInput(1,1)*AuInput(2,0)*AuInput(2,1)*TCD0(msqInput(0),msqInput(1),msqInput
      (1),msqInput(2))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(1,0)*
      AuInput(1,1)*AuInput(2,1)*TCD0(msqInput(0),msqInput(1),msqInput(2),msqInput(
      1))*Yu(1,0)*Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(1,0)*AuInput(2,0)*
      AuInput(2,1)*TCD0(msqInput(1),msqInput(0),msqInput(1),msqInput(2))*Yu(1,0)*
      Yu(1,1)*Yu(2,0)*Yu(2,1) + 3*MuInput*AuInput(1,0)*AuInput(1,1)*AuInput(2,0)*
      TCD0(msqInput(1),msqInput(0),msqInput(2),msqInput(1))*Yu(1,0)*Yu(1,1)*Yu(2,0
      )*Yu(2,1) + 3*MuInput*AuInput(0,0)*TCC0(msqInput(0),msqInput(1),msuInput(0))
      *Yu(0,0)*Yu(0,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) + 3*
      MuInput*AuInput(0,1)*TCC0(msqInput(1),msqInput(0),msuInput(0))*Yu(0,0)*Yu(0,
      1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) + 3*MuInput*AuInput
      (1,0)*TCC0(msqInput(0),msqInput(1),msuInput(1))*Yu(1,0)*Yu(1,1)*(Yu(0,0)*Yu(
      0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) + 3*MuInput*AuInput(1,1)*TCC0(
      msqInput(1),msqInput(0),msuInput(1))*Yu(1,0)*Yu(1,1)*(Yu(0,0)*Yu(0,1) + Yu(1
      ,0)*Yu(1,1) + Yu(2,0)*Yu(2,1)) + 3*MuInput*AuInput(2,0)*TCC0(msqInput(0),
      msqInput(1),msuInput(2))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1)
      + Yu(2,0)*Yu(2,1)) + 3*MuInput*AuInput(2,1)*TCC0(msqInput(1),msqInput(0),
      msuInput(2))*Yu(2,0)*Yu(2,1)*(Yu(0,0)*Yu(0,1) + Yu(1,0)*Yu(1,1) + Yu(2,0)*Yu
      (2,1)) + 3*MuInput*AuInput(0,2)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(0),
      msqInput(2),msqInput(0),msqInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*
      MuInput*AuInput(0,0)*AuInput(0,2)*AuInput(2,2)*TCD0(msqInput(0),msqInput(2),
      msqInput(2),msqInput(0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput
      (0,0)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(2),msqInput(0),msqInput(0),
      msqInput(2))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput(0,0)*
      AuInput(0,2)*AuInput(2,0)*TCD0(msqInput(2),msqInput(0),msqInput(2),msqInput(
      0))*Yu(0,0)*Yu(0,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput(1,2)*AuInput(2,0)*
      AuInput(2,2)*TCD0(msqInput(0),msqInput(2),msqInput(1),msqInput(2))*Yu(1,0)*
      Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput(1,0)*AuInput(1,2)*AuInput(2,2)*
      TCD0(msqInput(0),msqInput(2),msqInput(2),msqInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0
      )*Yu(2,2) + 3*MuInput*AuInput(1,0)*AuInput(2,0)*AuInput(2,2)*TCD0(msqInput(2
      ),msqInput(0),msqInput(1),msqInput(2))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*
      MuInput*AuInput(1,0)*AuInput(1,2)*AuInput(2,0)*TCD0(msqInput(2),msqInput(0),
      msqInput(2),msqInput(1))*Yu(1,0)*Yu(1,2)*Yu(2,0)*Yu(2,2) + 3*MuInput*AuInput
      (0,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(1),msqInput(2),msqInput(0),
      msqInput(2))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(0,1)*
      AuInput(0,2)*AuInput(2,2)*TCD0(msqInput(1),msqInput(2),msqInput(2),msqInput(
      0))*Yu(0,1)*Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(0,1)*AuInput(2,1)*
      AuInput(2,2)*TCD0(msqInput(2),msqInput(1),msqInput(0),msqInput(2))*Yu(0,1)*
      Yu(0,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(0,1)*AuInput(0,2)*AuInput(2,1)*
      TCD0(msqInput(2),msqInput(1),msqInput(2),msqInput(0))*Yu(0,1)*Yu(0,2)*Yu(2,1
      )*Yu(2,2) + 3*MuInput*AuInput(1,2)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(1
      ),msqInput(2),msqInput(1),msqInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*
      MuInput*AuInput(1,1)*AuInput(1,2)*AuInput(2,2)*TCD0(msqInput(1),msqInput(2),
      msqInput(2),msqInput(1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput
      (1,1)*AuInput(2,1)*AuInput(2,2)*TCD0(msqInput(2),msqInput(1),msqInput(1),
      msqInput(2))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(1,1)*
      AuInput(1,2)*AuInput(2,1)*TCD0(msqInput(2),msqInput(1),msqInput(2),msqInput(
      1))*Yu(1,1)*Yu(1,2)*Yu(2,1)*Yu(2,2) + 3*MuInput*AuInput(0,0)*TCC0(msqInput(0
      ),msuInput(0),msuInput(2))*Yu(0,0)*Yu(2,0)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1
      ) + Yu(0,2)*Yu(2,2)) + 3*MuInput*AuInput(2,0)*TCC0(msqInput(0),msuInput(2),
      msuInput(0))*Yu(0,0)*Yu(2,0)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu
      (2,2)) + 3*MuInput*AuInput(0,1)*TCC0(msqInput(1),msuInput(0),msuInput(2))*Yu
      (0,1)*Yu(2,1)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) + 3*
      MuInput*AuInput(2,1)*TCC0(msqInput(1),msuInput(2),msuInput(0))*Yu(0,1)*Yu(2,
      1)*(Yu(0,0)*Yu(2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) + 3*MuInput*AuInput
      (0,2)*TCC0(msqInput(2),msuInput(0),msuInput(2))*Yu(0,2)*Yu(2,2)*(Yu(0,0)*Yu(
      2,0) + Yu(0,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) + 3*MuInput*AuInput(2,2)*TCC0(
      msqInput(2),msuInput(2),msuInput(0))*Yu(0,2)*Yu(2,2)*(Yu(0,0)*Yu(2,0) + Yu(0
      ,1)*Yu(2,1) + Yu(0,2)*Yu(2,2)) + 3*MuInput*AuInput(1,0)*TCC0(msqInput(0),
      msuInput(1),msuInput(2))*Yu(1,0)*Yu(2,0)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1)
      + Yu(1,2)*Yu(2,2)) + 3*MuInput*AuInput(2,0)*TCC0(msqInput(0),msuInput(2),
      msuInput(1))*Yu(1,0)*Yu(2,0)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu
      (2,2)) + 3*MuInput*AuInput(1,1)*TCC0(msqInput(1),msuInput(1),msuInput(2))*Yu
      (1,1)*Yu(2,1)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) + 3*
      MuInput*AuInput(2,1)*TCC0(msqInput(1),msuInput(2),msuInput(1))*Yu(1,1)*Yu(2,
      1)*(Yu(1,0)*Yu(2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) + 3*MuInput*AuInput
      (1,2)*TCC0(msqInput(2),msuInput(1),msuInput(2))*Yu(1,2)*Yu(2,2)*(Yu(1,0)*Yu(
      2,0) + Yu(1,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) + 3*MuInput*AuInput(2,2)*TCC0(
      msqInput(2),msuInput(2),msuInput(1))*Yu(1,2)*Yu(2,2)*(Yu(1,0)*Yu(2,0) + Yu(1
      ,1)*Yu(2,1) + Yu(1,2)*Yu(2,2)) + 3*MuInput*AuInput(0,0)*TCC0(msqInput(0),
      msqInput(2),msuInput(0))*Yu(0,0)*Yu(0,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2)
      + Yu(2,0)*Yu(2,2)) + 3*MuInput*AuInput(0,2)*TCC0(msqInput(2),msqInput(0),
      msuInput(0))*Yu(0,0)*Yu(0,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu
      (2,2)) + 3*MuInput*AuInput(1,0)*TCC0(msqInput(0),msqInput(2),msuInput(1))*Yu
      (1,0)*Yu(1,2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) + 3*
      MuInput*AuInput(1,2)*TCC0(msqInput(2),msqInput(0),msuInput(1))*Yu(1,0)*Yu(1,
      2)*(Yu(0,0)*Yu(0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) + 3*MuInput*AuInput
      (2,0)*TCC0(msqInput(0),msqInput(2),msuInput(2))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(
      0,2) + Yu(1,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) + 3*MuInput*AuInput(2,2)*TCC0(
      msqInput(2),msqInput(0),msuInput(2))*Yu(2,0)*Yu(2,2)*(Yu(0,0)*Yu(0,2) + Yu(1
      ,0)*Yu(1,2) + Yu(2,0)*Yu(2,2)) + 3*MuInput*AuInput(0,1)*TCC0(msqInput(1),
      msqInput(2),msuInput(0))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2)
      + Yu(2,1)*Yu(2,2)) + 3*MuInput*AuInput(0,2)*TCC0(msqInput(2),msqInput(1),
      msuInput(0))*Yu(0,1)*Yu(0,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu
      (2,2)) + 3*MuInput*AuInput(1,1)*TCC0(msqInput(1),msqInput(2),msuInput(1))*Yu
      (1,1)*Yu(1,2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) + 3*
      MuInput*AuInput(1,2)*TCC0(msqInput(2),msqInput(1),msuInput(1))*Yu(1,1)*Yu(1,
      2)*(Yu(0,1)*Yu(0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) + 3*MuInput*AuInput
      (2,1)*TCC0(msqInput(1),msqInput(2),msuInput(2))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(
      0,2) + Yu(1,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)) + 3*MuInput*AuInput(2,2)*TCC0(
      msqInput(2),msqInput(1),msuInput(2))*Yu(2,1)*Yu(2,2)*(Yu(0,1)*Yu(0,2) + Yu(1
      ,1)*Yu(1,2) + Yu(2,1)*Yu(2,2)))));


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
