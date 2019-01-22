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

// File generated at Tue 22 Jan 2019 16:37:55

#include "HSSUSY_two_scale_high_scale_constraint.hpp"
#include "HSSUSY_two_scale_model.hpp"
#include "HSSUSY_info.hpp"
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
#define MODELCLASSNAME HSSUSY<Two_scale>

#if defined(ENABLE_HIMALAYA) && Himalaya_VERSION_MAJOR >= 2
#define FSHimalayaMh23L [&] () {                                        \
      MODEL->calculate_DRbar_masses();                                  \
                                                                        \
      himalaya::Parameters pars;                                        \
      const auto ThreeLoopAtAsAs = INPUTPARAMETER(ThreeLoopAtAsAs); \
      const auto DeltaLambda3L = INPUTPARAMETER(DeltaLambda3L); \
      const auto MuInput = INPUTPARAMETER(MuInput); \
      const auto TanBeta = INPUTPARAMETER(TanBeta); \
      const auto msq2 = INPUTPARAMETER(msq2); \
      const auto msd2 = INPUTPARAMETER(msd2); \
      const auto msu2 = INPUTPARAMETER(msu2); \
      const auto msl2 = INPUTPARAMETER(msl2); \
      const auto mse2 = INPUTPARAMETER(mse2); \
      const auto AtInput = INPUTPARAMETER(AtInput); \
      const auto AbInput = INPUTPARAMETER(AbInput); \
      const auto AtauInput = INPUTPARAMETER(AtauInput); \
      const auto M1Input = INPUTPARAMETER(M1Input); \
      const auto M2Input = INPUTPARAMETER(M2Input); \
      const auto M3Input = INPUTPARAMETER(M3Input); \
      const auto mAInput = INPUTPARAMETER(mAInput); \
      const auto g1 = MODELPARAMETER(g1); \
      const auto g2 = MODELPARAMETER(g2); \
      const auto g3 = MODELPARAMETER(g3); \
      const auto v = MODELPARAMETER(v); \
      const auto Yu = MODELPARAMETER(Yu); \
      const auto Yd = MODELPARAMETER(Yd); \
      const auto Ye = MODELPARAMETER(Ye); \
       \
       \
      pars.scale = MODELPARAMETER(scale); \
      pars.mu = Re(MuInput); \
      pars.g1 = Re(g1); \
      pars.g2 = Re(g2); \
      pars.g3 = Re(g3); \
      pars.vd = Re(v/Sqrt(1 + Sqr(TanBeta))); \
      pars.vu = Re((TanBeta*v)/Sqrt(1 + Sqr(TanBeta))); \
      pars.mq2 = Re(msq2); \
      pars.md2 = Re(msd2); \
      pars.mu2 = Re(msu2); \
      pars.ml2 = Re(msl2); \
      pars.me2 = Re(mse2); \
      pars.Au(2,2) = Re(AtInput); \
      pars.Ad(2,2) = Re(AbInput); \
      pars.Ae(2,2) = Re(AtauInput); \
      pars.Yu = Re(Yu); \
      pars.Yd = Re(Yd); \
      pars.Ye = Re(Ye); \
      pars.M1 = M1Input; \
      pars.M2 = M2Input; \
      pars.MG = M3Input; \
      pars.MA = mAInput; \
       \
      const double msbar_scheme = 1; \
      const double lambda_3L_eft = ThreeLoopAtAsAs; \
      const double lambda_3L_uncertainty = DeltaLambda3L; \
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
         model->get_problems().flag_bad_mass(HSSUSY_info::hh); \
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

HSSUSY_high_scale_constraint<Two_scale>::HSSUSY_high_scale_constraint(
   HSSUSY<Two_scale>* model_)
   : model(model_)
{
   initialize();
}

void HSSUSY_high_scale_constraint<Two_scale>::apply()
{
   check_model_ptr();

   


   update_scale();

   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto LambdaLoopOrder = INPUTPARAMETER(LambdaLoopOrder);
   const auto DeltaYt = INPUTPARAMETER(DeltaYt);
   const auto mAInput = INPUTPARAMETER(mAInput);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto M3Input = INPUTPARAMETER(M3Input);
   const auto Qmatch = INPUTPARAMETER(Qmatch);
   const auto MSUSY = INPUTPARAMETER(MSUSY);
   const auto M2Input = INPUTPARAMETER(M2Input);
   const auto M1Input = INPUTPARAMETER(M1Input);
   const auto AtauInput = INPUTPARAMETER(AtauInput);
   const auto AbInput = INPUTPARAMETER(AbInput);
   const auto DeltaOS = INPUTPARAMETER(DeltaOS);
   const auto TwoLoopAtAs = INPUTPARAMETER(TwoLoopAtAs);
   const auto TwoLoopAtAt = INPUTPARAMETER(TwoLoopAtAt);
   const auto TwoLoopAbAs = INPUTPARAMETER(TwoLoopAbAs);
   const auto TwoLoopAtAb = INPUTPARAMETER(TwoLoopAtAb);
   const auto TwoLoopAtauAtau = INPUTPARAMETER(TwoLoopAtauAtau);
   const auto DeltaEFT = INPUTPARAMETER(DeltaEFT);
   const auto msu2 = INPUTPARAMETER(msu2);
   const auto msq2 = INPUTPARAMETER(msq2);
   const auto msd2 = INPUTPARAMETER(msd2);
   const auto mse2 = INPUTPARAMETER(mse2);
   const auto msl2 = INPUTPARAMETER(msl2);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto v = MODELPARAMETER(v);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   MODEL->set_Lambdax(Re(IF(LambdaLoopOrder < 2 && DeltaYt >= 1, -((Power6(Yu(2,2)
      )*(Log(Sqr(mAInput)/Sqr(SCALE))*((-0.03515625*Log(msu2(2,2)/Sqr(SCALE)))/Sqr
      (TanBeta) + (0.005859375*Sqr(MuInput - AtInput*TanBeta)*(-12*Sqrt(msq2(2,2))
      *Sqrt(msu2(2,2))*Sqr(TanBeta)*TCF(1)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))) + Sqr(
      MuInput - AtInput*TanBeta)*TCF(2)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2)))))/(msq2(2
      ,2)*msu2(2,2)*Power6(TanBeta))) + Log(Sqr(MuInput)/Sqr(SCALE))*((-0.03515625
      *Log(msu2(2,2)/Sqr(SCALE))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) + (0.005859375*(
      1 + Sqr(TanBeta))*Sqr(MuInput - AtInput*TanBeta)*(-12*Sqrt(msq2(2,2))*Sqrt(
      msu2(2,2))*Sqr(TanBeta)*TCF(1)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))) + Sqr(
      MuInput - AtInput*TanBeta)*TCF(2)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2)))))/(msq2(2
      ,2)*msu2(2,2)*Power6(TanBeta))) + 0.005859375*Log(msu2(2,2)/Sqr(SCALE))*(3/
      Sqr(TanBeta) + (2*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt
      (msu2(2,2))))/(Sqrt(msq2(2,2))*Sqrt(msu2(2,2))) - (8*(1 + Sqr(TanBeta))*TCF(
      6)(Sqrt(msq2(2,2))/Sqrt(Sqr(MuInput))))/Sqr(TanBeta) - (4*(1 + Sqr(TanBeta))
      *TCF(6)(Sqrt(msu2(2,2))/Sqrt(Sqr(MuInput))))/Sqr(TanBeta)) - (0.0009765625*
      Sqr(MuInput - AtInput*TanBeta)*(-12*Sqrt(msq2(2,2))*Sqrt(msu2(2,2))*Sqr(
      TanBeta)*TCF(1)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))) + Sqr(MuInput - AtInput*
      TanBeta)*TCF(2)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))*(2*Sqr(MuInput - AtInput*
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))) + Sqrt(msq2(2,2))*Sqrt(msu2
      (2,2))*(3 - 8*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/Sqrt(Sqr(MuInput)))
      - 4*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/Sqrt(Sqr(MuInput))))))/(Power3
      (Sqrt(msq2(2,2)))*Power3(Sqrt(msu2(2,2)))*Power6(TanBeta)) + Log(msq2(2,2)/
      Sqr(SCALE))*((-0.03515625*Log(Sqr(mAInput)/Sqr(SCALE)))/Sqr(TanBeta) - (
      0.03515625*Log(Sqr(MuInput)/Sqr(SCALE))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) +
      0.005859375*(3/Sqr(TanBeta) + (2*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(
      msq2(2,2))/Sqrt(msu2(2,2))))/(Sqrt(msq2(2,2))*Sqrt(msu2(2,2))) - (8*(1 + Sqr
      (TanBeta))*TCF(6)(Sqrt(msq2(2,2))/Sqrt(Sqr(MuInput))))/Sqr(TanBeta) - (4*(1
      + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/Sqrt(Sqr(MuInput))))/Sqr(TanBeta)))))
      /Quad(3.141592653589793)) - (Quad(Yu(2,2))*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(
      SCALE))*(-0.0625*Log(msu2(2,2)/Sqr(SCALE)) + (0.010416666666666666*Sqr(
      MuInput - AtInput*TanBeta)*(-12*Sqrt(msq2(2,2))*Sqrt(msu2(2,2))*Sqr(TanBeta)
      *TCF(1)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))) + Sqr(MuInput - AtInput*TanBeta)*
      TCF(2)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2)))))/(msq2(2,2)*msu2(2,2)*Quad(TanBeta)
      )) + (0.010416666666666666*Sqr(MuInput - AtInput*TanBeta)*(-12*Sqrt(msq2(2,2
      ))*Sqrt(msu2(2,2))*Sqr(TanBeta)*TCF(1)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))) +
      Sqr(MuInput - AtInput*TanBeta)*TCF(2)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))*(
      TanBeta*Sqrt(Sqr(M3Input)) + TanBeta*Sqrt(Sqr(M3Input))*TCF(6)(Sqrt(msq2(2,2
      ))/Sqrt(Sqr(M3Input))) + TanBeta*Sqrt(Sqr(M3Input))*TCF(6)(Sqrt(msu2(2,2))/
      Sqrt(Sqr(M3Input))) + MuInput*TCF(9)(Sqrt(msq2(2,2))/Sqrt(Sqr(M3Input)),Sqrt
      (msu2(2,2))/Sqrt(Sqr(M3Input))) - AtInput*TanBeta*TCF(9)(Sqrt(msq2(2,2))/
      Sqrt(Sqr(M3Input)),Sqrt(msu2(2,2))/Sqrt(Sqr(M3Input)))))/(msq2(2,2)*msu2(2,2
      )*Power5(TanBeta)*Sqrt(Sqr(M3Input))) - (0.0625*Log(msu2(2,2)/Sqr(SCALE))*(
      TanBeta + TanBeta*TCF(6)(Sqrt(msq2(2,2))/Sqrt(Sqr(M3Input))) + TanBeta*TCF(6
      )(Sqrt(msu2(2,2))/Sqrt(Sqr(M3Input))) + (MuInput*TCF(9)(Sqrt(msq2(2,2))/Sqrt
      (Sqr(M3Input)),Sqrt(msu2(2,2))/Sqrt(Sqr(M3Input))))/Sqrt(Sqr(M3Input)) - (
      AtInput*TanBeta*TCF(9)(Sqrt(msq2(2,2))/Sqrt(Sqr(M3Input)),Sqrt(msu2(2,2))/
      Sqrt(Sqr(M3Input))))/Sqrt(Sqr(M3Input))))/TanBeta + Log(msq2(2,2)/Sqr(SCALE)
      )*(-0.0625*Log(Sqr(M3Input)/Sqr(SCALE)) - (0.0625*(TanBeta + TanBeta*TCF(6)(
      Sqrt(msq2(2,2))/Sqrt(Sqr(M3Input))) + TanBeta*TCF(6)(Sqrt(msu2(2,2))/Sqrt(
      Sqr(M3Input))) + (MuInput*TCF(9)(Sqrt(msq2(2,2))/Sqrt(Sqr(M3Input)),Sqrt(
      msu2(2,2))/Sqrt(Sqr(M3Input))))/Sqrt(Sqr(M3Input)) - (AtInput*TanBeta*TCF(9)
      (Sqrt(msq2(2,2))/Sqrt(Sqr(M3Input)),Sqrt(msu2(2,2))/Sqrt(Sqr(M3Input))))/
      Sqrt(Sqr(M3Input))))/TanBeta)))/Quad(3.141592653589793), 0) + IF(
      LambdaLoopOrder < 3 && DeltaYt >= 1, WHICH(IsCloseRel(Sqr(SCALE),msq2(2,2),
      0.01) && IsCloseRel(Sqr(SCALE),msu2(2,2),0.01) && IsCloseRel(SCALE,Abs(
      M3Input),0.01), (-0.000244140625*Quad(g3)*Quad(Yu(2,2))*((2*(-
      38.425925925925924 + Log(msq2(2,2)/Sqr(SCALE))*(-24.814814814814813 - (
      7.703703703703703*(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2))) + (
      13.185185185185185*(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2)) + (
      1.7777777777777777*Sqr(AtInput - MuInput/TanBeta))/msq2(2,2) +
      0.2222222222222222*Sqr(Log(msq2(2,2)/Sqr(SCALE)))) + 3*Sqr(-
      1.3333333333333333 - 1.3333333333333333*Log(msq2(2,2)/Sqr(SCALE)) + (
      1.3333333333333333*(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2))))*(12*Log(
      msq2(2,2)/Sqr(SCALE)) + (12*Sqr(AtInput - MuInput/TanBeta))/msq2(2,2) - Quad
      (AtInput - MuInput/TanBeta)/Sqr(msq2(2,2))) + (0.5 - 2*Log(msq2(2,2)/Sqr(
      SCALE)) + 2*(-1.3333333333333333 - 1.3333333333333333*Log(msq2(2,2)/Sqr(
      SCALE)) + (1.3333333333333333*(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2))))
      *((74.66666666666667*Cube(AtInput - MuInput/TanBeta))/Power3(Sqrt(msq2(2,2))
      ) - (64*(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2)) - (5.333333333333333*
      Power5(AtInput - MuInput/TanBeta))/Power(msq2(2,2),2.5) - (32*Sqr(AtInput -
      MuInput/TanBeta))/msq2(2,2) - 96*Sqr(Log(msq2(2,2)/Sqr(SCALE))) + Log(msq2(2
      ,2)/Sqr(SCALE))*((-21.333333333333332*Cube(AtInput - MuInput/TanBeta))/
      Power3(Sqrt(msq2(2,2))) + (128*(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2))
      - (128*Sqr(AtInput - MuInput/TanBeta))/msq2(2,2) + (5.333333333333333*Quad(
      AtInput - MuInput/TanBeta))/Sqr(msq2(2,2))) + (2.6666666666666665*Quad(
      AtInput - MuInput/TanBeta))/Sqr(msq2(2,2)))))/Power6(3.141592653589793),
      True, (-0.000244140625*Quad(g3)*Quad(Yu(2,2))*((0.5 + 0.08333333333333333*(-
      Log((1.02*msq2(2,2))/Sqr(SCALE)) - 10*Log(Power(msd2(0,0)*msd2(1,1)*msq2(0,0
      )*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)/Sqr(SCALE)) - Log((
      0.9800000000000001*msu2(2,2))/Sqr(SCALE))) - Log(Sqr(M3Input)/Sqr(SCALE)) +
      2*((-0.6664000000000001*msq2(2,2)*msu2(2,2))/((-1.02*msq2(2,2) + Sqr(M3Input
      ))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (0.6666666666666666*
      Quad(M3Input))/((-1.02*msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2
      ,2) + Sqr(M3Input))) + Log((1.02*msq2(2,2))/Sqr(SCALE))*((-
      2.7199999999999998*M3Input*(AtInput - MuInput/TanBeta)*msq2(2,2))/((1.02*
      msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input)))
      + (1.3599999999999999*msq2(2,2)*Sqr(M3Input))/Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input)) - (0.6936*Sqr(msq2(2,2)))/Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) +
      Log((0.9800000000000001*msu2(2,2))/Sqr(SCALE))*((-2.6133333333333333*M3Input
      *(AtInput - MuInput/TanBeta)*msu2(2,2))/((-1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))
      ) + (1.3066666666666666*msu2(2,2)*Sqr(M3Input))/Sqr(-0.9800000000000001*msu2
      (2,2) + Sqr(M3Input)) - (0.6402666666666668*Sqr(msu2(2,2)))/Sqr(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + Log(Sqr(M3Input)/Sqr(SCALE))
      *((2.6666666666666665*(AtInput - MuInput/TanBeta)*Cube(M3Input))/((-1.02*
      msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) -
      (0.6666666666666666*(-2.04*msq2(2,2)*Power6(M3Input) - 1.9600000000000002*
      msu2(2,2)*Power6(M3Input) + 2*Power8(M3Input) + 1.0404*Quad(M3Input)*Sqr(
      msq2(2,2)) + 0.9604000000000001*Quad(M3Input)*Sqr(msu2(2,2))))/(Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))
      )))*((-64.02561024409763*Sqr(M3Input)*Sqr(AtInput - MuInput/TanBeta))/(msq2(
      2,2)*msu2(2,2)) - (512*M3Input*Cube(AtInput - MuInput/TanBeta))/Sqr(1.02*
      msq2(2,2) - 0.9800000000000001*msu2(2,2)) + (32.01280512204882*Quad(AtInput
      - MuInput/TanBeta)*(-1.02*msq2(2,2)*Power6(M3Input) - 0.9800000000000001*
      msu2(2,2)*Power6(M3Input) - 4.998*msq2(2,2)*msu2(2,2)*Quad(M3Input) + 1.0404
      *Quad(M3Input)*Sqr(msq2(2,2)) + 5.0979600000000005*msu2(2,2)*Sqr(M3Input)*
      Sqr(msq2(2,2)) + 0.9604000000000001*Quad(M3Input)*Sqr(msu2(2,2)) +
      4.898040000000001*msq2(2,2)*Sqr(M3Input)*Sqr(msu2(2,2)) - 4.996000800000001*
      Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(msq2(2,2)*msu2(2,2)*(1.02*msq2(2,2) - Sqr(
      M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))*Sqr(1.02*msq2(2,2)
      - 0.9800000000000001*msu2(2,2))) + (16.00640256102441*(-2.04*msq2(2,2)*
      Power6(M3Input) - 1.9600000000000002*msu2(2,2)*Power6(M3Input) +
      8.996400000000001*msq2(2,2)*msu2(2,2)*Quad(M3Input) + 2.0808*Quad(M3Input)*
      Sqr(msq2(2,2)) - 6.117552000000001*msu2(2,2)*Sqr(M3Input)*Sqr(msq2(2,2)) +
      1.9208000000000003*Quad(M3Input)*Sqr(msu2(2,2)) - 5.8776480000000015*msq2(2,
      2)*Sqr(M3Input)*Sqr(msu2(2,2)) + 2.9976004800000005*Sqr(msq2(2,2))*Sqr(msu2(
      2,2))))/(msq2(2,2)*msu2(2,2)*(1.02*msq2(2,2) - Sqr(M3Input))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + PolyLog(2,1 - (
      0.9803921568627451*Sqr(M3Input))/msq2(2,2))*((-32*Quad(AtInput - MuInput/
      TanBeta)*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2) - 2*Sqr(M3Input)))/
      Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) + (256*(AtInput -
      MuInput/TanBeta)*Cube(M3Input))/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2
      ,2))*(-1.02*msq2(2,2) + Sqr(M3Input))) + (128*M3Input*Cube(AtInput - MuInput
      /TanBeta)*(-1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2) + 2*Sqr(M3Input)))
      /Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) - (64*Quad(M3Input))/
      Sqr(1.02*msq2(2,2) - Sqr(M3Input))) + Sqr(Log((1.02*msq2(2,2))/Sqr(SCALE)))*
      ((64*(AtInput - MuInput/TanBeta)*(2*Cube(M3Input) - 1.02*M3Input*msq2(2,2)))
      /((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(
      M3Input))) + (65.28*M3Input*msq2(2,2)*(1.02*msq2(2,2) + 0.9800000000000001*
      msu2(2,2))*Power5(AtInput - MuInput/TanBeta))/(Quad(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))) + (64*
      M3Input*Cube(AtInput - MuInput/TanBeta)*(2.9988000000000006*msq2(2,2)*msu2(2
      ,2) + 2*Quad(M3Input) - (3.06*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Sqr(
      M3Input) - 1.0404*Sqr(msq2(2,2))))/(Cube(1.02*msq2(2,2) - 0.9800000000000001
      *msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))) - (16*(3*Quad(M3Input) - 4.08*
      msq2(2,2)*Sqr(M3Input) + 2.0808*Sqr(msq2(2,2))))/Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input)) - (32*Sqr(AtInput - MuInput/TanBeta)*(4.244832000000001*Cube(msq2(
      2,2)) + (3.06*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Quad(M3Input) + Sqr(
      M3Input)*(3.9984000000000006*msq2(2,2)*msu2(2,2) - 8.3232*Sqr(msq2(2,2))) -
      2.039184*msu2(2,2)*Sqr(msq2(2,2))))/(Sqr(1.02*msq2(2,2) - 0.9800000000000001
      *msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) + (16*Quad(AtInput -
      MuInput/TanBeta)*(4.159935360000001*Cube(msq2(2,2))*msu2(2,2) + 2*(1.02*msq2
      (2,2) - 0.9800000000000001*msu2(2,2))*Power6(M3Input) + 7.996800000000001*
      msq2(2,2)*msu2(2,2)*Quad(M3Input) + 5.4121608*Quad(msq2(2,2)) -
      0.9992001600000001*Sqr(msq2(2,2))*Sqr(msu2(2,2)) - 2*Sqr(M3Input)*(
      4.244832000000001*Cube(msq2(2,2)) + 5.0979600000000005*msu2(2,2)*Sqr(msq2(2,
      2)) - 0.9796080000000001*msq2(2,2)*Sqr(msu2(2,2)))))/(Quad(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input)))) - 4*(Log
      ((0.9800000000000001*msu2(2,2))/Sqr(SCALE))*(6 + ((6.12*msq2(2,2))/Cube(1.02
      *msq2(2,2) - 0.9800000000000001*msu2(2,2)) + (5.880000000000001*msu2(2,2))/
      Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)))*Quad(AtInput - MuInput/
      TanBeta) - (12*Sqr(AtInput - MuInput/TanBeta))/(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))) + Log((1.02*msq2(2,2))/Sqr(SCALE))*(6 + ((-
      6.12*msq2(2,2))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) - (
      5.880000000000001*msu2(2,2))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2
      ,2)))*Quad(AtInput - MuInput/TanBeta) + (12*Sqr(AtInput - MuInput/TanBeta))/
      (1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))) + (12*Quad(AtInput -
      MuInput/TanBeta))/Sqr(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)))*((-
      0.6664000000000001*msq2(2,2)*msu2(2,2))/((-1.02*msq2(2,2) + Sqr(M3Input))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (0.6666666666666666*Quad(
      M3Input))/((-1.02*msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input))) + Log((1.02*msq2(2,2))/Sqr(SCALE))*((-2.7199999999999998*
      M3Input*(AtInput - MuInput/TanBeta)*msq2(2,2))/((1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))) + (
      1.3599999999999999*msq2(2,2)*Sqr(M3Input))/Sqr(-1.02*msq2(2,2) + Sqr(M3Input
      )) - (0.6936*Sqr(msq2(2,2)))/Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) + Log((
      0.9800000000000001*msu2(2,2))/Sqr(SCALE))*((-2.6133333333333333*M3Input*(
      AtInput - MuInput/TanBeta)*msu2(2,2))/((-1.02*msq2(2,2) + 0.9800000000000001
      *msu2(2,2))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (
      1.3066666666666666*msu2(2,2)*Sqr(M3Input))/Sqr(-0.9800000000000001*msu2(2,2)
      + Sqr(M3Input)) - (0.6402666666666668*Sqr(msu2(2,2)))/Sqr(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + Log(Sqr(M3Input)/Sqr(SCALE))
      *((2.6666666666666665*(AtInput - MuInput/TanBeta)*Cube(M3Input))/((-1.02*
      msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) -
      (0.6666666666666666*(-2.04*msq2(2,2)*Power6(M3Input) - 1.9600000000000002*
      msu2(2,2)*Power6(M3Input) + 2*Power8(M3Input) + 1.0404*Quad(M3Input)*Sqr(
      msq2(2,2)) + 0.9604000000000001*Quad(M3Input)*Sqr(msu2(2,2))))/(Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))
      )) + Log((1.02*msq2(2,2))/Sqr(SCALE))*((-128*M3Input*(AtInput - MuInput/
      TanBeta))/(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) + (128*M3Input*
      Cube(AtInput - MuInput/TanBeta)*(3.06*msq2(2,2) + 0.9800000000000001*msu2(2,
      2)))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) - (130.56*M3Input*
      msq2(2,2)*Power5(AtInput - MuInput/TanBeta))/(Cube(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))) + (32*(
      2.9988000000000006*msq2(2,2)*msu2(2,2) + 7*Quad(M3Input) - 2*(2.04*msq2(2,2)
      + 2.9400000000000004*msu2(2,2))*Sqr(M3Input))*Sqr(AtInput - MuInput/TanBeta)
      )/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(
      M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (16*(7*Power6(
      M3Input) - (7.140000000000001*msq2(2,2) + 5.880000000000001*msu2(2,2))*Quad(
      M3Input) - 2.039184*msu2(2,2)*Sqr(msq2(2,2)) + Sqr(M3Input)*(
      4.998000000000001*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)))))/((-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))*Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input))) - (16*Quad(AtInput - MuInput/TanBeta)*((7.140000000000001*msq2(2,
      2) + 2.9400000000000004*msu2(2,2))*Power6(M3Input) + 4*Power8(M3Input) -
      3.0587760000000004*(5.1*msq2(2,2) + 0.9800000000000001*msu2(2,2))*msu2(2,2)*
      Sqr(msq2(2,2)) - Quad(M3Input)*(14.994000000000002*msq2(2,2)*msu2(2,2) +
      30.171599999999998*Sqr(msq2(2,2)) + 5.762400000000001*Sqr(msu2(2,2))) + Sqr(
      M3Input)*(16.979328000000002*Cube(msq2(2,2)) + 31.607352000000002*msu2(2,2)*
      Sqr(msq2(2,2)) + 6.857256000000001*msq2(2,2)*Sqr(msu2(2,2)))))/(Cube(1.02*
      msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) + Log((0.9800000000000001
      *msu2(2,2))/Sqr(SCALE))*((-64*(AtInput - MuInput/TanBeta)*Cube(M3Input))/((-
      1.02*msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)
      )) - (64*M3Input*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Power5(
      AtInput - MuInput/TanBeta)*(-1.9992000000000003*msq2(2,2)*msu2(2,2) + (1.02*
      msq2(2,2) + 0.9800000000000001*msu2(2,2))*Sqr(M3Input)))/(Quad(1.02*msq2(2,2
      ) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (128*Cube(AtInput - MuInput/
      TanBeta)*(Cube(M3Input)*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2)) -
      1.9992000000000003*M3Input*msq2(2,2)*msu2(2,2)))/((-1.02*msq2(2,2) + Sqr(
      M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))*Sqr(1.02*msq2(2,2)
      - 0.9800000000000001*msu2(2,2))) + (16*(2*(1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*Power6(M3Input) + 3.9984000000000006*msq2(2,2)
      *(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*msu2(2,2)*Sqr(M3Input) -
      Quad(M3Input)*(7.996800000000001*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2))
      + 0.9604000000000001*Sqr(msu2(2,2))) - 1.9984003200000002*Sqr(msq2(2,2))*Sqr
      (msu2(2,2))))/(Sqr(-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*
      msu2(2,2) + Sqr(M3Input))) + (32*Sqr(AtInput - MuInput/TanBeta)*(2*(1.02*
      msq2(2,2) + 0.9800000000000001*msu2(2,2))*Power8(M3Input) + 3*Cube(1.02*msq2
      (2,2) + 0.9800000000000001*msu2(2,2))*Quad(M3Input) - 3.9984000000000006*
      msq2(2,2)*msu2(2,2)*Sqr(M3Input)*Sqr(1.02*msq2(2,2) + 0.9800000000000001*
      msu2(2,2)) + 1.9984003200000002*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,
      2))*Sqr(msq2(2,2))*Sqr(msu2(2,2)) - 2*Power6(M3Input)*(1.9992000000000003*
      msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) + 2.8812000000000006*Sqr(msu2(2,
      2)))))/(Sqr(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,
      2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (16*
      (1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Quad(AtInput - MuInput/
      TanBeta)*(4*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Power8(M3Input)
      - 7.996800000000001*msq2(2,2)*msu2(2,2)*Sqr(M3Input)*Sqr(1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2)) + 3.9968006400000005*(1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*Sqr(msq2(2,2))*Sqr(msu2(2,2)) - 2*Power6(
      M3Input)*(5.997600000000001*msq2(2,2)*msu2(2,2) + 5.202*Sqr(msq2(2,2)) +
      4.8020000000000005*Sqr(msu2(2,2))) + Quad(M3Input)*(5.306040000000001*Cube(
      msq2(2,2)) + 4.705960000000001*Cube(msu2(2,2)) + 19.372248000000003*msu2(2,2
      )*Sqr(msq2(2,2)) + 18.612552000000004*msq2(2,2)*Sqr(msu2(2,2)))))/(Quad(1.02
      *msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input
      ))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))))) + PolyLog(2,1 - (
      1.020408163265306*Sqr(M3Input))/msu2(2,2))*((128*M3Input*Cube(AtInput -
      MuInput/TanBeta)*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2) - 2*Sqr(
      M3Input)))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) + (32*Quad(
      AtInput - MuInput/TanBeta)*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2) -
      2*Sqr(M3Input)))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) + (256*
      (AtInput - MuInput/TanBeta)*Cube(M3Input))/((-1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))
      ) - (64*Quad(M3Input))/Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) +
      Sqr(Log((0.9800000000000001*msu2(2,2))/Sqr(SCALE)))*((64*(AtInput - MuInput/
      TanBeta)*(2*Cube(M3Input) - 0.9800000000000001*M3Input*msu2(2,2)))/((-1.02*
      msq2(2,2) + 0.9800000000000001*msu2(2,2))*(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input))) + (62.720000000000006*M3Input*(1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*msu2(2,2)*Power5(AtInput - MuInput/TanBeta))/(
      Quad(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-0.9800000000000001*
      msu2(2,2) + Sqr(M3Input))) + (64*M3Input*Cube(AtInput - MuInput/TanBeta)*(
      2.9988000000000006*msq2(2,2)*msu2(2,2) + 2*Quad(M3Input) - (1.02*msq2(2,2) +
      2.9400000000000004*msu2(2,2))*Sqr(M3Input) - 0.9604000000000001*Sqr(msu2(2,2
      ))))/(Cube(-1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (16*(3*Quad(M3Input) -
      3.9200000000000004*msu2(2,2)*Sqr(M3Input) + 1.9208000000000003*Sqr(msu2(2,2)
      )))/Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)) + (16*Quad(AtInput -
      MuInput/TanBeta)*(3.840063360000001*Cube(msu2(2,2))*msq2(2,2) - 2*(1.02*msq2
      (2,2) - 0.9800000000000001*msu2(2,2))*Power6(M3Input) + 7.996800000000001*
      msq2(2,2)*msu2(2,2)*Quad(M3Input) + 4.611840800000001*Quad(msu2(2,2)) +
      1.9600000000000002*msu2(2,2)*Sqr(M3Input)*(-4.998*msq2(2,2)*msu2(2,2) +
      1.0404*Sqr(msq2(2,2)) - 3.8416000000000006*Sqr(msu2(2,2))) -
      0.9992001600000001*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Quad(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input))) + (32*Sqr(AtInput - MuInput/TanBeta)*((1.02*msq2(2,2) -
      2.9400000000000004*msu2(2,2))*Quad(M3Input) + 1.9208000000000003*(1.02*msq2(
      2,2) - 1.9600000000000002*msu2(2,2))*Sqr(msu2(2,2)) + Sqr(M3Input)*(-
      3.9984000000000006*msq2(2,2)*msu2(2,2) + 7.683200000000001*Sqr(msu2(2,2)))))
      /(Sqr(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001
      *msu2(2,2) + Sqr(M3Input)))) + Log((0.9800000000000001*msu2(2,2))/Sqr(SCALE)
      )*((128*M3Input*(AtInput - MuInput/TanBeta))/(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2)) - (128*M3Input*Cube(AtInput - MuInput/TanBeta)
      *(1.02*msq2(2,2) + 2.9400000000000004*msu2(2,2)))/Cube(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2)) - (125.44000000000001*M3Input*msu2(2,2)*Power5
      (AtInput - MuInput/TanBeta))/(Cube(-1.02*msq2(2,2) + 0.9800000000000001*msu2
      (2,2))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (32*(
      2.9988000000000006*msq2(2,2)*msu2(2,2) + 7*Quad(M3Input) - 2*(3.06*msq2(2,2)
      + 1.9600000000000002*msu2(2,2))*Sqr(M3Input))*Sqr(AtInput - MuInput/TanBeta)
      )/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(
      M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (16*(7*Power6(
      M3Input) - (6.12*msq2(2,2) + 6.86*msu2(2,2))*Quad(M3Input) -
      1.9592160000000003*msq2(2,2)*Sqr(msu2(2,2)) + Sqr(M3Input)*(
      4.998000000000001*msq2(2,2)*msu2(2,2) + 2.8812000000000006*Sqr(msu2(2,2)))))
      /((-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input))) + (16*Quad(AtInput - MuInput/TanBeta)*((3.06*msq2(2,2) + 6.86*
      msu2(2,2))*Power6(M3Input) + 4*Power8(M3Input) - 2.9388240000000008*msq2(2,2
      )*(1.02*msq2(2,2) + 4.9*msu2(2,2))*Sqr(msu2(2,2)) - Quad(M3Input)*(
      14.994000000000002*msq2(2,2)*msu2(2,2) + 6.2424*Sqr(msq2(2,2)) +
      27.851600000000005*Sqr(msu2(2,2))) + Sqr(M3Input)*(15.059072000000004*Cube(
      msu2(2,2)) + 7.137144*msu2(2,2)*Sqr(msq2(2,2)) + 30.367848000000006*msq2(2,2
      )*Sqr(msu2(2,2)))))/(Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-
      1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input)))) + Log(Sqr(M3Input)/Sqr(SCALE))*((64.02561024409763*Quad(M3Input)
      *(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2) - Sqr(M3Input))*Sqr(AtInput
      - MuInput/TanBeta))/(msq2(2,2)*msu2(2,2)*(1.02*msq2(2,2) - Sqr(M3Input))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (128*Cube(M3Input)*Power5(
      AtInput - MuInput/TanBeta))/((-1.02*msq2(2,2) + Sqr(M3Input))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))*Sqr(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))) + Log((1.02*msq2(2,2))/Sqr(SCALE))*((-64*Cube
      (M3Input)*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Power5(AtInput -
      MuInput/TanBeta))/(Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-
      1.02*msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)
      )) + (64*(AtInput - MuInput/TanBeta)*Cube(M3Input)*(1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2) - 2*Sqr(M3Input)))/((1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (128*Cube(M3Input)*Cube(
      AtInput - MuInput/TanBeta)*(3.9984000000000006*msq2(2,2)*msu2(2,2) + 2*Quad(
      M3Input) - 2*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Sqr(M3Input) -
      1.0404*Sqr(msq2(2,2)) - 0.9604000000000001*Sqr(msu2(2,2))))/(Cube(1.02*msq2(
      2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (16*(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Quad(M3Input)*(-1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2) + 2*Sqr(M3Input)))/(Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (16*(1.02*
      msq2(2,2) + 0.9800000000000001*msu2(2,2))*Quad(M3Input)*Quad(AtInput -
      MuInput/TanBeta)*(2*Quad(M3Input) - 2*(1.02*msq2(2,2) + 0.9800000000000001*
      msu2(2,2))*Sqr(M3Input) + 1.0404*Sqr(msq2(2,2)) + 0.9604000000000001*Sqr(
      msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))
      - (32*Quad(M3Input)*Sqr(AtInput - MuInput/TanBeta)*(2*Quad(M3Input) - 2*(
      1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Sqr(M3Input) + 1.0404*Sqr(
      msq2(2,2)) + 0.9604000000000001*Sqr(msu2(2,2))))/((1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + Log((0.9800000000000001*
      msu2(2,2))/Sqr(SCALE))*((64*Cube(M3Input)*(1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*Power5(AtInput - MuInput/TanBeta))/(Cube(1.02*
      msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))*(
      -0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (64*(AtInput - MuInput/
      TanBeta)*Cube(M3Input)*(-1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2) + 2*
      Sqr(M3Input)))/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(
      2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (128*
      Cube(M3Input)*Cube(AtInput - MuInput/TanBeta)*(3.9984000000000006*msq2(2,2)*
      msu2(2,2) + 2*Quad(M3Input) - 2*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,
      2))*Sqr(M3Input) - 1.0404*Sqr(msq2(2,2)) - 0.9604000000000001*Sqr(msu2(2,2))
      ))/(Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) +
      Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (16*(1.02*
      msq2(2,2) - 0.9800000000000001*msu2(2,2))*Quad(M3Input)*(-1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2) + 2*Sqr(M3Input)))/(Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (16*(1.02*
      msq2(2,2) + 0.9800000000000001*msu2(2,2))*Quad(M3Input)*Quad(AtInput -
      MuInput/TanBeta)*(2*Quad(M3Input) - 2*(1.02*msq2(2,2) + 0.9800000000000001*
      msu2(2,2))*Sqr(M3Input) + 1.0404*Sqr(msq2(2,2)) + 0.9604000000000001*Sqr(
      msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))
      + (32*Quad(M3Input)*Sqr(AtInput - MuInput/TanBeta)*(2*Quad(M3Input) - 2*(
      1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Sqr(M3Input) + 1.0404*Sqr(
      msq2(2,2)) + 0.9604000000000001*Sqr(msu2(2,2))))/((1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) - (16.00640256102441*Quad(
      M3Input)*(2*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Power6(M3Input)
      + 2*(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Sqr(M3Input)*Sqr(1.02*
      msq2(2,2) - 0.9800000000000001*msu2(2,2)) + Quad(M3Input)*(
      1.9992000000000003*msq2(2,2)*msu2(2,2) - 4.1616*Sqr(msq2(2,2)) -
      3.8416000000000006*Sqr(msu2(2,2))) + 0.9996000000000002*msq2(2,2)*msu2(2,2)*
      (1.0404*Sqr(msq2(2,2)) + 0.9604000000000001*Sqr(msu2(2,2)))))/(msq2(2,2)*
      msu2(2,2)*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2
      ,2) + Sqr(M3Input))) - (32.01280512204882*Quad(AtInput - MuInput/TanBeta)*
      Sqr(M3Input)*((1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Power8(M3Input
      ) + Cube(1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Quad(M3Input) - 2*
      Power6(M3Input)*(0.9996000000000002*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,
      2)) + 0.9604000000000001*Sqr(msu2(2,2))) - 0.9996000000000002*msq2(2,2)*msu2
      (2,2)*Sqr(M3Input)*(3.9984000000000006*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2
      (2,2)) + 0.9604000000000001*Sqr(msu2(2,2))) + 0.9992001600000001*(1.02*msq2(
      2,2) + 0.9800000000000001*msu2(2,2))*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(msq2(2
      ,2)*msu2(2,2)*Sqr(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))
      )) + (Log((0.9800000000000001*msu2(2,2))/Sqr(SCALE))*(6 + ((6.12*msq2(2,2))/
      Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) + (5.880000000000001*
      msu2(2,2))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)))*Quad(AtInput
       - MuInput/TanBeta) - (12*Sqr(AtInput - MuInput/TanBeta))/(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))) + Log((1.02*msq2(2,2))/Sqr(SCALE))*(6 + ((-
      6.12*msq2(2,2))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)) - (
      5.880000000000001*msu2(2,2))/Cube(1.02*msq2(2,2) - 0.9800000000000001*msu2(2
      ,2)))*Quad(AtInput - MuInput/TanBeta) + (12*Sqr(AtInput - MuInput/TanBeta))/
      (1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))) + (12*Quad(AtInput -
      MuInput/TanBeta))/Sqr(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)))*(2*((-
      0.8888888888888888*(AtInput - MuInput/TanBeta)*(14*Cube(M3Input) - 3.06*
      M3Input*msq2(2,2) - 30*M3Input*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)
      *msu2(0,0)*msu2(1,1),0.16666666666666666) - 2.9400000000000004*M3Input*msu2(
      2,2)))/((-1.02*msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input))) + PolyLog(2,1 - (0.9803921568627451*Power(msd2(0,0)*msd2(1,1)
      *msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666))/msq2(2,2))*((
      -1.6666666666666667*(2.04*msq2(2,2) - 2*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*
      msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666))*(1.02*msq2(2,2)*Power(
      msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666) - 2*Quad(M3Input) - 3.06*msq2(2,2)*Sqr(M3Input) + 3*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Sqr(M3Input) + 1.0404*Sqr(msq2(2,2))))/Cube(1.02*msq2(2
      ,2) - Sqr(M3Input)) + (26.666666666666668*M3Input*(AtInput - MuInput/TanBeta
      )*Sqr(1.02*msq2(2,2) - Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,
      0)*msu2(1,1),0.16666666666666666)))/((1.02*msq2(2,2) - 0.9800000000000001*
      msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input)))) + Sqr(Log(Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)/Sqr(
      SCALE)))*((-0.8333333333333334*(2.04*msq2(2,2) - 2*Power(msd2(0,0)*msd2(1,1)
      *msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666))*(1.02*msq2(2,
      2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666) - 2*Quad(M3Input) - 3.06*msq2(2,2)*Sqr(M3Input) + 3*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Sqr(M3Input) + 1.0404*Sqr(msq2(2,2))))/Cube(1.02*msq2(2
      ,2) - Sqr(M3Input)) - (0.8333333333333334*(-2*Power(msd2(0,0)*msd2(1,1)*msq2
      (0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666) +
      1.9600000000000002*msu2(2,2))*(0.9800000000000001*Power(msd2(0,0)*msd2(1,1)*
      msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2) - 2*
      Quad(M3Input) + 3*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*Sqr(M3Input) - 2.9400000000000004*msu2(2,2)*
      Sqr(M3Input) + 0.9604000000000001*Sqr(msu2(2,2))))/Cube(0.9800000000000001*
      msu2(2,2) - Sqr(M3Input)) + (AtInput - MuInput/TanBeta)*((13.333333333333334
      *M3Input*Sqr(1.02*msq2(2,2) - Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)))/((1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) + (
      13.333333333333334*M3Input*Sqr(-Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1
      )*msu2(0,0)*msu2(1,1),0.16666666666666666) + 0.9800000000000001*msu2(2,2)))/
      ((-1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001*
      msu2(2,2) + Sqr(M3Input))))) + Sqr(Log((1.02*msq2(2,2))/Sqr(SCALE)))*((
      0.05555555555555555*(128*Power(M3Input,16) - 567.12*Power(M3Input,14)*msq2(2
      ,2) - 29.96401439808001*Cube(msq2(2,2))*Cube(msu2(2,2))*Cbrt(msd2(0,0)*msd2(
      1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1)) - 333.20000000000005*Power(
      M3Input,14)*msu2(2,2) - 395.83058400000004*Cube(msq2(2,2))*Power10(M3Input)
      + 4.705960000000001*Cube(msu2(2,2))*Power10(M3Input) - 91.8*msq2(2,2)*Cbrt(
      msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Power10(M3Input
      ) - 299.88*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)
      *msu2(1,1),0.16666666666666666)*msu2(2,2)*Power10(M3Input) + 88.2*Cbrt(msd2(
      0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*msu2(2,2)*Power10(
      M3Input) + 61.2*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2
      (0,0)*msu2(1,1),0.16666666666666666)*Power12(M3Input) + 1711.3152000000002*
      msq2(2,2)*msu2(2,2)*Power12(M3Input) - 58.800000000000004*Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(
      2,2)*Power12(M3Input) - 17.665584328532052*Cube(msu2(2,2))*Power5(msq2(2,2))
      - 403.2066528000001*Cube(msu2(2,2))*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2
      (0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power6(M3Input) +
      31.836240000000004*Cube(msq2(2,2))*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1
      ,1)*msu2(0,0)*msu2(1,1))*Power6(M3Input) + 254.12184000000008*Cube(msu2(2,2)
      )*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Power6(
      M3Input) + 561.5912736000001*Cube(msq2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(
      0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Power6(
      M3Input) - 46.3713937344*Power5(msq2(2,2))*Power6(M3Input) -
      3.3771358628888852*msu2(2,2)*Power7(msq2(2,2)) + 694.0914523200003*Cube(msu2
      (2,2))*msq2(2,2)*Power8(M3Input) - 191.01744000000002*Cube(msq2(2,2))*Power(
      msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Power8(M3Input) - 169.41456000000005*Cube(msu2(2,2))*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Power8(M3Input) + 970.3049227200003*Cube(msq2(2,2))*
      msu2(2,2)*Power8(M3Input) + 209.91600000000003*msq2(2,2)*Cbrt(msd2(0,0)*msd2
      (1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*msu2(2,2)*Power8(M3Input) +
      609.2682927609602*Cube(msq2(2,2))*Cube(msu2(2,2))*Quad(M3Input) -
      86.40142560000002*Cube(msu2(2,2))*msq2(2,2)*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,
      0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Quad(M3Input) - 93.59854560000002*Cube(
      msq2(2,2))*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))
      *msu2(2,2)*Quad(M3Input) - 248.85981304128003*msu2(2,2)*Power5(msq2(2,2))*
      Quad(M3Input) + 42.794171932032015*Power6(msq2(2,2))*Quad(M3Input) +
      112.44305278080002*msu2(2,2)*Power6(M3Input)*Quad(msq2(2,2)) + 103.91348736*
      Power8(M3Input)*Quad(msq2(2,2)) - 56.44893139200003*msq2(2,2)*Power6(M3Input
      )*Quad(msu2(2,2)) + 55.34208960000002*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*
      msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power6(M3Input)*Quad(msu2
      (2,2)) - 39.66183088000001*Power8(M3Input)*Quad(msu2(2,2)) +
      112.89786278400004*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(M3Input)*Quad(msu2(2,2)) -
      83.01313440000003*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1))*Quad(M3Input)*Quad(msu2(2,2)) + 20.966420154624547*Quad(msq2(2,2)
      )*Quad(msu2(2,2)) + 179.78408638848006*Cube(msq2(2,2))*Cube(msu2(2,2))*Power
      (msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Sqr(M3Input) + 51.87104103129985*msu2(2,2)*Power6(msq2(
      2,2))*Sqr(M3Input) - 10.338171008843524*Power7(msq2(2,2))*Sqr(M3Input) -
      121.2344022546317*Cube(msu2(2,2))*Quad(msq2(2,2))*Sqr(M3Input) -
      82.22125550833157*Cube(msq2(2,2))*Quad(msu2(2,2))*Sqr(M3Input) +
      56.44893139200002*msq2(2,2)*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1))*Quad(msu2(2,2))*Sqr(M3Input) + 124.848*Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*
      Power10(M3Input)*Sqr(msq2(2,2)) - 2252.2787280000002*msu2(2,2)*Power10(
      M3Input)*Sqr(msq2(2,2)) + 742.8456*Power12(M3Input)*Sqr(msq2(2,2)) -
      1161.3503619648004*Cube(msu2(2,2))*Power6(M3Input)*Sqr(msq2(2,2)) -
      214.11432000000002*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1))*msu2(2,2)*Power6(M3Input)*Sqr(msq2(2,2)) + 62.424*Cbrt(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Power8(M3Input)*Sqr(msq2(
      2,2)) - 183.52656000000002*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Power8(M3Input)*Sqr(msq2(
      2,2)) + 411.27078585600015*Cube(msu2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,
      0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(M3Input)*Sqr(msq2
      (2,2)) + 28.78895500992001*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2
      (0,0)*msu2(1,1))*Quad(msu2(2,2))*Sqr(msq2(2,2)) + 155.46035705356806*Quad(
      M3Input)*Quad(msu2(2,2))*Sqr(msq2(2,2)) - 146.88242352000003*Cube(msu2(2,2))
      *Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Sqr(
      M3Input)*Sqr(msq2(2,2)) - 172.73373005952004*Power(msd2(0,0)*msd2(1,1)*msq2(
      0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(msu2(2,2))*Sqr(
      M3Input)*Sqr(msq2(2,2)) - 1864.1940240000004*msq2(2,2)*Power10(M3Input)*Sqr(
      msu2(2,2)) + 172.872*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)
      *msu2(1,1),0.16666666666666666)*Power10(M3Input)*Sqr(msu2(2,2)) +
      251.62480000000005*Power12(M3Input)*Sqr(msu2(2,2)) - 1534.8913497792005*Cube
      (msq2(2,2))*Power6(M3Input)*Sqr(msu2(2,2)) - 88.16472*msq2(2,2)*Cbrt(msd2(0,
      0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Power6(M3Input)*Sqr(
      msu2(2,2)) - 1.081566387461146*Power6(msq2(2,2))*Sqr(msu2(2,2)) + 528.98832*
      msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Power8(M3Input)*Sqr(msu2(2,2)) - 259.30800000000005*
      Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Power8(
      M3Input)*Sqr(msu2(2,2)) - 550.3594481280002*Cube(msq2(2,2))*Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(
      M3Input)*Sqr(msu2(2,2)) + 328.50343948262406*Quad(M3Input)*Quad(msq2(2,2))*
      Sqr(msu2(2,2)) + 91.72657468800003*Cube(msq2(2,2))*Cbrt(msd2(0,0)*msd2(1,1)*
      msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Sqr(M3Input)*Sqr(msu2(2,2)) +
      39.23329052555137*Power5(msq2(2,2))*Sqr(M3Input)*Sqr(msu2(2,2)) -
      179.85602880000002*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*Power6(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))
      + 2768.7836433600005*Power8(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2)) +
      269.78404320000004*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1))*Quad(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Cube(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))*(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Quad(-1.02*msq2(2,2) + Sqr(M3Input))) + (
      7.3984*Sqr(M3Input)*Sqr(AtInput - MuInput/TanBeta)*Sqr(msq2(2,2)))/(Sqr(1.02
      *msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input
      ))) + (AtInput - MuInput/TanBeta)*((13.333333333333334*M3Input*Sqr(1.02*msq2
      (2,2) - Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)))/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(
      -1.02*msq2(2,2) + Sqr(M3Input))) + (0.4444444444444444*(1.02*Cube(M3Input)*
      msq2(2,2) - 36.260000000000005*Cube(M3Input)*msu2(2,2) - 6.9972*M3Input*msq2
      (2,2)*msu2(2,2) + 18*Power5(M3Input) + 3.1212*M3Input*Sqr(msq2(2,2)) +
      21.128800000000002*M3Input*Sqr(msu2(2,2))))/((-1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input))) + (0.4444444444444444*(-11.985605759232007*M3Input*Cube(msq2(2,2)
      )*Cube(msu2(2,2)) - 21.120348480000008*Cube(msu2(2,2))*msq2(2,2)*Power5(
      M3Input) + 164.31744672000005*Cube(msq2(2,2))*msu2(2,2)*Power5(M3Input) +
      8.8326464256*Cube(M3Input)*Power5(msq2(2,2)) + 10.819991871360003*M3Input*
      msu2(2,2)*Power5(msq2(2,2)) - 2.7117623904000006*Cube(M3Input)*Power5(msu2(2
      ,2)) + 2.7659976382080007*M3Input*msq2(2,2)*Power5(msu2(2,2)) -
      3.378487257792001*M3Input*Power6(msq2(2,2)) - 79.59060000000001*Cube(msq2(2,
      2))*Power7(M3Input) - 4.705960000000001*Cube(msu2(2,2))*Power7(M3Input) -
      16.9932*msq2(2,2)*msu2(2,2)*Power9(M3Input) - 112.4430527808*Cube(M3Input)*
      msu2(2,2)*Quad(msq2(2,2)) + 35.72026128*Power5(M3Input)*Quad(msq2(2,2)) +
      4.704077616000002*Cube(M3Input)*msq2(2,2)*Quad(msu2(2,2)) +
      4.611840800000001*Power5(M3Input)*Quad(msu2(2,2)) + 42.10629474240001*Cube(
      M3Input)*Cube(msu2(2,2))*Sqr(msq2(2,2)) - 33.646536000000005*msu2(2,2)*
      Power7(M3Input)*Sqr(msq2(2,2)) + 34.3332*Power9(M3Input)*Sqr(msq2(2,2)) -
      9.596318336640003*M3Input*Quad(msu2(2,2))*Sqr(msq2(2,2)) -
      11.211025795200003*Cube(M3Input)*Cube(msq2(2,2))*Sqr(msu2(2,2)) +
      48.000792000000004*msq2(2,2)*Power7(M3Input)*Sqr(msu2(2,2)) +
      29.107899700992004*M3Input*Quad(msq2(2,2))*Sqr(msu2(2,2)) -
      77.93761248000001*Power5(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Cube(-
      1.02*msq2(2,2) + Sqr(M3Input))*Sqr(1.02*msq2(2,2) - 0.9800000000000001*msu2(
      2,2))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (
      0.4444444444444444*(3.1836240000000005*M3Input*Cube(msq2(2,2)) -
      2.823576000000001*M3Input*Cube(msu2(2,2)) + 5.1*msq2(2,2)*Power5(M3Input) -
      4.9*msu2(2,2)*Power5(M3Input) - 5.202*Cube(M3Input)*Sqr(msq2(2,2)) -
      4.078368*M3Input*msu2(2,2)*Sqr(msq2(2,2)) + 4.8020000000000005*Cube(M3Input)
      *Sqr(msu2(2,2)) + 3.9184320000000006*M3Input*msq2(2,2)*Sqr(msu2(2,2))))/(Sqr
      (-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input))))) + PolyLog(2,1 - (1.020408163265306*Power(msd2(0,0)*msd2(1,1)*
      msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666))/msu2(2,2))*((-
      1.6666666666666667*(-2*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,
      0)*msu2(1,1),0.16666666666666666) + 1.9600000000000002*msu2(2,2))*(
      0.9800000000000001*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*msu2(2,2) - 2*Quad(M3Input) + 3*Power(msd2(0,
      0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*
      Sqr(M3Input) - 2.9400000000000004*msu2(2,2)*Sqr(M3Input) +
      0.9604000000000001*Sqr(msu2(2,2))))/Cube(0.9800000000000001*msu2(2,2) - Sqr(
      M3Input)) + (26.666666666666668*M3Input*(AtInput - MuInput/TanBeta)*Sqr(-
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666) + 0.9800000000000001*msu2(2,2)))/((-1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input)))) + Sqr(Log((0.9800000000000001*msu2(2,2))/Sqr(SCALE)))*((-
      0.4444444444444444*(AtInput - MuInput/TanBeta)*(3.1836240000000005*Cube(
      M3Input)*Cube(msq2(2,2)) + 93.17800800000002*Cube(M3Input)*Cube(msu2(2,2)) +
      3.840063360000001*M3Input*Cube(msu2(2,2))*msq2(2,2) - 56.47152000000001*
      M3Input*Cube(msu2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0
      ,0)*msu2(1,1),0.16666666666666666) + 30.6*Cube(M3Input)*msq2(2,2)*Cbrt(msd2(
      0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1)) - 3.1199515200000008
      *M3Input*Cube(msq2(2,2))*msu2(2,2) - 59.976000000000006*Cube(M3Input)*msq2(2
      ,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2) - 29.400000000000002*Cube(M3Input)*Cbrt(msd2(
      0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*msu2(2,2) -
      29.988000000000003*M3Input*msq2(2,2)*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2
      (1,1)*msu2(0,0)*msu2(1,1))*msu2(2,2) + 68.97240000000001*msq2(2,2)*msu2(2,2)
      *Power5(M3Input) - 18.36*msq2(2,2)*Power7(M3Input) + 17.64*msu2(2,2)*Power7(
      M3Input) - 24.903940320000007*M3Input*Quad(msu2(2,2)) - 11.215512*Cube(
      M3Input)*msu2(2,2)*Sqr(msq2(2,2)) + 1.0404*Power5(M3Input)*Sqr(msq2(2,2)) -
      57.79687200000001*Cube(M3Input)*msq2(2,2)*Sqr(msu2(2,2)) + 57.62400000000001
      *Cube(M3Input)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(
      1,1),0.16666666666666666)*Sqr(msu2(2,2)) + 58.77648000000001*M3Input*msq2(2,
      2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Sqr(msu2(2,2)) + 28.812000000000005*M3Input*Cbrt(msd2(0
      ,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Sqr(msu2(2,2)) -
      82.59440000000001*Power5(M3Input)*Sqr(msu2(2,2)) + 9.992001600000002*M3Input
      *Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Cube(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input))*Sqr(-1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))) + (
      0.05555555555555555*(28.235760000000006*Cube(msu2(2,2))*Cbrt(msd2(0,0)*msd2(
      1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1)) + 128*Power10(M3Input) +
      18.982336732800004*Power5(msu2(2,2)) - 61.2*msq2(2,2)*Power(msd2(0,0)*msd2(1
      ,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power6(
      M3Input) + 67.9728*msq2(2,2)*msu2(2,2)*Power6(M3Input) + 58.800000000000004*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2)*Power6(M3Input) + 44.88*msq2(2,2)*Power8(
      M3Input) - 544.88*msu2(2,2)*Power8(M3Input) + 9.550872000000002*Cube(msq2(2,
      2))*Quad(M3Input) - 326.5936240000001*Cube(msu2(2,2))*Quad(M3Input) + 91.8*
      msq2(2,2)*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*
      Quad(M3Input) - 119.95200000000001*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(
      0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Quad(
      M3Input) - 88.2*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(
      1,1))*msu2(2,2)*Quad(M3Input) - 15.993863894400006*msq2(2,2)*Quad(msu2(2,2))
      + 48.00079200000001*Cube(msu2(2,2))*msq2(2,2)*Sqr(M3Input) -
      169.41456000000005*Cube(msu2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(
      1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Sqr(M3Input) -
      6.2399030400000015*Cube(msq2(2,2))*msu2(2,2)*Sqr(M3Input) -
      59.976000000000006*msq2(2,2)*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1))*msu2(2,2)*Sqr(M3Input) + 44.273671680000014*Quad(msu2(2
      ,2))*Sqr(M3Input) - 0.9792161568000003*Cube(msu2(2,2))*Sqr(msq2(2,2)) -
      2.0808*Power6(M3Input)*Sqr(msq2(2,2)) - 33.646536000000005*msu2(2,2)*Quad(
      M3Input)*Sqr(msq2(2,2)) - 3.0575524896000013*Cube(msq2(2,2))*Sqr(msu2(2,2))
      - 29.388240000000007*msq2(2,2)*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1))*Sqr(msu2(2,2)) + 674.2008000000001*Power6(M3Input)*Sqr(
      msu2(2,2)) - 138.124728*msq2(2,2)*Quad(M3Input)*Sqr(msu2(2,2)) +
      115.24800000000002*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*Quad(M3Input)*Sqr(msu2(2,2)) +
      176.32944000000003*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*Sqr(M3Input)*Sqr(msu2(2,2)) +
      57.62400000000001*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1))*Sqr(M3Input)*Sqr(msu2(2,2)) + 35.971205760000004*Sqr(M3Input)*Sqr
      (msq2(2,2))*Sqr(msu2(2,2))))/((-1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2
      ))*Quad(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (6.8295111111111115
      *Sqr(M3Input)*Sqr(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/(Sqr(-1.02*msq2
      (2,2) + 0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input)))) + PolyLog(2,1 - (0.9803921568627451*Sqr(M3Input))/msq2(2,2))
      *((-0.1111111111111111*(18.977209118784007*Cube(msq2(2,2))*Cube(msu2(2,2)) +
      299.88*msq2(2,2)*Power10(M3Input) + 339.08000000000004*msu2(2,2)*Power10(
      M3Input) - 128*Power12(M3Input) + 3.2459975614080006*msu2(2,2)*Power5(msq2(2
      ,2)) - 9.550872000000002*Cube(msq2(2,2))*Power6(M3Input) +
      12.235496000000003*Cube(msu2(2,2))*Power6(M3Input) - 901.6392000000001*msq2(
      2,2)*msu2(2,2)*Power8(M3Input) - 375.3661934400001*Cube(msu2(2,2))*msq2(2,2)
      *Quad(M3Input) + 157.03755984000006*Cube(msq2(2,2))*msu2(2,2)*Quad(M3Input)
      - 21.6486432*Quad(M3Input)*Quad(msq2(2,2)) + 34.12762192000001*Quad(M3Input)
      *Quad(msu2(2,2)) + 9.9367272288*Power5(msq2(2,2))*Sqr(M3Input) -
      43.492124188800005*msu2(2,2)*Quad(msq2(2,2))*Sqr(M3Input) +
      31.987727788800015*msq2(2,2)*Quad(msu2(2,2))*Sqr(M3Input) + 255.917592*msu2(
      2,2)*Power6(M3Input)*Sqr(msq2(2,2)) - 101.9592*Power8(M3Input)*Sqr(msq2(2,2)
      ) - 22.071532174272008*Quad(msu2(2,2))*Sqr(msq2(2,2)) + 163.52909818560005*
      Cube(msu2(2,2))*Sqr(M3Input)*Sqr(msq2(2,2)) + 1004.0982000000001*msq2(2,2)*
      Power6(M3Input)*Sqr(msu2(2,2)) - 268.91200000000003*Power8(M3Input)*Sqr(msu2
      (2,2)) + 1.0395678464640001*Quad(msq2(2,2))*Sqr(msu2(2,2)) -
      41.78655069120001*Cube(msq2(2,2))*Sqr(M3Input)*Sqr(msu2(2,2)) -
      416.6664667200001*Quad(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Cube(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))*(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) + (
      0.8888888888888888*(AtInput - MuInput/TanBeta)*(1.02*Cube(M3Input)*msq2(2,2)
      - 36.260000000000005*Cube(M3Input)*msu2(2,2) - 6.9972*M3Input*msq2(2,2)*msu2
      (2,2) + 18*Power5(M3Input) + 3.1212*M3Input*Sqr(msq2(2,2)) +
      21.128800000000002*M3Input*Sqr(msu2(2,2))))/((-1.02*msq2(2,2) +
      0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input)))) + Log((0.9800000000000001*msu2(2,2))/Sqr(SCALE))*((
      0.1111111111111111*(130.84286287161606*Cube(msq2(2,2))*Cube(msu2(2,2)) +
      53.04*msq2(2,2)*Power10(M3Input) - 301.84000000000003*msu2(2,2)*Power10(
      M3Input) - 2.7659976382080007*msq2(2,2)*Power5(msu2(2,2)) +
      56.24402400000001*Cube(msq2(2,2))*Power6(M3Input) + 16.941456000000006*Cube(
      msu2(2,2))*Power6(M3Input) - 89.964*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2
      (0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Power6(
      M3Input) + 763.6944000000001*msq2(2,2)*msu2(2,2)*Power8(M3Input) +
      24.000396000000006*Cube(msu2(2,2))*msq2(2,2)*Quad(M3Input) +
      28.235760000000006*Cube(msu2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(
      1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(M3Input) +
      169.51736592000006*Cube(msq2(2,2))*msu2(2,2)*Quad(M3Input) -
      97.77102496000003*Quad(M3Input)*Quad(msu2(2,2)) - 57.600950400000016*Cube(
      msu2(2,2))*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)
      *msu2(1,1),0.16666666666666666)*Sqr(M3Input) - 93.59854560000002*Cube(msq2(2
      ,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2)*Sqr(M3Input) - 8.135287171200002*Power5(msu2(
      2,2))*Sqr(M3Input) - 9.5470516512*msu2(2,2)*Quad(msq2(2,2))*Sqr(M3Input) +
      262.48753097280013*msq2(2,2)*Quad(msu2(2,2))*Sqr(M3Input) +
      29.37648470400001*Cube(msu2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1
      ,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Sqr(msq2(2,2)) - 634.186224*
      msu2(2,2)*Power6(M3Input)*Sqr(msq2(2,2)) - 109.242*Power8(M3Input)*Sqr(msq2(
      2,2)) + 183.52656000000002*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Quad(M3Input)*Sqr(msq2(2,
      2)) - 119.95397920800004*Quad(msu2(2,2))*Sqr(msq2(2,2)) - 240.88717457280006
      *Cube(msu2(2,2))*Sqr(M3Input)*Sqr(msq2(2,2)) - 30.575524896000008*Cube(msq2(
      2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Sqr(msu2(2,2)) - 964.9138800000002*msq2(2,2)*Power6(
      M3Input)*Sqr(msu2(2,2)) + 86.436*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,
      1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power6(M3Input)*Sqr(msu2(2,2)) +
      350.54600000000005*Power8(M3Input)*Sqr(msu2(2,2)) - 205.71768000000003*msq2(
      2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Quad(M3Input)*Sqr(msu2(2,2)) - 3.118703539392001*Quad(
      msq2(2,2))*Sqr(msu2(2,2)) - 276.1989082272001*Cube(msq2(2,2))*Sqr(M3Input)*
      Sqr(msu2(2,2)) + 941.2465507200001*Quad(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2
      )) + 149.88002400000002*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0
      ,0)*msu2(1,1),0.16666666666666666)*Sqr(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2)
      )))/(Cube(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))*(1.02*msq2(2,2) -
      0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) + Log(
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)/Sqr(SCALE))*((-1.1111111111111112*(0.9411920000000003*
      Cube(msu2(2,2)) - 2.9400000000000004*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2
      (1,1)*msu2(0,0)*msu2(1,1))*msu2(2,2) + 6*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)
      *msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(M3Input) +
      1.9600000000000002*msu2(2,2)*Quad(M3Input) - 9*Cbrt(msd2(0,0)*msd2(1,1)*msq2
      (0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Sqr(M3Input) + 17.64*Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(
      2,2)*Sqr(M3Input) - 2.8812000000000006*Sqr(M3Input)*Sqr(msu2(2,2))))/Cube(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input)) + (4.444444444444445*(AtInput -
      MuInput/TanBeta)*(-6*M3Input*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1)) + 0.9800000000000001*Cube(M3Input)*msu2(2,2) +
      11.760000000000002*M3Input*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2) - 0.9604000000000001*
      M3Input*Sqr(msu2(2,2))))/((-1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))*
      Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + (0.8888888888888888*(
      AtInput - MuInput/TanBeta)*(2.823576000000001*M3Input*Cube(msu2(2,2)) -
      52.97880000000001*Cube(M3Input)*msq2(2,2)*msu2(2,2) - 29.400000000000002*
      Cube(M3Input)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1
      ,1),0.16666666666666666)*msu2(2,2) + 29.988000000000003*M3Input*msq2(2,2)*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2) - 16.32*msq2(2,2)*Power5(M3Input) +
      53.900000000000006*msu2(2,2)*Power5(M3Input) + 16*Power7(M3Input) + 3.058776
      *M3Input*msu2(2,2)*Sqr(msq2(2,2)) - 58.58440000000001*Cube(M3Input)*Sqr(msu2
      (2,2)) + 51.91922400000001*M3Input*msq2(2,2)*Sqr(msu2(2,2))))/((1.02*msq2(2,
      2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + PolyLog(2,1 - (
      1.020408163265306*Sqr(M3Input))/msu2(2,2))*((0.8888888888888888*(AtInput -
      MuInput/TanBeta)*(-37.74*Cube(M3Input)*msq2(2,2) + 0.9800000000000001*Cube(
      M3Input)*msu2(2,2) - 6.9972*M3Input*msq2(2,2)*msu2(2,2) + 18*Power5(M3Input)
      + 22.8888*M3Input*Sqr(msq2(2,2)) + 2.8812000000000006*M3Input*Sqr(msu2(2,2))
      ))/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) +
      Sqr(M3Input))) + (0.1111111111111111*(18.977209118784007*Cube(msq2(2,2))*
      Cube(msu2(2,2)) + 352.92*msq2(2,2)*Power10(M3Input) + 288.12*msu2(2,2)*
      Power10(M3Input) - 128*Power12(M3Input) + 2.7659976382080007*msq2(2,2)*
      Power5(msu2(2,2)) + 13.795704000000002*Cube(msq2(2,2))*Power6(M3Input) -
      8.470728000000003*Cube(msu2(2,2))*Power6(M3Input) - 901.6392000000001*msq2(2
      ,2)*msu2(2,2)*Power8(M3Input) + 144.96239184000004*Cube(msu2(2,2))*msq2(2,2)
      *Quad(M3Input) - 406.63368144000015*Cube(msq2(2,2))*msu2(2,2)*Quad(M3Input)
      + 40.04998992*Quad(M3Input)*Quad(msq2(2,2)) - 18.447363200000005*Quad(
      M3Input)*Quad(msu2(2,2)) + 8.135287171200002*Power5(msu2(2,2))*Sqr(M3Input)
      + 36.06663957120001*msu2(2,2)*Quad(msq2(2,2))*Sqr(M3Input) -
      38.57343645120002*msq2(2,2)*Quad(msu2(2,2))*Sqr(M3Input) +
      1045.0818000000002*msu2(2,2)*Power6(M3Input)*Sqr(msq2(2,2)) - 291.312*Power8
      (M3Input)*Sqr(msq2(2,2)) + 0.9596318336640003*Quad(msu2(2,2))*Sqr(msq2(2,2))
      - 40.14786242880001*Cube(msu2(2,2))*Sqr(M3Input)*Sqr(msq2(2,2)) +
      245.88160800000003*msq2(2,2)*Power6(M3Input)*Sqr(msu2(2,2)) -
      94.11920000000002*Power8(M3Input)*Sqr(msu2(2,2)) - 23.910060468672*Quad(msq2
      (2,2))*Sqr(msu2(2,2)) + 170.20375525440005*Cube(msq2(2,2))*Sqr(M3Input)*Sqr(
      msu2(2,2)) - 416.6664667200001*Quad(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))))
      /(Cube(-1.02*msq2(2,2) + Sqr(M3Input))*(1.02*msq2(2,2) - 0.9800000000000001*
      msu2(2,2))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + Log((1.02*
      msq2(2,2))/Sqr(SCALE))*(Log(Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)/Sqr(SCALE))*((-1.1111111111111112*(
      1.0612080000000002*Cube(msq2(2,2)) - 3.06*msq2(2,2)*Cbrt(msd2(0,0)*msd2(1,1)
      *msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1)) + 2.04*msq2(2,2)*Quad(M3Input) + 6
      *Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Quad(M3Input) + 18.36*msq2(2,2)*Power(msd2(0,0)*msd2(1,
      1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Sqr(M3Input)
      - 9*Cbrt(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1))*Sqr(
      M3Input) - 3.1212*Sqr(M3Input)*Sqr(msq2(2,2))))/Cube(-1.02*msq2(2,2) + Sqr(
      M3Input)) + (4.444444444444445*(AtInput - MuInput/TanBeta)*(1.02*Cube(
      M3Input)*msq2(2,2) + 12.24*M3Input*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(
      0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666) - 6*M3Input*Cbrt(
      msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1)) - 1.0404*
      M3Input*Sqr(msq2(2,2))))/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*
      Sqr(-1.02*msq2(2,2) + Sqr(M3Input)))) - (0.8888888888888888*(AtInput -
      MuInput/TanBeta)*(3.1836240000000005*M3Input*Cube(msq2(2,2)) - 30.6*Cube(
      M3Input)*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666) - 52.97880000000001*Cube(M3Input)*msq2(2,2)*
      msu2(2,2) + 29.988000000000003*M3Input*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*
      msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2) +
      56.1*msq2(2,2)*Power5(M3Input) - 15.680000000000001*msu2(2,2)*Power5(M3Input
      ) + 16*Power7(M3Input) - 63.4644*Cube(M3Input)*Sqr(msq2(2,2)) + 54.038376*
      M3Input*msu2(2,2)*Sqr(msq2(2,2)) + 2.9388240000000003*M3Input*msq2(2,2)*Sqr(
      msu2(2,2))))/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))*Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input))) + Log((0.9800000000000001*msu2(2,2))/Sqr(SCALE))*((-
      14.216533333333334*msq2(2,2)*msu2(2,2)*Sqr(M3Input)*Sqr(AtInput - MuInput/
      TanBeta))/((-1.02*msq2(2,2) + Sqr(M3Input))*(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input))*Sqr(-1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))) + (
      0.1111111111111111*(2.996401439808001*Cube(msq2(2,2))*Cube(msu2(2,2)) -
      14.280000000000001*msq2(2,2)*Power10(M3Input) - 1.9600000000000002*msu2(2,2)
      *Power10(M3Input) - 3.2459975614080006*msu2(2,2)*Power5(msq2(2,2)) -
      37.14228000000001*Cube(msq2(2,2))*Power6(M3Input) - 0.9411920000000003*Cube(
      msu2(2,2))*Power6(M3Input) + 27.988800000000005*msq2(2,2)*msu2(2,2)*Power8(
      M3Input) + 18.240300960000006*Cube(msu2(2,2))*msq2(2,2)*Quad(M3Input) -
      77.998788*Cube(msq2(2,2))*msu2(2,2)*Quad(M3Input) + 31.39053264*Quad(M3Input
      )*Quad(msq2(2,2)) - 9.9367272288*Power5(msq2(2,2))*Sqr(M3Input) +
      37.12742308800001*msu2(2,2)*Quad(msq2(2,2))*Sqr(M3Input) + 54.038376*msu2(2,
      2)*Power6(M3Input)*Sqr(msq2(2,2)) + 17.686799999999998*Power8(M3Input)*Sqr(
      msq2(2,2)) - 8.812945411200003*Cube(msu2(2,2))*Sqr(M3Input)*Sqr(msq2(2,2)) -
      55.83765600000001*msq2(2,2)*Power6(M3Input)*Sqr(msu2(2,2)) +
      2.8812000000000006*Power8(M3Input)*Sqr(msu2(2,2)) - 4.1582713858560005*Quad(
      msq2(2,2))*Sqr(msu2(2,2)) - 1.0191841632000003*Cube(msq2(2,2))*Sqr(M3Input)*
      Sqr(msu2(2,2)) + 22.981603680000003*Quad(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,
      2))))/(Cube(-1.02*msq2(2,2) + Sqr(M3Input))*Cube(-0.9800000000000001*msu2(2,
      2) + Sqr(M3Input))) + (0.4444444444444444*(AtInput - MuInput/TanBeta)*(-
      8.640142560000003*Cube(M3Input)*Cube(msu2(2,2))*msq2(2,2) +
      11.439822240000003*Cube(M3Input)*Cube(msq2(2,2))*msu2(2,2) +
      9.550872000000002*Cube(msq2(2,2))*Power5(M3Input) - 0.9411920000000003*Cube(
      msu2(2,2))*Power5(M3Input) + 6.6244848192*M3Input*Power5(msq2(2,2)) -
      33.9864*msq2(2,2)*msu2(2,2)*Power7(M3Input) - 10.8243216*Cube(M3Input)*Quad(
      msq2(2,2)) - 21.215670336000002*M3Input*msu2(2,2)*Quad(msq2(2,2)) +
      7.833729254400002*M3Input*Cube(msu2(2,2))*Sqr(msq2(2,2)) + 29.568168*msu2(2,
      2)*Power5(M3Input)*Sqr(msq2(2,2)) + 1.0404*Power7(M3Input)*Sqr(msq2(2,2)) +
      38.72899820160001*M3Input*Cube(msq2(2,2))*Sqr(msu2(2,2)) + 57.79687200000001
      *msq2(2,2)*Power5(M3Input)*Sqr(msu2(2,2)) + 0.9604000000000001*Power7(
      M3Input)*Sqr(msu2(2,2)) - 87.92961408000001*Cube(M3Input)*Sqr(msq2(2,2))*Sqr
      (msu2(2,2))))/(Sqr(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))
      ) + (0.1111111111111111*(-130.84286287161606*Cube(msq2(2,2))*Cube(msu2(2,2))
      + 314.16*msq2(2,2)*Power10(M3Input) - 50.96000000000001*msu2(2,2)*Power10(
      M3Input) + 3.2459975614080006*msu2(2,2)*Power5(msq2(2,2)) -
      19.101744000000004*Cube(msq2(2,2))*Power6(M3Input) - 49.88317600000001*Cube(
      msu2(2,2))*Power6(M3Input) + 89.96400000000001*msq2(2,2)*Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(
      2,2)*Power6(M3Input) - 763.6944000000001*msq2(2,2)*msu2(2,2)*Power8(M3Input)
      - 156.48258192000003*Cube(msu2(2,2))*msq2(2,2)*Quad(M3Input) -
      31.836240000000004*Cube(msq2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(
      1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(M3Input) -
      25.999596000000007*Cube(msq2(2,2))*msu2(2,2)*Quad(M3Input) + 114.73780896*
      Quad(M3Input)*Quad(msq2(2,2)) + 86.40142560000002*Cube(msu2(2,2))*msq2(2,2)*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Sqr(M3Input) + 62.39903040000002*Cube(msq2(2,2))*Power(
      msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2)*Sqr(M3Input) + 9.9367272288*Power5(msq2(2,2))
      *Sqr(M3Input) - 295.95860118720003*msu2(2,2)*Quad(msq2(2,2))*Sqr(M3Input) +
      8.467339708800004*msq2(2,2)*Quad(msu2(2,2))*Sqr(M3Input) + 29.37648470400001
      *Cube(msu2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*Sqr(msq2(2,2)) - 93.636*Power(msd2(0,0)*msd2(
      1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power6(
      M3Input)*Sqr(msq2(2,2)) + 1004.29812*msu2(2,2)*Power6(M3Input)*Sqr(msq2(2,2)
      ) - 379.746*Power8(M3Input)*Sqr(msq2(2,2)) + 214.11432000000002*Power(msd2(0
      ,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*
      msu2(2,2)*Quad(M3Input)*Sqr(msq2(2,2)) + 2.8788955009920008*Quad(msu2(2,2))*
      Sqr(msq2(2,2)) + 265.3675784928001*Cube(msu2(2,2))*Sqr(M3Input)*Sqr(msq2(2,2
      )) - 30.575524896000008*Cube(msq2(2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*
      msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Sqr(msu2(2,2)) +
      609.316176*msq2(2,2)*Power6(M3Input)*Sqr(msu2(2,2)) + 100.84200000000001*
      Power8(M3Input)*Sqr(msu2(2,2)) - 176.32944*msq2(2,2)*Power(msd2(0,0)*msd2(1,
      1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Quad(M3Input
      )*Sqr(msu2(2,2)) + 129.94598080800003*Quad(msq2(2,2))*Sqr(msu2(2,2)) +
      250.71930414720006*Cube(msq2(2,2))*Sqr(M3Input)*Sqr(msu2(2,2)) -
      941.2465507200002*Quad(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2)) -
      149.88002400000002*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*Sqr(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/
      (Cube(-1.02*msq2(2,2) + Sqr(M3Input))*(1.02*msq2(2,2) - 0.9800000000000001*
      msu2(2,2))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + PolyLog(2,1
       - Sqr(M3Input)/Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2
      (1,1),0.16666666666666666))*(0.5*(-2*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*
      msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666) + 2*Sqr(M3Input))*((-
      13.599999999999998*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666))/Cube(1.02*msq2(2,2) - Sqr(M3Input)
      ) - (13.066666666666666*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0
      ,0)*msu2(1,1),0.16666666666666666)*msu2(2,2))/Cube(0.9800000000000001*msu2(2
      ,2) - Sqr(M3Input)) - 3.3333333333333335/(1.02*msq2(2,2) - Sqr(M3Input)) -
      3.3333333333333335/(0.9800000000000001*msu2(2,2) - Sqr(M3Input)) + (13.872*
      Sqr(msq2(2,2)))/Cube(1.02*msq2(2,2) - Sqr(M3Input)) + (12.805333333333333*
      Sqr(msu2(2,2)))/Cube(0.9800000000000001*msu2(2,2) - Sqr(M3Input)) - (10.2*
      msq2(2,2))/Sqr(1.02*msq2(2,2) - Sqr(M3Input)) + (10*Power(msd2(0,0)*msd2(1,1
      )*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666))/Sqr(1.02*
      msq2(2,2) - Sqr(M3Input)) + (10*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1
      )*msu2(0,0)*msu2(1,1),0.16666666666666666))/Sqr(0.9800000000000001*msu2(2,2)
      - Sqr(M3Input)) - (9.8*msu2(2,2))/Sqr(0.9800000000000001*msu2(2,2) - Sqr(
      M3Input))) + (26.666666666666668*(AtInput - MuInput/TanBeta)*(1.02*Cube(
      M3Input)*msq2(2,2) - 2*Cube(M3Input)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*
      msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666) + 1.02*M3Input*msq2(2,2)*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666) + 0.9800000000000001*Cube(M3Input)*msu2(2,2) -
      1.9992000000000003*M3Input*msq2(2,2)*msu2(2,2) + 0.9800000000000001*M3Input*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2))*(-Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1
      ,1)*msu2(0,0)*msu2(1,1),0.16666666666666666) + Sqr(M3Input)))/(Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))
      ) + Sqr(Log(Sqr(M3Input)/Sqr(SCALE)))*((-0.8888888888888888*(AtInput -
      MuInput/TanBeta)*(27*Power11(M3Input) + 219.91200000000003*msq2(2,2)*msu2(2,
      2)*Power7(M3Input) - 83.64*msq2(2,2)*Power9(M3Input) - 80.36000000000001*
      msu2(2,2)*Power9(M3Input) - 144.78206400000002*msu2(2,2)*Power5(M3Input)*Sqr
      (msq2(2,2)) + 59.3028*Power7(M3Input)*Sqr(msq2(2,2)) - 139.10433600000005*
      msq2(2,2)*Power5(M3Input)*Sqr(msu2(2,2)) + 54.74280000000001*Power7(M3Input)
      *Sqr(msu2(2,2)) + 86.93041392*Cube(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/
      (Cube(-1.02*msq2(2,2) + Sqr(M3Input))*Cube(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input))) + (0.1111111111111111*(262*Power(M3Input,16) -
      895.5600000000001*Power(M3Input,14)*msq2(2,2) - 860.44*Power(M3Input,14)*
      msu2(2,2) - 339.5865600000001*Cube(msq2(2,2))*Power10(M3Input) -
      301.18144000000007*Cube(msu2(2,2))*Power10(M3Input) + 2970.8112000000006*
      msq2(2,2)*msu2(2,2)*Power12(M3Input) + 963.8559033600003*Cube(msu2(2,2))*
      msq2(2,2)*Power8(M3Input) + 1044.1437753600003*Cube(msq2(2,2))*msu2(2,2)*
      Power8(M3Input) - 451.4578169310722*Cube(msq2(2,2))*Cube(msu2(2,2))*Quad(
      M3Input) - 167.6037956544*msu2(2,2)*Power6(M3Input)*Quad(msq2(2,2)) +
      53.03917584*Power8(M3Input)*Quad(msq2(2,2)) - 148.64885266560006*msq2(2,2)*
      Power6(M3Input)*Quad(msu2(2,2)) + 45.19603984000002*Power8(M3Input)*Quad(
      msu2(2,2)) - 35.94243455078493*Quad(msq2(2,2))*Quad(msu2(2,2)) +
      146.7038144929997*Cube(msu2(2,2))*Quad(msq2(2,2))*Sqr(M3Input) +
      140.9507237285684*Cube(msq2(2,2))*Quad(msu2(2,2))*Sqr(M3Input) - 2989.443744
      *msu2(2,2)*Power10(M3Input)*Sqr(msq2(2,2)) + 920.754*Power12(M3Input)*Sqr(
      msq2(2,2)) - 360.35154570240013*Cube(msu2(2,2))*Power6(M3Input)*Sqr(msq2(2,2
      )) - 0.9596318336640003*Quad(M3Input)*Quad(msu2(2,2))*Sqr(msq2(2,2)) -
      2872.2106560000007*msq2(2,2)*Power10(M3Input)*Sqr(msu2(2,2)) +
      849.9540000000002*Power12(M3Input)*Sqr(msu2(2,2)) - 375.05977205760007*Cube(
      msq2(2,2))*Power6(M3Input)*Sqr(msu2(2,2)) - 1.0395678464640001*Quad(M3Input)
      *Quad(msq2(2,2))*Sqr(msu2(2,2)) + 2402.0771846400003*Power8(M3Input)*Sqr(
      msq2(2,2))*Sqr(msu2(2,2))))/(Quad(-1.02*msq2(2,2) + Sqr(M3Input))*Quad(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + (7.111111111111111*Power6(
      M3Input)*Sqr(AtInput - MuInput/TanBeta))/(Sqr(-1.02*msq2(2,2) + Sqr(M3Input)
      )*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + Log(Sqr(M3Input)/Sqr
      (SCALE))*((0.1111111111111111*(-47.94242303692803*Cube(msq2(2,2))*Cube(msu2(
      2,2)) + 1980.8400000000001*msq2(2,2)*Power10(M3Input) + 180*Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*
      Power10(M3Input) + 1903.16*msu2(2,2)*Power10(M3Input) - 864*Power12(M3Input)
      + 372.4840080000001*Cube(msq2(2,2))*Power6(M3Input) + 330.3583920000001*Cube
      (msu2(2,2))*Power6(M3Input) - 179.928*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*
      msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*
      Power6(M3Input) - 244.8*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1
      ,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power8(M3Input) -
      4086.3648000000003*msq2(2,2)*msu2(2,2)*Power8(M3Input) - 235.20000000000002*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2)*Power8(M3Input) - 549.1290604800002*Cube(msu2
      (2,2))*msq2(2,2)*Quad(M3Input) - 95.50872000000001*Cube(msq2(2,2))*Power(
      msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Quad(M3Input) - 84.70728000000003*Cube(msu2(2,2))*Power
      (msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Quad(M3Input) - 594.8707564800002*Cube(msq2(2,2))*msu2(
      2,2)*Quad(M3Input) - 9.74188944*Quad(M3Input)*Quad(msq2(2,2)) -
      8.301313440000003*Quad(M3Input)*Quad(msu2(2,2)) - 28.800475200000008*Cube(
      msu2(2,2))*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)
      *msu2(1,1),0.16666666666666666)*Sqr(M3Input) - 31.199515200000008*Cube(msq2(
      2,2))*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2)*Sqr(M3Input) - 3.1823505504000003*msu2(2,2)*
      Quad(msq2(2,2))*Sqr(M3Input) - 2.8224465696000007*msq2(2,2)*Quad(msu2(2,2))*
      Sqr(M3Input) + 280.908*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,
      0)*msu2(1,1),0.16666666666666666)*Power6(M3Input)*Sqr(msq2(2,2)) +
      2908.8959760000002*msu2(2,2)*Power6(M3Input)*Sqr(msq2(2,2)) - 1518.984*
      Power8(M3Input)*Sqr(msq2(2,2)) + 91.76328000000001*Power(msd2(0,0)*msd2(1,1)
      *msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Quad
      (M3Input)*Sqr(msq2(2,2)) + 337.8295740960001*Cube(msu2(2,2))*Sqr(M3Input)*
      Sqr(msq2(2,2)) + 2794.8216240000006*msq2(2,2)*Power6(M3Input)*Sqr(msu2(2,2))
      + 259.30800000000005*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)
      *msu2(1,1),0.16666666666666666)*Power6(M3Input)*Sqr(msu2(2,2)) -
      1402.1840000000002*Power8(M3Input)*Sqr(msu2(2,2)) + 88.16472000000002*msq2(2
      ,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Quad(M3Input)*Sqr(msu2(2,2)) + 351.6185363040001*Cube(
      msq2(2,2))*Sqr(M3Input)*Sqr(msu2(2,2)) - 1892.4851030400005*Quad(M3Input)*
      Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Cube(-1.02*msq2(2,2) + Sqr(M3Input))*Cube(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))) + Log((0.9800000000000001*msu2
      (2,2))/Sqr(SCALE))*((0.2222222222222222*(64*Power(M3Input,16) - 208.08*Power
      (M3Input,14)*msq2(2,2) - 176.4*Power(M3Input,14)*msu2(2,2) -
      107.18200800000001*Cube(msq2(2,2))*Power10(M3Input) - 29.176952000000007*
      Cube(msu2(2,2))*Power10(M3Input) + 603.7584*msq2(2,2)*msu2(2,2)*Power12(
      M3Input) + 5.755487885583209*Cube(msq2(2,2))*Power5(msu2(2,2)) -
      6.327445577600002*Power5(msu2(2,2))*Power6(M3Input) + 12.480205920000003*
      Cube(msu2(2,2))*msq2(2,2)*Power8(M3Input) + 417.0335198400001*Cube(msq2(2,2)
      )*msu2(2,2)*Power8(M3Input) - 52.93642543660802*Cube(msq2(2,2))*Cube(msu2(2,
      2))*Quad(M3Input) + 11.985989765568004*msq2(2,2)*Power5(msu2(2,2))*Quad(
      M3Input) - 80.6195472768*msu2(2,2)*Power6(M3Input)*Quad(msq2(2,2)) +
      14.07161808*Power8(M3Input)*Quad(msq2(2,2)) - 46.09996063680002*msq2(2,2)*
      Power6(M3Input)*Quad(msu2(2,2)) + 27.67104480000001*Power8(M3Input)*Quad(
      msu2(2,2)) - 5.9904057584641555*Quad(msq2(2,2))*Quad(msu2(2,2)) +
      24.450635748833285*Cube(msu2(2,2))*Quad(msq2(2,2))*Sqr(M3Input) -
      5.872946822023683*Cube(msq2(2,2))*Quad(msu2(2,2))*Sqr(M3Input) -
      794.2621680000001*msu2(2,2)*Power10(M3Input)*Sqr(msq2(2,2)) + 243.4536*
      Power12(M3Input)*Sqr(msq2(2,2)) - 12.729810038400004*Cube(msu2(2,2))*Power6(
      M3Input)*Sqr(msq2(2,2)) + 58.53754185350402*Quad(M3Input)*Quad(msu2(2,2))*
      Sqr(msq2(2,2)) - 16.927905545832964*Power5(msu2(2,2))*Sqr(M3Input)*Sqr(msq2(
      2,2)) - 361.4753520000001*msq2(2,2)*Power10(M3Input)*Sqr(msu2(2,2)) +
      117.16880000000002*Power12(M3Input)*Sqr(msu2(2,2)) - 243.58501500480008*Cube
      (msq2(2,2))*Power6(M3Input)*Sqr(msu2(2,2)) + 44.70141739795201*Quad(M3Input)
      *Quad(msq2(2,2))*Sqr(msu2(2,2)) + 502.59768048000007*Power8(M3Input)*Sqr(
      msq2(2,2))*Sqr(msu2(2,2))))/(Cube(-1.02*msq2(2,2) + Sqr(M3Input))*(1.02*msq2
      (2,2) - 0.9800000000000001*msu2(2,2))*Quad(-0.9800000000000001*msu2(2,2) +
      Sqr(M3Input))) - (0.4444444444444444*(AtInput - MuInput/TanBeta)*(-
      16.320269280000005*Cube(M3Input)*Cube(msu2(2,2))*msq2(2,2) + 68*Power11(
      M3Input) + 2.823576000000001*Cube(msu2(2,2))*Power5(M3Input) +
      352.85880000000003*msq2(2,2)*msu2(2,2)*Power7(M3Input) - 139.74*msq2(2,2)*
      Power9(M3Input) - 179.34*msu2(2,2)*Power9(M3Input) + 11.750593881600004*
      M3Input*Cube(msu2(2,2))*Sqr(msq2(2,2)) - 187.60492800000003*msu2(2,2)*Power5
      (M3Input)*Sqr(msq2(2,2)) + 78.03*Power7(M3Input)*Sqr(msq2(2,2)) - 210.61572*
      msq2(2,2)*Power5(M3Input)*Sqr(msu2(2,2)) + 115.24800000000002*Power7(M3Input
      )*Sqr(msu2(2,2)) + 104.91601680000001*Cube(M3Input)*Sqr(msq2(2,2))*Sqr(msu2(
      2,2))))/(Cube(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))*(1.02*msq2(2,2)
      - 0.9800000000000001*msu2(2,2))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) + (
      13.937777777777779*msu2(2,2)*Quad(M3Input)*Sqr(AtInput - MuInput/TanBeta))/(
      (1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-1.02*msq2(2,2) + Sqr(
      M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + Log((1.02*
      msq2(2,2))/Sqr(SCALE))*((-0.2222222222222222*(64*Power(M3Input,16) - 183.6*
      Power(M3Input,14)*msq2(2,2) - 199.92000000000002*Power(M3Input,14)*msu2(2,2)
      - 32.897448000000004*Cube(msq2(2,2))*Power10(M3Input) - 95.06039200000002*
      Cube(msu2(2,2))*Power10(M3Input) + 603.7584*msq2(2,2)*msu2(2,2)*Power12(
      M3Input) + 6.234912115952488*Cube(msu2(2,2))*Power5(msq2(2,2)) -
      7.7285656224*Power5(msq2(2,2))*Power6(M3Input) + 384.96635184000013*Cube(
      msu2(2,2))*msq2(2,2)*Power8(M3Input) + 13.519789920000004*Cube(msq2(2,2))*
      msu2(2,2)*Power8(M3Input) - 52.93642543660802*Cube(msq2(2,2))*Cube(msu2(2,2)
      )*Quad(M3Input) + 14.065989432768003*msu2(2,2)*Power5(msq2(2,2))*Quad(
      M3Input) - 51.978392323200005*msu2(2,2)*Power6(M3Input)*Quad(msq2(2,2)) +
      32.4729648*Power8(M3Input)*Quad(msq2(2,2)) - 71.50197976320003*msq2(2,2)*
      Power6(M3Input)*Quad(msu2(2,2)) + 11.990786080000005*Power8(M3Input)*Quad(
      msu2(2,2)) - 5.9904057584641555*Quad(msq2(2,2))*Quad(msu2(2,2)) -
      6.112658937208322*Cube(msu2(2,2))*Quad(msq2(2,2))*Sqr(M3Input) +
      23.49178728809473*Cube(msq2(2,2))*Quad(msu2(2,2))*Sqr(M3Input) -
      376.22944800000005*msu2(2,2)*Power10(M3Input)*Sqr(msq2(2,2)) + 126.9288*
      Power12(M3Input)*Sqr(msq2(2,2)) - 234.03266147520006*Cube(msu2(2,2))*Power6(
      M3Input)*Sqr(msq2(2,2)) + 41.264168847552014*Quad(M3Input)*Quad(msu2(2,2))*
      Sqr(msq2(2,2)) - 763.1146320000001*msq2(2,2)*Power10(M3Input)*Sqr(msu2(2,2))
      + 224.73360000000002*Power12(M3Input)*Sqr(msu2(2,2)) - 13.249394121600005*
      Cube(msq2(2,2))*Power6(M3Input)*Sqr(msu2(2,2)) + 63.41363863430401*Quad(
      M3Input)*Quad(msq2(2,2))*Sqr(msu2(2,2)) - 19.086465661079043*Power5(msq2(2,2
      ))*Sqr(M3Input)*Sqr(msu2(2,2)) + 502.59768048000007*Power8(M3Input)*Sqr(msq2
      (2,2))*Sqr(msu2(2,2))))/(Cube(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))*
      (1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Quad(-1.02*msq2(2,2) + Sqr(
      M3Input))) - (14.506666666666666*msq2(2,2)*Quad(M3Input)*Sqr(AtInput -
      MuInput/TanBeta))/((1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input))*Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input))) + (0.4444444444444444*(AtInput - MuInput/TanBeta)*(-
      17.679725280000003*Cube(M3Input)*Cube(msq2(2,2))*msu2(2,2) + 68*Power11(
      M3Input) + 3.1836240000000005*Cube(msq2(2,2))*Power5(M3Input) +
      352.85880000000003*msq2(2,2)*msu2(2,2)*Power7(M3Input) - 186.66*msq2(2,2)*
      Power9(M3Input) - 134.26000000000002*msu2(2,2)*Power9(M3Input) -
      219.21228000000002*msu2(2,2)*Power5(M3Input)*Sqr(msq2(2,2)) + 124.848*Power7
      (M3Input)*Sqr(msq2(2,2)) + 12.230209958400003*M3Input*Cube(msq2(2,2))*Sqr(
      msu2(2,2)) - 180.24787200000003*msq2(2,2)*Power5(M3Input)*Sqr(msu2(2,2)) +
      72.03000000000002*Power7(M3Input)*Sqr(msu2(2,2)) + 104.91601680000001*Cube(
      M3Input)*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Cube(-1.02*msq2(2,2) + Sqr(M3Input
      ))*(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2))*Sqr(-0.9800000000000001*
      msu2(2,2) + Sqr(M3Input)))) + Log(Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1
      ,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)/Sqr(SCALE))*((
      0.1111111111111111*(204.*msq2(2,2)*Power10(M3Input) + 196.00000000000003*
      msu2(2,2)*Power10(M3Input) - 40*Power12(M3Input) + 21.224160000000005*Cube(
      msq2(2,2))*Power6(M3Input) + 18.823840000000004*Cube(msu2(2,2))*Power6(
      M3Input) - 839.6640000000001*msq2(2,2)*msu2(2,2)*Power8(M3Input) -
      134.40221760000006*Cube(msu2(2,2))*msq2(2,2)*Quad(M3Input) -
      145.59773760000004*Cube(msq2(2,2))*msu2(2,2)*Quad(M3Input) +
      428.22864000000004*msu2(2,2)*Power6(M3Input)*Sqr(msq2(2,2)) - 62.424*Power8(
      M3Input)*Sqr(msq2(2,2)) + 411.43536000000006*msq2(2,2)*Power6(M3Input)*Sqr(
      msu2(2,2)) - 57.62400000000001*Power8(M3Input)*Sqr(msu2(2,2))))/(Cube(-1.02*
      msq2(2,2) + Sqr(M3Input))*Cube(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))
      ) - (4.444444444444445*(AtInput - MuInput/TanBeta)*(-10.995600000000001*Cube
      (M3Input)*msq2(2,2)*msu2(2,2) + 5.1*msq2(2,2)*Power5(M3Input) + 4.9*msu2(2,2
      )*Power5(M3Input) + Power7(M3Input)))/(Sqr(-1.02*msq2(2,2) + Sqr(M3Input))*
      Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + (0.8888888888888888*(
      AtInput - MuInput/TanBeta)*(30.6*Cube(M3Input)*msq2(2,2)*Power(msd2(0,0)*
      msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666) +
      58.97640000000001*Cube(M3Input)*msq2(2,2)*msu2(2,2) + 29.400000000000002*
      Cube(M3Input)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1
      ,1),0.16666666666666666)*msu2(2,2) - 76.5*msq2(2,2)*Power5(M3Input) - 60*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Power5(M3Input) - 73.5*msu2(2,2)*Power5(M3Input) + 85*
      Power7(M3Input) + 3.1212*Cube(M3Input)*Sqr(msq2(2,2)) + 2.8812000000000006*
      Cube(M3Input)*Sqr(msu2(2,2))))/(Sqr(-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + PolyLog(2,1 - (
      1.040816326530612*msq2(2,2))/msu2(2,2))*((0.05555555555555555*(-24.48*msq2(2
      ,2)*Power10(M3Input) + 23.520000000000003*msu2(2,2)*Power10(M3Input) -
      6.491995122816001*msu2(2,2)*Power5(msq2(2,2)) + 5.531995276416001*msq2(2,2)*
      Power5(msu2(2,2)) - 72.16214400000001*Cube(msq2(2,2))*Power6(M3Input) +
      64.00105600000002*Cube(msu2(2,2))*Power6(M3Input) + 180.48297792000005*Cube(
      msu2(2,2))*msq2(2,2)*Quad(M3Input) - 195.51696192000003*Cube(msq2(2,2))*msu2
      (2,2)*Quad(M3Input) + 62.78106528*Quad(M3Input)*Quad(msq2(2,2)) -
      53.49735328000002*Quad(M3Input)*Quad(msu2(2,2)) - 19.8734544576*Power5(msq2(
      2,2))*Sqr(M3Input) + 16.270574342400003*Power5(msu2(2,2))*Sqr(M3Input) +
      74.25484617600002*msu2(2,2)*Quad(msq2(2,2))*Sqr(M3Input) - 65.85708662400003
      *msq2(2,2)*Quad(msu2(2,2))*Sqr(M3Input) + 224.31024000000002*msu2(2,2)*
      Power6(M3Input)*Sqr(msq2(2,2)) + 29.1312*Power8(M3Input)*Sqr(msq2(2,2)) +
      7.677054669312002*Quad(msu2(2,2))*Sqr(msq2(2,2)) - 15.667458508800005*Cube(
      msu2(2,2))*Sqr(M3Input)*Sqr(msq2(2,2)) - 215.51376000000005*msq2(2,2)*Power6
      (M3Input)*Sqr(msu2(2,2)) - 26.891200000000005*Power8(M3Input)*Sqr(msu2(2,2))
      - 8.316542771712001*Quad(msq2(2,2))*Sqr(msu2(2,2)) + 16.306946611200004*Cube
      (msq2(2,2))*Sqr(M3Input)*Sqr(msu2(2,2))))/(Cube(1.02*msq2(2,2) - Sqr(M3Input
      ))*Cube(-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (0.4444444444444444
      *(AtInput - MuInput/TanBeta)*(6.367248000000001*M3Input*Cube(msq2(2,2)) -
      5.647152000000002*M3Input*Cube(msu2(2,2)) + 10.2*msq2(2,2)*Power5(M3Input) -
      9.8*msu2(2,2)*Power5(M3Input) - 10.404*Cube(M3Input)*Sqr(msq2(2,2)) -
      8.156736*M3Input*msu2(2,2)*Sqr(msq2(2,2)) + 9.604000000000001*Cube(M3Input)*
      Sqr(msu2(2,2)) + 7.836864000000001*M3Input*msq2(2,2)*Sqr(msu2(2,2))))/(Sqr(-
      1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(
      M3Input)))) + Log(Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)/Sqr(SCALE))*((26.666666666666668*M3Input*(
      AtInput - MuInput/TanBeta)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666))/((-1.02*msq2(2,2) + Sqr(M3Input))*
      (-0.9800000000000001*msu2(2,2) + Sqr(M3Input))) - (1.1111111111111112*(-
      23.46*msq2(2,2)*Power6(M3Input) + 18*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*
      msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power6(M3Input) -
      22.540000000000003*msu2(2,2)*Power6(M3Input) + 12*Power8(M3Input) - 15.3*
      msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Quad(M3Input) + 43.982400000000005*msq2(2,2)*msu2(2,2)*
      Quad(M3Input) - 14.700000000000001*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(
      1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Quad(M3Input) -
      11.995200000000002*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Sqr(M3Input) + 3.058776*
      Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*msu2(2,2)*Sqr(msq2(2,2)) + 11.4444*Quad(M3Input)*Sqr(
      msq2(2,2)) + 9.3636*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*Sqr(M3Input)*Sqr(msq2(2,2)) - 21.411432*msu2(
      2,2)*Sqr(M3Input)*Sqr(msq2(2,2)) + 2.9388240000000003*msq2(2,2)*Power(msd2(0
      ,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*
      Sqr(msu2(2,2)) + 10.564400000000001*Quad(M3Input)*Sqr(msu2(2,2)) -
      20.571768000000002*msq2(2,2)*Sqr(M3Input)*Sqr(msu2(2,2)) + 8.643600000000001
      *Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Sqr(M3Input)*Sqr(msu2(2,2)) + 9.992001600000002*Sqr(
      msq2(2,2))*Sqr(msu2(2,2))))/(Sqr(-1.02*msq2(2,2) + Sqr(M3Input))*Sqr(-
      0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + (0.05555555555555555*(-
      5.760095040000002*Cube(msu2(2,2))*msq2(2,2) - 6.2399030400000015*Cube(msq2(2
      ,2))*msu2(2,2) - 287.64*msq2(2,2)*Power6(M3Input) - 360*Power(msd2(0,0)*msd2
      (1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Power6(
      M3Input) - 276.36*msu2(2,2)*Power6(M3Input) + 291*Power8(M3Input) + 306.*
      msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2(1,1),
      0.16666666666666666)*Quad(M3Input) - 111.9552*msq2(2,2)*msu2(2,2)*Quad(
      M3Input) + 294.*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*msu2
      (1,1),0.16666666666666666)*msu2(2,2)*Quad(M3Input) - 19.101744000000004*Cube
      (msq2(2,2))*Sqr(M3Input) - 16.941456000000006*Cube(msu2(2,2))*Sqr(M3Input) +
      239.90400000000002*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*msu2(2,2)*Sqr(M3Input) -
      61.175520000000006*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0)*
      msu2(1,1),0.16666666666666666)*msu2(2,2)*Sqr(msq2(2,2)) + 96.7572*Quad(
      M3Input)*Sqr(msq2(2,2)) - 187.272*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1
      ,1)*msu2(0,0)*msu2(1,1),0.16666666666666666)*Sqr(M3Input)*Sqr(msq2(2,2)) +
      212.07513600000001*msu2(2,2)*Sqr(M3Input)*Sqr(msq2(2,2)) -
      58.776480000000014*msq2(2,2)*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*
      msu2(0,0)*msu2(1,1),0.16666666666666666)*Sqr(msu2(2,2)) + 89.31720000000001*
      Quad(M3Input)*Sqr(msu2(2,2)) + 203.75846400000003*msq2(2,2)*Sqr(M3Input)*Sqr
      (msu2(2,2)) - 172.872*Power(msd2(0,0)*msd2(1,1)*msq2(0,0)*msq2(1,1)*msu2(0,0
      )*msu2(1,1),0.16666666666666666)*Sqr(M3Input)*Sqr(msu2(2,2)) -
      168.86482704000002*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Sqr(-1.02*msq2(2,2) +
      Sqr(M3Input))*Sqr(-0.9800000000000001*msu2(2,2) + Sqr(M3Input)))) + 3*Sqr((-
      0.6664000000000001*msq2(2,2)*msu2(2,2))/((Sqr(M3Input) - 1.02*msq2(2,2))*(
      Sqr(M3Input) - 0.9800000000000001*msu2(2,2))) + (0.6666666666666666*Quad(
      M3Input))/((Sqr(M3Input) - 1.02*msq2(2,2))*(Sqr(M3Input) -
      0.9800000000000001*msu2(2,2))) + Log((1.02*msq2(2,2))/Sqr(SCALE))*((-
      2.7199999999999998*M3Input*(AtInput - MuInput/TanBeta)*msq2(2,2))/((Sqr(
      M3Input) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.9800000000000001*msu2(2,2)))
      + (1.3599999999999999*Sqr(M3Input)*msq2(2,2))/Sqr(Sqr(M3Input) - 1.02*msq2(2
      ,2)) - (0.6936*Sqr(msq2(2,2)))/Sqr(Sqr(M3Input) - 1.02*msq2(2,2))) + Log((
      0.9800000000000001*msu2(2,2))/Sqr(SCALE))*((-2.6133333333333333*M3Input*(
      AtInput - MuInput/TanBeta)*msu2(2,2))/((Sqr(M3Input) - 0.9800000000000001*
      msu2(2,2))*(-1.02*msq2(2,2) + 0.9800000000000001*msu2(2,2))) + (
      1.3066666666666666*Sqr(M3Input)*msu2(2,2))/Sqr(Sqr(M3Input) -
      0.9800000000000001*msu2(2,2)) - (0.6402666666666668*Sqr(msu2(2,2)))/Sqr(Sqr(
      M3Input) - 0.9800000000000001*msu2(2,2))) + Log(Sqr(M3Input)/Sqr(SCALE))*((
      2.6666666666666665*(AtInput - MuInput/TanBeta)*Cube(M3Input))/((Sqr(M3Input)
      - 1.02*msq2(2,2))*(Sqr(M3Input) - 0.9800000000000001*msu2(2,2))) - (
      0.6666666666666666*(-2.04*msq2(2,2)*Power6(M3Input) - 1.9600000000000002*
      msu2(2,2)*Power6(M3Input) + 2*Power8(M3Input) + 1.0404*Sqr(msq2(2,2))*Quad(
      M3Input) + 0.9604000000000001*Sqr(msu2(2,2))*Quad(M3Input)))/(Sqr(Sqr(
      M3Input) - 1.02*msq2(2,2))*Sqr(Sqr(M3Input) - 0.9800000000000001*msu2(2,2)))
      )))))/Power6(3.141592653589793)), 0) + IF(Abs(Qmatch) > 1, WHICH(
      LambdaLoopOrder < 1, (0.0625*(0.01*(27*Log(MSUSY/Qmatch)*Quad(g1) + 225*Log(
      MSUSY/Qmatch)*Quad(g2) + 27*Log(MSUSY/Qmatch)*Quad(g1)*Quad(Cos(2*ArcTan(
      TanBeta))) + 75*Log(MSUSY/Qmatch)*Quad(g2)*Quad(Cos(2*ArcTan(TanBeta))) -
      1200*Log(MSUSY/Qmatch)*Quad(Yd(2,2)) - 400*Log(MSUSY/Qmatch)*Quad(Ye(2,2)) -
      1200*Log(MSUSY/Qmatch)*Quad(Yu(2,2)) + 90*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)
      + 90*Log(MSUSY/Qmatch)*Quad(Cos(2*ArcTan(TanBeta)))*Sqr(g1)*Sqr(g2) - 27*Log
      (MSUSY/Qmatch)*Quad(g1)*Sqr(Cos(2*ArcTan(TanBeta))) - 225*Log(MSUSY/Qmatch)*
      Quad(g2)*Sqr(Cos(2*ArcTan(TanBeta))) - 180*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)
      *Sqr(Cos(2*ArcTan(TanBeta))) + 180*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Yd(2,2)) + 300*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Yd(2,2)) + 60*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Ye(2,2)) + 100*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Ye(2,2)) + 180*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Yu(2,2)) + 300*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Yu(2,2))) + 0.25*(-((4.92*Log(MSUSY/Qmatch)*Quad(g1) -
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta)))) -
      (0.6*Sqr(g1) + Sqr(g2))*(-12*TanBeta*Cos(2*ArcTan(TanBeta))*Log(MSUSY/Qmatch
      )*Sin(2*ArcTan(TanBeta))*Sqr(Yd(2,2)) - 4*TanBeta*Cos(2*ArcTan(TanBeta))*Log
      (MSUSY/Qmatch)*Sin(2*ArcTan(TanBeta))*Sqr(Ye(2,2)) + (12*Cos(2*ArcTan(
      TanBeta))*Log(MSUSY/Qmatch)*Sin(2*ArcTan(TanBeta))*Sqr(Yu(2,2)))/TanBeta))))
      /Sqr(3.141592653589793), LambdaLoopOrder < 2, (0.00390625*(1.971*Power6(g1)*
      Sqr(Log(MSUSY/Qmatch)) - 24.375*Power6(g2)*Sqr(Log(MSUSY/Qmatch)) + 2.25*
      Cube(0.6*Sqr(g1) + Sqr(g2))*Power6(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/
      Qmatch)) - 180*Power6(Yd(2,2))*Sqr(Log(MSUSY/Qmatch)) - 28*Power6(Ye(2,2))*
      Sqr(Log(MSUSY/Qmatch)) - 180*Power6(Yu(2,2))*Sqr(Log(MSUSY/Qmatch)) - 8.925*
      Quad(g2)*Sqr(g1)*Sqr(Log(MSUSY/Qmatch)) + 16.8*Quad(Yd(2,2))*Sqr(g1)*Sqr(Log
      (MSUSY/Qmatch)) + 21.6*Quad(Ye(2,2))*Sqr(g1)*Sqr(Log(MSUSY/Qmatch)) + 31.2*
      Quad(Yu(2,2))*Sqr(g1)*Sqr(Log(MSUSY/Qmatch)) + 1.665*Quad(g1)*Sqr(g2)*Sqr(
      Log(MSUSY/Qmatch)) + 108*Quad(Yd(2,2))*Sqr(g2)*Sqr(Log(MSUSY/Qmatch)) + 36*
      Quad(Ye(2,2))*Sqr(g2)*Sqr(Log(MSUSY/Qmatch)) + 108*Quad(Yu(2,2))*Sqr(g2)*Sqr
      (Log(MSUSY/Qmatch)) - 2.025*Quad(Cos(2*ArcTan(TanBeta)))*Sqr(g1)*Sqr(0.6*Sqr
      (g1) + Sqr(g2))*Sqr(Log(MSUSY/Qmatch)) - 10.125*Quad(Cos(2*ArcTan(TanBeta)))
      *Sqr(g2)*Sqr(0.6*Sqr(g1) + Sqr(g2))*Sqr(Log(MSUSY/Qmatch)) + 192*Quad(Yd(2,2
      ))*Sqr(g3)*Sqr(Log(MSUSY/Qmatch)) + 192*Quad(Yu(2,2))*Sqr(g3)*Sqr(Log(MSUSY/
      Qmatch)) - 0.63*Quad(g1)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))
      *Sqr(Log(MSUSY/Qmatch)) + 24*Quad(g2)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch)) - 4.5*Quad(Yd(2,2))*(0.6*Sqr(g1) +
      Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch)) - 7.5*Quad(Ye(2,
      2))*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch
      )) - 4.5*Quad(Yu(2,2))*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*
      Sqr(Log(MSUSY/Qmatch)) + 6.75*Sqr(g1)*Sqr(g2)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(
      Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch)) + 1.62*Quad(g1)*Sqr(Log(MSUSY
      /Qmatch))*Sqr(Yd(2,2)) + 13.5*Quad(g2)*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) -
      48*Quad(Ye(2,2))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) - 108*Quad(Yu(2,2))*Sqr
      (Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) + 5.4*Sqr(g1)*Sqr(g2)*Sqr(Log(MSUSY/Qmatch)
      )*Sqr(Yd(2,2)) + 13.5*Quad(Cos(2*ArcTan(TanBeta)))*Sqr(0.6*Sqr(g1) + Sqr(g2)
      )*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) - 6.15*Sqr(g1)*(0.6*Sqr(g1) + Sqr(g2))
      *Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) - 33.75*Sqr
      (g2)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yd(2,2)) - 24*(0.6*Sqr(g1) + Sqr(g2))*Sqr(g3)*Sqr(Cos(2*ArcTan(
      TanBeta)))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) + 0.54*Quad(g1)*Sqr(Log(MSUSY
      /Qmatch))*Sqr(Ye(2,2)) + 4.5*Quad(g2)*Sqr(Log(MSUSY/Qmatch))*Sqr(Ye(2,2)) -
      48*Quad(Yd(2,2))*Sqr(Log(MSUSY/Qmatch))*Sqr(Ye(2,2)) - 48*Quad(Yu(2,2))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Ye(2,2)) + 1.8*Sqr(g1)*Sqr(g2)*Sqr(Log(MSUSY/Qmatch))
      *Sqr(Ye(2,2)) + 4.5*Quad(Cos(2*ArcTan(TanBeta)))*Sqr(0.6*Sqr(g1) + Sqr(g2))*
      Sqr(Log(MSUSY/Qmatch))*Sqr(Ye(2,2)) - 4.05*Sqr(g1)*(0.6*Sqr(g1) + Sqr(g2))*
      Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch))*Sqr(Ye(2,2)) - 11.25*Sqr(
      g2)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch
      ))*Sqr(Ye(2,2)) + 18*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr
      (Log(MSUSY/Qmatch))*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + 1.62*Quad(g1)*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)) + 13.5*Quad(g2)*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)) -
      108*Quad(Yd(2,2))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)) - 48*Quad(Ye(2,2))*Sqr
      (Log(MSUSY/Qmatch))*Sqr(Yu(2,2)) + 5.4*Sqr(g1)*Sqr(g2)*Sqr(Log(MSUSY/Qmatch)
      )*Sqr(Yu(2,2)) + 13.5*Quad(Cos(2*ArcTan(TanBeta)))*Sqr(0.6*Sqr(g1) + Sqr(g2)
      )*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)) - 7.95*Sqr(g1)*(0.6*Sqr(g1) + Sqr(g2))
      *Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)) - 33.75*Sqr
      (g2)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)) - 24*(0.6*Sqr(g1) + Sqr(g2))*Sqr(g3)*Sqr(Cos(2*ArcTan(
      TanBeta)))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)) + 45*(0.6*Sqr(g1) + Sqr(g2))*
      Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2))*Sqr(Yu(2,2))
      + 18*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Ye(2,2))*Sqr(Yu(2,2)) - 16*Sqr(3.141592653589793)*((0.0625*(-
      1.476*Log(MSUSY/Qmatch)*Power6(g1) + 1.9*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(g1)
      - 2.46*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(g2) + 12.666666666666666*Log(MSUSY/
      Qmatch)*Power6(g2)*(0.75 - 0.16666666666666666*Sqr(Cos(2*ArcTan(TanBeta))))
      - (0.6666666666666666*Cos(2*ArcTan(TanBeta))*Quad(g2)*Sin(2*ArcTan(TanBeta))
      *(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta))/Sqr(
      3.141592653589793) + (0.010416666666666666*((4*Cos(2*ArcTan(TanBeta))*(2*Log
      (Sqr(M2Input)/Sqr(Qmatch))*Quad(g2) + Log(Sqr(MuInput)/Sqr(Qmatch))*(0.36*
      Quad(g1) + Quad(g2)))*Sin(2*ArcTan(TanBeta))*(3*Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta - Sqr(Cos(2*ArcTan(TanBeta)))*(-
      25.333333333333332*Log(MSUSY/Qmatch)*Log(Sqr(M2Input)/Sqr(Qmatch))*Power6(g2
      ) + Log(Sqr(MuInput)/Sqr(Qmatch))*(5.904*Log(MSUSY/Qmatch)*Power6(g1) -
      12.666666666666666*Log(MSUSY/Qmatch)*Power6(g2)) + (0.36*Quad(g1) + Quad(g2)
      )*(-1.2*Log(MSUSY/Qmatch)*Sqr(g1) - 6*Log(MSUSY/Qmatch)*Sqr(g2) + 6*Log(
      MSUSY/Qmatch)*Sqr(Yd(2,2)) + 6*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) +
      2*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) + 2*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2
      ,2)) + 6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + (6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))
      /Sqr(TanBeta)))))/Sqr(3.141592653589793) + (0.0625*(0.5*((0.3*(0.6*Sqr(g1) +
      Sqr(g2))*(Sqr(g1) + 5*Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta))) - (0.36*Quad(g1))
      /Sqr(1 + Sqr(TanBeta)) - (5*Quad(g2))/Sqr(1 + Sqr(TanBeta)) - (0.36*Quad(g1)
      *Quad(TanBeta))/Sqr(1 + Sqr(TanBeta)) - (5*Quad(g2)*Quad(TanBeta))/Sqr(1 +
      Sqr(TanBeta)) - (2.4*Sqr(g1)*Sqr(g2)*Sqr(TanBeta))/Sqr(1 + Sqr(TanBeta)) - (
      0.08*(5*Sqr(g2) + 3*Sqr(g1)*Sqr(TanBeta))*(3*Sqr(g1) + 5*Sqr(g2)*Sqr(TanBeta
      )))/Sqr(1 + Sqr(TanBeta)))*(-1.2*Log(MSUSY/Qmatch)*Sqr(g1) - 6*Log(MSUSY/
      Qmatch)*Sqr(g2) + 6*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + 6*Log(MSUSY/Qmatch)*Sqr
      (TanBeta)*Sqr(Yd(2,2)) + 2*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) + 2*Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) + 6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + (6*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta)) + Log(Sqr(MuInput)/Sqr(Qmatch)
      )*((-0.144*(41*Log(MSUSY/Qmatch)*Power6(g1) - 30*Log(MSUSY/Qmatch)*Quad(g1)*
      Sqr(TanBeta)*Sqr(Yd(2,2)) - 10*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(TanBeta)*Sqr(
      Ye(2,2)) + 30*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(Yu(2,2))))/Sqr(1 + Sqr(TanBeta)
      ) + (3.3333333333333335*(19*Log(MSUSY/Qmatch)*Power6(g2) + 18*Log(MSUSY/
      Qmatch)*Quad(g2)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 6*Log(MSUSY/Qmatch)*Quad(g2)*
      Sqr(TanBeta)*Sqr(Ye(2,2)) - 18*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(Yu(2,2))))/Sqr
      (1 + Sqr(TanBeta)) - (0.144*(41*Log(MSUSY/Qmatch)*Power6(g1)*Quad(TanBeta) +
      30*Log(MSUSY/Qmatch)*Quad(g1)*Quad(TanBeta)*Sqr(Yd(2,2)) + 10*Log(MSUSY/
      Qmatch)*Quad(g1)*Quad(TanBeta)*Sqr(Ye(2,2)) - 30*Log(MSUSY/Qmatch)*Quad(g1)*
      Sqr(TanBeta)*Sqr(Yu(2,2))))/Sqr(1 + Sqr(TanBeta)) + (3.3333333333333335*(19*
      Log(MSUSY/Qmatch)*Power6(g2)*Quad(TanBeta) - 18*Log(MSUSY/Qmatch)*Quad(g2)*
      Quad(TanBeta)*Sqr(Yd(2,2)) - 6*Log(MSUSY/Qmatch)*Quad(g2)*Quad(TanBeta)*Sqr(
      Ye(2,2)) + 18*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Sqr(1 +
      Sqr(TanBeta)) - (0.16*(-95*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(g1)*Sqr(TanBeta) +
      123*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(g2)*Sqr(TanBeta) - 90*Log(MSUSY/Qmatch)*
      Quad(TanBeta)*Sqr(g1)*Sqr(g2)*Sqr(Yd(2,2)) + 90*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(g2)*Sqr(TanBeta)*Sqr(Yd(2,2)) - 30*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(
      g1)*Sqr(g2)*Sqr(Ye(2,2)) + 30*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)*Sqr(TanBeta)
      *Sqr(Ye(2,2)) - 90*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)*Sqr(Yu(2,2)) + 90*Log(
      MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Sqr(1 + Sqr(
      TanBeta)) - 2*((0.0026666666666666666*(5*Sqr(g2) + 3*Sqr(g1)*Sqr(TanBeta))*(
      369*Log(MSUSY/Qmatch)*Quad(g1) - 475*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(TanBeta)
      - 270*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 450*Log(MSUSY/
      Qmatch)*Sqr(g2)*Sqr(TanBeta)*Sqr(Yd(2,2)) - 90*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (TanBeta)*Sqr(Ye(2,2)) + 150*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta)*Sqr(Ye(2
      ,2)) + 270*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(Yu(2,2)) - 450*Log(MSUSY/Qmatch)*
      Sqr(g2)*Sqr(Yu(2,2))))/Sqr(1 + Sqr(TanBeta)) + (0.0026666666666666666*(3*Sqr
      (g1) + 5*Sqr(g2)*Sqr(TanBeta))*(-475*Log(MSUSY/Qmatch)*Quad(g2) + 369*Log(
      MSUSY/Qmatch)*Quad(g1)*Sqr(TanBeta) + 270*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(
      TanBeta)*Sqr(Yd(2,2)) - 450*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta)*Sqr(Yd(2,
      2)) + 90*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 150*Log(MSUSY
      /Qmatch)*Sqr(g2)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 270*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(Yu(2,2)) + 450*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(Yu(2,2))))/Sqr(1 + Sqr(
      TanBeta))) + 0.5*(0.04*(123*Log(MSUSY/Qmatch)*Quad(g1) - 475*Log(MSUSY/
      Qmatch)*Quad(g2))*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta))) + 0.6*
      (Sqr(g1) + 5*Sqr(g2))*((4.92*Log(MSUSY/Qmatch)*Quad(g1) - 6.333333333333333*
      Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta))) - (4*Cos(2*ArcTan(
      TanBeta))*Sin(2*ArcTan(TanBeta))*(0.6*Sqr(g1) + Sqr(g2))*(3*Log(MSUSY/Qmatch
      )*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) -
      3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta)))) + 0.4*((TanBeta*Sqr(g1)*(-
      9.84*Log(MSUSY/Qmatch)*Quad(g1) + 0.25*((4.92*Log(MSUSY/Qmatch)*Quad(g1) -
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta))) -
      (4*Cos(2*ArcTan(TanBeta))*Sin(2*ArcTan(TanBeta))*(0.6*Sqr(g1) + Sqr(g2))*(3*
      Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)
      *Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta))*TCf0(M1Input/
      MuInput))/(1 + Sqr(TanBeta)) + (-1.2*Sqr(g1) + 0.25*(0.6*Sqr(g1) + Sqr(g2))*
      Sqr(Cos(2*ArcTan(TanBeta))))*((-2*TanBeta*Sqr(g1)*(3*Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2)))*TCf0(M1Input/MuInput))/(1 + Sqr(TanBeta)) + (Sqr
      (g1)*(3*TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) + TanBeta*
      Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Ye(2,2)) - (3*Log(MSUSY/Qmatch)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta)*TCf0(M1Input/MuInput) + TanBeta*((
      M1Input*Sqr(g1)*(0.6*Log(MSUSY/Qmatch)*Sqr(g1) + 3*Log(MSUSY/Qmatch)*Sqr(g2)
      - 3*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd
      (2,2)) - Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) - Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr
      (Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2
      ,2)))/Sqr(TanBeta))*TCD1f0(M1Input/MuInput))/MuInput + 8.2*Log(MSUSY/Qmatch)
      *Quad(g1)*TCf0(M1Input/MuInput)))/(1 + Sqr(TanBeta)))) + 2*((TanBeta*Sqr(g2)
      *(12.666666666666666*Log(MSUSY/Qmatch)*Quad(g2) + 0.25*((4.92*Log(MSUSY/
      Qmatch)*Quad(g1) - 6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*
      ArcTan(TanBeta))) - (4*Cos(2*ArcTan(TanBeta))*Sin(2*ArcTan(TanBeta))*(0.6*
      Sqr(g1) + Sqr(g2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))
      /TanBeta))*TCf0(M2Input/MuInput))/(1 + Sqr(TanBeta)) + (-2*Sqr(g2) + 0.25*(
      0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta))))*((-2*TanBeta*Sqr(g2)*(3*
      Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)
      *Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))*TCf0(M2Input/MuInput))/(1
      + Sqr(TanBeta)) + (Sqr(g2)*(3*TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*
      Sqr(Yd(2,2)) + TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Ye(2,2)) - (
      3*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta)*TCf0(M2Input/
      MuInput) + TanBeta*((M2Input*Sqr(g2)*(0.6*Log(MSUSY/Qmatch)*Sqr(g1) + 3*Log(
      MSUSY/Qmatch)*Sqr(g2) - 3*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) - 3*Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) - Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) - Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) -
      (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta))*TCD1f0(M2Input/MuInput))/
      MuInput - 6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2)*TCf0(M2Input/MuInput)
      ))/(1 + Sqr(TanBeta)))) + 0.08333333333333333*(4.92*Log(MSUSY/Qmatch)*Quad(
      g1)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*TCg0(M1Input/MuInput
      ) + 0.6*Sqr(g1)*((M1Input*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta))
      )*(0.6*Log(MSUSY/Qmatch)*Sqr(g1) + 3*Log(MSUSY/Qmatch)*Sqr(g2) - 3*Log(MSUSY
      /Qmatch)*Sqr(Yd(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) - Log(
      MSUSY/Qmatch)*Sqr(Ye(2,2)) - Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3
      *Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(
      TanBeta))*TCD1g0(M1Input/MuInput))/MuInput + ((4.92*Log(MSUSY/Qmatch)*Quad(
      g1) - 6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta
      ))) - (4*Cos(2*ArcTan(TanBeta))*Sin(2*ArcTan(TanBeta))*(0.6*Sqr(g1) + Sqr(g2
      ))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta)*TCg0(
      M1Input/MuInput))) + 0.25*(-6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2)*(
      0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*TCg0(M2Input/MuInput) +
      Sqr(g2)*((M2Input*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*(0.6*
      Log(MSUSY/Qmatch)*Sqr(g1) + 3*Log(MSUSY/Qmatch)*Sqr(g2) - 3*Log(MSUSY/Qmatch
      )*Sqr(Yd(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) - Log(MSUSY/
      Qmatch)*Sqr(Ye(2,2)) - Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      )*TCD1g0(M2Input/MuInput))/MuInput + ((4.92*Log(MSUSY/Qmatch)*Quad(g1) -
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta))) -
      (4*Cos(2*ArcTan(TanBeta))*Sin(2*ArcTan(TanBeta))*(0.6*Sqr(g1) + Sqr(g2))*(3*
      Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)
      *Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta)*TCg0(M2Input/
      MuInput))) - 0.5833333333333334*((0.36*M1Input*(Quad(g1) + Quad(g1)*Quad(
      TanBeta))*(0.6*Log(MSUSY/Qmatch)*Sqr(g1) + 3*Log(MSUSY/Qmatch)*Sqr(g2) - 3*
      Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2
      )) - Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) - Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(
      2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))
      )/Sqr(TanBeta))*TCD1f(1)(M1Input/MuInput))/(MuInput*Sqr(1 + Sqr(TanBeta))) +
      (0.144*(41*Log(MSUSY/Qmatch)*Power6(g1) + 41*Log(MSUSY/Qmatch)*Power6(g1)*
      Quad(TanBeta) + 30*Log(MSUSY/Qmatch)*Quad(g1)*Quad(TanBeta)*Sqr(Yd(2,2)) -
      30*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*Log(MSUSY/
      Qmatch)*Quad(g1)*Quad(TanBeta)*Sqr(Ye(2,2)) - 10*Log(MSUSY/Qmatch)*Quad(g1)*
      Sqr(TanBeta)*Sqr(Ye(2,2)) + 30*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(Yu(2,2)) - 30*
      Log(MSUSY/Qmatch)*Quad(g1)*Sqr(TanBeta)*Sqr(Yu(2,2)))*TCf(1)(M1Input/MuInput
      ))/Sqr(1 + Sqr(TanBeta))) - 2.25*((M2Input*(Quad(g2) + Quad(g2)*Quad(TanBeta
      ))*(0.6*Log(MSUSY/Qmatch)*Sqr(g1) + 3*Log(MSUSY/Qmatch)*Sqr(g2) - 3*Log(
      MSUSY/Qmatch)*Sqr(Yd(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) -
      Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) - Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2))
      - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(
      TanBeta))*TCD1f(2)(M2Input/MuInput))/(MuInput*Sqr(1 + Sqr(TanBeta))) - (
      0.6666666666666666*(19*Log(MSUSY/Qmatch)*Power6(g2) + 19*Log(MSUSY/Qmatch)*
      Power6(g2)*Quad(TanBeta) - 18*Log(MSUSY/Qmatch)*Quad(g2)*Quad(TanBeta)*Sqr(
      Yd(2,2)) + 18*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(TanBeta)*Sqr(Yd(2,2)) - 6*Log(
      MSUSY/Qmatch)*Quad(g2)*Quad(TanBeta)*Sqr(Ye(2,2)) + 6*Log(MSUSY/Qmatch)*Quad
      (g2)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 18*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(Yu(2,2))
      + 18*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(TanBeta)*Sqr(Yu(2,2)))*TCf(2)(M2Input/
      MuInput))/Sqr(1 + Sqr(TanBeta))) - 0.54*((-4*Quad(g1)*Sqr(TanBeta)*(3*Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr
      (Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))*TCf(3)(M1Input/MuInput))/Sqr(1
       + Sqr(TanBeta)) + (2*TanBeta*Quad(g1)*(3*TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr
      (TanBeta))*Sqr(Yd(2,2)) + TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(
      Ye(2,2)) - (3*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta)*
      TCf(3)(M1Input/MuInput) + Sqr(TanBeta)*((M1Input*Quad(g1)*(0.6*Log(MSUSY/
      Qmatch)*Sqr(g1) + 3*Log(MSUSY/Qmatch)*Sqr(g2) - 3*Log(MSUSY/Qmatch)*Sqr(Yd(2
      ,2)) - 3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) - Log(MSUSY/Qmatch)*Sqr
      (Ye(2,2)) - Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch
      )*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta))*TCD1f(3)(
      M1Input/MuInput))/MuInput + 16.4*Log(MSUSY/Qmatch)*Power6(g1)*TCf(3)(M1Input
      /MuInput)))/Sqr(1 + Sqr(TanBeta))) - 3.5*((-4*Quad(g2)*Sqr(TanBeta)*(3*Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr
      (Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))*TCf(4)(M2Input/MuInput))/Sqr(1
       + Sqr(TanBeta)) + (2*TanBeta*Quad(g2)*(3*TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr
      (TanBeta))*Sqr(Yd(2,2)) + TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(
      Ye(2,2)) - (3*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta)*
      TCf(4)(M2Input/MuInput) + Sqr(TanBeta)*((M2Input*Quad(g2)*(0.6*Log(MSUSY/
      Qmatch)*Sqr(g1) + 3*Log(MSUSY/Qmatch)*Sqr(g2) - 3*Log(MSUSY/Qmatch)*Sqr(Yd(2
      ,2)) - 3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) - Log(MSUSY/Qmatch)*Sqr
      (Ye(2,2)) - Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch
      )*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta))*TCD1f(4)(
      M2Input/MuInput))/MuInput - 12.666666666666666*Log(MSUSY/Qmatch)*Power6(g2)*
      TCf(4)(M2Input/MuInput)))/Sqr(1 + Sqr(TanBeta))) - 1.6*((-4*Sqr(g1)*Sqr(g2)*
      Sqr(TanBeta)*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))*TCf(5)
      (M1Input/MuInput,M2Input/MuInput))/Sqr(1 + Sqr(TanBeta)) + (2*TanBeta*Sqr(g1
      )*Sqr(g2)*(3*TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) +
      TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Ye(2,2)) - (3*Log(MSUSY/
      Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta)*TCf(5)(M1Input/MuInput,
      M2Input/MuInput) + Sqr(TanBeta)*((0.2*Sqr(g1)*Sqr(g2)*(3*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(TanBeta) + 15*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta) - 15*Log(
      MSUSY/Qmatch)*Quad(TanBeta)*Sqr(Yd(2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(TanBeta)
      *Sqr(Yd(2,2)) - 5*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(Ye(2,2)) - 5*Log(MSUSY
      /Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - 15*
      Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2)))*(M2Input*TCD01f(5)(M1Input/
      MuInput,M2Input/MuInput) + M1Input*TCD10f(5)(M1Input/MuInput,M2Input/MuInput
      )))/(MuInput*Sqr(TanBeta)) + (-6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2)*
      Sqr(g1) + 8.2*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(g2))*TCf(5)(M1Input/MuInput,
      M2Input/MuInput)))/Sqr(1 + Sqr(TanBeta))) - 1.1666666666666667*((0.12*(Sqr(
      g1)*Sqr(g2) + Quad(TanBeta)*Sqr(g1)*Sqr(g2))*(3*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(TanBeta) + 15*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta) - 15*Log(MSUSY/
      Qmatch)*Quad(TanBeta)*Sqr(Yd(2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(
      Yd(2,2)) - 5*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(Ye(2,2)) - 5*Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - 15*
      Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2)))*(M2Input*TCD01f(6)(M1Input/
      MuInput,M2Input/MuInput) + M1Input*TCD10f(6)(M1Input/MuInput,M2Input/MuInput
      )))/(MuInput*Sqr(TanBeta)*Sqr(1 + Sqr(TanBeta))) + (0.04*(-95*Log(MSUSY/
      Qmatch)*Quad(g2)*Sqr(g1) - 95*Log(MSUSY/Qmatch)*Quad(g2)*Quad(TanBeta)*Sqr(
      g1) + 123*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(g2) + 123*Log(MSUSY/Qmatch)*Quad(g1
      )*Quad(TanBeta)*Sqr(g2) + 180*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(g1)*Sqr(g2
      )*Sqr(Yd(2,2)) - 180*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)*Sqr(TanBeta)*Sqr(Yd(2
      ,2)) + 60*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(g1)*Sqr(g2)*Sqr(Ye(2,2)) - 60*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)*Sqr(TanBeta)*Sqr(Ye(2,2)) + 180*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(g2)*Sqr(Yu(2,2)) - 180*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(g2)
      *Sqr(TanBeta)*Sqr(Yu(2,2)))*TCf(6)(M1Input/MuInput,M2Input/MuInput))/Sqr(1 +
      Sqr(TanBeta))) + 0.2*((4*Sqr(g1)*Sqr(g2)*Sqr(TanBeta)*(3*Log(MSUSY/Qmatch)*
      Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))*TCf(7)(M1Input/MuInput,M2Input/MuInput))/Sqr
      (1 + Sqr(TanBeta)) - (2*TanBeta*Sqr(g1)*Sqr(g2)*(3*TanBeta*Log(MSUSY/Qmatch)
      *(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) + TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(
      TanBeta))*Sqr(Ye(2,2)) - (3*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)
      ))/TanBeta)*TCf(7)(M1Input/MuInput,M2Input/MuInput) + Sqr(TanBeta)*((0.2*Sqr
      (g1)*Sqr(g2)*(3*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta) + 15*Log(MSUSY/Qmatch
      )*Sqr(g2)*Sqr(TanBeta) - 15*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(Yd(2,2)) -
      15*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) - 5*Log(MSUSY/Qmatch)*Quad(
      TanBeta)*Sqr(Ye(2,2)) - 5*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 15*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,
      2)))*(M2Input*TCD01f(7)(M1Input/MuInput,M2Input/MuInput) + M1Input*TCD10f(7)
      (M1Input/MuInput,M2Input/MuInput)))/(MuInput*Sqr(TanBeta)) + (-
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(g1) + 8.2*Log(MSUSY/Qmatch)
      *Quad(g1)*Sqr(g2))*TCf(7)(M1Input/MuInput,M2Input/MuInput)))/Sqr(1 + Sqr(
      TanBeta))) - 2.0655911179772892*((0.025819888974716116*g1*g2*TanBeta*(123*g2
      *Cube(g1)*Log(MSUSY/Qmatch) - 95*g1*Cube(g2)*Log(MSUSY/Qmatch))*TCf(8)(
      M1Input/MuInput,M2Input/MuInput))/(1 + Sqr(TanBeta)) + 0.7745966692414834*g1
      *g2*((-2*g1*g2*TanBeta*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))*
      TCf(8)(M1Input/MuInput,M2Input/MuInput))/(1 + Sqr(TanBeta)) + (g1*g2*(3*
      TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) + TanBeta*Log(
      MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Ye(2,2)) - (3*Log(MSUSY/Qmatch)*(1 +
      Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta)*TCf(8)(M1Input/MuInput,M2Input/MuInput)
      + TanBeta*((0.2*g1*g2*(3*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta) + 15*Log(
      MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta) - 15*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(
      Yd(2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) - 5*Log(MSUSY/
      Qmatch)*Quad(TanBeta)*Sqr(Ye(2,2)) - 5*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye
      (2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - 15*Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Yu(2,2)))*(M2Input*TCD01f(8)(M1Input/MuInput,M2Input/MuInput) +
      M1Input*TCD10f(8)(M1Input/MuInput,M2Input/MuInput)))/(MuInput*Sqr(TanBeta))
      + (4.1*g2*Cube(g1)*Log(MSUSY/Qmatch) - 3.1666666666666665*g1*Cube(g2)*Log(
      MSUSY/Qmatch))*TCf(8)(M1Input/MuInput,M2Input/MuInput)))/(1 + Sqr(TanBeta)))
      )))/Sqr(3.141592653589793) + (0.0625*(Log(mse2(2,2)/Sqr(Qmatch))*(0.5*(-0.6*
      Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(Ye(2,2)))*Ye(2,2)*(10*Cube(Ye(2,2))*Log
      (MSUSY/Qmatch) - 9*Log(MSUSY/Qmatch)*Sqr(g1)*Ye(2,2) - 9*Log(MSUSY/Qmatch)*
      Sqr(g2)*Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Ye(2,2) + 12*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2))*Ye(2,2)) + Sqr(Ye(2,2))*(-0.6*(8.2*Cos(2*ArcTan(TanBeta
      ))*Log(MSUSY/Qmatch)*Quad(g1) - (2*Sin(2*ArcTan(TanBeta))*Sqr(g1)*(3*Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr
      (Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta) + 0.5*Ye(2,2)*(10*
      Cube(Ye(2,2))*Log(MSUSY/Qmatch) - 9*Log(MSUSY/Qmatch)*Sqr(g1)*Ye(2,2) - 9*
      Log(MSUSY/Qmatch)*Sqr(g2)*Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Ye(2,2
      ) + 12*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Ye(2,2)))) + Log(msl2(2,2)/Sqr(Qmatch)
      )*(0.5*(-0.5*Cos(2*ArcTan(TanBeta))*(-0.6*Sqr(g1) + Sqr(g2)) + Sqr(Ye(2,2)))
      *Ye(2,2)*(10*Cube(Ye(2,2))*Log(MSUSY/Qmatch) - 9*Log(MSUSY/Qmatch)*Sqr(g1)*
      Ye(2,2) - 9*Log(MSUSY/Qmatch)*Sqr(g2)*Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yd(
      2,2))*Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Ye(2,2)) + Sqr(Ye(2,2))*(
      0.5*(-(Cos(2*ArcTan(TanBeta))*(-4.92*Log(MSUSY/Qmatch)*Quad(g1) -
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))) + (2*Sin(2*ArcTan(TanBeta))*(
      -0.6*Sqr(g1) + Sqr(g2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log
      (MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))
      )/TanBeta) + 0.5*Ye(2,2)*(10*Cube(Ye(2,2))*Log(MSUSY/Qmatch) - 9*Log(MSUSY/
      Qmatch)*Sqr(g1)*Ye(2,2) - 9*Log(MSUSY/Qmatch)*Sqr(g2)*Ye(2,2) + 12*Log(MSUSY
      /Qmatch)*Sqr(Yd(2,2))*Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Ye(2,2))))
      + (2*((0.03333333333333333*TanBeta*Cube(AtauInput - MuInput*TanBeta)*Quad(Ye
      (2,2))*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*MuInput*Log(MSUSY/Qmatch)*
      Sqr(g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + 30*MuInput*Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)
      ) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)))*TCF(2)(Sqrt(msl2
      (2,2)/mse2(2,2))))/Sqrt(mse2(2,2)*msl2(2,2)) + (-0.4*TanBeta*(AtauInput -
      MuInput*TanBeta)*Quad(Ye(2,2))*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*
      MuInput*Log(MSUSY/Qmatch)*Sqr(g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)
      ) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(
      MSUSY/Qmatch)*Sqr(Ye(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(
      Ye(2,2))) + Cube(Ye(2,2))*Sqr(AtauInput - MuInput*TanBeta)*(10*Cube(Ye(2,2))
      *Log(MSUSY/Qmatch) - 9*Log(MSUSY/Qmatch)*Sqr(g1)*Ye(2,2) - 9*Log(MSUSY/
      Qmatch)*Sqr(g2)*Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Ye(2,2) + 12*Log
      (MSUSY/Qmatch)*Sqr(Yu(2,2))*Ye(2,2)))*(TCF(1)(Sqrt(msl2(2,2)/mse2(2,2))) - (
      0.08333333333333333*Sqr(AtauInput - MuInput*TanBeta)*TCF(2)(Sqrt(msl2(2,2)/
      mse2(2,2))))/Sqrt(mse2(2,2)*msl2(2,2)))))/Sqrt(mse2(2,2)*msl2(2,2)) + (0.25*
      (-0.4*TanBeta*(AtauInput - MuInput*TanBeta)*Cos(2*ArcTan(TanBeta))*Sqr(Ye(2,
      2))*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*MuInput*Log(MSUSY/Qmatch)*Sqr
      (g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + 30*MuInput*Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)
      ) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)))*(-0.9*Sqr(g1)*
      TCF(3)(Sqrt(msl2(2,2)/mse2(2,2))) + 0.3*Sqr(g1)*TCF(4)(Sqrt(msl2(2,2)/mse2(2
      ,2))) - Sqr(g2)*TCF(4)(Sqrt(msl2(2,2)/mse2(2,2)))) + Sqr(AtauInput - MuInput
      *TanBeta)*(0.006666666666666667*Cos(2*ArcTan(TanBeta))*Sqr(Ye(2,2))*(-1107*
      Log(MSUSY/Qmatch)*Quad(g1)*TCF(3)(Sqrt(msl2(2,2)/mse2(2,2))) + 369*Log(MSUSY
      /Qmatch)*Quad(g1)*TCF(4)(Sqrt(msl2(2,2)/mse2(2,2))) + 950*Log(MSUSY/Qmatch)*
      Quad(g2)*TCF(4)(Sqrt(msl2(2,2)/mse2(2,2)))) + ((-2*Sin(2*ArcTan(TanBeta))*
      Sqr(Ye(2,2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/
      TanBeta + 0.5*Cos(2*ArcTan(TanBeta))*Ye(2,2)*(10*Cube(Ye(2,2))*Log(MSUSY/
      Qmatch) - 9*Log(MSUSY/Qmatch)*Sqr(g1)*Ye(2,2) - 9*Log(MSUSY/Qmatch)*Sqr(g2)*
      Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Ye(2,2) + 12*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2))*Ye(2,2)))*(-0.9*Sqr(g1)*TCF(3)(Sqrt(msl2(2,2)/mse2(2,2))) + 0.3
      *Sqr(g1)*TCF(4)(Sqrt(msl2(2,2)/mse2(2,2))) - Sqr(g2)*TCF(4)(Sqrt(msl2(2,2)/
      mse2(2,2)))))))/Sqrt(mse2(2,2)*msl2(2,2)) + (0.08333333333333333*(0.4*
      TanBeta*(AtauInput - MuInput*TanBeta)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Ye(2,2))*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*
      MuInput*Log(MSUSY/Qmatch)*Sqr(g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)
      ) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(
      MSUSY/Qmatch)*Sqr(Ye(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(
      Ye(2,2))) - Sqr(AtauInput - MuInput*TanBeta)*((4.92*Log(MSUSY/Qmatch)*Quad(
      g1) - 6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta
      )))*Sqr(Ye(2,2)) + (0.6*Sqr(g1) + Sqr(g2))*((-4*Cos(2*ArcTan(TanBeta))*Sin(2
      *ArcTan(TanBeta))*Sqr(Ye(2,2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)
      ) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu
      (2,2))))/TanBeta + 0.5*Sqr(Cos(2*ArcTan(TanBeta)))*Ye(2,2)*(10*Cube(Ye(2,2))
      *Log(MSUSY/Qmatch) - 9*Log(MSUSY/Qmatch)*Sqr(g1)*Ye(2,2) - 9*Log(MSUSY/
      Qmatch)*Sqr(g2)*Ye(2,2) + 12*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Ye(2,2) + 12*Log
      (MSUSY/Qmatch)*Sqr(Yu(2,2))*Ye(2,2)))))*TCF(5)(Sqrt(msl2(2,2)/mse2(2,2))))/
      Sqrt(mse2(2,2)*msl2(2,2))))/Sqr(3.141592653589793) + (0.0625*(3*((Sqr(Yd(2,2
      ))*(-0.2*Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(Yd(2,2)))*(0.4*Log(MSUSY/
      Qmatch)*msd2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.4*
      Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1
      ) + 0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(2,2
      )*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msl2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/
      Qmatch)*msq2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.4*
      Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1
      ) - 0.8*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(2,2
      )*Sqr(g1) - 0.5333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.4*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.4*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 4*Log(MSUSY/Qmatch)*
      Sqr(mAInput)*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))
      *Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)
      ) + 4*Log(MSUSY/Qmatch)*Sqr(AbInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) - 4*Log
      (MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))))/msd2(2,2) +
      Log(msd2(2,2)/Sqr(Qmatch))*(0.5*(-0.2*Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(
      Yd(2,2)))*Yd(2,2)*(18*Cube(Yd(2,2))*Log(MSUSY/Qmatch) - Log(MSUSY/Qmatch)*
      Sqr(g1)*Yd(2,2) - 9*Log(MSUSY/Qmatch)*Sqr(g2)*Yd(2,2) - 32*Log(MSUSY/Qmatch)
      *Sqr(g3)*Yd(2,2) + 4*Log(MSUSY/Qmatch)*Sqr(Ye(2,2))*Yd(2,2) + 6*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2))*Yd(2,2)) + Sqr(Yd(2,2))*(0.2*(-8.2*Cos(2*ArcTan(TanBeta
      ))*Log(MSUSY/Qmatch)*Quad(g1) + (2*Sin(2*ArcTan(TanBeta))*Sqr(g1)*(3*Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr
      (Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta) + 0.5*Yd(2,2)*(18*
      Cube(Yd(2,2))*Log(MSUSY/Qmatch) - Log(MSUSY/Qmatch)*Sqr(g1)*Yd(2,2) - 9*Log(
      MSUSY/Qmatch)*Sqr(g2)*Yd(2,2) - 32*Log(MSUSY/Qmatch)*Sqr(g3)*Yd(2,2) + 4*Log
      (MSUSY/Qmatch)*Sqr(Ye(2,2))*Yd(2,2) + 6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Yd(2,
      2))))) + 3*((Sqr(Yd(2,2))*(-0.1*Cos(2*ArcTan(TanBeta))*(Sqr(g1) + 5*Sqr(g2))
      + Sqr(Yd(2,2)))*(0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0
      )*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/msq2(2,2) +
      Log(msq2(2,2)/Sqr(Qmatch))*(0.5*(-0.1*Cos(2*ArcTan(TanBeta))*(Sqr(g1) + 5*
      Sqr(g2)) + Sqr(Yd(2,2)))*Yd(2,2)*(18*Cube(Yd(2,2))*Log(MSUSY/Qmatch) - Log(
      MSUSY/Qmatch)*Sqr(g1)*Yd(2,2) - 9*Log(MSUSY/Qmatch)*Sqr(g2)*Yd(2,2) - 32*Log
      (MSUSY/Qmatch)*Sqr(g3)*Yd(2,2) + 4*Log(MSUSY/Qmatch)*Sqr(Ye(2,2))*Yd(2,2) +
      6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Yd(2,2)) + Sqr(Yd(2,2))*(0.5*(-(Cos(2*
      ArcTan(TanBeta))*(1.64*Log(MSUSY/Qmatch)*Quad(g1) - 6.333333333333333*Log(
      MSUSY/Qmatch)*Quad(g2))) + (2*Sin(2*ArcTan(TanBeta))*(0.2*Sqr(g1) + Sqr(g2))
      *(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta) + 0.5*Yd
      (2,2)*(18*Cube(Yd(2,2))*Log(MSUSY/Qmatch) - Log(MSUSY/Qmatch)*Sqr(g1)*Yd(2,2
      ) - 9*Log(MSUSY/Qmatch)*Sqr(g2)*Yd(2,2) - 32*Log(MSUSY/Qmatch)*Sqr(g3)*Yd(2,
      2) + 4*Log(MSUSY/Qmatch)*Sqr(Ye(2,2))*Yd(2,2) + 6*Log(MSUSY/Qmatch)*Sqr(Yu(2
      ,2))*Yd(2,2))))) + 6*(((-0.4*TanBeta*(AbInput - MuInput*TanBeta)*Quad(Yd(2,2
      ))*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*MuInput*Log(MSUSY/Qmatch)*Sqr(
      g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + 30*MuInput*Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)
      ) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2))))/Sqrt(msd2(2,2)*
      msq2(2,2)) + Sqr(AbInput - MuInput*TanBeta)*((-0.5*Quad(Yd(2,2))*(msq2(2,2)*
      (0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(1,1)*
      Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.4*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.5333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (M1Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.4
      *Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.4*Log(MSUSY
      /Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 4*Log(MSUSY
      /Qmatch)*Sqr(mAInput)*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(
      TanBeta))*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*Sqr(AbInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,
      2)) - 4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))) +
      msd2(2,2)*(0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1
      ) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msd2(2,2)*
      msq2(2,2)*Sqrt(msd2(2,2)*msq2(2,2))) + (Cube(Yd(2,2))*(18*Cube(Yd(2,2))*Log(
      MSUSY/Qmatch) - Log(MSUSY/Qmatch)*Sqr(g1)*Yd(2,2) - 9*Log(MSUSY/Qmatch)*Sqr(
      g2)*Yd(2,2) - 32*Log(MSUSY/Qmatch)*Sqr(g3)*Yd(2,2) + 4*Log(MSUSY/Qmatch)*Sqr
      (Ye(2,2))*Yd(2,2) + 6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Yd(2,2)))/Sqrt(msd2(2,2
      )*msq2(2,2))))*(TCF(1)(Sqrt(msq2(2,2)/msd2(2,2))) - (0.08333333333333333*Sqr
      (AbInput - MuInput*TanBeta)*TCF(2)(Sqrt(msq2(2,2)/msd2(2,2))))/Sqrt(msd2(2,2
      )*msq2(2,2))) + (Quad(Yd(2,2))*Sqr(AbInput - MuInput*TanBeta)*((0.5*msd2(2,2
      )*Sqrt(msq2(2,2)/msd2(2,2))*((msq2(2,2)*(-0.4*Log(MSUSY/Qmatch)*Sqr(g1) - (
      0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*
      msd2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(
      MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(0,0
      )*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msd2(2,2) +
      (0.4*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)
      *msq2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1))/msd2(2,2) + (0.8*Log(
      MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(1,1
      )*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1))/msd2(2,2) +
      (0.5333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input))/msd2(2,2) + (
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msd2(2,2) - (0.4*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msd2(2,2)*(1 + Sqr(TanBeta))) + (
      0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(msd2(2,2)*(1 + Sqr
      (TanBeta))) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*
      Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*
      Sqr(AbInput)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr
      (Yd(2,2)))/msd2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yd(2,2)))/msd2(
      2,2) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) -
      (4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log
      (MSUSY/Qmatch)*Sqr(AbInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) + (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2)))/msd2(2,2)
      + (0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)
      *Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))/msd2(2,2))*TCD1F(1)(Sqrt(msq2(2,
      2)/msd2(2,2))))/msq2(2,2) + 0.08333333333333333*((0.4*TanBeta*(AbInput -
      MuInput*TanBeta)*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*MuInput*Log(
      MSUSY/Qmatch)*Sqr(g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + 30*
      MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(MSUSY/
      Qmatch)*Sqr(Ye(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)
      ))*TCF(2)(Sqrt(msq2(2,2)/msd2(2,2))))/Sqrt(msd2(2,2)*msq2(2,2)) - Sqr(
      AbInput - MuInput*TanBeta)*((0.5*msd2(2,2)*Sqrt(msq2(2,2)/msd2(2,2))*((msq2(
      2,2)*(-0.4*Log(MSUSY/Qmatch)*Sqr(g1) - (0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(
      g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*
      Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2
      (1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/msd2(2,
      2) + (0.4*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/
      Qmatch)*msl2(1,1)*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(
      g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*
      Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2
      (2,2)*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/msd2(2,
      2) + (0.8*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1))/msd2(2,2) + (0.5333333333333333*Log(MSUSY/Qmatch)
      *Sqr(g1)*Sqr(M1Input))/msd2(2,2) + (10.666666666666666*Log(MSUSY/Qmatch)*Sqr
      (g3)*Sqr(M3Input))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/
      (msd2(2,2)*(1 + Sqr(TanBeta))) + (0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)
      *Sqr(TanBeta))/(msd2(2,2)*(1 + Sqr(TanBeta))) - (4*Log(MSUSY/Qmatch)*msq2(2,
      2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yd(2,2)))/
      msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AbInput)*Sqr(Yd(2,2)))/msd2(2,2) - (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yd(2,2)))/msd2(2,2) + (4*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(
      TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(
      TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AbInput)*Sqr(
      TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(
      TanBeta)*Sqr(Yd(2,2)))/msd2(2,2)))/msd2(2,2) + (0.2*Log(MSUSY/Qmatch)*msd2(0
      ,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch
      )*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(
      MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) -
      0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*
      Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1
      ) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY
      /Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)
      *Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(
      TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 +
      Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta))/msd2(2,2))*TCD1F(2)(Sqrt(msq2(2,2)/msd2(2,2))))/
      (msq2(2,2)*Sqrt(msd2(2,2)*msq2(2,2))) - (0.5*(msq2(2,2)*(0.4*Log(MSUSY/
      Qmatch)*msd2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.4*
      Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1
      ) + 0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(2,2
      )*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msl2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/
      Qmatch)*msq2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.4*
      Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1
      ) - 0.8*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(2,2
      )*Sqr(g1) - 0.5333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.4*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.4*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 4*Log(MSUSY/Qmatch)*
      Sqr(mAInput)*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))
      *Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)
      ) + 4*Log(MSUSY/Qmatch)*Sqr(AbInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) - 4*Log
      (MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))) + msd2(2,2)*(
      0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*TCF(2)(Sqrt(msq2(2,2)/msd2(2,2)
      )))/(msd2(2,2)*msq2(2,2)*Sqrt(msd2(2,2)*msq2(2,2)))))))/Sqrt(msd2(2,2)*msq2(
      2,2))) + 0.75*((-0.4*TanBeta*(AbInput - MuInput*TanBeta)*Cos(2*ArcTan(
      TanBeta))*Sqr(Yd(2,2))*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*MuInput*
      Log(MSUSY/Qmatch)*Sqr(g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + 30*
      MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(MSUSY/
      Qmatch)*Sqr(Ye(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)
      ))*(-0.3*Sqr(g1)*TCF(3)(Sqrt(msq2(2,2)/msd2(2,2))) + (-0.3*Sqr(g1) - Sqr(g2)
      )*TCF(4)(Sqrt(msq2(2,2)/msd2(2,2)))))/Sqrt(msd2(2,2)*msq2(2,2)) + Sqr(
      AbInput - MuInput*TanBeta)*((-0.5*Cos(2*ArcTan(TanBeta))*Sqr(Yd(2,2))*(msq2(
      2,2)*(0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(1
      ,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch
      )*mse2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.4*Log(
      MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) -
      0.4*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(2,2)*
      Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*
      msq2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/
      Qmatch)*msu2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.8*
      Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.5333333333333333*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(M1Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) -
      (0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta))
      + 4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msq2(2
      ,2)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr
      (TanBeta))*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*Sqr(AbInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yd(2,2)) - 4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yd(
      2,2))) + msd2(2,2)*(0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0
      )*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*(-0.3*Sqr(g1)
      *TCF(3)(Sqrt(msq2(2,2)/msd2(2,2))) + (-0.3*Sqr(g1) - Sqr(g2))*TCF(4)(Sqrt(
      msq2(2,2)/msd2(2,2)))))/(msd2(2,2)*msq2(2,2)*Sqrt(msd2(2,2)*msq2(2,2))) + (
      Cos(2*ArcTan(TanBeta))*Sqr(Yd(2,2))*((0.5*msd2(2,2)*Sqrt(msq2(2,2)/msd2(2,2)
      )*(-0.3*Sqr(g1) - Sqr(g2))*((msq2(2,2)*(-0.4*Log(MSUSY/Qmatch)*Sqr(g1) - (
      0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*
      msd2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(
      MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(0,0
      )*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msd2(2,2) +
      (0.4*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)
      *msq2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1))/msd2(2,2) + (0.8*Log(
      MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(1,1
      )*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1))/msd2(2,2) +
      (0.5333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input))/msd2(2,2) + (
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msd2(2,2) - (0.4*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msd2(2,2)*(1 + Sqr(TanBeta))) + (
      0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(msd2(2,2)*(1 + Sqr
      (TanBeta))) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*
      Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*
      Sqr(AbInput)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr
      (Yd(2,2)))/msd2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yd(2,2)))/msd2(
      2,2) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) -
      (4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log
      (MSUSY/Qmatch)*Sqr(AbInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) + (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2)))/msd2(2,2)
      + (0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)
      *Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))/msd2(2,2))*TCD1F(4)(Sqrt(msq2(2,
      2)/msd2(2,2))))/msq2(2,2) - 0.3*((0.5*msd2(2,2)*Sqrt(msq2(2,2)/msd2(2,2))*
      Sqr(g1)*((msq2(2,2)*(-0.4*Log(MSUSY/Qmatch)*Sqr(g1) - (0.4*Log(MSUSY/Qmatch)
      *msd2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(
      MSUSY/Qmatch)*mse2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(2,2
      )*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1))/msd2(2,2) +
      (0.4*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)
      *msl2(2,2)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(
      MSUSY/Qmatch)*msq2(2,2)*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(0,0
      )*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1))/msd2(2,2) +
      (0.8*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1))/msd2(2,2) + (0.5333333333333333*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input))/msd2(2,2) + (10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*Sqr(
      g1)*Sqr(mAInput))/(msd2(2,2)*(1 + Sqr(TanBeta))) + (0.4*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(msd2(2,2)*(1 + Sqr(TanBeta))) - (4*Log(
      MSUSY/Qmatch)*msq2(2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*msu2(
      2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AbInput)*Sqr(Yd(2,2)
      ))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yd(2,2)))/msd2(2,2) + (
      4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch
      )*msu2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(
      AbInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2)))/msd2(2,2) + (0.2*Log(MSUSY/
      Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2
      )*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2
      )*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*
      Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666*Log(MSUSY/Qmatch
      )*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 +
      Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1
       + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta))/msd2(2,2))*TCD1F(3)(Sqrt(msq2(2,2)/msd2(2,2)))
      )/msq2(2,2) + 8.2*Log(MSUSY/Qmatch)*Quad(g1)*TCF(3)(Sqrt(msq2(2,2)/msd2(2,2)
      ))) + (-2.46*Log(MSUSY/Qmatch)*Quad(g1) + 6.333333333333333*Log(MSUSY/Qmatch
      )*Quad(g2))*TCF(4)(Sqrt(msq2(2,2)/msd2(2,2)))) + ((-2*Sin(2*ArcTan(TanBeta))
      *Sqr(Yd(2,2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/
      TanBeta + 0.5*Cos(2*ArcTan(TanBeta))*Yd(2,2)*(18*Cube(Yd(2,2))*Log(MSUSY/
      Qmatch) - Log(MSUSY/Qmatch)*Sqr(g1)*Yd(2,2) - 9*Log(MSUSY/Qmatch)*Sqr(g2)*Yd
      (2,2) - 32*Log(MSUSY/Qmatch)*Sqr(g3)*Yd(2,2) + 4*Log(MSUSY/Qmatch)*Sqr(Ye(2,
      2))*Yd(2,2) + 6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Yd(2,2)))*(-0.3*Sqr(g1)*TCF(3
      )(Sqrt(msq2(2,2)/msd2(2,2))) + (-0.3*Sqr(g1) - Sqr(g2))*TCF(4)(Sqrt(msq2(2,2
      )/msd2(2,2)))))/Sqrt(msd2(2,2)*msq2(2,2)))) + 0.05*((0.4*TanBeta*(AbInput -
      MuInput*TanBeta)*(3*Sqr(g1) + 5*Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Yd(
      2,2))*(-3*MuInput*Log(MSUSY/Qmatch)*Sqr(g1) - 15*MuInput*Log(MSUSY/Qmatch)*
      Sqr(g2) + 30*MuInput*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + 30*MuInput*Log(MSUSY/
      Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)
      ) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)))*TCF(5)(Sqrt(msq2
      (2,2)/msd2(2,2))))/Sqrt(msd2(2,2)*msq2(2,2)) - Sqr(AbInput - MuInput*TanBeta
      )*((-0.5*(3*Sqr(g1) + 5*Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Yd(2,2))*(
      msq2(2,2)*(0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*
      msd2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/
      Qmatch)*mse2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.4*
      Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(2,2
      )*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*
      msq2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/
      Qmatch)*msu2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.8*
      Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.5333333333333333*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(M1Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) -
      (0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta))
      + 4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msq2(2
      ,2)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr
      (TanBeta))*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*Sqr(AbInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yd(2,2)) - 4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yd(
      2,2))) + msd2(2,2)*(0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0
      )*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*TCF(5)(Sqrt(
      msq2(2,2)/msd2(2,2))))/(msd2(2,2)*msq2(2,2)*Sqrt(msd2(2,2)*msq2(2,2))) + ((
      0.5*msd2(2,2)*Sqrt(msq2(2,2)/msd2(2,2))*(3*Sqr(g1) + 5*Sqr(g2))*Sqr(Cos(2*
      ArcTan(TanBeta)))*Sqr(Yd(2,2))*((msq2(2,2)*(-0.4*Log(MSUSY/Qmatch)*Sqr(g1) -
      (0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)
      *msd2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1))/msd2(2,2) - (0.4*Log(
      MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(0,0
      )*Sqr(g1))/msd2(2,2) + (0.4*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msd2(2,2) +
      (0.4*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)
      *msq2(0,0)*Sqr(g1))/msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/
      msd2(2,2) - (0.4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1))/msd2(2,2) + (0.8*Log(
      MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(1,1
      )*Sqr(g1))/msd2(2,2) + (0.8*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1))/msd2(2,2) +
      (0.5333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input))/msd2(2,2) + (
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msd2(2,2) - (0.4*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msd2(2,2)*(1 + Sqr(TanBeta))) + (
      0.4*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(msd2(2,2)*(1 + Sqr
      (TanBeta))) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*
      Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*
      Sqr(AbInput)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr
      (Yd(2,2)))/msd2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yd(2,2)))/msd2(
      2,2) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) -
      (4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) - (4*Log
      (MSUSY/Qmatch)*Sqr(AbInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2) + (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*Sqr(TanBeta)*Sqr(Yd(2,2)))/msd2(2,2)))/msd2(2,2)
      + (0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)
      *Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))/msd2(2,2))*TCD1F(5)(Sqrt(msq2(2,
      2)/msd2(2,2))))/msq2(2,2) + ((24.6*Log(MSUSY/Qmatch)*Quad(g1) -
      31.666666666666668*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*
      Sqr(Yd(2,2)) + (3*Sqr(g1) + 5*Sqr(g2))*((-4*Cos(2*ArcTan(TanBeta))*Sin(2*
      ArcTan(TanBeta))*Sqr(Yd(2,2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2))
      + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2
      ,2))))/TanBeta + 0.5*Sqr(Cos(2*ArcTan(TanBeta)))*Yd(2,2)*(18*Cube(Yd(2,2))*
      Log(MSUSY/Qmatch) - Log(MSUSY/Qmatch)*Sqr(g1)*Yd(2,2) - 9*Log(MSUSY/Qmatch)*
      Sqr(g2)*Yd(2,2) - 32*Log(MSUSY/Qmatch)*Sqr(g3)*Yd(2,2) + 4*Log(MSUSY/Qmatch)
      *Sqr(Ye(2,2))*Yd(2,2) + 6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))*Yd(2,2))))*TCF(5)(
      Sqrt(msq2(2,2)/msd2(2,2))))/Sqrt(msd2(2,2)*msq2(2,2))))))/Sqr(
      3.141592653589793) + (0.0625*(-0.1875*(0.005333333333333333*(1107*Log(MSUSY/
      Qmatch)*Power6(g1) - 2375*Log(MSUSY/Qmatch)*Power6(g2) - 1425*Log(MSUSY/
      Qmatch)*Quad(g2)*Sqr(g1) + 1845*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(g2))*Sqr(Sin(
      4*ArcTan(TanBeta))) + (8*Cos(4*ArcTan(TanBeta))*Sin(4*ArcTan(TanBeta))*(0.36
      *Quad(g1) + Quad(g2) + 1.2*Sqr(g1)*Sqr(g2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta
      )*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2))))/TanBeta) + 0.00020833333333333335*Log(Sqr(mAInput)/
      Sqr(Qmatch))*(4280.4*Log(MSUSY/Qmatch)*Power6(g1) - 16783.333333333332*Log(
      MSUSY/Qmatch)*Power6(g2) - 3990*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(g1) + 5166*
      Log(MSUSY/Qmatch)*Quad(g1)*Sqr(g2) - 9*(Cos(8*ArcTan(TanBeta))*(147.6*Log(
      MSUSY/Qmatch)*Power6(g1) - 316.6666666666667*Log(MSUSY/Qmatch)*Power6(g2) -
      190*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(g1) + 246*Log(MSUSY/Qmatch)*Quad(g1)*Sqr(
      g2)) - (8*Sin(8*ArcTan(TanBeta))*(9*Quad(g1) + 25*Quad(g2) + 30*Sqr(g1)*Sqr(
      g2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta) - 4*(Cos
      (4*ArcTan(TanBeta))*(147.6*Log(MSUSY/Qmatch)*Power6(g1) - 2216.6666666666665
      *Log(MSUSY/Qmatch)*Power6(g2) - 570*Log(MSUSY/Qmatch)*Quad(g2)*Sqr(g1) + 738
      *Log(MSUSY/Qmatch)*Quad(g1)*Sqr(g2)) - (4*Sin(4*ArcTan(TanBeta))*(9*Quad(g1)
      + 175*Quad(g2) + 90*Sqr(g1)*Sqr(g2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(
      Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)
      *Sqr(Yu(2,2))))/TanBeta)) + 0.0033333333333333335*((-4*Cos(2*ArcTan(TanBeta)
      )*(6*(Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2
      )/Sqr(Qmatch)))*Quad(g1) + 18*(Log(mse2(0,0)/Sqr(Qmatch)) + Log(mse2(1,1)/
      Sqr(Qmatch)) + Log(mse2(2,2)/Sqr(Qmatch)))*Quad(g1) + 24*(Log(msu2(0,0)/Sqr(
      Qmatch)) + Log(msu2(1,1)/Sqr(Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)))*Quad(g1)
      + 3*(Log(msq2(0,0)/Sqr(Qmatch)) + Log(msq2(1,1)/Sqr(Qmatch)) + Log(msq2(2,2)
      /Sqr(Qmatch)))*(Quad(g1) + 25*Quad(g2)) + (Log(msl2(0,0)/Sqr(Qmatch)) + Log(
      msl2(1,1)/Sqr(Qmatch)) + Log(msl2(2,2)/Sqr(Qmatch)))*(9*Quad(g1) + 25*Quad(
      g2)))*Sin(2*ArcTan(TanBeta))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2))
      + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2
      ,2))))/TanBeta + Sqr(Cos(2*ArcTan(TanBeta)))*(295.2*Log(MSUSY/Qmatch)*(Log(
      mse2(0,0)/Sqr(Qmatch)) + Log(mse2(1,1)/Sqr(Qmatch)) + Log(mse2(2,2)/Sqr(
      Qmatch)))*Power6(g1) + 0.13333333333333333*(Log(msl2(0,0)/Sqr(Qmatch)) + Log
      (msl2(1,1)/Sqr(Qmatch)) + Log(msl2(2,2)/Sqr(Qmatch)))*(1107*Log(MSUSY/Qmatch
      )*Power6(g1) - 2375*Log(MSUSY/Qmatch)*Power6(g2)) + 6*(16.4*Log(MSUSY/Qmatch
      )*(Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2)/
      Sqr(Qmatch)))*Power6(g1) + (Quad(g1)*(0.4*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1
      ) + 0.4*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msd2(2,2
      )*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*
      mse2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msl2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1
      ) + 0.4*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.4*Log(MSUSY/Qmatch)*msq2(2,2
      )*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      msu2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) -
      0.5333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.4*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.4*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 4*Log(MSUSY/Qmatch)*
      Sqr(mAInput)*Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))
      *Sqr(Yd(2,2)) + 4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)
      ) + 4*Log(MSUSY/Qmatch)*Sqr(AbInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) - 4*Log
      (MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))))/msd2(2,2)) +
      24*(16.4*Log(MSUSY/Qmatch)*(Log(msu2(0,0)/Sqr(Qmatch)) + Log(msu2(1,1)/Sqr(
      Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)))*Power6(g1) + (Quad(g1)*(-0.8*Log(
      MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) -
      0.8*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(0,0)*
      Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      mse2(2,2)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) + 0.8*Log(MSUSY/
      Qmatch)*msl2(1,1)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) - 0.8*
      Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1
      ) - 0.8*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(0,0
      )*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*
      msu2(2,2)*Sqr(g1) - 2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input
      ) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) - (0.8*Log(
      MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) + (0.8*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (4*Log(MSUSY
      /Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2
      (2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/
      msu2(2,2)) + 3*((Log(msq2(0,0)/Sqr(Qmatch)) + Log(msq2(1,1)/Sqr(Qmatch)) +
      Log(msq2(2,2)/Sqr(Qmatch)))*(16.4*Log(MSUSY/Qmatch)*Power6(g1) -
      316.6666666666667*Log(MSUSY/Qmatch)*Power6(g2)) + ((Quad(g1) + 25*Quad(g2))*
      (0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/msq2(2,2)))) + 3*((Sqr(Yu(2,2))
      *(0.4*Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(Yu(2,2)))*(-0.8*Log(MSUSY/Qmatch)
      *msd2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) - 0.8*Log(MSUSY
      /Qmatch)*msd2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) - 0.8*
      Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1
      ) + 0.8*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(1,1
      )*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      msq2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/
      Qmatch)*msq2(2,2)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) + 1.6*
      Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1
      ) - 2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) - (0.8*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) + (0.8*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/msu2(2,2) +
      Log(msu2(2,2)/Sqr(Qmatch))*(0.1*(0.4*Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(Yu
      (2,2)))*Yu(2,2)*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 17*Log(MSUSY/Qmatch)*
      Sqr(g1)*Yu(2,2) - 45*Log(MSUSY/Qmatch)*Sqr(g2)*Yu(2,2) - 160*Log(MSUSY/
      Qmatch)*Sqr(g3)*Yu(2,2) + 30*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Yu(2,2) + 20*Log
      (MSUSY/Qmatch)*Sqr(Ye(2,2))*Yu(2,2)) + Sqr(Yu(2,2))*(0.4*(8.2*Cos(2*ArcTan(
      TanBeta))*Log(MSUSY/Qmatch)*Quad(g1) - (2*Sin(2*ArcTan(TanBeta))*Sqr(g1)*(3*
      Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)
      *Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta) + 0.1*Yu(2,2)*(
      90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 17*Log(MSUSY/Qmatch)*Sqr(g1)*Yu(2,2) -
      45*Log(MSUSY/Qmatch)*Sqr(g2)*Yu(2,2) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(2,2)
      + 30*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Yu(2,2) + 20*Log(MSUSY/Qmatch)*Sqr(Ye(2,
      2))*Yu(2,2))))) + 3*((Sqr(Yu(2,2))*(-0.1*Cos(2*ArcTan(TanBeta))*(Sqr(g1) - 5
      *Sqr(g2)) + Sqr(Yu(2,2)))*(0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log
      (MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) +
      0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/
      Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1
      )*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(
      M2Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY
      /Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2
      (2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/
      msq2(2,2) + Log(msq2(2,2)/Sqr(Qmatch))*(0.1*(-0.1*Cos(2*ArcTan(TanBeta))*(
      Sqr(g1) - 5*Sqr(g2)) + Sqr(Yu(2,2)))*Yu(2,2)*(90*Cube(Yu(2,2))*Log(MSUSY/
      Qmatch) - 17*Log(MSUSY/Qmatch)*Sqr(g1)*Yu(2,2) - 45*Log(MSUSY/Qmatch)*Sqr(g2
      )*Yu(2,2) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(2,2) + 30*Log(MSUSY/Qmatch)*Sqr
      (Yd(2,2))*Yu(2,2) + 20*Log(MSUSY/Qmatch)*Sqr(Ye(2,2))*Yu(2,2)) + Sqr(Yu(2,2)
      )*(0.5*(Cos(2*ArcTan(TanBeta))*(-1.64*Log(MSUSY/Qmatch)*Quad(g1) -
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2)) - (2*Sin(2*ArcTan(TanBeta))*(-
      0.2*Sqr(g1) + Sqr(g2))*(3*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))
      /TanBeta) + 0.1*Yu(2,2)*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 17*Log(MSUSY/
      Qmatch)*Sqr(g1)*Yu(2,2) - 45*Log(MSUSY/Qmatch)*Sqr(g2)*Yu(2,2) - 160*Log(
      MSUSY/Qmatch)*Sqr(g3)*Yu(2,2) + 30*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Yu(2,2) +
      20*Log(MSUSY/Qmatch)*Sqr(Ye(2,2))*Yu(2,2))))) + 6*(((2*(AtInput - MuInput/
      TanBeta)*Quad(Yu(2,2))*(1.7333333333333334*M1Input*Log(MSUSY/Qmatch)*Sqr(g2)
      + 6*M2Input*Log(MSUSY/Qmatch)*Sqr(g2) + 10.666666666666666*M3Input*Log(MSUSY
      /Qmatch)*Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (0.6*(-(MuInput*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta))
      - 5*MuInput*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta) + 10*MuInput*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2)
      )))/Cube(TanBeta)))/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AtInput - MuInput/
      TanBeta)*((-0.5*Quad(Yu(2,2))*(msq2(2,2)*(-0.8*Log(MSUSY/Qmatch)*msd2(0,0)*
      Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      msd2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/
      Qmatch)*mse2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) + 0.8*
      Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1
      ) + 0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(0,0
      )*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      msq2(2,2)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) + 1.6*Log(MSUSY/
      Qmatch)*msu2(1,1)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) -
      2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) - (0.8*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) + (0.8*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(
      0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2
      (2,2)*msu2(2,2))) + (0.2*Cube(Yu(2,2))*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) -
      17*Log(MSUSY/Qmatch)*Sqr(g1)*Yu(2,2) - 45*Log(MSUSY/Qmatch)*Sqr(g2)*Yu(2,2)
      - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(2,2) + 30*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*
      Yu(2,2) + 20*Log(MSUSY/Qmatch)*Sqr(Ye(2,2))*Yu(2,2)))/Sqrt(msq2(2,2)*msu2(2,
      2))))*(TCF(1)(Sqrt(msq2(2,2)/msu2(2,2))) - (0.08333333333333333*Sqr(AtInput
      - MuInput/TanBeta)*TCF(2)(Sqrt(msq2(2,2)/msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,
      2))) + (Quad(Yu(2,2))*Sqr(AtInput - MuInput/TanBeta)*((0.5*Sqrt(msq2(2,2)/
      msu2(2,2))*msu2(2,2)*((msq2(2,2)*(-1.6*Log(MSUSY/Qmatch)*Sqr(g1) + (0.8*Log(
      MSUSY/Qmatch)*msd2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msd2(1,1
      )*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1))/msu2(2,2) +
      (0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)
      *mse2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/
      msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1))/msu2(2,2) - (0.8*Log(
      MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2(2,2
      )*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1))/msu2(2,2) +
      (0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)
      *msq2(2,2)*Sqr(g1))/msu2(2,2) - (1.6*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/
      msu2(2,2) - (1.6*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1))/msu2(2,2) + (
      2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input))/msu2(2,2) + (
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) + (0.8*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msu2(2,2)*(1 + Sqr(TanBeta))) - (
      0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(msu2(2,2)*(1 + Sqr
      (TanBeta))) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)*msq2(2
      ,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2))
      )/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4
      *Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,
      2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput
      )*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)
      *Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,2) + (0.2*Log(MSUSY/Qmatch)*
      msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1
      ) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1
      )*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1
      ) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY
      /Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)
      *Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(
      TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 +
      Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2,2))*TCD1F(1)(Sqrt(msq2(2,2)/msu2(2,2))))/
      msq2(2,2) + 0.08333333333333333*((-2*(AtInput - MuInput/TanBeta)*(
      1.7333333333333334*M1Input*Log(MSUSY/Qmatch)*Sqr(g2) + 6*M2Input*Log(MSUSY/
      Qmatch)*Sqr(g2) + 10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12
      *AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      0.6*(-(MuInput*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta)) - 5*MuInput*Log(MSUSY
      /Qmatch)*Sqr(g2)*Sqr(TanBeta) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) +
      10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta))*TCF(
      2)(Sqrt(msq2(2,2)/msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) - Sqr(AtInput -
      MuInput/TanBeta)*((0.5*Sqrt(msq2(2,2)/msu2(2,2))*msu2(2,2)*((msq2(2,2)*(-1.6
      *Log(MSUSY/Qmatch)*Sqr(g1) + (0.8*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1))/msu2(
      2,2) + (0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/
      Qmatch)*msd2(2,2)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(
      g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*
      Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2
      (0,0)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msu2(2,
      2) - (0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/
      Qmatch)*msq2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(
      g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1))/msu2(2,2) - (1.6*
      Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/msu2(2,2) - (1.6*Log(MSUSY/Qmatch)*msu2
      (1,1)*Sqr(g1))/msu2(2,2) + (2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (M1Input))/msu2(2,2) + (10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msu2(2,2
      )*(1 + Sqr(TanBeta))) - (0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(
      TanBeta))/(msu2(2,2)*(1 + Sqr(TanBeta))) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))
      - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta))))/msu2(2,2) + (0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(
      MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) +
      0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/
      Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1
      )*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(
      M2Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY
      /Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2
      (2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))/
      msu2(2,2))*TCD1F(2)(Sqrt(msq2(2,2)/msu2(2,2))))/(msq2(2,2)*Sqrt(msq2(2,2)*
      msu2(2,2))) - (0.5*(msq2(2,2)*(-0.8*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) -
      0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msd2(2,2)*
      Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      mse2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) + 0.8*Log(MSUSY/
      Qmatch)*msl2(0,0)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) + 0.8*
      Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1
      ) - 0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(2,2
      )*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*
      msu2(1,1)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) -
      2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) - (0.8*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) + (0.8*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(
      0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*TCF(2)(Sqrt(msq2(2,2)/msu2(2,2)
      )))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2)))))))/Sqrt(msq2(2,2)*msu2(
      2,2))) + 0.75*((2*(AtInput - MuInput/TanBeta)*Cos(2*ArcTan(TanBeta))*Sqr(Yu(
      2,2))*(1.7333333333333334*M1Input*Log(MSUSY/Qmatch)*Sqr(g2) + 6*M2Input*Log(
      MSUSY/Qmatch)*Sqr(g2) + 10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3)
      + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) - (0.6*(-(MuInput*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta)) - 5*MuInput*Log(
      MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2
      )) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta))
      *(0.6*Sqr(g1)*TCF(3)(Sqrt(msq2(2,2)/msu2(2,2))) + Sqr(g2)*TCF(4)(Sqrt(msq2(2
      ,2)/msu2(2,2)))))/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AtInput - MuInput/TanBeta)
      *((-0.5*Cos(2*ArcTan(TanBeta))*Sqr(Yu(2,2))*(msq2(2,2)*(-0.8*Log(MSUSY/
      Qmatch)*msd2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) - 0.8*
      Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1
      ) - 0.8*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(2,2
      )*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*
      msl2(1,1)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/
      Qmatch)*msq2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) - 0.8*
      Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1
      ) + 1.6*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(2,2
      )*Sqr(g1) - 2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) - (0.8*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) + (0.8*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(
      0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*(0.6*Sqr(g1)*TCF(3)(Sqrt(msq2(2
      ,2)/msu2(2,2))) + Sqr(g2)*TCF(4)(Sqrt(msq2(2,2)/msu2(2,2)))))/(msq2(2,2)*
      msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) + (Cos(2*ArcTan(TanBeta))*Sqr(Yu(2,2))*
      ((0.5*Sqrt(msq2(2,2)/msu2(2,2))*msu2(2,2)*Sqr(g2)*((msq2(2,2)*(-1.6*Log(
      MSUSY/Qmatch)*Sqr(g1) + (0.8*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1))/msu2(2,2)
      + (0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/
      Qmatch)*msd2(2,2)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(
      g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*
      Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2
      (0,0)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msu2(2,
      2) - (0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/
      Qmatch)*msq2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(
      g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1))/msu2(2,2) - (1.6*
      Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/msu2(2,2) - (1.6*Log(MSUSY/Qmatch)*msu2
      (1,1)*Sqr(g1))/msu2(2,2) + (2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (M1Input))/msu2(2,2) + (10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msu2(2,2
      )*(1 + Sqr(TanBeta))) - (0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(
      TanBeta))/(msu2(2,2)*(1 + Sqr(TanBeta))) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))
      - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(
      TanBeta))))/msu2(2,2) + (0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(
      MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) +
      0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/
      Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1
      )*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(
      M2Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY
      /Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2
      (2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))/
      msu2(2,2))*TCD1F(4)(Sqrt(msq2(2,2)/msu2(2,2))))/msq2(2,2) + 0.6*((0.5*Sqrt(
      msq2(2,2)/msu2(2,2))*msu2(2,2)*Sqr(g1)*((msq2(2,2)*(-1.6*Log(MSUSY/Qmatch)*
      Sqr(g1) + (0.8*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*Log(
      MSUSY/Qmatch)*msd2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msd2(2,2
      )*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/msu2(2,2) +
      (0.8*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)
      *mse2(2,2)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1))/
      msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1))/msu2(2,2) - (0.8*Log(
      MSUSY/Qmatch)*msl2(2,2)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(0,0
      )*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/msu2(2,2) +
      (0.8*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1))/msu2(2,2) - (1.6*Log(MSUSY/Qmatch)
      *msu2(0,0)*Sqr(g1))/msu2(2,2) - (1.6*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1))/
      msu2(2,2) + (2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input))/msu2
      (2,2) + (10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2
      ) + (0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msu2(2,2)*(1 + Sqr(TanBeta
      ))) - (0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(msu2(2,2)*(
      1 + Sqr(TanBeta))) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)
      *msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(
      Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2
      ,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)
      *msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr
      (AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(
      mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,2) + (0.2*Log(MSUSY
      /Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2
      )*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*
      msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1
      ) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2
      )*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*
      Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666*Log(MSUSY/Qmatch
      )*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 +
      Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1
       + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2,2))*TCD1F(3)(Sqrt(msq2(2,2)/msu2(2,2)))
      )/msq2(2,2) + 8.2*Log(MSUSY/Qmatch)*Quad(g1)*TCF(3)(Sqrt(msq2(2,2)/msu2(2,2)
      ))) - 6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2)*TCF(4)(Sqrt(msq2(2,2)/
      msu2(2,2)))) + ((-2*Sin(2*ArcTan(TanBeta))*Sqr(Yu(2,2))*(3*Log(MSUSY/Qmatch)
      *Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)) - 3
      *Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta + 0.1*Cos(2*ArcTan(TanBeta))*Yu(2,
      2)*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 17*Log(MSUSY/Qmatch)*Sqr(g1)*Yu(2,2
      ) - 45*Log(MSUSY/Qmatch)*Sqr(g2)*Yu(2,2) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(
      2,2) + 30*Log(MSUSY/Qmatch)*Sqr(Yd(2,2))*Yu(2,2) + 20*Log(MSUSY/Qmatch)*Sqr(
      Ye(2,2))*Yu(2,2)))*(0.6*Sqr(g1)*TCF(3)(Sqrt(msq2(2,2)/msu2(2,2))) + Sqr(g2)*
      TCF(4)(Sqrt(msq2(2,2)/msu2(2,2)))))/Sqrt(msq2(2,2)*msu2(2,2)))) + 0.25*((-2*
      (AtInput - MuInput/TanBeta)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta
      )))*Sqr(Yu(2,2))*(1.7333333333333334*M1Input*Log(MSUSY/Qmatch)*Sqr(g2) + 6*
      M2Input*Log(MSUSY/Qmatch)*Sqr(g2) + 10.666666666666666*M3Input*Log(MSUSY/
      Qmatch)*Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,
      2)))/Sqr(TanBeta) - (0.6*(-(MuInput*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(TanBeta))
      - 5*MuInput*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(TanBeta) + 10*MuInput*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2)
      )))/Cube(TanBeta))*TCF(5)(Sqrt(msq2(2,2)/msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,
      2)) - Sqr(AtInput - MuInput/TanBeta)*((-0.5*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(
      2*ArcTan(TanBeta)))*Sqr(Yu(2,2))*(msq2(2,2)*(-0.8*Log(MSUSY/Qmatch)*msd2(0,0
      )*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      msd2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) - 0.8*Log(MSUSY/
      Qmatch)*mse2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1) + 0.8*
      Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) + 0.8*Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1
      ) + 0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(0,0
      )*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) - 0.8*Log(MSUSY/Qmatch)*
      msq2(2,2)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) + 1.6*Log(MSUSY/
      Qmatch)*msu2(1,1)*Sqr(g1) + 1.6*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1) -
      2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) -
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) - (0.8*Log(MSUSY/
      Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(TanBeta)) + (0.8*Log(MSUSY/Qmatch)*
      Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)
      *Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)
      *(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch
      )*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(
      0.2*Log(MSUSY/Qmatch)*msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*
      Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      mse2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*mse2(2,2)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*
      Log(MSUSY/Qmatch)*msl2(1,1)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1
      ) + 0.2*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1
      )*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*
      msu2(0,0)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(g1) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*
      Sqr(M1Input) - 6*Log(MSUSY/Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr
      (mAInput))/(1 + Sqr(TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*TCF(5)(Sqrt(msq2(2,2)/msu2(2,2)
      )))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) + ((0.5*Sqrt(msq2(2,2)/
      msu2(2,2))*msu2(2,2)*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr
      (Yu(2,2))*((msq2(2,2)*(-1.6*Log(MSUSY/Qmatch)*Sqr(g1) + (0.8*Log(MSUSY/
      Qmatch)*msd2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(
      g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msd2(2,2)*Sqr(g1))/msu2(2,2) + (0.8*
      Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2
      (1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1))/msu2(2,
      2) - (0.8*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/
      Qmatch)*msl2(1,1)*Sqr(g1))/msu2(2,2) - (0.8*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(
      g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2(0,0)*Sqr(g1))/msu2(2,2) + (0.8*
      Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1))/msu2(2,2) + (0.8*Log(MSUSY/Qmatch)*msq2
      (2,2)*Sqr(g1))/msu2(2,2) - (1.6*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1))/msu2(2,
      2) - (1.6*Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1))/msu2(2,2) + (
      2.1333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input))/msu2(2,2) + (
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) + (0.8*
      Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(msu2(2,2)*(1 + Sqr(TanBeta))) - (
      0.8*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(msu2(2,2)*(1 + Sqr
      (TanBeta))) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)*msq2(2
      ,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2))
      )/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4
      *Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,
      2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput
      )*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)
      *Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,2) + (0.2*Log(MSUSY/Qmatch)*
      msd2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msd2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msd2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(0,0)*Sqr(g1) + 0.2*
      Log(MSUSY/Qmatch)*mse2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*mse2(2,2)*Sqr(g1
      ) - 0.2*Log(MSUSY/Qmatch)*msl2(0,0)*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(1,1
      )*Sqr(g1) - 0.2*Log(MSUSY/Qmatch)*msl2(2,2)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*
      msq2(0,0)*Sqr(g1) + 0.2*Log(MSUSY/Qmatch)*msq2(1,1)*Sqr(g1) + 0.2*Log(MSUSY/
      Qmatch)*msq2(2,2)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(0,0)*Sqr(g1) - 0.4*
      Log(MSUSY/Qmatch)*msu2(1,1)*Sqr(g1) - 0.4*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(g1
      ) - 0.13333333333333333*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(M1Input) - 6*Log(MSUSY
      /Qmatch)*Sqr(g2)*Sqr(M2Input) - 10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)
      *Sqr(M3Input) + (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput))/(1 + Sqr(
      TanBeta)) - (0.2*Log(MSUSY/Qmatch)*Sqr(g1)*Sqr(mAInput)*Sqr(TanBeta))/(1 +
      Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2,2))*TCD1F(5)(Sqrt(msq2(2,2)/msu2(2,2))))/
      msq2(2,2) + (Sqr(Yu(2,2))*((4.92*Log(MSUSY/Qmatch)*Quad(g1) -
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sqr(Cos(2*ArcTan(TanBeta))) -
      (4*Cos(2*ArcTan(TanBeta))*Sin(2*ArcTan(TanBeta))*(0.6*Sqr(g1) + Sqr(g2))*(3*
      Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)
      *Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta) + 0.1*(0.6*Sqr(
      g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta)))*Yu(2,2)*(90*Cube(Yu(2,2))*Log(
      MSUSY/Qmatch) - 17*Log(MSUSY/Qmatch)*Sqr(g1)*Yu(2,2) - 45*Log(MSUSY/Qmatch)*
      Sqr(g2)*Yu(2,2) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(2,2) + 30*Log(MSUSY/
      Qmatch)*Sqr(Yd(2,2))*Yu(2,2) + 20*Log(MSUSY/Qmatch)*Sqr(Ye(2,2))*Yu(2,2)))*
      TCF(5)(Sqrt(msq2(2,2)/msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))))))/Sqr(
      3.141592653589793)) + 0.25*(-(Sqr(Cos(2*ArcTan(TanBeta)))*(40.344*Power6(g1)
      *Sqr(Log(MSUSY/Qmatch)) + 40.111111111111114*Power6(g2)*Sqr(Log(MSUSY/Qmatch
      )))) + (4*Cos(2*ArcTan(TanBeta))*(4.92*Log(MSUSY/Qmatch)*Quad(g1) -
      6.333333333333333*Log(MSUSY/Qmatch)*Quad(g2))*Sin(2*ArcTan(TanBeta))*(3*Log(
      MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yd(2,2)) + Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr
      (Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/TanBeta - (0.6*Sqr(g1) + Sqr(
      g2))*((4*Sqr(Sin(2*ArcTan(TanBeta)))*Sqr(3*Sqr(TanBeta)*Log(MSUSY/Qmatch)*
      Sqr(Yd(2,2)) + Sqr(TanBeta)*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) - 3*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2))))/Sqr(TanBeta) + 2*Cos(2*ArcTan(TanBeta))*((-2*Cos(2*
      ArcTan(TanBeta))*Sqr(3*Sqr(TanBeta)*Log(MSUSY/Qmatch)*Sqr(Yd(2,2)) + Sqr(
      TanBeta)*Log(MSUSY/Qmatch)*Sqr(Ye(2,2)) - 3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))
      /Sqr(TanBeta) + Sin(2*ArcTan(TanBeta))*(-((((-6*Log(MSUSY/Qmatch)*Quad(
      TanBeta)*Sqr(Yd(2,2)))/(1 + Sqr(TanBeta)) - (6*Log(MSUSY/Qmatch)*Sqr(TanBeta
      )*Sqr(Yd(2,2)))/(1 + Sqr(TanBeta)) - (2*Log(MSUSY/Qmatch)*Quad(TanBeta)*Sqr(
      Ye(2,2)))/(1 + Sqr(TanBeta)) - (2*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Ye(2,2)
      ))/(1 + Sqr(TanBeta)) + (6*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/(1 + Sqr(TanBeta)
      ) + (6*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2)))/(1 + Sqr(TanBeta)))*(3*
      TanBeta*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yd(2,2)) + TanBeta*Log(
      MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Ye(2,2)) - (3*Log(MSUSY/Qmatch)*(1 +
      Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta))/(1 + Sqr(TanBeta))) - (2*(22.5*TanBeta
      *Quad(Yd(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)) + 4.5*TanBeta*
      Quad(Ye(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)) - (13.5*Quad(Yu(2
      ,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Cube(TanBeta) - 1.4*
      TanBeta*Sqr(g1)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) - 9*
      TanBeta*Sqr(g2)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) - 16*
      TanBeta*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2)) - 1.8
      *TanBeta*Sqr(g1)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Ye(2,2)) - 3*
      TanBeta*Sqr(g2)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Ye(2,2)) + 9*
      TanBeta*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2))*Sqr(Ye(2,2
      )) + (2.6*Sqr(g1)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/
      TanBeta + (9*Sqr(g2)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))
      /TanBeta + (16*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)
      ))/TanBeta - (9*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yd(2,2))*
      Sqr(Yu(2,2)))/TanBeta - (3*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(
      Ye(2,2))*Sqr(Yu(2,2)))/TanBeta - (6*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr
      (Yu(2,2))*(0.008333333333333333*Sqr(g1) + 0.375*Sqr(g2) - 1.3333333333333333
      *Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(
      TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(
      TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(
      msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(
      msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2
      (2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input
      )/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/
      M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(
      msu2(2,2))/M3Input))/M3Input)))/TanBeta))/(1 + Sqr(TanBeta))))))))/Quad(
      3.141592653589793), LambdaLoopOrder < 3, (0.000244140625*(-1548*Cube(Log(
      MSUSY/Qmatch))*Power8(Yu(2,2)) - 2944*Cube(Log(MSUSY/Qmatch))*Quad(g3)*Quad(
      Yu(2,2)) + 4416*Cube(Log(MSUSY/Qmatch))*Power6(Yu(2,2))*Sqr(g3) - 256*Quad(
      3.141592653589793)*((0.005208333333333333*(-(Quad(Yu(2,2))*Sqr(g3)*((
      0.13333333333333333*(-20*AtInput*TanBeta*Cube(MuInput) + 2*M3Input*TanBeta*
      Cube(MuInput) - 20*MuInput*Cube(AtInput)*Cube(TanBeta) - 12*MuInput*Cube(
      M3Input)*Cube(TanBeta) + 5*Quad(MuInput) - 2*M3Input*Cube(AtInput)*Quad(
      TanBeta) + 12*AtInput*Cube(M3Input)*Quad(TanBeta) + 5*Quad(AtInput)*Quad(
      TanBeta) + 12*Quad(M3Input)*Quad(TanBeta) + 6*M3Input*MuInput*Cube(TanBeta)*
      Sqr(AtInput) + 84*AtInput*MuInput*Cube(TanBeta)*Sqr(M3Input) - 42*Quad(
      TanBeta)*Sqr(AtInput)*Sqr(M3Input) - 6*AtInput*M3Input*Sqr(MuInput)*Sqr(
      TanBeta) + 30*Sqr(AtInput)*Sqr(MuInput)*Sqr(TanBeta) - 42*Sqr(M3Input)*Sqr(
      MuInput)*Sqr(TanBeta))*(90*AtInput*Cube(TanBeta)*Log(MSUSY/Qmatch)*Sqr(g3) +
      160*M3Input*Cube(TanBeta)*Log(MSUSY/Qmatch)*Sqr(g3) - 90*MuInput*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(TanBeta) - 90*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) +
      180*AtInput*TanBeta*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 180*AtInput*Cube(
      TanBeta)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - 90*MuInput*Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Yu(2,2))))/(Power5(M3Input)*Power7(TanBeta)) + (36*Log(Sqrt(
      msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)) - 2*((-
      0.06666666666666667*(AtInput - MuInput/TanBeta)*Log(Sqrt(msq2(2,2)*msu2(2,2)
      )/Sqr(Qmatch))*(6*AtInput*MuInput*TanBeta - 8*M3Input*MuInput*TanBeta - 3*
      Sqr(MuInput) + 8*AtInput*M3Input*Sqr(TanBeta) - 3*Sqr(AtInput)*Sqr(TanBeta)
      + 24*Sqr(M3Input)*Sqr(TanBeta))*(90*AtInput*Cube(TanBeta)*Log(MSUSY/Qmatch)*
      Sqr(g3) + 160*M3Input*Cube(TanBeta)*Log(MSUSY/Qmatch)*Sqr(g3) - 90*MuInput*
      Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(TanBeta) - 90*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu
      (2,2)) + 180*AtInput*TanBeta*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 180*AtInput*
      Cube(TanBeta)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - 90*MuInput*Log(MSUSY/Qmatch)*
      Sqr(TanBeta)*Sqr(Yu(2,2))))/(Power5(TanBeta)*Quad(M3Input)) + (24 - (24*(
      AtInput - MuInput/TanBeta))/M3Input + Cube(AtInput - MuInput/TanBeta)/Cube(
      M3Input) - (4*Sqr(AtInput - MuInput/TanBeta))/Sqr(M3Input))*((Log(Sqrt(msq2(
      2,2)*msu2(2,2))/Sqr(Qmatch))*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*
      Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr
      (TanBeta) - (0.6*(10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log
      (MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta)))/M3Input + (
      AtInput - MuInput/TanBeta)*((6*Log(MSUSY/Qmatch)*Log(Sqrt(msq2(2,2)*msu2(2,2
      ))/Sqr(Qmatch))*Sqr(g3))/M3Input + (0.5*(msq2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(M3Input*msq2(2,2)
      *msu2(2,2))))))) - (((AtInput - MuInput/TanBeta)*(24 + (12*(AtInput -
      MuInput/TanBeta))/M3Input - Cube(AtInput - MuInput/TanBeta)/Cube(M3Input) +
      (2*Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) - (28*Sqr(AtInput -
      MuInput/TanBeta))/Sqr(M3Input)))/M3Input - (2*(AtInput - MuInput/TanBeta)*
      Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*(24 - (24*(AtInput - MuInput/
      TanBeta))/M3Input + Cube(AtInput - MuInput/TanBeta)/Cube(M3Input) - (4*Sqr(
      AtInput - MuInput/TanBeta))/Sqr(M3Input)))/M3Input + 36*Sqr(Log(Sqrt(msq2(2,
      2)*msu2(2,2))/Sqr(Qmatch))))*(-14*Log(MSUSY/Qmatch)*Quad(g3)*Quad(Yu(2,2)) +
      0.2*Cube(Yu(2,2))*Sqr(g3)*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 160*Log(
      MSUSY/Qmatch)*Sqr(g3)*Yu(2,2)))))/Quad(3.141592653589793) + (0.01171875*((
      Power6(Yu(2,2))*(1 + Sqr(TanBeta))*((-50.09958959346384*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/(1 + Sqr(TanBeta)) - (2*(msq2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,
      2)) + (3*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)) + 9.5*((6*Log(MSUSY/Qmatch)*Sqr(MuInput)*(
      1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)) - (
      0.5*Sqr(MuInput)*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*
      Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*
      msu2(2,2)))) + (-13 + 27*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)) - (3*Sqr
      (MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) + (6*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(
      Qmatch))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) + 19.6878144/(1 + Sqr(
      TanBeta)) - (24*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))/(1 + Sqr(TanBeta
      )) - (3*(0.007049244444444251 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*
      Sqr(-MuInput + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta))))/
      (Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))))*((2*(AtInput -
      MuInput/TanBeta)*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12
      *AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      0.6*(10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log(MSUSY/Qmatch
      )*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta)))/Sqrt(msq2(2,2)*msu2(2,2)) - (
      0.5*Sqr(AtInput - MuInput/TanBeta)*(msq2(2,2)*(-10.666666666666666*Log(MSUSY
      /Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2
      )))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2
      (2,2)*msu2(2,2)))) + 13*((6*Log(MSUSY/Qmatch)*Log(Sqrt(msq2(2,2)*msu2(2,2))/
      Sqr(Qmatch))*Sqr(Yu(2,2)))/(1 + Sqr(TanBeta)) + (0.5*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*(1 + Sqr(TanBeta)))) - 3*((6*Log(MSUSY/
      Qmatch)*Sqr(Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*Sqr(Yu(2,2)))/(1 +
      Sqr(TanBeta)) + (Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*(1 + Sqr(TanBeta)))) - 6*((-0.5*Log(Sqrt(
      msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(MuInput)*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) + ((6*Log(MSUSY/
      Qmatch)*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (0.5*Sqr(MuInput)*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)))/Sqrt(msq2(2,2)*msu2(2,2))) + 0.5*((6*Log(
      MSUSY/Qmatch)*(-1 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*Power6(
      MuInput - AtInput*TanBeta)*Sqr(Yu(2,2)))/(Power3(Sqrt(msq2(2,2)*msu2(2,2)))*
      Quad(TanBeta)*(1 + Sqr(TanBeta))) + ((6*(-1 + Log(Sqrt(msq2(2,2)*msu2(2,2))/
      Sqr(Qmatch)))*Power5(MuInput - AtInput*TanBeta)*((3*AtInput*Log(MSUSY/Qmatch
      )*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta + (3*MuInput*Log(MSUSY/Qmatch)*(1
       + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - TanBeta*(10.666666666666666*
      M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(Power3(Sqrt(msq2(2,2)*msu2(2,2)))*
      Quad(TanBeta)) + Power6(MuInput - AtInput*TanBeta)*((-1.5*(-1 + Log(Sqrt(
      msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(msq2(2,2)*(-10.666666666666666*Log(MSUSY
      /Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2
      )))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Power3(
      Sqrt(msq2(2,2)*msu2(2,2)))*Quad(TanBeta)) + ((12*Log(MSUSY/Qmatch)*(-1 + Log
      (Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Power6(TanBeta) + (0.5*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr
      (g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Quad(TanBeta)))
      /Power3(Sqrt(msq2(2,2)*msu2(2,2)))))/(1 + Sqr(TanBeta))) + 3*((2*(-MuInput +
      AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*(
      0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(12*
      MuInput*Cot(2*ArcTan(TanBeta))*Csc(2*ArcTan(TanBeta))*Log(MSUSY/Qmatch)*Sqr(
      Yu(2,2)) - (3*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      TanBeta - (3*MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + TanBeta*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (
      12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))
      )/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) + Sqr(-MuInput
       + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*((6*Log(MSUSY
      /Qmatch)*(0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*
      Sqr(Yu(2,2)))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) +
      ((-0.5*(0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(
      msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta
      )) + ((6*Log(MSUSY/Qmatch)*(0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,
      2))/Sqr(Qmatch)))*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Quad(TanBeta) + (0.5*(
      msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqr(TanBeta)))/Sqrt(msq2(2,2)*msu2(2,
      2)))/(1 + Sqr(TanBeta)))) + (Sqr(AtInput - MuInput/TanBeta)*((118.1268864*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/(1 + Sqr(TanBeta)) + (13.5*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)) - 3*((6*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)) - (
      0.5*Sqr(MuInput)*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*
      Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*
      msu2(2,2)))) - 24*((6*Log(MSUSY/Qmatch)*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(
      Qmatch))*Sqr(Yu(2,2)))/(1 + Sqr(TanBeta)) + (0.5*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*(1 + Sqr(TanBeta)))) + 6*((-0.5*Log(Sqrt(
      msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(MuInput)*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) + ((6*Log(MSUSY/
      Qmatch)*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (0.5*Sqr(MuInput)*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)))/Sqrt(msq2(2,2)*msu2(2,2))) - 3*((2*(-
      MuInput + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*(
      0.007049244444444251 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(12*
      MuInput*Cot(2*ArcTan(TanBeta))*Csc(2*ArcTan(TanBeta))*Log(MSUSY/Qmatch)*Sqr(
      Yu(2,2)) - (3*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      TanBeta - (3*MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + TanBeta*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (
      12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))
      )/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) + Sqr(-MuInput
       + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*((6*Log(MSUSY
      /Qmatch)*(0.007049244444444251 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))
      *Sqr(Yu(2,2)))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) +
      ((-0.5*(0.007049244444444251 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(
      msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta
      )) + ((6*Log(MSUSY/Qmatch)*(0.007049244444444251 + Log(Sqrt(msq2(2,2)*msu2(2
      ,2))/Sqr(Qmatch)))*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Quad(TanBeta) + (0.5*(
      msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqr(TanBeta)))/Sqrt(msq2(2,2)*msu2(2,
      2)))/(1 + Sqr(TanBeta))))))/Sqrt(msq2(2,2)*msu2(2,2)) + 12*(((-MuInput +
      AtInput*TanBeta)*(0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(
      Qmatch)))*(12*MuInput*Cot(2*ArcTan(TanBeta))*Csc(2*ArcTan(TanBeta))*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))
      *Sqr(Yu(2,2)))/TanBeta - (3*MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr
      (Yu(2,2)))/Sqr(TanBeta) + TanBeta*(10.666666666666666*M3Input*Log(MSUSY/
      Qmatch)*Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,
      2)))/Sqr(TanBeta))))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(
      TanBeta))) + (-MuInput + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(
      TanBeta)))*((6*(-MuInput + AtInput*TanBeta)*Log(MSUSY/Qmatch)*(
      0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*Sqr(Yu(2,2
      )))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) + (((
      0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*((-3*
      AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta - (3*
      MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      TanBeta*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12*AtInput*
      Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(Sqrt(
      msq2(2,2)*msu2(2,2))*Sqr(TanBeta)) + (-MuInput + AtInput*TanBeta)*((-0.5*(
      0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(msq2(2,2)
      *(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)) + (
      (6*Log(MSUSY/Qmatch)*(0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/
      Sqr(Qmatch)))*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Quad(TanBeta) + (0.5*(msq2(2,
      2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(
      MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)
      *msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqr(TanBeta)))/Sqrt(msq2(2,2)*msu2(2,
      2))))/(1 + Sqr(TanBeta)))) - 2*((Cube(-MuInput + AtInput*TanBeta)*(-
      0.4373952 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(12*MuInput*Cot(2*
      ArcTan(TanBeta))*Csc(2*ArcTan(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*
      AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/TanBeta - (3*
      MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      TanBeta*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12*AtInput*
      Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2
      )*msu2(2,2)*Quad(TanBeta)*(1 + Sqr(TanBeta))) + (-MuInput + AtInput*TanBeta
      + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*((6*Cube(-MuInput + AtInput*
      TanBeta)*Log(MSUSY/Qmatch)*(-0.4373952 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(
      Qmatch)))*Sqr(Yu(2,2)))/(msq2(2,2)*msu2(2,2)*Quad(TanBeta)*(1 + Sqr(TanBeta)
      )) + ((3*(-0.4373952 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*Sqr(-
      MuInput + AtInput*TanBeta)*((-3*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))
      *Sqr(Yu(2,2)))/TanBeta - (3*MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr
      (Yu(2,2)))/Sqr(TanBeta) + TanBeta*(10.666666666666666*M3Input*Log(MSUSY/
      Qmatch)*Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,
      2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Quad(TanBeta)) + Cube(-MuInput +
      AtInput*TanBeta)*(((-0.4373952 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))
      *((10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) - 4*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2))
      )/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4
      *Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2
      ,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,
      2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2
      )))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)
      ))/(msu2(2,2)*Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Quad(TanBeta)) + (((-
      0.4373952 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*((10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msq2(2,2) - 2*Log(MSUSY/Qmatch)*Sqr
      (Yu(2,2)) - (2*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/msq2(2,2) - (2*Log(
      MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/msq2(2,2) + (2*Log(MSUSY/Qmatch)*
      Sqr(MuInput)*Sqr(Yu(2,2)))/msq2(2,2) - (2*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/
      Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta)) - (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta)) - (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta))))/(msq2(2,2)*Quad(TanBeta)) + ((12*Log(MSUSY/Qmatch)*(-0.4373952 +
      Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))
      /Power6(TanBeta) + (0.5*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*
      Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Quad(TanBeta)))
      /msq2(2,2))/msu2(2,2)))/(1 + Sqr(TanBeta)))) + 0.25*((25 - 26*Log(Sqrt(msq2(
      2,2)*msu2(2,2))/Sqr(Qmatch)) + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)) - (4*
      Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2
      (2,2)) - 25/(1 + Sqr(TanBeta)) + (24*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(
      Qmatch)))/(1 + Sqr(TanBeta)) + (2*(-0.04145706666666667 + Log(Sqrt(msq2(2,2)
      *msu2(2,2))/Sqr(Qmatch)))*Sqr(-MuInput + AtInput*TanBeta + 2*MuInput*TanBeta
      *Csc(2*ArcTan(TanBeta))))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(
      TanBeta))))*(Quad(AtInput - MuInput/TanBeta)*(((10.666666666666666*Log(MSUSY
      /Qmatch)*Sqr(g3)*Sqr(M3Input))/msq2(2,2) - 2*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))
      - (2*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/msq2(2,2) - (2*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/msq2(2,2) + (2*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*Sqr(Yu(2,2)))/msq2(2,2) - (2*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta)) - (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta)) - (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(
      TanBeta)))/(msq2(2,2)*msu2(2,2)) + ((10.666666666666666*Log(MSUSY/Qmatch)*
      Sqr(g3)*Sqr(M3Input))/msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(
      MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(
      2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(
      MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(
      MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)))/(msq2(2,2
      )*msu2(2,2))) + (4*Cube(AtInput - MuInput/TanBeta)*(10.666666666666666*
      M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (0.6*(10*MuInput*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube
      (TanBeta)))/(msq2(2,2)*msu2(2,2))) + (Quad(AtInput - MuInput/TanBeta)*((-150
      *Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/(1 + Sqr(TanBeta)) + (6*Log(MSUSY/Qmatch)*
      Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr
      (TanBeta)) - (13*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*
      Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)) - (0.5*Sqr(
      MuInput)*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*
      Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))
      ) + 24*((6*Log(MSUSY/Qmatch)*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(
      Yu(2,2)))/(1 + Sqr(TanBeta)) + (0.5*(msq2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,
      2)*(1 + Sqr(TanBeta)))) - 4*((-0.5*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)
      )*Sqr(MuInput)*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr
      (M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (
      4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu
      (2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr
      (g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2
      ))) + ((6*Log(MSUSY/Qmatch)*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(
      MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (0.5*Sqr(MuInput)*(
      msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)))/Sqrt(msq2(2,2)*msu2(2,2))) + 2*((2*
      (-MuInput + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*(-
      0.04145706666666667 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(12*
      MuInput*Cot(2*ArcTan(TanBeta))*Csc(2*ArcTan(TanBeta))*Log(MSUSY/Qmatch)*Sqr(
      Yu(2,2)) - (3*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      TanBeta - (3*MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + TanBeta*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (
      12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))
      )/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) + Sqr(-MuInput
       + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*((6*Log(MSUSY
      /Qmatch)*(-0.04145706666666667 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))
      *Sqr(Yu(2,2)))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) +
      ((-0.5*(-0.04145706666666667 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*(
      msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta
      )) + ((6*Log(MSUSY/Qmatch)*(-0.04145706666666667 + Log(Sqrt(msq2(2,2)*msu2(2
      ,2))/Sqr(Qmatch)))*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Quad(TanBeta) + (0.5*(
      msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(
      M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2
      *Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqr(TanBeta)))/Sqrt(msq2(2,2)*msu2(2,
      2)))/(1 + Sqr(TanBeta))))))/(msq2(2,2)*msu2(2,2)))))/Sqr(TanBeta) + (-0.5 -
      4*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)) + (9.5*Sqr(MuInput))/Sqrt(msq2(
      2,2)*msu2(2,2)) - (6*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(MuInput)
      )/Sqrt(msq2(2,2)*msu2(2,2)) - 8.34993159891064/(1 + Sqr(TanBeta)) + (13*Log(
      Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))/(1 + Sqr(TanBeta)) - (2*(-MuInput +
      AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta)))*Cube(-MuInput +
      AtInput*TanBeta)*(-0.4373952 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))))/
      (msq2(2,2)*msu2(2,2)*Quad(TanBeta)*(1 + Sqr(TanBeta))) + (0.5*(-1 + Log(Sqrt
      (msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*Power6(MuInput - AtInput*TanBeta))/(
      Power3(Sqrt(msq2(2,2)*msu2(2,2)))*Quad(TanBeta)*(1 + Sqr(TanBeta))) + (12*(-
      MuInput + AtInput*TanBeta)*(-MuInput + AtInput*TanBeta + 2*MuInput*TanBeta*
      Csc(2*ArcTan(TanBeta)))*(0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))
      /Sqr(Qmatch))))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta)))
      + (3*(0.04173653333333333 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*Sqr(
      -MuInput + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta))))/(
      Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta))) + (0.25*Quad(
      AtInput - MuInput/TanBeta)*(25 - 26*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch
      )) + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)) - (4*Log(Sqrt(msq2(2,2)*msu2(2,2
      ))/Sqr(Qmatch))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) - 25/(1 + Sqr(
      TanBeta)) + (24*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))/(1 + Sqr(TanBeta
      )) + (2*(-0.04145706666666667 + Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch)))*
      Sqr(-MuInput + AtInput*TanBeta + 2*MuInput*TanBeta*Csc(2*ArcTan(TanBeta))))/
      (Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1 + Sqr(TanBeta)))))/(msq2(2,2)*
      msu2(2,2)) + (Sqr(AtInput - MuInput/TanBeta)*(-13 + 27*Log(Sqrt(msq2(2,2)*
      msu2(2,2))/Sqr(Qmatch)) - (3*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) + (6*
      Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2
      (2,2)) + 19.6878144/(1 + Sqr(TanBeta)) - (24*Log(Sqrt(msq2(2,2)*msu2(2,2))/
      Sqr(Qmatch)))/(1 + Sqr(TanBeta)) - (3*(0.007049244444444251 + Log(Sqrt(msq2(
      2,2)*msu2(2,2))/Sqr(Qmatch)))*Sqr(-MuInput + AtInput*TanBeta + 2*MuInput*
      TanBeta*Csc(2*ArcTan(TanBeta))))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(TanBeta)*(1
      + Sqr(TanBeta)))))/Sqrt(msq2(2,2)*msu2(2,2)) + 3*Sqr(Log(Sqrt(msq2(2,2)*msu2
      (2,2))/Sqr(Qmatch))) - (3*Sqr(Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(Qmatch))))/(
      1 + Sqr(TanBeta)))*((-6*Log(MSUSY/Qmatch)*Power8(Yu(2,2))*(1 + Sqr(TanBeta))
      )/Sqr(TanBeta) + (1 + Sqr(TanBeta))*((6*Log(MSUSY/Qmatch)*Power8(Yu(2,2))*(1
       + Sqr(TanBeta)))/Quad(TanBeta) + (0.3*Power5(Yu(2,2))*(90*Cube(Yu(2,2))*Log
      (MSUSY/Qmatch) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(2,2)))/Sqr(TanBeta)))))/
      Quad(3.141592653589793)) - 3*((0.2*Cube(Yu(2,2))*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))*(90*Cube(Yu(2,2))*Log(MSUSY/
      Qmatch) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(2,2)))/msu2(2,2) + Log(msu2(2,2)/
      Sqr(Qmatch))*(0.01*Sqr(Yu(2,2))*Sqr(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 160
      *Sqr(g3)*Log(MSUSY/Qmatch)*Yu(2,2)) + 2*Sqr(Yu(2,2))*(0.0025*Sqr(90*Cube(Yu(
      2,2))*Log(MSUSY/Qmatch) - 160*Sqr(g3)*Log(MSUSY/Qmatch)*Yu(2,2)) + 2*Yu(2,2)
      *(30.375*Power5(Yu(2,2))*Sqr(Log(MSUSY/Qmatch)) - 72*Cube(Yu(2,2))*Sqr(g3)*
      Sqr(Log(MSUSY/Qmatch)) + 88*Quad(g3)*Sqr(Log(MSUSY/Qmatch))*Yu(2,2)))) + 0.5
      *Quad(Yu(2,2))*((((10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input)
      )/msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)*msq2(2
      ,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2))
      )/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4
      *Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,
      2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput
      )*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)
      *Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)))*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/msu2(2,2) + (2*((0.1111111111111111*
      Log(MSUSY/Qmatch)*(-6 + Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(
      Qmatch)) + Log(msd2(2,2)/Sqr(Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log
      (msq2(1,1)/Sqr(Qmatch)) + 2*Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(
      Qmatch)) + Log(msu2(1,1)/Sqr(Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(
      Sqr(M3Input)/Sqr(Qmatch)))*Quad(g3)*Sqr(M3Input))/Sqr(3.141592653589793) +
      96*Quad(g3)*Sqr(M3Input)*Sqr(Log(MSUSY/Qmatch)) + (12*Quad(Yu(2,2))*Sqr(
      mAInput)*Sqr(Log(MSUSY/Qmatch)))/Sqr(TanBeta) + (36*Quad(Yu(2,2))*Sqr(
      mAInput)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (36*msq2
      (2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(
      TanBeta) + (36*msu2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch)))/Quad(TanBeta) + (84*Quad(Yu(2,2))*Sqr(AtInput)*Sqr(1 + Sqr(TanBeta
      ))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) - (48*Quad(Yu(2,2))*Sqr(MuInput)*
      Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) - (
      21.333333333333332*Sqr(g3)*Sqr(mAInput)*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))
      /Sqr(TanBeta) + (42.666666666666664*AtInput*M3Input*Sqr(g3)*(1 + Sqr(TanBeta
      ))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*
      msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) - (21.333333333333332*msu2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*Sqr(
      AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr
      (TanBeta) - (42.666666666666664*Sqr(g3)*Sqr(M3Input)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (21.333333333333332*Sqr(g3)*
      Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2))*(-
      1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/
      Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(
      TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(
      msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))
      *TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF
      (6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(
      Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(
      Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2
      ))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))/Sqr(TanBeta) + (8*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(
      g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(
      TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(
      TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(
      msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(
      msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2
      (2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input
      )/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/
      M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(
      msu2(2,2))/M3Input))/M3Input)))/Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msu2(2,2
      )*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2)
      )*((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(
      Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(
      AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(
      msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/
      Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(
      TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(
      6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput -
      MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/
      M3Input)))/Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*
      Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(
      Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (8*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta)))/msu2(2,2))) - 3*((0.2*Cube(Yu(2,2))*(-10.666666666666666*Log
      (MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))*(90*Cube(Yu(2,2))*Log(MSUSY/
      Qmatch) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(2,2)))/msq2(2,2) + Log(msq2(2,2)/
      Sqr(Qmatch))*(0.01*Sqr(Yu(2,2))*Sqr(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 160
      *Sqr(g3)*Log(MSUSY/Qmatch)*Yu(2,2)) + 2*Sqr(Yu(2,2))*(0.0025*Sqr(90*Cube(Yu(
      2,2))*Log(MSUSY/Qmatch) - 160*Sqr(g3)*Log(MSUSY/Qmatch)*Yu(2,2)) + 2*Yu(2,2)
      *(30.375*Power5(Yu(2,2))*Sqr(Log(MSUSY/Qmatch)) - 72*Cube(Yu(2,2))*Sqr(g3)*
      Sqr(Log(MSUSY/Qmatch)) + 88*Quad(g3)*Sqr(Log(MSUSY/Qmatch))*Yu(2,2)))) + 0.5
      *Quad(Yu(2,2))*((((10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input)
      )/msq2(2,2) - 2*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (2*Log(MSUSY/Qmatch)*msu2(2
      ,2)*Sqr(Yu(2,2)))/msq2(2,2) - (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2))
      )/msq2(2,2) + (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msq2(2,2) - (2
      *Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*msu2(2,
      2)*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(TanBeta)) - (2*Log(MSUSY/Qmatch)*Sqr(AtInput
      )*Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(TanBeta)) - (2*Log(MSUSY/Qmatch)*Sqr(mAInput)
      *Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(TanBeta)) + (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/(msq2(2,2)*Sqr(TanBeta)))*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/msq2(2,2) + (2*((0.1111111111111111*
      Log(MSUSY/Qmatch)*(-6 + Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(
      Qmatch)) + Log(msd2(2,2)/Sqr(Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log
      (msq2(1,1)/Sqr(Qmatch)) + 2*Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(
      Qmatch)) + Log(msu2(1,1)/Sqr(Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(
      Sqr(M3Input)/Sqr(Qmatch)))*Quad(g3)*Sqr(M3Input))/Sqr(3.141592653589793) +
      96*Quad(g3)*Sqr(M3Input)*Sqr(Log(MSUSY/Qmatch)) + (6*Quad(Yu(2,2))*Sqr(
      mAInput)*Sqr(Log(MSUSY/Qmatch)))/Sqr(TanBeta) + (18*Quad(Yu(2,2))*Sqr(
      mAInput)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (18*msq2
      (2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(
      TanBeta) + (18*msu2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch)))/Quad(TanBeta) + (42*Quad(Yu(2,2))*Sqr(AtInput)*Sqr(1 + Sqr(TanBeta
      ))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) - (24*Quad(Yu(2,2))*Sqr(MuInput)*
      Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) - (
      10.666666666666666*Sqr(g3)*Sqr(mAInput)*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))
      /Sqr(TanBeta) + (21.333333333333332*AtInput*M3Input*Sqr(g3)*(1 + Sqr(TanBeta
      ))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (10.666666666666666*
      msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) - (10.666666666666666*msu2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (10.666666666666666*Sqr(
      AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr
      (TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(M3Input)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (10.666666666666666*Sqr(g3)*
      Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2))*(-
      1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/
      Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(
      TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(
      msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))
      *TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF
      (6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(
      Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(
      Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2
      ))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(
      g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(
      TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(
      TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(
      msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(
      msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2
      (2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input
      )/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/
      M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(
      msu2(2,2))/M3Input))/M3Input)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2
      )*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2)
      )*((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(
      Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(
      AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(
      msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/
      Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(
      TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(
      6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput -
      MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/
      M3Input)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*
      Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(
      Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta)))/msq2(2,2))) - 6*(((2*(AtInput - MuInput/TanBeta)*Quad(Yu(2,2
      ))*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12*AtInput*Log(
      MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (0.6*(10*
      MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log(MSUSY/Qmatch)*Sqr(
      TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta)))/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(
      AtInput - MuInput/TanBeta)*((-0.5*Quad(Yu(2,2))*(msq2(2,2)*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) + (0.2*Cube(Yu(2
      ,2))*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(
      2,2)))/Sqrt(msq2(2,2)*msu2(2,2))))*((0.5*Sqrt(msq2(2,2)/msu2(2,2))*msu2(2,2)
      *((msq2(2,2)*((10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/
      msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)*msq2(2,2
      )*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/
      msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,2
      )*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)
      *Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,2) + (-10.666666666666666*
      Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2,2))*TCD1F(1)(
      Sqrt(msq2(2,2)/msu2(2,2))))/msq2(2,2) + 0.08333333333333333*((-2*(AtInput -
      MuInput/TanBeta)*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12
      *AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      0.6*(10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log(MSUSY/Qmatch
      )*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta))*TCF(2)(Sqrt(msq2(2,2)/msu2(2,2)
      )))/Sqrt(msq2(2,2)*msu2(2,2)) - Sqr(AtInput - MuInput/TanBeta)*((0.5*Sqrt(
      msq2(2,2)/msu2(2,2))*msu2(2,2)*((msq2(2,2)*((10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) -
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch
      )*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) -
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*
      Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2
      (2,2) + (-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log
      (MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch
      )*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))/msu2(2,2))*TCD1F(2)(Sqrt(msq2(2,2)/msu2(2,2))))/(msq2(2,2)*Sqrt(
      msq2(2,2)*msu2(2,2))) - (0.5*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*TCF(2)(Sqrt(msq2(2,2)/msu2(2,2)
      )))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2)))))) + (TCF(1)(Sqrt(msq2(2
      ,2)/msu2(2,2))) - (0.08333333333333333*Sqr(AtInput - MuInput/TanBeta)*TCF(2)
      (Sqrt(msq2(2,2)/msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)))*(2*(AtInput -
      MuInput/TanBeta)*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*Sqr(g3) + (12
      *AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      0.6*(10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log(MSUSY/Qmatch
      )*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta))*((-0.5*Quad(Yu(2,2))*(msq2(2,2)
      *(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) +
      msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) + (0.2*Cube(Yu(2
      ,2))*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 160*Log(MSUSY/Qmatch)*Sqr(g3)*Yu(
      2,2)))/Sqrt(msq2(2,2)*msu2(2,2))) + (Quad(Yu(2,2))*(Sqr(10.666666666666666*
      Sqr(g3)*M3Input*Log(MSUSY/Qmatch) + (12*AtInput*(1 + Sqr(TanBeta))*Log(MSUSY
      /Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (0.6*(10*MuInput*Log(MSUSY/Qmatch)*Sqr
      (Yu(2,2)) + 10*MuInput*Sqr(TanBeta)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/Cube(
      TanBeta)) + 2*(AtInput - MuInput/TanBeta)*((-0.1111111111111111*M3Input*Log(
      MSUSY/Qmatch)*(-6 + Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch))
      + Log(msd2(2,2)/Sqr(Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1
      )/Sqr(Qmatch)) + 2*Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) +
      Log(msu2(1,1)/Sqr(Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input
      )/Sqr(Qmatch)))*Quad(g3))/Sqr(3.141592653589793) - 64*M3Input*Quad(g3)*Sqr(
      Log(MSUSY/Qmatch)) + (144*AtInput*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch)))/Quad(TanBeta) - (64*AtInput*Sqr(g3)*(1 + Sqr(TanBeta))*
      Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (64*M3Input*Sqr(g3)*(1 +
      Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (3*MuInput
      *Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(3*Log(MSUSY/Qmatch)*Sqr(
      Yu(2,2)) + (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta)))/Cube(TanBeta) +
      (24*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-
      1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/
      Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(
      TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(
      msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))
      *TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF
      (6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(
      Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(
      Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2
      ))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))/Sqr(TanBeta) - ((-0.5*
      MuInput*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*(-45*Quad(Yu(2,2)) - 45*
      Quad(Yu(2,2))*Sqr(TanBeta) + 32*Sqr(g3)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Quad(
      TanBeta) + (6*MuInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-
      1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/
      Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(
      TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(
      msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))
      *TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF
      (6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(
      Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(
      Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2
      ))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))/Sqr(TanBeta))/TanBeta - (
      MuInput*((-3*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-3*Log(MSUSY
      /Qmatch)*Sqr(Yu(2,2)) - (3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta)))/
      Sqr(TanBeta) - ((0.5*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*(-27*Quad(Yu(
      2,2)) - 27*Quad(Yu(2,2))*Sqr(TanBeta) + 32*Sqr(g3)*Sqr(TanBeta)*Sqr(Yu(2,2))
      ))/Cube(TanBeta) - (6*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-
      1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/
      Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(
      TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(
      msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))
      *TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF
      (6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(
      Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(
      Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2
      ))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))/TanBeta)/TanBeta))/TanBeta))
      )/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AtInput - MuInput/TanBeta)*((-0.1*Cube(Yu(
      2,2))*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input)
      + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(
      MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*
      Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      )/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*
      Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta)))*(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 160*Log(
      MSUSY/Qmatch)*Sqr(g3)*Yu(2,2)))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2
      ))) + (Quad(Yu(2,2))*((0.015*Sqr(90*Cube(Yu(2,2))*Log(MSUSY/Qmatch) - 160*
      Sqr(g3)*Log(MSUSY/Qmatch)*Yu(2,2)))/Sqr(Yu(2,2)) + (4*(30.375*Power5(Yu(2,2)
      )*Sqr(Log(MSUSY/Qmatch)) - 72*Cube(Yu(2,2))*Sqr(g3)*Sqr(Log(MSUSY/Qmatch)) +
      88*Quad(g3)*Sqr(Log(MSUSY/Qmatch))*Yu(2,2)))/Yu(2,2)))/Sqrt(msq2(2,2)*msu2(2
      ,2)) + (0.5*Quad(Yu(2,2))*((0.75*Sqr(msu2(2,2)*(-10.666666666666666*Sqr(g3)*
      Sqr(M3Input)*Log(MSUSY/Qmatch) + (2*Sqr(mAInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,
      2)))/Sqr(TanBeta) + (2*Sqr(AtInput)*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr
      (Yu(2,2)))/Sqr(TanBeta) - (2*Sqr(MuInput)*(1 + Sqr(TanBeta))*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)
      *msq2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msq2(2,2)*(-
      10.666666666666666*Sqr(g3)*Sqr(M3Input)*Log(MSUSY/Qmatch) + (4*Sqr(mAInput)*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*(1 + Sqr(
      TanBeta))*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*(1 +
      Sqr(TanBeta))*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(Sqr
      (msq2(2,2))*Sqr(msu2(2,2))) - ((-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3
      )*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta))*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*
      Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta)) + msq2(2,2)*((0.1111111111111111*Log(MSUSY/
      Qmatch)*(-6 + Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(
      msd2(2,2)/Sqr(Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(
      Qmatch)) + 2*Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(
      msu2(1,1)/Sqr(Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/
      Sqr(Qmatch)))*Quad(g3)*Sqr(M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*
      Sqr(M3Input)*Sqr(Log(MSUSY/Qmatch)) + (12*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log
      (MSUSY/Qmatch)))/Sqr(TanBeta) + (36*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(
      TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (36*msq2(2,2)*Quad(Yu(2,2)
      )*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (36*msu2(2,2
      )*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta)
      + (84*Quad(Yu(2,2))*Sqr(AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)
      ))/Quad(TanBeta) - (48*Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch)))/Quad(TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(mAInput)
      *Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (42.666666666666664*
      AtInput*M3Input*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2
      )))/Sqr(TanBeta) - (21.333333333333332*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*
      Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*msu2
      (2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (21.333333333333332*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (42.666666666666664*Sqr(g3)*
      Sqr(M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (21.333333333333332*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(
      mAInput)*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-
      1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/
      Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (8*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta)) + msu2(2,2)*((0.1111111111111111*Log(MSUSY/Qmatch)*(-6 + Log(
      msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2)/Sqr(
      Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(Qmatch)) + 2*
      Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(msu2(1,1)/Sqr(
      Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/Sqr(Qmatch)))*
      Quad(g3)*Sqr(M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*Sqr(M3Input)*Sqr
      (Log(MSUSY/Qmatch)) + (6*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log(MSUSY/Qmatch)))/
      Sqr(TanBeta) + (18*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch)))/Quad(TanBeta) + (18*msq2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(
      TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (18*msu2(2,2)*Quad(Yu(2,2)
      )*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (42*Quad(Yu(
      2,2))*Sqr(AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(
      TanBeta) - (24*Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch)))/Quad(TanBeta) - (10.666666666666666*Sqr(g3)*Sqr(mAInput)*Sqr
      (Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (21.333333333333332*AtInput
      *M3Input*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr
      (TanBeta) - (10.666666666666666*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log
      (MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (10.666666666666666*msu2(2,2)*
      Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (10.666666666666666*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(
      M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) + (10.666666666666666*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput
      )*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*
      Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(
      Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta)))/(msq2(2,2)*msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)))) + (Quad(
      Yu(2,2))*Sqr(AtInput - MuInput/TanBeta)*((0.125*msu2(2,2)*Sqr((msq2(2,2)*((
      10.666666666666666*Sqr(g3)*Sqr(M3Input)*Log(MSUSY/Qmatch))/msu2(2,2) - 4*Log
      (MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Sqr(AtInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))
      )/msu2(2,2) + (4*Sqr(MuInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/msu2(2,2) - (4
      *Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Sqr(AtInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))
      /(msu2(2,2)*Sqr(TanBeta)) - (4*Sqr(mAInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/
      (msu2(2,2)*Sqr(TanBeta)) + (4*Sqr(MuInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/(
      msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2
      (2,2)*Sqr(TanBeta))))/msu2(2,2) + (-10.666666666666666*Sqr(g3)*Sqr(M3Input)*
      Log(MSUSY/Qmatch) + (2*Sqr(mAInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Sqr(AtInput)*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2))
      )/Sqr(TanBeta) - (2*Sqr(MuInput)*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu
      (2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr
      (Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*msu2(2,2)*
      Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2,2))*TCD2F(1)(Sqrt(msq2(2,2)/msu2(2,2))))/
      msq2(2,2) + 0.5*Sqrt(msq2(2,2)/msu2(2,2))*TCD1F(1)(Sqrt(msq2(2,2)/msu2(2,2))
      )*((-0.25*Sqr(msu2(2,2))*Sqr((msq2(2,2)*((10.666666666666666*Sqr(g3)*Sqr(
      M3Input)*Log(MSUSY/Qmatch))/msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (
      4*Sqr(AtInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Sqr(MuInput)*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*
      Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) -
      (4*Sqr(AtInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (
      4*Sqr(mAInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4
      *Sqr(MuInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*
      Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,
      2) + (-10.666666666666666*Sqr(g3)*Sqr(M3Input)*Log(MSUSY/Qmatch) + (2*Sqr(
      mAInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Sqr(AtInput)*(1 +
      Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Sqr(MuInput)
      *(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 +
      Sqr(TanBeta))*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1
       + Sqr(TanBeta))*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta))/
      msu2(2,2)))/Sqr(msq2(2,2)) + (msu2(2,2)*((((10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) -
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch
      )*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) -
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*
      Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)))*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/
      msu2(2,2) + ((0.1111111111111111*Log(MSUSY/Qmatch)*(-6 + Log(msd2(0,0)/Sqr(
      Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2)/Sqr(Qmatch)) + 2*Log(
      msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(Qmatch)) + 2*Log(msq2(2,2)/Sqr(
      Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(msu2(1,1)/Sqr(Qmatch)) + Log(
      msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/Sqr(Qmatch)))*Quad(g3)*Sqr(
      M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*Sqr(M3Input)*Sqr(Log(MSUSY/
      Qmatch)) + (6*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log(MSUSY/Qmatch)))/Sqr(TanBeta
      ) + (18*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))
      )/Quad(TanBeta) + (18*msq2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch)))/Quad(TanBeta) + (18*msu2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(
      TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (42*Quad(Yu(2,2))*Sqr(
      AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) - (24*
      Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/
      Quad(TanBeta) - (10.666666666666666*Sqr(g3)*Sqr(mAInput)*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (21.333333333333332*AtInput*M3Input*
      Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (10.666666666666666*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (10.666666666666666*msu2(2,2)*Sqr(g3)*
      (1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      10.666666666666666*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(
      M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) + (10.666666666666666*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput
      )*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*
      Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(
      Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta))/msu2(2,2) + (msq2(2,2)*((((-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) + 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) +
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch
      )*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*
      Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)))*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/
      msu2(2,2) - ((0.1111111111111111*Log(MSUSY/Qmatch)*(-6 + Log(msd2(0,0)/Sqr(
      Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2)/Sqr(Qmatch)) + 2*Log(
      msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(Qmatch)) + 2*Log(msq2(2,2)/Sqr(
      Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(msu2(1,1)/Sqr(Qmatch)) + Log(
      msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/Sqr(Qmatch)))*Quad(g3)*Sqr(
      M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*Sqr(M3Input)*Sqr(Log(MSUSY/
      Qmatch)) + (12*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log(MSUSY/Qmatch)))/Sqr(
      TanBeta) + (36*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch)))/Quad(TanBeta) + (36*msq2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*
      Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (36*msu2(2,2)*Quad(Yu(2,2))*Sqr(1 +
      Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (84*Quad(Yu(2,2))*Sqr(
      AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) - (48*
      Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/
      Quad(TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(mAInput)*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (42.666666666666664*AtInput*M3Input*
      Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (21.333333333333332*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*msu2(2,2)*Sqr(g3)*
      (1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      21.333333333333332*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (42.666666666666664*Sqr(g3)*Sqr(
      M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) + (21.333333333333332*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(mAInput
      )*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*
      Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(
      Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (8*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta))/msu2(2,2)))/msu2(2,2)))/msq2(2,2)) + 0.08333333333333333*(-2*
      (AtInput - MuInput/TanBeta)*(10.666666666666666*M3Input*Log(MSUSY/Qmatch)*
      Sqr(g3) + (12*AtInput*Log(MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr
      (TanBeta) - (0.6*(10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Log
      (MSUSY/Qmatch)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta))*((0.5*Sqrt(msq2(2,
      2)/msu2(2,2))*msu2(2,2)*((msq2(2,2)*((10.666666666666666*Log(MSUSY/Qmatch)*
      Sqr(g3)*Sqr(M3Input))/msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(
      MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(
      AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(
      2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(
      MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(
      MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,2
      ) + (-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(
      MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)
      *msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/
      Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(
      MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta))/msu2(2,2))*TCD1F(2)(Sqrt(msq2(2,2)/msu2(2,2))))/(msq2(2,2)*Sqrt(
      msq2(2,2)*msu2(2,2))) - (0.5*(msq2(2,2)*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666*Log(
      MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*TCF(2)(Sqrt(msq2(2,2)/msu2(2,2)
      )))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2)))) - (TCF(2)(Sqrt(msq2(2,2
      )/msu2(2,2)))*(Sqr(10.666666666666666*Sqr(g3)*M3Input*Log(MSUSY/Qmatch) + (
      12*AtInput*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) -
      (0.6*(10*MuInput*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + 10*MuInput*Sqr(TanBeta)*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2))))/Cube(TanBeta)) + 2*(AtInput - MuInput/
      TanBeta)*((-0.1111111111111111*M3Input*Log(MSUSY/Qmatch)*(-6 + Log(msd2(0,0)
      /Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2)/Sqr(Qmatch)) + 2*
      Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(Qmatch)) + 2*Log(msq2(2,2)/
      Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(msu2(1,1)/Sqr(Qmatch)) + Log
      (msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/Sqr(Qmatch)))*Quad(g3))/Sqr(
      3.141592653589793) - 64*M3Input*Quad(g3)*Sqr(Log(MSUSY/Qmatch)) + (144*
      AtInput*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(
      TanBeta) - (64*AtInput*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr
      (Yu(2,2)))/Sqr(TanBeta) + (64*M3Input*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (3*MuInput*Log(MSUSY/Qmatch)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2))*(3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) + (3*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta)))/Cube(TanBeta) + (24*AtInput*Log(
      MSUSY/Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) -
      Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) +
      (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25
      *Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/
      Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/
      MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/
      MuInput))/Sqr(TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(
      Qmatch)) + TCF(6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input)
      - ((AtInput - MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2)
      )/M3Input))/M3Input)))/Sqr(TanBeta) - ((-0.5*MuInput*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*(-45*Quad(Yu(2,2)) - 45*Quad(Yu(2,2))*Sqr(TanBeta) + 32*
      Sqr(g3)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Quad(TanBeta) + (6*MuInput*Log(MSUSY/
      Qmatch)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(
      Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (
      0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*
      Sqr(AtInput - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt
      (msq2(2,2)*msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))
      /Sqr(TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr
      (TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF
      (6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput -
      MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/
      M3Input)))/Sqr(TanBeta))/TanBeta - (MuInput*((-3*Log(MSUSY/Qmatch)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2))*(-3*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (3*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta)))/Sqr(TanBeta) - ((0.5*(1 + Sqr(TanBeta))
      *Sqr(Log(MSUSY/Qmatch))*(-27*Quad(Yu(2,2)) - 27*Quad(Yu(2,2))*Sqr(TanBeta) +
      32*Sqr(g3)*Sqr(TanBeta)*Sqr(Yu(2,2))))/Cube(TanBeta) - (6*Log(MSUSY/Qmatch)*
      (1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*
      ((0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr
      (MuInput)/Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput
      - MuInput/TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*
      msu2(2,2)) + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(
      TanBeta) + (0.5*(1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(
      TanBeta)) - 1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(
      6)(Sqrt(msq2(2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput -
      MuInput/TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/
      M3Input)))/TanBeta)/TanBeta))/TanBeta)))/Sqrt(msq2(2,2)*msu2(2,2)) - Sqr(
      AtInput - MuInput/TanBeta)*((-0.25*Sqrt(msq2(2,2)/msu2(2,2))*((msq2(2,2)*((
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) - 4*Log
      (MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/
      msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) + (4*
      Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/
      (msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(
      msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(
      msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(
      msu2(2,2)*Sqr(TanBeta))))/msu2(2,2) + (-10.666666666666666*Log(MSUSY/Qmatch)
      *Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2,2))*(msq2(2,2)*(-10.666666666666666*Log
      (MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msu2(2,2)*(-10.666666666666666
      *Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta
      ))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1
      + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(
      MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))*TCD1F(2)(Sqrt(msq2(
      2,2)/msu2(2,2))))/(Sqrt(msq2(2,2)*msu2(2,2))*Sqr(msq2(2,2))) + (0.5*TCF(2)(
      Sqrt(msq2(2,2)/msu2(2,2)))*((0.75*Sqr(msu2(2,2)*(-10.666666666666666*Sqr(g3)
      *Sqr(M3Input)*Log(MSUSY/Qmatch) + (2*Sqr(mAInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Sqr(AtInput)*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Sqr(MuInput)*(1 + Sqr(TanBeta))*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)
      *msq2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta)) + msq2(2,2)*(-
      10.666666666666666*Sqr(g3)*Sqr(M3Input)*Log(MSUSY/Qmatch) + (4*Sqr(mAInput)*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*(1 + Sqr(
      TanBeta))*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*(1 +
      Sqr(TanBeta))*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta))))/(Sqr
      (msq2(2,2))*Sqr(msu2(2,2))) - ((-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3
      )*Sqr(M3Input) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta))*(-10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*
      Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta)
      + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/
      Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2)))/Sqr(TanBeta)) + msq2(2,2)*((0.1111111111111111*Log(MSUSY/
      Qmatch)*(-6 + Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(
      msd2(2,2)/Sqr(Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(
      Qmatch)) + 2*Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(
      msu2(1,1)/Sqr(Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/
      Sqr(Qmatch)))*Quad(g3)*Sqr(M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*
      Sqr(M3Input)*Sqr(Log(MSUSY/Qmatch)) + (12*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log
      (MSUSY/Qmatch)))/Sqr(TanBeta) + (36*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(
      TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (36*msq2(2,2)*Quad(Yu(2,2)
      )*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (36*msu2(2,2
      )*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta)
      + (84*Quad(Yu(2,2))*Sqr(AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)
      ))/Quad(TanBeta) - (48*Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch)))/Quad(TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(mAInput)
      *Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (42.666666666666664*
      AtInput*M3Input*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2
      )))/Sqr(TanBeta) - (21.333333333333332*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*
      Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*msu2
      (2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (21.333333333333332*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (42.666666666666664*Sqr(g3)*
      Sqr(M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (21.333333333333332*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(
      mAInput)*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-
      1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/
      Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (8*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta)) + msu2(2,2)*((0.1111111111111111*Log(MSUSY/Qmatch)*(-6 + Log(
      msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2)/Sqr(
      Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(Qmatch)) + 2*
      Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(msu2(1,1)/Sqr(
      Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/Sqr(Qmatch)))*
      Quad(g3)*Sqr(M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*Sqr(M3Input)*Sqr
      (Log(MSUSY/Qmatch)) + (6*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log(MSUSY/Qmatch)))/
      Sqr(TanBeta) + (18*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch)))/Quad(TanBeta) + (18*msq2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(
      TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (18*msu2(2,2)*Quad(Yu(2,2)
      )*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (42*Quad(Yu(
      2,2))*Sqr(AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(
      TanBeta) - (24*Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch)))/Quad(TanBeta) - (10.666666666666666*Sqr(g3)*Sqr(mAInput)*Sqr
      (Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (21.333333333333332*AtInput
      *M3Input*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr
      (TanBeta) - (10.666666666666666*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log
      (MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (10.666666666666666*msu2(2,2)*
      Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (10.666666666666666*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(
      M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) + (10.666666666666666*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(mAInput
      )*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*
      Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(
      Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta)))/(msq2(2,2)*msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((0.125*
      msu2(2,2)*Sqr((msq2(2,2)*((10.666666666666666*Sqr(g3)*Sqr(M3Input)*Log(MSUSY
      /Qmatch))/msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Sqr(AtInput)*Log
      (MSUSY/Qmatch)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Sqr(MuInput)*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(
      2,2) - (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Sqr(AtInput)*Log
      (MSUSY/Qmatch)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Sqr(mAInput)*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Sqr(MuInput)*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*
      msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,2) + (-
      10.666666666666666*Sqr(g3)*Sqr(M3Input)*Log(MSUSY/Qmatch) + (2*Sqr(mAInput)*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*Sqr(AtInput)*(1 + Sqr(
      TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Sqr(MuInput)*(1
      + Sqr(TanBeta))*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(
      TanBeta))*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 +
      Sqr(TanBeta))*Log(MSUSY/Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2
      ,2))*TCD2F(2)(Sqrt(msq2(2,2)/msu2(2,2))))/msq2(2,2) + 0.5*Sqrt(msq2(2,2)/
      msu2(2,2))*TCD1F(2)(Sqrt(msq2(2,2)/msu2(2,2)))*((-0.25*Sqr(msu2(2,2))*Sqr((
      msq2(2,2)*((10.666666666666666*Sqr(g3)*Sqr(M3Input)*Log(MSUSY/Qmatch))/msu2(
      2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Sqr(AtInput)*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/msu2(2,2) + (4*Sqr(MuInput)*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/
      msu2(2,2) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(
      MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Sqr(AtInput)*Log(MSUSY/Qmatch)
      *Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Sqr(mAInput)*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Sqr(MuInput)*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(
      Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta))))/msu2(2,2) + (-10.666666666666666*Sqr(g3
      )*Sqr(M3Input)*Log(MSUSY/Qmatch) + (2*Sqr(mAInput)*Log(MSUSY/Qmatch)*Sqr(Yu(
      2,2)))/Sqr(TanBeta) + (2*Sqr(AtInput)*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)*
      Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Sqr(MuInput)*(1 + Sqr(TanBeta))*Log(MSUSY/
      Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/Qmatch)
      *msq2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (2*(1 + Sqr(TanBeta))*Log(MSUSY/
      Qmatch)*msu2(2,2)*Sqr(Yu(2,2)))/Sqr(TanBeta))/msu2(2,2)))/Sqr(msq2(2,2)) + (
      msu2(2,2)*((((10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input))/
      msu2(2,2) - 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) - (4*Log(MSUSY/Qmatch)*msq2(2,2
      )*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/
      msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*
      Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*msq2(2,2
      )*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(AtInput)
      *Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*Log(MSUSY/Qmatch)*Sqr(mAInput)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)))*(-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input) + (2*Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)
      ))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2
      ,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + (2*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2)))/Sqr(TanBeta) - (2*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/msu2(2,2) + ((0.1111111111111111*Log(
      MSUSY/Qmatch)*(-6 + Log(msd2(0,0)/Sqr(Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch))
      + Log(msd2(2,2)/Sqr(Qmatch)) + 2*Log(msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1
      )/Sqr(Qmatch)) + 2*Log(msq2(2,2)/Sqr(Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) +
      Log(msu2(1,1)/Sqr(Qmatch)) + Log(msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input
      )/Sqr(Qmatch)))*Quad(g3)*Sqr(M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*
      Sqr(M3Input)*Sqr(Log(MSUSY/Qmatch)) + (6*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log(
      MSUSY/Qmatch)))/Sqr(TanBeta) + (18*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(
      TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (18*msq2(2,2)*Quad(Yu(2,2)
      )*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (18*msu2(2,2
      )*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta)
      + (42*Quad(Yu(2,2))*Sqr(AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)
      ))/Quad(TanBeta) - (24*Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch)))/Quad(TanBeta) - (10.666666666666666*Sqr(g3)*Sqr(mAInput)
      *Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (21.333333333333332*
      AtInput*M3Input*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2
      )))/Sqr(TanBeta) - (10.666666666666666*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*
      Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (10.666666666666666*msu2
      (2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) - (10.666666666666666*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*Sqr(g3)*
      Sqr(M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(
      TanBeta) + (10.666666666666666*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(
      Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(
      mAInput)*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-
      1 + 2*Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/
      Sqr(Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta))/msu2(2,2) + (msq2(2,2)*((((-10.666666666666666*Log(MSUSY/
      Qmatch)*Sqr(g3)*Sqr(M3Input))/msu2(2,2) + 4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)) +
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch
      )*Sqr(AtInput)*Sqr(Yu(2,2)))/msu2(2,2) - (4*Log(MSUSY/Qmatch)*Sqr(MuInput)*
      Sqr(Yu(2,2)))/msu2(2,2) + (4*Log(MSUSY/Qmatch)*Sqr(Yu(2,2)))/Sqr(TanBeta) +
      (4*Log(MSUSY/Qmatch)*msq2(2,2)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*
      Log(MSUSY/Qmatch)*Sqr(AtInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) + (4*
      Log(MSUSY/Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)) - (4*
      Log(MSUSY/Qmatch)*Sqr(MuInput)*Sqr(Yu(2,2)))/(msu2(2,2)*Sqr(TanBeta)))*(-
      10.666666666666666*Log(MSUSY/Qmatch)*Sqr(g3)*Sqr(M3Input) + (4*Log(MSUSY/
      Qmatch)*Sqr(mAInput)*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*msq2(
      2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/Qmatch)*
      msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (4*Log(MSUSY/
      Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (4*Log(
      MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2)))/Sqr(TanBeta)))/
      msu2(2,2) - ((0.1111111111111111*Log(MSUSY/Qmatch)*(-6 + Log(msd2(0,0)/Sqr(
      Qmatch)) + Log(msd2(1,1)/Sqr(Qmatch)) + Log(msd2(2,2)/Sqr(Qmatch)) + 2*Log(
      msq2(0,0)/Sqr(Qmatch)) + 2*Log(msq2(1,1)/Sqr(Qmatch)) + 2*Log(msq2(2,2)/Sqr(
      Qmatch)) + Log(msu2(0,0)/Sqr(Qmatch)) + Log(msu2(1,1)/Sqr(Qmatch)) + Log(
      msu2(2,2)/Sqr(Qmatch)) + 12*Log(Sqr(M3Input)/Sqr(Qmatch)))*Quad(g3)*Sqr(
      M3Input))/Sqr(3.141592653589793) + 96*Quad(g3)*Sqr(M3Input)*Sqr(Log(MSUSY/
      Qmatch)) + (12*Quad(Yu(2,2))*Sqr(mAInput)*Sqr(Log(MSUSY/Qmatch)))/Sqr(
      TanBeta) + (36*Quad(Yu(2,2))*Sqr(mAInput)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch)))/Quad(TanBeta) + (36*msq2(2,2)*Quad(Yu(2,2))*Sqr(1 + Sqr(TanBeta))*
      Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (36*msu2(2,2)*Quad(Yu(2,2))*Sqr(1 +
      Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) + (84*Quad(Yu(2,2))*Sqr(
      AtInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/Quad(TanBeta) - (48*
      Quad(Yu(2,2))*Sqr(MuInput)*Sqr(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch)))/
      Quad(TanBeta) - (21.333333333333332*Sqr(g3)*Sqr(mAInput)*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (42.666666666666664*AtInput*M3Input*
      Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta)
      - (21.333333333333332*msq2(2,2)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (21.333333333333332*msu2(2,2)*Sqr(g3)*
      (1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (
      21.333333333333332*Sqr(AtInput)*Sqr(g3)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/
      Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) - (42.666666666666664*Sqr(g3)*Sqr(
      M3Input)*(1 + Sqr(TanBeta))*Sqr(Log(MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta
      ) + (21.333333333333332*Sqr(g3)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Log(
      MSUSY/Qmatch))*Sqr(Yu(2,2)))/Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(mAInput
      )*Sqr(Yu(2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*
      Log(Sqr(mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(
      Qmatch))*(1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/
      TanBeta)*TCF(5)(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2))
      + ((1 + Sqr(TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(
      1 + Sqr(TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msq2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*msu2(2,2)*(1 + Sqr(TanBeta))*Sqr(Yu(2,2
      ))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) + (8*Log(MSUSY/Qmatch)*Sqr(AtInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta) - (8*Log(MSUSY/Qmatch)*Sqr(MuInput)*(1 + Sqr(TanBeta))*Sqr(Yu(
      2,2))*(-1.3333333333333333*Sqr(g3) - Sqr(Yu(2,2))*((0.375*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(Qmatch))))/Sqr(TanBeta) + (0.75*Log(Sqr(MuInput)/Sqr(Qmatch))*(
      1 + Sqr(TanBeta)))/Sqr(TanBeta) - (0.25*Sqr(AtInput - MuInput/TanBeta)*TCF(5
      )(Sqrt(msq2(2,2))/Sqrt(msu2(2,2))))/Sqrt(msq2(2,2)*msu2(2,2)) + ((1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msq2(2,2))/MuInput))/Sqr(TanBeta) + (0.5*(1 + Sqr(
      TanBeta))*TCF(6)(Sqrt(msu2(2,2))/MuInput))/Sqr(TanBeta)) -
      1.3333333333333333*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(Qmatch)) + TCF(6)(Sqrt(msq2
      (2,2))/M3Input) + TCF(6)(Sqrt(msu2(2,2))/M3Input) - ((AtInput - MuInput/
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msu2(2,2))/M3Input))/M3Input)))
      /Sqr(TanBeta))/msu2(2,2)))/msu2(2,2)))/msq2(2,2)))/Sqrt(msq2(2,2)*msu2(2,2))
      ))))/Sqrt(msq2(2,2)*msu2(2,2)))))/Power6(3.141592653589793), True, 0), 0) +
      0.25*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta))) + IF(
      LambdaLoopOrder >= 3, FSHimalayaMh23L, 0)*UnitStep(-3 + LambdaLoopOrder) + (
      IF(DeltaOS >= 1, IF(TwoLoopAtAs >= 1, Quad(Yu(2,2))*Sqr(g3)*WHICH(Abs(-1 +
      Sqrt(Abs(msq2(2,2)))/Sqrt(Abs(msu2(2,2)))) < 0.01 && Abs(-1 + M3Input/Sqrt(
      Abs(msu2(2,2)))) < 0.01 && Abs(-1 + M3Input/Sqrt(Abs(msq2(2,2)))) < 0.01, (
      0.010416666666666666*Re((-6 + (18*(AtInput - MuInput/TanBeta))/Sqrt(Abs(msq2
      (2,2))) + Cube(AtInput - MuInput/TanBeta)/Power3(Sqrt(Abs(msq2(2,2)))) - (9*
      Sqr(AtInput - MuInput/TanBeta))/Abs(msq2(2,2)))*Sqr(2 + (AtInput - MuInput/
      TanBeta)/Sqrt(Abs(msq2(2,2)))) + Log(Sqr(SCALE)/Abs(msq2(2,2)))*(-12 + (24*(
      AtInput - MuInput/TanBeta))/Sqrt(Abs(msq2(2,2))) - (4*Cube(AtInput - MuInput
      /TanBeta))/Power3(Sqrt(Abs(msq2(2,2)))) - (6*Sqr(AtInput - MuInput/TanBeta))
      /Abs(msq2(2,2)) + Quad(AtInput - MuInput/TanBeta)/Sqr(Abs(msq2(2,2))))))/
      Quad(3.141592653589793), Abs(-1 + Sqrt(Abs(msq2(2,2)))/Sqrt(Abs(msu2(2,2))))
      < 0.085, (0.015625*(1 - (5.88235294117647*(Sqrt(Abs(msq2(2,2))) - 0.915*Sqrt
      (Abs(msu2(2,2)))))/Sqrt(Abs(msu2(2,2))))*Re((-2.8532880624229446*Quad(
      M3Input)*((0.8372250000000001*Abs(msu2(2,2))*(3 + 2*Log((1.194422049031025*
      Sqr(SCALE))/Abs(msu2(2,2))) - Log((0.8372250000000001*Abs(msu2(2,2)))/Sqr(
      M3Input))*(-2 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)) + (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(
      1 - (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(-1 + (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) - (2*
      Quad(M3Input)*((Abs(msu2(2,2))*(3 + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))) -
      Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)) + Abs(
      msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) + (
      1.4266440312114723*(AtInput - MuInput/TanBeta)*Power7(M3Input)*((-
      162.52540772058205*Abs(msu2(2,2))*((-0.3255499999999998*Abs(msu2(2,2))*Sqr(
      AtInput - MuInput/TanBeta))/Quad(M3Input) - 0.17766242741323135*((Abs(msu2(2
      ,2))*((2*Abs(msu2(2,2)))/Sqr(M3Input) - Sqr(AtInput - MuInput/TanBeta)/Sqr(
      M3Input)))/Sqr(M3Input) - (0.8372250000000001*Abs(msu2(2,2))*((4*Abs(msu2(2,
      2)))/Sqr(M3Input) + Sqr(AtInput - MuInput/TanBeta)/Sqr(M3Input)))/Sqr(
      M3Input) + (1.4018914012500003*Sqr(Abs(msu2(2,2))))/Quad(M3Input)))*((
      24.5737981876824*(-((Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-1 + (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(M3Input)) + (
      0.8372250000000001*Abs(msu2(2,2))*Log((0.8372250000000001*Abs(msu2(2,2)))/
      Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input))*Sqr(AtInput
      - MuInput/TanBeta))/(Abs(msu2(2,2))*(-1 + (0.8372250000000001*Abs(msu2(2,2))
      )/Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input))) - (4.7776881961241*Quad(
      M3Input)*((Abs(msu2(2,2))*ComplexLog(1 - (0.8372250000000001*Abs(msu2(2,2)))
      /Sqr(M3Input))*(-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(
      M3Input) + (0.8372250000000001*Abs(msu2(2,2))*(ComplexLog(1 - Abs(msu2(2,2))
      /Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)) - (2*Abs(msu2(2,2))*(2 +
      Log((1.*Sqr(SCALE))/Sqr(M3Input))))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(Abs(
      msu2(2,2))) + (1.4266440312114723*(AtInput - MuInput/TanBeta)*Cube(M3Input)*
      ((0.8372250000000001*Abs(msu2(2,2))*(3 + 2*Log((1.194422049031025*Sqr(SCALE)
      )/Abs(msu2(2,2))) - Log((0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*(-
      2 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)) + (0.8372250000000001*
      Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(-1 + (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) - (
      6.1434495469206*M3Input*(AtInput - MuInput/TanBeta)*(-8*ComplexLog(1 - Abs(
      msu2(2,2))/Sqr(M3Input)) + (1.1394249999999992*Abs(msu2(2,2)))/Sqr(M3Input)
      + (4*Abs(msu2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input
      ) + (6.535025000000001*Abs(msu2(2,2))*Log((0.8372250000000001*Abs(msu2(2,2))
      )/Sqr(SCALE)))/Sqr(M3Input) - (8.162775*Abs(msu2(2,2))*Log((1.*Abs(msu2(2,2)
      ))/Sqr(SCALE)))/Sqr(M3Input) - (0.6510999999999996*Abs(msu2(2,2))*Log((1.*
      Sqr(SCALE))/Sqr(M3Input)))/Sqr(M3Input) + Abs(msu2(2,2))/((1 - Abs(msu2(2,2)
      )/Sqr(M3Input))*Sqr(M3Input)) + Abs(msu2(2,2))/((1 - (0.8372250000000001*Abs
      (msu2(2,2)))/Sqr(M3Input))*Sqr(M3Input)) + (0.8372250000000001*Abs(msu2(2,2)
      ))/((-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(M3Input)) +
      (0.8372250000000001*Abs(msu2(2,2)))/((-1 + Abs(msu2(2,2))/Sqr(M3Input))*Sqr(
      M3Input)) + (4*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))/Abs
      (msu2(2,2)) + (0.1627749999999999*Abs(msu2(2,2))*Log((0.8372250000000001*Abs
      (msu2(2,2)))/Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + (0.8372250000000001*Abs(
      msu2(2,2)))/Sqr(M3Input))) + (0.1627749999999999*Abs(msu2(2,2))*Log(Abs(msu2
      (2,2))/Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))) -
      (4.7776881961241*ComplexLog(1 - (0.8372250000000001*Abs(msu2(2,2)))/Sqr(
      M3Input))*Sqr(M3Input)*Sqr(-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(
      M3Input)))/Abs(msu2(2,2))))/Abs(msu2(2,2)) + ((AtInput - MuInput/TanBeta)*
      Cube(M3Input)*((Abs(msu2(2,2))*(3 + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))) -
      Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)) + Abs(
      msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2)))))/Sqr(
      M3Input) + ((AtInput - MuInput/TanBeta)*((-0.06460317591804718*Sqr(M3Input)*
      ((ComplexLog(1 - (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(
      msu2(2,2)))*Sqr(-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Quad
      (M3Input) - (0.8372250000000001*Abs(msu2(2,2))*((Abs(msu2(2,2))*((-
      0.4883249999999997*Abs(msu2(2,2)))/Sqr(M3Input) + (1.6744500000000002*Abs(
      msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) - (2*Abs(msu2(2
      ,2))*Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) - (
      0.8372250000000001*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs
      (msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (Log((
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*(-2 + (0.8372250000000001*
      Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2))))/Quad(M3Input) + (
      0.8372250000000001*Abs(msu2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)
      )*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input)))/Abs(
      msu2(2,2)) + (0.7009457006250002*(AtInput - MuInput/TanBeta)*Sqr(Abs(msu2(2,
      2)))*((-9.283011849255823*(AtInput - MuInput/TanBeta)*Power7(M3Input)*((
      ComplexLog(1 - (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(
      msu2(2,2)))*Sqr(-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Quad
      (M3Input) - (0.8372250000000001*Abs(msu2(2,2))*((Abs(msu2(2,2))*((-
      0.4883249999999997*Abs(msu2(2,2)))/Sqr(M3Input) + (1.6744500000000002*Abs(
      msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) - (2*Abs(msu2(2
      ,2))*Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) - (
      0.8372250000000001*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs
      (msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (Log((
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*(-2 + (0.8372250000000001*
      Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2))))/Quad(M3Input) + (
      0.8372250000000001*Abs(msu2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)
      )*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input)))/Quad(
      Abs(msu2(2,2))) + 1.2191256391440917*(-0.1627749999999999*((24.5737981876824
      *(-((Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-1 + (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(M3Input)) + (
      0.8372250000000001*Abs(msu2(2,2))*Log((0.8372250000000001*Abs(msu2(2,2)))/
      Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input))*Sqr(AtInput
      - MuInput/TanBeta))/(Abs(msu2(2,2))*(-1 + (0.8372250000000001*Abs(msu2(2,2))
      )/Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input))) - (4.7776881961241*Quad(
      M3Input)*((Abs(msu2(2,2))*ComplexLog(1 - (0.8372250000000001*Abs(msu2(2,2)))
      /Sqr(M3Input))*(-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(
      M3Input) + (0.8372250000000001*Abs(msu2(2,2))*(ComplexLog(1 - Abs(msu2(2,2))
      /Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)) - (2*Abs(msu2(2,2))*(2 +
      Log((1.*Sqr(SCALE))/Sqr(M3Input))))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(Abs(
      msu2(2,2))) - (1.4266440312114723*(AtInput - MuInput/TanBeta)*Cube(M3Input)*
      ((0.8372250000000001*Abs(msu2(2,2))*(3 + 2*Log((1.194422049031025*Sqr(SCALE)
      )/Abs(msu2(2,2))) - Log((0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*(-
      2 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)) + (0.8372250000000001*
      Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(-1 + (
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) - (
      6.1434495469206*M3Input*(AtInput - MuInput/TanBeta)*(-8*ComplexLog(1 - Abs(
      msu2(2,2))/Sqr(M3Input)) + (1.1394249999999992*Abs(msu2(2,2)))/Sqr(M3Input)
      + (4*Abs(msu2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input
      ) + (6.535025000000001*Abs(msu2(2,2))*Log((0.8372250000000001*Abs(msu2(2,2))
      )/Sqr(SCALE)))/Sqr(M3Input) - (8.162775*Abs(msu2(2,2))*Log((1.*Abs(msu2(2,2)
      ))/Sqr(SCALE)))/Sqr(M3Input) - (0.6510999999999996*Abs(msu2(2,2))*Log((1.*
      Sqr(SCALE))/Sqr(M3Input)))/Sqr(M3Input) + Abs(msu2(2,2))/((1 - Abs(msu2(2,2)
      )/Sqr(M3Input))*Sqr(M3Input)) + Abs(msu2(2,2))/((1 - (0.8372250000000001*Abs
      (msu2(2,2)))/Sqr(M3Input))*Sqr(M3Input)) + (0.8372250000000001*Abs(msu2(2,2)
      ))/((-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(M3Input)) +
      (0.8372250000000001*Abs(msu2(2,2)))/((-1 + Abs(msu2(2,2))/Sqr(M3Input))*Sqr(
      M3Input)) + (4*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))/Abs
      (msu2(2,2)) + (0.1627749999999999*Abs(msu2(2,2))*Log((0.8372250000000001*Abs
      (msu2(2,2)))/Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + (0.8372250000000001*Abs(
      msu2(2,2)))/Sqr(M3Input))) + (0.1627749999999999*Abs(msu2(2,2))*Log(Abs(msu2
      (2,2))/Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))) -
      (4.7776881961241*ComplexLog(1 - (0.8372250000000001*Abs(msu2(2,2)))/Sqr(
      M3Input))*Sqr(M3Input)*Sqr(-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(
      M3Input)))/Abs(msu2(2,2))))/Abs(msu2(2,2)) + (3*(AtInput - MuInput/TanBeta)*
      Cube(M3Input)*((Abs(msu2(2,2))*(3 + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))) -
      Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)) + Abs(
      msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2)))) + (
      7.166532294186149*(AtInput - MuInput/TanBeta)*Power7(M3Input)*((ComplexLog(1
       - (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2)))*Sqr
      (-1 + (0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input)))/Quad(M3Input) - (
      0.8372250000000001*Abs(msu2(2,2))*((Abs(msu2(2,2))*((-0.4883249999999997*Abs
      (msu2(2,2)))/Sqr(M3Input) + (1.6744500000000002*Abs(msu2(2,2))*Log((1.*Sqr(
      SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) - (2*Abs(msu2(2,2))*Log((
      1.194422049031025*Sqr(SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) - (
      0.8372250000000001*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs
      (msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (Log((
      0.8372250000000001*Abs(msu2(2,2)))/Sqr(M3Input))*(-2 + (0.8372250000000001*
      Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2))))/Quad(M3Input) + (
      0.8372250000000001*Abs(msu2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)
      )*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input)))/Quad(
      Abs(msu2(2,2))))))/Power5(M3Input)))/M3Input))/Quad(Abs(msu2(2,2)))))/Quad(
      3.141592653589793) + (0.0009435645454673102*(Sqrt(Abs(msq2(2,2))) - 0.915*
      Sqrt(Abs(msu2(2,2))))*Re((-2*Quad(M3Input)*((Abs(msu2(2,2))*(3 + 2*Log((1.*
      Sqr(SCALE))/Abs(msu2(2,2))) - Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(
      msu2(2,2))/Sqr(M3Input)) + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) +
      ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(
      M3Input))))/Sqr(Abs(msu2(2,2))) - (1.4431485685359067*Quad(M3Input)*((
      1.177225*Abs(msu2(2,2))*(3 + 2*Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(
      2,2))) - Log((1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*(-2 + (1.177225*Abs(
      msu2(2,2)))/Sqr(M3Input)) + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(
      M3Input) + ComplexLog(1 - (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(-1 + (
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) + (
      0.7215742842679533*(AtInput - MuInput/TanBeta)*Power7(M3Input)*((
      248.96826777593753*Abs(msu2(2,2))*((0.35444999999999993*Abs(msu2(2,2))*Sqr(
      AtInput - MuInput/TanBeta))/Quad(M3Input) + 0.16315997398484572*((Abs(msu2(2
      ,2))*((2*Abs(msu2(2,2)))/Sqr(M3Input) - Sqr(AtInput - MuInput/TanBeta)/Sqr(
      M3Input)))/Sqr(M3Input) - (1.177225*Abs(msu2(2,2))*((4*Abs(msu2(2,2)))/Sqr(
      M3Input) + Sqr(AtInput - MuInput/TanBeta)/Sqr(M3Input)))/Sqr(M3Input) + (
      2.7717174012499997*Sqr(Abs(msu2(2,2))))/Quad(M3Input)))*((-
      22.570179150797014*((1.177225*Abs(msu2(2,2))*Log((1.177225*Abs(msu2(2,2)))/
      Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) - (Abs(msu2(2
      ,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(
      M3Input)))/Sqr(M3Input))*Sqr(AtInput - MuInput/TanBeta))/(Abs(msu2(2,2))*(-1
       + Abs(msu2(2,2))/Sqr(M3Input))*(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)
      )) - (3.397821147189365*Quad(M3Input)*((Abs(msu2(2,2))*ComplexLog(1 - (
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(
      M3Input)))/Sqr(M3Input) + (1.177225*Abs(msu2(2,2))*(ComplexLog(1 - Abs(msu2(
      2,2))/Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)) - (2*Abs(msu2(2,2))*(
      2 + Log((1.*Sqr(SCALE))/Sqr(M3Input))))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(
      Abs(msu2(2,2))) + ((AtInput - MuInput/TanBeta)*Cube(M3Input)*((Abs(msu2(2,2)
      )*(3 + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))) - Log(Abs(msu2(2,2))/Sqr(
      M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)) + Abs(msu2(2,2))/Sqr(M3Input)))
      /Sqr(M3Input) + ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(
      msu2(2,2))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) + (0.7215742842679533*(
      AtInput - MuInput/TanBeta)*Cube(M3Input)*((1.177225*Abs(msu2(2,2))*(3 + 2*
      Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(2,2))) - Log((1.177225*Abs(msu2
      (2,2)))/Sqr(M3Input))*(-2 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)) + (
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - (
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/
      Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) + (5.642544787699253*M3Input*(AtInput -
      MuInput/TanBeta)*(-8*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)) - (
      1.2405749999999998*Abs(msu2(2,2)))/Sqr(M3Input) + (4*Abs(msu2(2,2))*
      ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) - (7.822775*Abs(
      msu2(2,2))*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)))/Sqr(M3Input) + (9.595025*Abs
      (msu2(2,2))*Log((1.177225*Abs(msu2(2,2)))/Sqr(SCALE)))/Sqr(M3Input) + (
      0.7088999999999999*Abs(msu2(2,2))*Log((1.*Sqr(SCALE))/Sqr(M3Input)))/Sqr(
      M3Input) + Abs(msu2(2,2))/((1 - (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(
      M3Input)) + Abs(msu2(2,2))/((1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))
      + (1.177225*Abs(msu2(2,2)))/((-1 + Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input)
      ) + (1.177225*Abs(msu2(2,2)))/((-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))
      *Sqr(M3Input)) + (4*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input)
      )/Abs(msu2(2,2)) - (0.17722499999999997*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/
      Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))) - (
      0.17722499999999997*Abs(msu2(2,2))*Log((1.177225*Abs(msu2(2,2)))/Sqr(M3Input
      )))/(Sqr(M3Input)*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))) - (
      3.397821147189365*ComplexLog(1 - (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr
      (M3Input)*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)))/Abs(msu2(2,2))))
      /Abs(msu2(2,2))))/Sqr(M3Input) + ((AtInput - MuInput/TanBeta)*((
      0.0500481932404575*Sqr(M3Input)*((-1.177225*Abs(msu2(2,2))*((Abs(msu2(2,2))*
      ((0.5316749999999999*Abs(msu2(2,2)))/Sqr(M3Input) - (2*Abs(msu2(2,2))*Log((
      0.8494552867973413*Sqr(SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) + (2.35445*Abs(
      msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))))/Sqr(M3Input) - (1.177225*Abs
      (msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(msu2(2,2))/Sqr(
      M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (Log((1.177225*Abs(msu2(2,2)))/Sqr(
      M3Input))*(-2 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2))))
      /Quad(M3Input) + (1.177225*Abs(msu2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input)
      + (ComplexLog(1 - (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2))
      )*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)))/Quad(M3Input)))/Abs(msu2
      (2,2)) + (1.3858587006249998*(AtInput - MuInput/TanBeta)*Sqr(Abs(msu2(2,2)))
      *((3.8502754206973924*(AtInput - MuInput/TanBeta)*Power7(M3Input)*((-
      1.177225*Abs(msu2(2,2))*((Abs(msu2(2,2))*((0.5316749999999999*Abs(msu2(2,2))
      )/Sqr(M3Input) - (2*Abs(msu2(2,2))*Log((0.8494552867973413*Sqr(SCALE))/Abs(
      msu2(2,2))))/Sqr(M3Input) + (2.35445*Abs(msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(
      msu2(2,2))))/Sqr(M3Input) - (1.177225*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(
      M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) +
      (Log((1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*(-2 + (1.177225*Abs(msu2(2,2)))
      /Sqr(M3Input))*Sqr(Abs(msu2(2,2))))/Quad(M3Input) + (1.177225*Abs(msu2(2,2))
      *ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(
      M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (ComplexLog(1 - (1.177225*Abs(msu2(
      2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2)))*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/
      Sqr(M3Input)))/Quad(M3Input)))/Quad(Abs(msu2(2,2))) - 0.7967249877351067*((
      5.096731720784048*(AtInput - MuInput/TanBeta)*Power7(M3Input)*((-1.177225*
      Abs(msu2(2,2))*((Abs(msu2(2,2))*((0.5316749999999999*Abs(msu2(2,2)))/Sqr(
      M3Input) - (2*Abs(msu2(2,2))*Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(2,
      2))))/Sqr(M3Input) + (2.35445*Abs(msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(msu2(2,
      2))))/Sqr(M3Input) - (1.177225*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input
      ))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (Log((
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*(-2 + (1.177225*Abs(msu2(2,2)))/Sqr(
      M3Input))*Sqr(Abs(msu2(2,2))))/Quad(M3Input) + (1.177225*Abs(msu2(2,2))*
      ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(
      M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (ComplexLog(1 - (1.177225*Abs(msu2(
      2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2)))*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/
      Sqr(M3Input)))/Quad(M3Input)))/Quad(Abs(msu2(2,2))) + 0.17722499999999997*((
      -22.570179150797014*((1.177225*Abs(msu2(2,2))*Log((1.177225*Abs(msu2(2,2)))/
      Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) - (Abs(msu2(2
      ,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(
      M3Input)))/Sqr(M3Input))*Sqr(AtInput - MuInput/TanBeta))/(Abs(msu2(2,2))*(-1
       + Abs(msu2(2,2))/Sqr(M3Input))*(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)
      )) - (3.397821147189365*Quad(M3Input)*((Abs(msu2(2,2))*ComplexLog(1 - (
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(
      M3Input)))/Sqr(M3Input) + (1.177225*Abs(msu2(2,2))*(ComplexLog(1 - Abs(msu2(
      2,2))/Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)) - (2*Abs(msu2(2,2))*(
      2 + Log((1.*Sqr(SCALE))/Sqr(M3Input))))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(
      Abs(msu2(2,2))) + (3*(AtInput - MuInput/TanBeta)*Cube(M3Input)*((Abs(msu2(2,
      2))*(3 + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))) - Log(Abs(msu2(2,2))/Sqr(
      M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)) + Abs(msu2(2,2))/Sqr(M3Input)))
      /Sqr(M3Input) + ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(
      msu2(2,2))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) - (0.7215742842679533*(
      AtInput - MuInput/TanBeta)*Cube(M3Input)*((1.177225*Abs(msu2(2,2))*(3 + 2*
      Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(2,2))) - Log((1.177225*Abs(msu2
      (2,2)))/Sqr(M3Input))*(-2 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)) + (
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - (
      1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/
      Sqr(M3Input))))/Sqr(Abs(msu2(2,2))) + (5.642544787699253*M3Input*(AtInput -
      MuInput/TanBeta)*(-8*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)) - (
      1.2405749999999998*Abs(msu2(2,2)))/Sqr(M3Input) + (4*Abs(msu2(2,2))*
      ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) - (7.822775*Abs(
      msu2(2,2))*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)))/Sqr(M3Input) + (9.595025*Abs
      (msu2(2,2))*Log((1.177225*Abs(msu2(2,2)))/Sqr(SCALE)))/Sqr(M3Input) + (
      0.7088999999999999*Abs(msu2(2,2))*Log((1.*Sqr(SCALE))/Sqr(M3Input)))/Sqr(
      M3Input) + Abs(msu2(2,2))/((1 - (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(
      M3Input)) + Abs(msu2(2,2))/((1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))
      + (1.177225*Abs(msu2(2,2)))/((-1 + Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input)
      ) + (1.177225*Abs(msu2(2,2)))/((-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))
      *Sqr(M3Input)) + (4*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input)
      )/Abs(msu2(2,2)) - (0.17722499999999997*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/
      Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))) - (
      0.17722499999999997*Abs(msu2(2,2))*Log((1.177225*Abs(msu2(2,2)))/Sqr(M3Input
      )))/(Sqr(M3Input)*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))) - (
      3.397821147189365*ComplexLog(1 - (1.177225*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr
      (M3Input)*Sqr(-1 + (1.177225*Abs(msu2(2,2)))/Sqr(M3Input)))/Abs(msu2(2,2))))
      /Abs(msu2(2,2))))))/Power5(M3Input)))/M3Input))/Quad(Abs(msu2(2,2)))))/Sqrt(
      Abs(msu2(2,2))), Abs(-1 + Sqrt(Abs(msq2(2,2)))/M3Input) < 0.085, (
      0.0009435645454673102*(M3Input - 0.915*Sqrt(Abs(msq2(2,2))))*Re(-
      2.7717174012499997*(0.022663710673270766*ComplexLog(0.15054471320265872) +
      0.8494552867973413*(3.6617324413227936 + 2*Log((1.*Sqr(SCALE))/Abs(msq2(2,2)
      )))) - (2.7717174012499997*Sqr(Abs(msq2(2,2)))*((0.8494552867973413*Abs(msu2
      (2,2))*(3 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)) - (-2 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((0.8494552867973413*
      Abs(msu2(2,2)))/Abs(msq2(2,2))) + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2)))))/
      Abs(msq2(2,2)) + ComplexLog(1 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2
      (2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))))/Sqr(
      Abs(msu2(2,2))) + (1.7701422470949426*(AtInput - MuInput/TanBeta)*Power3(
      Sqrt(Abs(msq2(2,2))))*((0.5206694477168091*((1.6989105735946823*(
      0.8494552867973413 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr
      (AtInput - MuInput/TanBeta))/Abs(msq2(2,2)) + Log((1.*Abs(msq2(2,2)))/Abs(
      msu2(2,2)))*(1.4431485685359067 + (0.8494552867973413*Abs(msu2(2,2))*((
      1.6989105735946826*Abs(msu2(2,2)))/Abs(msq2(2,2)) - (0.8494552867973412*Sqr(
      AtInput - MuInput/TanBeta))/Abs(msq2(2,2))))/Abs(msq2(2,2)) -
      0.8494552867973413*((3.397821147189365*Abs(msu2(2,2)))/Abs(msq2(2,2)) + (
      0.8494552867973412*Sqr(AtInput - MuInput/TanBeta))/Abs(msq2(2,2)))))*Sqr(Abs
      (msu2(2,2)))*((-5.543434802499999*Abs(msq2(2,2))*((-0.12788100252938794*Abs(
      msu2(2,2))*ComplexLog(0.15054471320265872))/Abs(msq2(2,2)) +
      0.8494552867973413*((-1 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))
      )*ComplexLog(1 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) - (
      1.6989105735946826*Abs(msu2(2,2))*(2 + Log((0.8494552867973413*Sqr(SCALE))/
      Abs(msq2(2,2)))))/Abs(msq2(2,2)))))/Abs(msu2(2,2)) + (1.2772891249999998*(
      AtInput - MuInput/TanBeta)*(0.022663710673270766*ComplexLog(
      0.15054471320265872) + 0.8494552867973413*(3.6617324413227936 + 2*Log((1.*
      Sqr(SCALE))/Abs(msq2(2,2))))))/Sqrt(Abs(msq2(2,2))) + (22.57017915079702*(-
      0.13859710249514381*(-1 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))
      ) + (0.12788100252938794*Abs(msu2(2,2))*Log((0.8494552867973413*Abs(msu2(2,2
      )))/Abs(msq2(2,2))))/Abs(msq2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/(Abs(
      msq2(2,2))*(0.8494552867973413 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(
      msq2(2,2)))*(-1 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))) + (
      0.9216589861751152*(AtInput - MuInput/TanBeta)*(-11.588731795280644 + (
      11.588731795280644*Abs(msu2(2,2)))/Abs(msq2(2,2)) + (0.8494552867973413*Abs(
      msu2(2,2)))/(Abs(msq2(2,2))*(1 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(
      msq2(2,2)))) + 0.8494552867973413/(-1 + (0.8494552867973413*Abs(msu2(2,2)))/
      Abs(msq2(2,2))) - 7.199173001148221*(-0.8494552867973413 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) - 0.1067211471893647*
      ComplexLog(0.15054471320265872) - 8*ComplexLog(1 - (0.8494552867973413*Abs(
      msu2(2,2)))/Abs(msq2(2,2))) + (4.7089*Abs(msq2(2,2))*ComplexLog(1 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))/Abs(msu2(2,2)) + (
      3.397821147189365*Abs(msu2(2,2))*ComplexLog(1 - (0.8494552867973413*Abs(msu2
      (2,2)))/Abs(msq2(2,2))))/Abs(msq2(2,2)) + 7.6450975811760715*Log((1.*Abs(
      msq2(2,2)))/Sqr(SCALE)) - (0.8494552867973413*Abs(msu2(2,2))*Log((1.*Abs(
      msq2(2,2)))/Sqr(SCALE)))/Abs(msq2(2,2)) + 0.8494552867973413*Log((1.*Abs(
      msu2(2,2)))/Sqr(SCALE)) - (7.6450975811760715*Abs(msu2(2,2))*Log((1.*Abs(
      msu2(2,2)))/Sqr(SCALE)))/Abs(msq2(2,2)) + 3.397821147189365*Log((
      0.8494552867973413*Sqr(SCALE))/Abs(msq2(2,2))) - (3.397821147189365*Abs(msu2
      (2,2))*Log((0.8494552867973413*Sqr(SCALE))/Abs(msq2(2,2))))/Abs(msq2(2,2)) -
      (0.8494552867973413*Log((0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))
      /Sqr(-1 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) + (
      0.8494552867973413*Abs(msu2(2,2))*Log((0.8494552867973413*Abs(msu2(2,2)))/
      Abs(msq2(2,2))))/(Abs(msq2(2,2))*Sqr(-1 + (0.8494552867973413*Abs(msu2(2,2))
      )/Abs(msq2(2,2))))))/(Sqrt(Abs(msq2(2,2)))*(0.8494552867973413 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))) + (1.2772891249999998*(
      AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs(msq2(2,2))))*((0.8494552867973413
      *Abs(msu2(2,2))*(3 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)) - (-
      2 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) + 2*Log((1.*Sqr(SCALE))/
      Abs(msu2(2,2)))))/Abs(msq2(2,2)) + ComplexLog(1 - (0.8494552867973413*Abs(
      msu2(2,2)))/Abs(msq2(2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msu2(2,2)))/Abs
      (msq2(2,2)))))/Sqr(Abs(msu2(2,2)))))/(Cube(0.8494552867973413 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr(Abs(msq2(2,2)))) + (
      0.9216589861751152*(AtInput - MuInput/TanBeta)*((2*(-1.6989105735946826 + (
      1.6989105735946826*Abs(msu2(2,2)))/Abs(msq2(2,2)) + (0.8494552867973413 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((1.*Abs(msq2(2,2)))/
      Abs(msu2(2,2))))*((0.01635355080792133*ComplexLog(0.15054471320265872)*Sqr(
      Abs(msu2(2,2))))/Sqr(Abs(msq2(2,2))) - 0.8494552867973413*((
      0.8494552867973413*Abs(msu2(2,2))*(2.548365860392024 - (2.548365860392024*
      Abs(msu2(2,2)))/Abs(msq2(2,2)) - 0.8494552867973413*(-2 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((0.8494552867973413*
      Abs(msu2(2,2)))/Abs(msq2(2,2))) - (1.6989105735946826*Abs(msu2(2,2))*Log((1.
      *Sqr(SCALE))/Abs(msq2(2,2))))/Abs(msq2(2,2)) + 1.6989105735946826*Log((1.*
      Sqr(SCALE))/Abs(msu2(2,2)))))/Abs(msq2(2,2)) + (0.13545597786404023*Sqr(Abs(
      msu2(2,2))))/Sqr(Abs(msq2(2,2))) + 0.8494552867973413*ComplexLog(1 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr(-1 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))))/Sqr(0.8494552867973413
       - (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) + (0.6650454232884363
      *(AtInput - MuInput/TanBeta)*((2.5545782499999996*(AtInput - MuInput/TanBeta
      )*(-0.8494552867973413 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))
      + 0.8494552867973413*Log((1.*Abs(msq2(2,2)))/Abs(msu2(2,2))))*((
      0.01635355080792133*ComplexLog(0.15054471320265872)*Sqr(Abs(msu2(2,2))))/Sqr
      (Abs(msq2(2,2))) - 0.8494552867973413*((0.8494552867973413*Abs(msu2(2,2))*(
      2.548365860392024 - (2.548365860392024*Abs(msu2(2,2)))/Abs(msq2(2,2)) -
      0.8494552867973413*(-2 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))
      *Log((0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) - (
      1.6989105735946826*Abs(msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(msq2(2,2))))/Abs(
      msq2(2,2)) + 1.6989105735946826*Log((1.*Sqr(SCALE))/Abs(msu2(2,2)))))/Abs(
      msq2(2,2)) + (0.13545597786404023*Sqr(Abs(msu2(2,2))))/Sqr(Abs(msq2(2,2))) +
      0.8494552867973413*ComplexLog(1 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(
      msq2(2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))))/
      (Sqrt(Abs(msq2(2,2)))*Cube(0.8494552867973413 - (0.8494552867973413*Abs(msu2
      (2,2)))/Abs(msq2(2,2)))) - (1.177225*Abs(msq2(2,2))*(-1.6989105735946826 + (
      1.6989105735946826*Abs(msu2(2,2)))/Abs(msq2(2,2)) + (0.8494552867973413 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((1.*Abs(msq2(2,2)))/
      Abs(msu2(2,2))))*((10.620853482569656*(AtInput - MuInput/TanBeta)*Power(Abs(
      msq2(2,2)),2.5)*((0.01635355080792133*ComplexLog(0.15054471320265872)*Sqr(
      Abs(msu2(2,2))))/Sqr(Abs(msq2(2,2))) - 0.8494552867973413*((
      0.8494552867973413*Abs(msu2(2,2))*(2.548365860392024 - (2.548365860392024*
      Abs(msu2(2,2)))/Abs(msq2(2,2)) - 0.8494552867973413*(-2 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((0.8494552867973413*
      Abs(msu2(2,2)))/Abs(msq2(2,2))) - (1.6989105735946826*Abs(msu2(2,2))*Log((1.
      *Sqr(SCALE))/Abs(msq2(2,2))))/Abs(msq2(2,2)) + 1.6989105735946826*Log((1.*
      Sqr(SCALE))/Abs(msu2(2,2)))))/Abs(msq2(2,2)) + (0.13545597786404023*Sqr(Abs(
      msu2(2,2))))/Sqr(Abs(msq2(2,2))) + 0.8494552867973413*ComplexLog(1 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr(-1 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))))/Cube(Abs(msu2(2,2))) +
      (-1 + (1.*Abs(msq2(2,2)))/Abs(msu2(2,2)))*((-5.543434802499999*Abs(msq2(2,2)
      )*((-0.12788100252938794*Abs(msu2(2,2))*ComplexLog(0.15054471320265872))/Abs
      (msq2(2,2)) + 0.8494552867973413*((-1 + (0.8494552867973413*Abs(msu2(2,2)))/
      Abs(msq2(2,2)))*ComplexLog(1 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(
      2,2))) - (1.6989105735946826*Abs(msu2(2,2))*(2 + Log((0.8494552867973413*Sqr
      (SCALE))/Abs(msq2(2,2)))))/Abs(msq2(2,2)))))/Abs(msu2(2,2)) - (
      1.2772891249999998*(AtInput - MuInput/TanBeta)*(0.022663710673270766*
      ComplexLog(0.15054471320265872) + 0.8494552867973413*(3.6617324413227936 + 2
      *Log((1.*Sqr(SCALE))/Abs(msq2(2,2))))))/Sqrt(Abs(msq2(2,2))) + (
      22.57017915079702*(-0.13859710249514381*(-1 + (0.8494552867973413*Abs(msu2(2
      ,2)))/Abs(msq2(2,2))) + (0.12788100252938794*Abs(msu2(2,2))*Log((
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))/Abs(msq2(2,2)))*Sqr(
      AtInput - MuInput/TanBeta))/(Abs(msq2(2,2))*(0.8494552867973413 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*(-1 + (0.8494552867973413
      *Abs(msu2(2,2)))/Abs(msq2(2,2)))) + (0.9216589861751152*(AtInput - MuInput/
      TanBeta)*(-11.588731795280644 + (11.588731795280644*Abs(msu2(2,2)))/Abs(msq2
      (2,2)) + (0.8494552867973413*Abs(msu2(2,2)))/(Abs(msq2(2,2))*(1 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))) + 0.8494552867973413/(-1
       + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) - 7.199173001148221*(
      -0.8494552867973413 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) -
      0.1067211471893647*ComplexLog(0.15054471320265872) - 8*ComplexLog(1 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) + (4.7089*Abs(msq2(2,2))*
      ComplexLog(1 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))/Abs(msu2
      (2,2)) + (3.397821147189365*Abs(msu2(2,2))*ComplexLog(1 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))/Abs(msq2(2,2)) +
      7.6450975811760715*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)) - (0.8494552867973413
      *Abs(msu2(2,2))*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)))/Abs(msq2(2,2)) +
      0.8494552867973413*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)) - (7.6450975811760715
      *Abs(msu2(2,2))*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)))/Abs(msq2(2,2)) +
      3.397821147189365*Log((0.8494552867973413*Sqr(SCALE))/Abs(msq2(2,2))) - (
      3.397821147189365*Abs(msu2(2,2))*Log((0.8494552867973413*Sqr(SCALE))/Abs(
      msq2(2,2))))/Abs(msq2(2,2)) - (0.8494552867973413*Log((0.8494552867973413*
      Abs(msu2(2,2)))/Abs(msq2(2,2))))/Sqr(-1 + (0.8494552867973413*Abs(msu2(2,2))
      )/Abs(msq2(2,2))) + (0.8494552867973413*Abs(msu2(2,2))*Log((
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))/(Abs(msq2(2,2))*Sqr(-1 +
      (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))))/(Sqrt(Abs(msq2(2,2)))
      *(0.8494552867973413 - (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))))
      + (3.8318673749999994*(AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs(msq2(2,2))
      ))*((0.8494552867973413*Abs(msu2(2,2))*(3 + (0.8494552867973413*Abs(msu2(2,2
      )))/Abs(msq2(2,2)) - (-2 + (0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)
      ))*Log((0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2))) + 2*Log((1.*Sqr(
      SCALE))/Abs(msu2(2,2)))))/Abs(msq2(2,2)) + ComplexLog(1 - (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr(-1 + (
      0.8494552867973413*Abs(msu2(2,2)))/Abs(msq2(2,2)))))/Sqr(Abs(msu2(2,2))))))/
      (Abs(msu2(2,2))*Quad(-1 + (1.*Abs(msq2(2,2)))/Abs(msu2(2,2))))))/Sqrt(Abs(
      msq2(2,2)))))/Sqrt(Abs(msq2(2,2)))))/Sqr(Abs(msu2(2,2)))))/Sqrt(Abs(msq2(2,2
      ))) + (0.015625*(1 - (5.88235294117647*(M3Input - 0.915*Sqrt(Abs(msq2(2,2)))
      ))/Sqrt(Abs(msq2(2,2))))*Re(-1.4018914012500006*(0.037799933149422246*
      ComplexLog(-0.1944220490310249) + 1.194422049031025*(4.33754298327075 + 2*
      Log((0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2))))) - (1.4018914012500006*
      Sqr(Abs(msq2(2,2)))*((1.194422049031025*Abs(msu2(2,2))*(3 + (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)) - (-2 + (1.194422049031025*
      Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((1.194422049031025*Abs(msu2(2,2)))/Abs(
      msq2(2,2))) + 2*Log((0.9999999999999999*Sqr(SCALE))/Abs(msu2(2,2)))))/Abs(
      msq2(2,2)) + ComplexLog(1 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)
      ))*Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))))/Sqr(Abs(
      msu2(2,2))) + (0.5369670767482759*(AtInput - MuInput/TanBeta)*Power3(Sqrt(
      Abs(msq2(2,2))))*((2.0353131917913196*((2.38884409806205*(1.194422049031025
      - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr(AtInput - MuInput/
      TanBeta))/Abs(msq2(2,2)) + Log((0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2(
      2,2)))*(2.853288062422944 + (1.194422049031025*Abs(msu2(2,2))*((
      2.38884409806205*Abs(msu2(2,2)))/Abs(msq2(2,2)) - (1.194422049031025*Sqr(
      AtInput - MuInput/TanBeta))/Abs(msq2(2,2))))/Abs(msq2(2,2)) -
      1.194422049031025*((4.7776881961241*Abs(msu2(2,2)))/Abs(msq2(2,2)) + (
      1.194422049031025*Sqr(AtInput - MuInput/TanBeta))/Abs(msq2(2,2)))))*Sqr(Abs(
      msu2(2,2)))*((0.7660608750000002*(AtInput - MuInput/TanBeta)*(
      0.037799933149422246*ComplexLog(-0.1944220490310249) + 1.194422049031025*(
      4.33754298327075 + 2*Log((0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2))))))/
      Sqrt(Abs(msq2(2,2))) - (2.803782802500001*Abs(msq2(2,2))*((
      0.23222198218044715*Abs(msu2(2,2))*ComplexLog(-0.1944220490310249))/Abs(msq2
      (2,2)) + 1.194422049031025*((-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(
      msq2(2,2)))*ComplexLog(1 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))
      ) - (2.38884409806205*Abs(msu2(2,2))*(2 + Log((1.194422049031025*Sqr(SCALE))
      /Abs(msq2(2,2)))))/Abs(msq2(2,2)))))/Abs(msu2(2,2)) - (24.573798187682407*(
      0.21220392058673745*(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))
      - (0.23222198218044715*Abs(msu2(2,2))*Log((1.194422049031025*Abs(msu2(2,2)))
      /Abs(msq2(2,2))))/Abs(msq2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/(Abs(msq2(
      2,2))*(1.194422049031025 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))
      )*(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))) + (
      1.0928961748633879*(AtInput - MuInput/TanBeta)*(-2.2175047962965717 + (
      2.2175047962965717*Abs(msu2(2,2)))/Abs(msq2(2,2)) + (1.194422049031025*Abs(
      msu2(2,2)))/(Abs(msq2(2,2))*(1 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2
      (2,2)))) + 4.700072529518397*(-1.194422049031025 + (1.194422049031025*Abs(
      msu2(2,2)))/Abs(msq2(2,2))) + 1.194422049031025/(-1 + (1.194422049031025*Abs
      (msu2(2,2)))/Abs(msq2(2,2))) - 0.1265881961241002*ComplexLog(-
      0.1944220490310249) - 8*ComplexLog(1 - (1.194422049031025*Abs(msu2(2,2)))/
      Abs(msq2(2,2))) + (3.3489000000000004*Abs(msq2(2,2))*ComplexLog(1 - (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))))/Abs(msu2(2,2)) + (
      4.7776881961241*Abs(msu2(2,2))*ComplexLog(1 - (1.194422049031025*Abs(msu2(2,
      2)))/Abs(msq2(2,2))))/Abs(msq2(2,2)) + 10.749798441279225*Log((
      0.9999999999999999*Abs(msq2(2,2)))/Sqr(SCALE)) - (1.194422049031025*Abs(msu2
      (2,2))*Log((0.9999999999999999*Abs(msq2(2,2)))/Sqr(SCALE)))/Abs(msq2(2,2)) +
      1.194422049031025*Log((0.9999999999999999*Abs(msu2(2,2)))/Sqr(SCALE)) - (
      10.749798441279225*Abs(msu2(2,2))*Log((0.9999999999999999*Abs(msu2(2,2)))/
      Sqr(SCALE)))/Abs(msq2(2,2)) + 4.7776881961241*Log((1.194422049031025*Sqr(
      SCALE))/Abs(msq2(2,2))) - (4.7776881961241*Abs(msu2(2,2))*Log((
      1.194422049031025*Sqr(SCALE))/Abs(msq2(2,2))))/Abs(msq2(2,2)) - (
      1.194422049031025*Log((1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))))/
      Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))) + (
      1.194422049031025*Abs(msu2(2,2))*Log((1.194422049031025*Abs(msu2(2,2)))/Abs(
      msq2(2,2))))/(Abs(msq2(2,2))*Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs
      (msq2(2,2))))))/(Sqrt(Abs(msq2(2,2)))*(1.194422049031025 - (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))) + (0.7660608750000002*(
      AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs(msq2(2,2))))*((1.194422049031025*
      Abs(msu2(2,2))*(3 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)) - (-2
      + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((1.194422049031025*
      Abs(msu2(2,2)))/Abs(msq2(2,2))) + 2*Log((0.9999999999999999*Sqr(SCALE))/Abs(
      msu2(2,2)))))/Abs(msq2(2,2)) + ComplexLog(1 - (1.194422049031025*Abs(msu2(2,
      2)))/Abs(msq2(2,2)))*Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,
      2)))))/Sqr(Abs(msu2(2,2)))))/(Cube(1.194422049031025 - (1.194422049031025*
      Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr(Abs(msq2(2,2)))) + (1.0928961748633879*(
      AtInput - MuInput/TanBeta)*((2*(-2.38884409806205 + (2.38884409806205*Abs(
      msu2(2,2)))/Abs(msq2(2,2)) + (1.194422049031025 + (1.194422049031025*Abs(
      msu2(2,2)))/Abs(msq2(2,2)))*Log((0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2
      (2,2))))*((0.05392704900781591*ComplexLog(-0.1944220490310249)*Sqr(Abs(msu2(
      2,2))))/Sqr(Abs(msq2(2,2))) - 1.194422049031025*((1.194422049031025*Abs(msu2
      (2,2))*(3.5832661470930747 - (3.5832661470930747*Abs(msu2(2,2)))/Abs(msq2(2,
      2)) - 1.194422049031025*(-2 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,
      2)))*Log((1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))) - (
      2.38884409806205*Abs(msu2(2,2))*Log((0.9999999999999999*Sqr(SCALE))/Abs(msq2
      (2,2))))/Abs(msq2(2,2)) + 2.38884409806205*Log((0.9999999999999999*Sqr(SCALE
      ))/Abs(msu2(2,2)))))/Abs(msq2(2,2)) - (0.20418262657451347*Sqr(Abs(msu2(2,2)
      )))/Sqr(Abs(msq2(2,2))) + 1.194422049031025*ComplexLog(1 - (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Sqr(-1 + (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))))))/Sqr(1.194422049031025 -
      (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))) + (1.5591738046027015*(
      AtInput - MuInput/TanBeta)*((1.5321217500000004*(AtInput - MuInput/TanBeta)*
      (-1.194422049031025 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)) +
      1.194422049031025*Log((0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2(2,2))))*(
      (0.05392704900781591*ComplexLog(-0.1944220490310249)*Sqr(Abs(msu2(2,2))))/
      Sqr(Abs(msq2(2,2))) - 1.194422049031025*((1.194422049031025*Abs(msu2(2,2))*(
      3.5832661470930747 - (3.5832661470930747*Abs(msu2(2,2)))/Abs(msq2(2,2)) -
      1.194422049031025*(-2 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*
      Log((1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))) - (2.38884409806205*
      Abs(msu2(2,2))*Log((0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2))))/Abs(msq2
      (2,2)) + 2.38884409806205*Log((0.9999999999999999*Sqr(SCALE))/Abs(msu2(2,2))
      )))/Abs(msq2(2,2)) - (0.20418262657451347*Sqr(Abs(msu2(2,2))))/Sqr(Abs(msq2(
      2,2))) + 1.194422049031025*ComplexLog(1 - (1.194422049031025*Abs(msu2(2,2)))
      /Abs(msq2(2,2)))*Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))
      )))/(Sqrt(Abs(msq2(2,2)))*Cube(1.194422049031025 - (1.194422049031025*Abs(
      msu2(2,2)))/Abs(msq2(2,2)))) - (0.8372250000000001*Abs(msq2(2,2))*(-
      2.38884409806205 + (2.38884409806205*Abs(msu2(2,2)))/Abs(msq2(2,2)) + (
      1.194422049031025 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((
      0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2(2,2))))*((-1 + (
      0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2(2,2)))*((-0.7660608750000002*(
      AtInput - MuInput/TanBeta)*(0.037799933149422246*ComplexLog(-
      0.1944220490310249) + 1.194422049031025*(4.33754298327075 + 2*Log((
      0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2))))))/Sqrt(Abs(msq2(2,2))) - (
      2.803782802500001*Abs(msq2(2,2))*((0.23222198218044715*Abs(msu2(2,2))*
      ComplexLog(-0.1944220490310249))/Abs(msq2(2,2)) + 1.194422049031025*((-1 + (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*ComplexLog(1 - (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))) - (2.38884409806205*Abs(
      msu2(2,2))*(2 + Log((1.194422049031025*Sqr(SCALE))/Abs(msq2(2,2)))))/Abs(
      msq2(2,2)))))/Abs(msu2(2,2)) - (24.573798187682407*(0.21220392058673745*(-1
      + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))) - (0.23222198218044715*
      Abs(msu2(2,2))*Log((1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))))/Abs(
      msq2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/(Abs(msq2(2,2))*(
      1.194422049031025 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*(-1 +
      (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))) + (1.0928961748633879*(
      AtInput - MuInput/TanBeta)*(-2.2175047962965717 + (2.2175047962965717*Abs(
      msu2(2,2)))/Abs(msq2(2,2)) + (1.194422049031025*Abs(msu2(2,2)))/(Abs(msq2(2,
      2))*(1 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))) +
      4.700072529518397*(-1.194422049031025 + (1.194422049031025*Abs(msu2(2,2)))/
      Abs(msq2(2,2))) + 1.194422049031025/(-1 + (1.194422049031025*Abs(msu2(2,2)))
      /Abs(msq2(2,2))) - 0.1265881961241002*ComplexLog(-0.1944220490310249) - 8*
      ComplexLog(1 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))) + (
      3.3489000000000004*Abs(msq2(2,2))*ComplexLog(1 - (1.194422049031025*Abs(msu2
      (2,2)))/Abs(msq2(2,2))))/Abs(msu2(2,2)) + (4.7776881961241*Abs(msu2(2,2))*
      ComplexLog(1 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))))/Abs(msq2(
      2,2)) + 10.749798441279225*Log((0.9999999999999999*Abs(msq2(2,2)))/Sqr(SCALE
      )) - (1.194422049031025*Abs(msu2(2,2))*Log((0.9999999999999999*Abs(msq2(2,2)
      ))/Sqr(SCALE)))/Abs(msq2(2,2)) + 1.194422049031025*Log((0.9999999999999999*
      Abs(msu2(2,2)))/Sqr(SCALE)) - (10.749798441279225*Abs(msu2(2,2))*Log((
      0.9999999999999999*Abs(msu2(2,2)))/Sqr(SCALE)))/Abs(msq2(2,2)) +
      4.7776881961241*Log((1.194422049031025*Sqr(SCALE))/Abs(msq2(2,2))) - (
      4.7776881961241*Abs(msu2(2,2))*Log((1.194422049031025*Sqr(SCALE))/Abs(msq2(2
      ,2))))/Abs(msq2(2,2)) - (1.194422049031025*Log((1.194422049031025*Abs(msu2(2
      ,2)))/Abs(msq2(2,2))))/Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(
      2,2))) + (1.194422049031025*Abs(msu2(2,2))*Log((1.194422049031025*Abs(msu2(2
      ,2)))/Abs(msq2(2,2))))/(Abs(msq2(2,2))*Sqr(-1 + (1.194422049031025*Abs(msu2(
      2,2)))/Abs(msq2(2,2))))))/(Sqrt(Abs(msq2(2,2)))*(1.194422049031025 - (
      1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))) + (2.298182625000001*(
      AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs(msq2(2,2))))*((1.194422049031025*
      Abs(msu2(2,2))*(3 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)) - (-2
      + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((1.194422049031025*
      Abs(msu2(2,2)))/Abs(msq2(2,2))) + 2*Log((0.9999999999999999*Sqr(SCALE))/Abs(
      msu2(2,2)))))/Abs(msq2(2,2)) + ComplexLog(1 - (1.194422049031025*Abs(msu2(2,
      2)))/Abs(msq2(2,2)))*Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,
      2)))))/Sqr(Abs(msu2(2,2)))) + (3.2218024604896556*(AtInput - MuInput/TanBeta
      )*Power(Abs(msq2(2,2)),2.5)*((0.05392704900781591*ComplexLog(-
      0.1944220490310249)*Sqr(Abs(msu2(2,2))))/Sqr(Abs(msq2(2,2))) -
      1.194422049031025*((1.194422049031025*Abs(msu2(2,2))*(3.5832661470930747 - (
      3.5832661470930747*Abs(msu2(2,2)))/Abs(msq2(2,2)) - 1.194422049031025*(-2 +
      (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2)))*Log((1.194422049031025*
      Abs(msu2(2,2)))/Abs(msq2(2,2))) - (2.38884409806205*Abs(msu2(2,2))*Log((
      0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2))))/Abs(msq2(2,2)) +
      2.38884409806205*Log((0.9999999999999999*Sqr(SCALE))/Abs(msu2(2,2)))))/Abs(
      msq2(2,2)) - (0.20418262657451347*Sqr(Abs(msu2(2,2))))/Sqr(Abs(msq2(2,2))) +
      1.194422049031025*ComplexLog(1 - (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2
      (2,2)))*Sqr(-1 + (1.194422049031025*Abs(msu2(2,2)))/Abs(msq2(2,2))))))/Cube(
      Abs(msu2(2,2)))))/(Abs(msu2(2,2))*Quad(-1 + (0.9999999999999999*Abs(msq2(2,2
      )))/Abs(msu2(2,2))))))/Sqrt(Abs(msq2(2,2)))))/Sqrt(Abs(msq2(2,2)))))/Sqr(Abs
      (msu2(2,2)))))/Quad(3.141592653589793), Abs(-1 + Sqrt(Abs(msu2(2,2)))/
      M3Input) < 0.085, (0.0009435645454673102*(M3Input - 0.915*Sqrt(Abs(msu2(2,2)
      )))*Re(-2.7717174012499997*(0.022663710673270766*ComplexLog(
      0.15054471320265872) + 0.8494552867973413*(3.6617324413227936 + 2*Log((1.*
      Sqr(SCALE))/Abs(msu2(2,2))))) + (1.7701422470949426*(AtInput - MuInput/
      TanBeta)*Power3(Sqrt(Abs(msu2(2,2))))*((0.9216589861751152*(AtInput -
      MuInput/TanBeta)*((0.6650454232884363*(AtInput - MuInput/TanBeta)*Sqr(Abs(
      msq2(2,2)))*((2.5545782499999996*(AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs
      (msu2(2,2))))*(0.8494552867973413 - (0.8494552867973413*Abs(msq2(2,2)))/Abs(
      msu2(2,2)) + (0.8494552867973413*Abs(msq2(2,2))*Log((1.*Abs(msq2(2,2)))/Abs(
      msu2(2,2))))/Abs(msu2(2,2)))*((-0.8494552867973413*Abs(msq2(2,2))*((
      0.01925180884985518*Abs(msq2(2,2))*ComplexLog(0.15054471320265872))/Abs(msu2
      (2,2)) + 0.7215742842679533*(-2 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(
      msu2(2,2)))*Log((0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))) +
      0.8494552867973413*(-2.548365860392024 + (2.3889036968510293*Abs(msq2(2,2)))
      /Abs(msu2(2,2)) - 1.6989105735946826*Log((1.*Sqr(SCALE))/Abs(msq2(2,2))) + (
      1.6989105735946826*Abs(msq2(2,2))*Log((1.*Sqr(SCALE))/Abs(msu2(2,2))))/Abs(
      msu2(2,2)))))/Abs(msu2(2,2)) + 0.7215742842679533*ComplexLog(1 - (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/(Cube(-
      0.8494552867973413 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr
      (Abs(msq2(2,2)))) - (1.177225*(1.6989105735946826 - (1.6989105735946826*Abs(
      msq2(2,2)))/Abs(msu2(2,2)) + (0.8494552867973413 + (0.8494552867973413*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*Log((1.*Abs(msq2(2,2)))/Abs(msu2(2,2))))*((
      10.620853482569656*(AtInput - MuInput/TanBeta)*Sqrt(Abs(msu2(2,2)))*((-
      0.8494552867973413*Abs(msq2(2,2))*((0.01925180884985518*Abs(msq2(2,2))*
      ComplexLog(0.15054471320265872))/Abs(msu2(2,2)) + 0.7215742842679533*(-2 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.8494552867973413*
      Abs(msq2(2,2)))/Abs(msu2(2,2))) + 0.8494552867973413*(-2.548365860392024 + (
      2.3889036968510293*Abs(msq2(2,2)))/Abs(msu2(2,2)) - 1.6989105735946826*Log((
      1.*Sqr(SCALE))/Abs(msq2(2,2))) + (1.6989105735946826*Abs(msq2(2,2))*Log((1.*
      Sqr(SCALE))/Abs(msu2(2,2))))/Abs(msu2(2,2)))))/Abs(msu2(2,2)) +
      0.7215742842679533*ComplexLog(1 - (0.8494552867973413*Abs(msq2(2,2)))/Abs(
      msu2(2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/
      Abs(msq2(2,2)) + (-1 + (1.*Abs(msq2(2,2)))/Abs(msu2(2,2)))*((-
      5.543434802499999*Abs(msu2(2,2))*(0.8494552867973413*(-1 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*ComplexLog(1 - (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))) + (0.8494552867973413*Abs
      (msq2(2,2))*(-0.15054471320265872*ComplexLog(0.15054471320265872) -
      1.6989105735946826*(2 + Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(2,2))))
      ))/Abs(msu2(2,2))))/Abs(msq2(2,2)) + (3.8318673749999994*(AtInput - MuInput/
      TanBeta)*(0.022663710673270766*ComplexLog(0.15054471320265872) +
      0.8494552867973413*(3.6617324413227936 + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2)
      )))))/Sqrt(Abs(msu2(2,2))) + (22.57017915079702*(0.13859710249514381*(-1 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))) - (0.12788100252938794*
      Abs(msq2(2,2))*Log((0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))))/Abs(
      msu2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/((-1 + (0.8494552867973413*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*(-0.8494552867973413 + (0.8494552867973413*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*Abs(msu2(2,2))) - (1.2772891249999998*(AtInput -
      MuInput/TanBeta)*Power3(Sqrt(Abs(msu2(2,2))))*((0.8494552867973413*Abs(msq2(
      2,2))*(3 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)) - (-2 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.8494552867973413*
      Abs(msq2(2,2)))/Abs(msu2(2,2))) + 2*Log((1.*Sqr(SCALE))/Abs(msq2(2,2)))))/
      Abs(msu2(2,2)) + ComplexLog(1 - (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2
      (2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/Sqr(
      Abs(msq2(2,2))) + (0.9216589861751152*(AtInput - MuInput/TanBeta)*(
      5.473356228886606 + 0.8494552867973413/(1 - (0.8494552867973413*Abs(msq2(2,2
      )))/Abs(msu2(2,2))) - (5.473356228886606*Abs(msq2(2,2)))/Abs(msu2(2,2)) + (
      0.8494552867973413*Abs(msq2(2,2)))/((-1 + (0.8494552867973413*Abs(msq2(2,2))
      )/Abs(msu2(2,2)))*Abs(msu2(2,2))) + 0.10672114718936498*ComplexLog(
      0.15054471320265872) - 0.8494552867973413*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)
      ) + (7.6450975811760715*Abs(msq2(2,2))*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)))/
      Abs(msu2(2,2)) - 7.6450975811760715*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)) + (
      0.8494552867973413*Abs(msq2(2,2))*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)))/Abs(
      msu2(2,2)) - 3.397821147189365*Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(
      2,2))) + (3.397821147189365*Abs(msq2(2,2))*Log((0.8494552867973413*Sqr(SCALE
      ))/Abs(msu2(2,2))))/Abs(msu2(2,2)) + ((0.8494552867973413 - (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.8494552867973413*
      Abs(msq2(2,2)))/Abs(msu2(2,2))))/Sqr(-1 + (0.8494552867973413*Abs(msq2(2,2))
      )/Abs(msu2(2,2))) - (4.7089*Abs(msu2(2,2))*ComplexLog(1 - (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))))/Abs(msq2(2,2))))/((-
      0.8494552867973413 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*
      Sqrt(Abs(msu2(2,2)))))))/Quad(-1 + (1.*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/
      Power(Abs(msu2(2,2)),2.5) + (2*(1.6989105735946826 - (1.6989105735946826*Abs
      (msq2(2,2)))/Abs(msu2(2,2)) + (0.8494552867973413 + (0.8494552867973413*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*Log((1.*Abs(msq2(2,2)))/Abs(msu2(2,2))))*((-
      0.8494552867973413*Abs(msq2(2,2))*((0.01925180884985518*Abs(msq2(2,2))*
      ComplexLog(0.15054471320265872))/Abs(msu2(2,2)) + 0.7215742842679533*(-2 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.8494552867973413*
      Abs(msq2(2,2)))/Abs(msu2(2,2))) + 0.8494552867973413*(-2.548365860392024 + (
      2.3889036968510293*Abs(msq2(2,2)))/Abs(msu2(2,2)) - 1.6989105735946826*Log((
      1.*Sqr(SCALE))/Abs(msq2(2,2))) + (1.6989105735946826*Abs(msq2(2,2))*Log((1.*
      Sqr(SCALE))/Abs(msu2(2,2))))/Abs(msu2(2,2)))))/Abs(msu2(2,2)) +
      0.7215742842679533*ComplexLog(1 - (0.8494552867973413*Abs(msq2(2,2)))/Abs(
      msu2(2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/
      Sqr(-0.8494552867973413 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))
      )))/Sqrt(Abs(msu2(2,2))) + (0.5206694477168091*Sqr(Abs(msq2(2,2)))*((-
      5.543434802499999*Abs(msu2(2,2))*(0.8494552867973413*(-1 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*ComplexLog(1 - (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))) + (0.8494552867973413*Abs
      (msq2(2,2))*(-0.15054471320265872*ComplexLog(0.15054471320265872) -
      1.6989105735946826*(2 + Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(2,2))))
      ))/Abs(msu2(2,2))))/Abs(msq2(2,2)) + (1.2772891249999998*(AtInput - MuInput/
      TanBeta)*(0.022663710673270766*ComplexLog(0.15054471320265872) +
      0.8494552867973413*(3.6617324413227936 + 2*Log((1.*Sqr(SCALE))/Abs(msu2(2,2)
      )))))/Sqrt(Abs(msu2(2,2))) + (22.57017915079702*(0.13859710249514381*(-1 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))) - (0.12788100252938794*
      Abs(msq2(2,2))*Log((0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))))/Abs(
      msu2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/((-1 + (0.8494552867973413*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*(-0.8494552867973413 + (0.8494552867973413*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*Abs(msu2(2,2))) + (1.2772891249999998*(AtInput -
      MuInput/TanBeta)*Power3(Sqrt(Abs(msu2(2,2))))*((0.8494552867973413*Abs(msq2(
      2,2))*(3 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)) - (-2 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.8494552867973413*
      Abs(msq2(2,2)))/Abs(msu2(2,2))) + 2*Log((1.*Sqr(SCALE))/Abs(msq2(2,2)))))/
      Abs(msu2(2,2)) + ComplexLog(1 - (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2
      (2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/Sqr(
      Abs(msq2(2,2))) + (0.9216589861751152*(AtInput - MuInput/TanBeta)*(
      5.473356228886606 + 0.8494552867973413/(1 - (0.8494552867973413*Abs(msq2(2,2
      )))/Abs(msu2(2,2))) - (5.473356228886606*Abs(msq2(2,2)))/Abs(msu2(2,2)) + (
      0.8494552867973413*Abs(msq2(2,2)))/((-1 + (0.8494552867973413*Abs(msq2(2,2))
      )/Abs(msu2(2,2)))*Abs(msu2(2,2))) + 0.10672114718936498*ComplexLog(
      0.15054471320265872) - 0.8494552867973413*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)
      ) + (7.6450975811760715*Abs(msq2(2,2))*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)))/
      Abs(msu2(2,2)) - 7.6450975811760715*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)) + (
      0.8494552867973413*Abs(msq2(2,2))*Log((1.*Abs(msu2(2,2)))/Sqr(SCALE)))/Abs(
      msu2(2,2)) - 3.397821147189365*Log((0.8494552867973413*Sqr(SCALE))/Abs(msu2(
      2,2))) + (3.397821147189365*Abs(msq2(2,2))*Log((0.8494552867973413*Sqr(SCALE
      ))/Abs(msu2(2,2))))/Abs(msu2(2,2)) + ((0.8494552867973413 - (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.8494552867973413*
      Abs(msq2(2,2)))/Abs(msu2(2,2))))/Sqr(-1 + (0.8494552867973413*Abs(msq2(2,2))
      )/Abs(msu2(2,2))) - (4.7089*Abs(msu2(2,2))*ComplexLog(1 - (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))))/Abs(msq2(2,2))))/((-
      0.8494552867973413 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*
      Sqrt(Abs(msu2(2,2)))))*((1.6989105735946823*(-0.8494552867973413 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(AtInput - MuInput/
      TanBeta))/Abs(msu2(2,2)) + Log((1.*Abs(msq2(2,2)))/Abs(msu2(2,2)))*(
      0.8494552867973413*(1.6989105735946826 - (0.8494552867973412*Sqr(AtInput -
      MuInput/TanBeta))/Abs(msu2(2,2))) - (0.8494552867973413*Abs(msq2(2,2))*(
      3.397821147189365 + (0.8494552867973412*Sqr(AtInput - MuInput/TanBeta))/Abs(
      msu2(2,2))))/Abs(msu2(2,2)) + (1.4431485685359067*Sqr(Abs(msq2(2,2))))/Sqr(
      Abs(msu2(2,2))))))/(Cube(-0.8494552867973413 + (0.8494552867973413*Abs(msq2(
      2,2)))/Abs(msu2(2,2)))*Sqr(Abs(msu2(2,2))))))/Sqr(Abs(msq2(2,2))) - (
      2.7717174012499997*((0.8494552867973413*Abs(msq2(2,2))*(3 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)) - (-2 + (
      0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.8494552867973413*
      Abs(msq2(2,2)))/Abs(msu2(2,2))) + 2*Log((1.*Sqr(SCALE))/Abs(msq2(2,2)))))/
      Abs(msu2(2,2)) + ComplexLog(1 - (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2
      (2,2)))*Sqr(-1 + (0.8494552867973413*Abs(msq2(2,2)))/Abs(msu2(2,2))))*Sqr(
      Abs(msu2(2,2))))/Sqr(Abs(msq2(2,2)))))/Sqrt(Abs(msu2(2,2))) + (0.015625*(1 -
      (5.88235294117647*(M3Input - 0.915*Sqrt(Abs(msu2(2,2)))))/Sqrt(Abs(msu2(2,2)
      )))*Re(-1.4018914012500006*(0.037799933149422246*ComplexLog(-
      0.1944220490310249) + 1.194422049031025*(4.33754298327075 + 2*Log((
      0.9999999999999999*Sqr(SCALE))/Abs(msu2(2,2))))) + (0.5369670767482759*(
      AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs(msu2(2,2))))*((1.0928961748633879
      *(AtInput - MuInput/TanBeta)*((2*(2.38884409806205 - (2.38884409806205*Abs(
      msq2(2,2)))/Abs(msu2(2,2)) + (1.194422049031025 + (1.194422049031025*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*Log((0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2
      (2,2))))*((-1.194422049031025*Abs(msq2(2,2))*((0.04514907360556868*Abs(msq2(
      2,2))*ComplexLog(-0.1944220490310249))/Abs(msu2(2,2)) + 1.426644031211472*(-
      2 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + 1.194422049031025*(-
      3.5832661470930747 + (3.7542129466269216*Abs(msq2(2,2)))/Abs(msu2(2,2)) -
      2.38884409806205*Log((0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2))) + (
      2.38884409806205*Abs(msq2(2,2))*Log((0.9999999999999999*Sqr(SCALE))/Abs(msu2
      (2,2))))/Abs(msu2(2,2)))))/Abs(msu2(2,2)) + 1.426644031211472*ComplexLog(1 -
      (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/Sqr(-1.194422049031025 +
      (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + (1.5591738046027015*(
      AtInput - MuInput/TanBeta)*Sqr(Abs(msq2(2,2)))*((1.5321217500000004*(AtInput
       - MuInput/TanBeta)*Power3(Sqrt(Abs(msu2(2,2))))*(1.194422049031025 - (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)) + (1.194422049031025*Abs(
      msq2(2,2))*Log((0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2(2,2))))/Abs(msu2
      (2,2)))*((-1.194422049031025*Abs(msq2(2,2))*((0.04514907360556868*Abs(msq2(2
      ,2))*ComplexLog(-0.1944220490310249))/Abs(msu2(2,2)) + 1.426644031211472*(-2
       + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((1.194422049031025
      *Abs(msq2(2,2)))/Abs(msu2(2,2))) + 1.194422049031025*(-3.5832661470930747 +
      (3.7542129466269216*Abs(msq2(2,2)))/Abs(msu2(2,2)) - 2.38884409806205*Log((
      0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2))) + (2.38884409806205*Abs(msq2(
      2,2))*Log((0.9999999999999999*Sqr(SCALE))/Abs(msu2(2,2))))/Abs(msu2(2,2)))))
      /Abs(msu2(2,2)) + 1.426644031211472*ComplexLog(1 - (1.194422049031025*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(
      msu2(2,2)))))/(Cube(-1.194422049031025 + (1.194422049031025*Abs(msq2(2,2)))/
      Abs(msu2(2,2)))*Sqr(Abs(msq2(2,2)))) - (0.8372250000000001*(2.38884409806205
       - (2.38884409806205*Abs(msq2(2,2)))/Abs(msu2(2,2)) + (1.194422049031025 + (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((0.9999999999999999*
      Abs(msq2(2,2)))/Abs(msu2(2,2))))*((3.2218024604896556*(AtInput - MuInput/
      TanBeta)*Sqrt(Abs(msu2(2,2)))*((-1.194422049031025*Abs(msq2(2,2))*((
      0.04514907360556868*Abs(msq2(2,2))*ComplexLog(-0.1944220490310249))/Abs(msu2
      (2,2)) + 1.426644031211472*(-2 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2
      (2,2)))*Log((1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) +
      1.194422049031025*(-3.5832661470930747 + (3.7542129466269216*Abs(msq2(2,2)))
      /Abs(msu2(2,2)) - 2.38884409806205*Log((0.9999999999999999*Sqr(SCALE))/Abs(
      msq2(2,2))) + (2.38884409806205*Abs(msq2(2,2))*Log((0.9999999999999999*Sqr(
      SCALE))/Abs(msu2(2,2))))/Abs(msu2(2,2)))))/Abs(msu2(2,2)) +
      1.426644031211472*ComplexLog(1 - (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2
      (2,2)))*Sqr(-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/Abs(
      msq2(2,2)) + (-1 + (0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2(2,2)))*((
      2.298182625000001*(AtInput - MuInput/TanBeta)*(0.037799933149422246*
      ComplexLog(-0.1944220490310249) + 1.194422049031025*(4.33754298327075 + 2*
      Log((0.9999999999999999*Sqr(SCALE))/Abs(msu2(2,2))))))/Sqrt(Abs(msu2(2,2)))
      - (2.803782802500001*Abs(msu2(2,2))*(1.194422049031025*(-1 + (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*ComplexLog(1 - (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + (1.194422049031025*Abs(
      msq2(2,2))*(0.1944220490310249*ComplexLog(-0.1944220490310249) -
      2.38884409806205*(2 + Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2,2))))))/
      Abs(msu2(2,2))))/Abs(msq2(2,2)) - (24.573798187682407*(-0.21220392058673745*
      (-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + (
      0.23222198218044715*Abs(msq2(2,2))*Log((1.194422049031025*Abs(msq2(2,2)))/
      Abs(msu2(2,2))))/Abs(msu2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/((-
      1.194422049031025 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*(-1 +
      (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Abs(msu2(2,2))) - (
      0.7660608750000002*(AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs(msu2(2,2))))*
      ((1.194422049031025*Abs(msq2(2,2))*(3 + (1.194422049031025*Abs(msq2(2,2)))/
      Abs(msu2(2,2)) - (-2 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*
      Log((1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + 2*Log((
      0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2)))))/Abs(msu2(2,2)) + ComplexLog
      (1 - (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/Sqr(Abs(msq2(2,2))) + (
      1.0928961748633879*(AtInput - MuInput/TanBeta)*(7.831375057598367 +
      1.194422049031025/(1 - (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) -
      (7.8313750575983665*Abs(msq2(2,2)))/Abs(msu2(2,2)) + (1.194422049031025*Abs(
      msq2(2,2)))/((-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Abs(
      msu2(2,2))) + 0.1265881961241*ComplexLog(-0.1944220490310249) -
      1.194422049031025*Log((0.9999999999999999*Abs(msq2(2,2)))/Sqr(SCALE)) + (
      10.749798441279225*Abs(msq2(2,2))*Log((0.9999999999999999*Abs(msq2(2,2)))/
      Sqr(SCALE)))/Abs(msu2(2,2)) - 10.749798441279225*Log((0.9999999999999999*Abs
      (msu2(2,2)))/Sqr(SCALE)) + (1.194422049031025*Abs(msq2(2,2))*Log((
      0.9999999999999999*Abs(msu2(2,2)))/Sqr(SCALE)))/Abs(msu2(2,2)) -
      4.7776881961241*Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2,2))) + (
      4.7776881961241*Abs(msq2(2,2))*Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2
      ,2))))/Abs(msu2(2,2)) + ((1.194422049031025 - (1.194422049031025*Abs(msq2(2,
      2)))/Abs(msu2(2,2)))*Log((1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))))
      /Sqr(-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) - (
      3.3489000000000004*Abs(msu2(2,2))*ComplexLog(1 - (1.194422049031025*Abs(msq2
      (2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2
      (2,2))))/Abs(msq2(2,2))))/((-1.194422049031025 + (1.194422049031025*Abs(msq2
      (2,2)))/Abs(msu2(2,2)))*Sqrt(Abs(msu2(2,2)))))))/Quad(-1 + (
      0.9999999999999999*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/Power(Abs(msu2(2,2)),
      2.5)))/Sqrt(Abs(msu2(2,2))) + (2.0353131917913196*Sqr(Abs(msq2(2,2)))*((
      0.7660608750000002*(AtInput - MuInput/TanBeta)*(0.037799933149422246*
      ComplexLog(-0.1944220490310249) + 1.194422049031025*(4.33754298327075 + 2*
      Log((0.9999999999999999*Sqr(SCALE))/Abs(msu2(2,2))))))/Sqrt(Abs(msu2(2,2)))
      - (2.803782802500001*Abs(msu2(2,2))*(1.194422049031025*(-1 + (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*ComplexLog(1 - (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + (1.194422049031025*Abs(
      msq2(2,2))*(0.1944220490310249*ComplexLog(-0.1944220490310249) -
      2.38884409806205*(2 + Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2,2))))))/
      Abs(msu2(2,2))))/Abs(msq2(2,2)) - (24.573798187682407*(-0.21220392058673745*
      (-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + (
      0.23222198218044715*Abs(msq2(2,2))*Log((1.194422049031025*Abs(msq2(2,2)))/
      Abs(msu2(2,2))))/Abs(msu2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/((-
      1.194422049031025 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*(-1 +
      (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Abs(msu2(2,2))) + (
      0.7660608750000002*(AtInput - MuInput/TanBeta)*Power3(Sqrt(Abs(msu2(2,2))))*
      ((1.194422049031025*Abs(msq2(2,2))*(3 + (1.194422049031025*Abs(msq2(2,2)))/
      Abs(msu2(2,2)) - (-2 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*
      Log((1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + 2*Log((
      0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2)))))/Abs(msu2(2,2)) + ComplexLog
      (1 - (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))))/Sqr(Abs(msq2(2,2))) + (
      1.0928961748633879*(AtInput - MuInput/TanBeta)*(7.831375057598367 +
      1.194422049031025/(1 - (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) -
      (7.8313750575983665*Abs(msq2(2,2)))/Abs(msu2(2,2)) + (1.194422049031025*Abs(
      msq2(2,2)))/((-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Abs(
      msu2(2,2))) + 0.1265881961241*ComplexLog(-0.1944220490310249) -
      1.194422049031025*Log((0.9999999999999999*Abs(msq2(2,2)))/Sqr(SCALE)) + (
      10.749798441279225*Abs(msq2(2,2))*Log((0.9999999999999999*Abs(msq2(2,2)))/
      Sqr(SCALE)))/Abs(msu2(2,2)) - 10.749798441279225*Log((0.9999999999999999*Abs
      (msu2(2,2)))/Sqr(SCALE)) + (1.194422049031025*Abs(msq2(2,2))*Log((
      0.9999999999999999*Abs(msu2(2,2)))/Sqr(SCALE)))/Abs(msu2(2,2)) -
      4.7776881961241*Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2,2))) + (
      4.7776881961241*Abs(msq2(2,2))*Log((1.194422049031025*Sqr(SCALE))/Abs(msu2(2
      ,2))))/Abs(msu2(2,2)) + ((1.194422049031025 - (1.194422049031025*Abs(msq2(2,
      2)))/Abs(msu2(2,2)))*Log((1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))))
      /Sqr(-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) - (
      3.3489000000000004*Abs(msu2(2,2))*ComplexLog(1 - (1.194422049031025*Abs(msq2
      (2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2
      (2,2))))/Abs(msq2(2,2))))/((-1.194422049031025 + (1.194422049031025*Abs(msq2
      (2,2)))/Abs(msu2(2,2)))*Sqrt(Abs(msu2(2,2)))))*((2.38884409806205*(-
      1.194422049031025 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(
      AtInput - MuInput/TanBeta))/Abs(msu2(2,2)) + Log((0.9999999999999999*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*(1.194422049031025*(2.38884409806205 - (
      1.194422049031025*Sqr(AtInput - MuInput/TanBeta))/Abs(msu2(2,2))) - (
      1.194422049031025*Abs(msq2(2,2))*(4.7776881961241 + (1.194422049031025*Sqr(
      AtInput - MuInput/TanBeta))/Abs(msu2(2,2))))/Abs(msu2(2,2)) + (
      2.853288062422944*Sqr(Abs(msq2(2,2))))/Sqr(Abs(msu2(2,2))))))/(Cube(-
      1.194422049031025 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(
      Abs(msu2(2,2))))))/Sqr(Abs(msq2(2,2))) - (1.4018914012500006*((
      1.194422049031025*Abs(msq2(2,2))*(3 + (1.194422049031025*Abs(msq2(2,2)))/Abs
      (msu2(2,2)) - (-2 + (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Log((
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))) + 2*Log((
      0.9999999999999999*Sqr(SCALE))/Abs(msq2(2,2)))))/Abs(msu2(2,2)) + ComplexLog
      (1 - (1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2)))*Sqr(-1 + (
      1.194422049031025*Abs(msq2(2,2)))/Abs(msu2(2,2))))*Sqr(Abs(msu2(2,2))))/Sqr(
      Abs(msq2(2,2)))))/Quad(3.141592653589793), True, (0.015625*Re((-2.*Quad(
      M3Input)*((1.*Abs(msq2(2,2))*(3 + 2*Log((1.*Sqr(SCALE))/Abs(msq2(2,2))) -
      Log((1.*Abs(msq2(2,2)))/Sqr(M3Input))*(-2 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)
      ) + (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - (1.*Abs
      (msq2(2,2)))/Sqr(M3Input))*Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input))))/Sqr(
      Abs(msq2(2,2))) - (2*Quad(M3Input)*((Abs(msu2(2,2))*(3 + 2*Log(Sqr(SCALE)/
      Abs(msu2(2,2))) - Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(msu2(2,2))/Sqr(
      M3Input)) + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - Abs(
      msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))))/Sqr(Abs(
      msu2(2,2))) + (1.*(AtInput - MuInput/TanBeta)*Power7(M3Input)*((1.*Sqr(Abs(
      msq2(2,2)))*((2*((1.*Abs(msq2(2,2)))/Sqr(M3Input) - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(AtInput - MuInput/TanBeta))/Sqr(M3Input) + Log((1.*Abs(msq2(2,
      2)))/Abs(msu2(2,2)))*((Abs(msu2(2,2))*((2*Abs(msu2(2,2)))/Sqr(M3Input) - Sqr
      (AtInput - MuInput/TanBeta)/Sqr(M3Input)))/Sqr(M3Input) - (1.*Abs(msq2(2,2))
      *((4*Abs(msu2(2,2)))/Sqr(M3Input) + Sqr(AtInput - MuInput/TanBeta)/Sqr(
      M3Input)))/Sqr(M3Input) + (2.*Sqr(Abs(msq2(2,2))))/Quad(M3Input)))*Sqr(Abs(
      msu2(2,2)))*((-4.*Quad(M3Input)*((Abs(msu2(2,2))*ComplexLog(1 - (1.*Abs(msq2
      (2,2)))/Sqr(M3Input))*(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/Sqr(M3Input)
      + (1.*Abs(msq2(2,2))*(ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*(-1 + Abs(
      msu2(2,2))/Sqr(M3Input)) - (2*Abs(msu2(2,2))*(2 + Log(Sqr(SCALE)/Sqr(M3Input
      ))))/Sqr(M3Input)))/Sqr(M3Input)))/(Abs(msq2(2,2))*Abs(msu2(2,2))) - (4*(-((
      Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-1 + (1.*Abs(msq2(2,2)))/
      Sqr(M3Input)))/Sqr(M3Input)) + (1.*Abs(msq2(2,2))*Log((1.*Abs(msq2(2,2)))/
      Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input))*Sqr(AtInput
      - MuInput/TanBeta))/((-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input))*((1.*Abs(msq2(2
      ,2)))/Sqr(M3Input) - Abs(msu2(2,2))/Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(M3Input)) + (1.*(AtInput - MuInput/TanBeta)*Cube(M3Input)*((1.
      *Abs(msq2(2,2))*(3 + 2*Log((1.*Sqr(SCALE))/Abs(msq2(2,2))) - Log((1.*Abs(
      msq2(2,2)))/Sqr(M3Input))*(-2 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)) + (1.*Abs(
      msq2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1 - (1.*Abs(msq2(2,2)))
      /Sqr(M3Input))*Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input))))/Sqr(Abs(msq2(2,2
      ))) + ((AtInput - MuInput/TanBeta)*(-8*ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input)) - (7.*Abs(msq2(2,2)))/Sqr(M3Input) + (7*Abs(msu2(2,2)))/Sqr(
      M3Input) + (4*Abs(msu2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)))/
      Sqr(M3Input) + (9.*Abs(msq2(2,2))*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)))/Sqr(
      M3Input) - (Abs(msu2(2,2))*Log((1.*Abs(msq2(2,2)))/Sqr(SCALE)))/Sqr(M3Input)
      + (1.*Abs(msq2(2,2))*Log(Abs(msu2(2,2))/Sqr(SCALE)))/Sqr(M3Input) - (9*Abs(
      msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(SCALE)))/Sqr(M3Input) + (4.*Abs(msq2(2,2))
      *Log(Sqr(SCALE)/Sqr(M3Input)))/Sqr(M3Input) - (4*Abs(msu2(2,2))*Log(Sqr(
      SCALE)/Sqr(M3Input)))/Sqr(M3Input) + Abs(msu2(2,2))/((1 - (1.*Abs(msq2(2,2))
      )/Sqr(M3Input))*Sqr(M3Input)) + (1.*Abs(msq2(2,2)))/((-1 + (1.*Abs(msq2(2,2)
      ))/Sqr(M3Input))*Sqr(M3Input)) + Abs(msu2(2,2))/((1 - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(M3Input)) + (1.*Abs(msq2(2,2)))/((-1 + Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(M3Input)) + (4*ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr
      (M3Input))/Abs(msu2(2,2)) + (Log((1.*Abs(msq2(2,2)))/Sqr(M3Input))*((-1.*Abs
      (msq2(2,2)))/Sqr(M3Input) + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(-1 + (1.*Abs(
      msq2(2,2)))/Sqr(M3Input)) - (1.*Abs(msq2(2,2))*Log(Abs(msu2(2,2))/Sqr(
      M3Input)))/(Sqr(M3Input)*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))) + (Abs(msu2(
      2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + Abs(msu2(2,2)
      )/Sqr(M3Input))) - (4.*ComplexLog(1 - (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(
      M3Input)*Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/Abs(msq2(2,2))))/(
      M3Input*((1.*Abs(msq2(2,2)))/Sqr(M3Input) - Abs(msu2(2,2))/Sqr(M3Input))) +
      ((AtInput - MuInput/TanBeta)*Cube(M3Input)*((Abs(msu2(2,2))*(3 + 2*Log(Sqr(
      SCALE)/Abs(msu2(2,2))) - Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(msu2(2,2
      ))/Sqr(M3Input)) + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + ComplexLog(1
       - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))))/Sqr(
      Abs(msu2(2,2)))))/(Cube((1.*Abs(msq2(2,2)))/Sqr(M3Input) - Abs(msu2(2,2))/
      Sqr(M3Input))*Power8(M3Input)) + ((AtInput - MuInput/TanBeta)*((2*(Log((1.*
      Abs(msq2(2,2)))/Abs(msu2(2,2)))*((1.*Abs(msq2(2,2)))/Sqr(M3Input) + Abs(msu2
      (2,2))/Sqr(M3Input)) - (2.*Abs(msq2(2,2)))/Sqr(M3Input) + (2*Abs(msu2(2,2)))
      /Sqr(M3Input))*((ComplexLog(1 - (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(Abs(
      msu2(2,2)))*Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/Quad(M3Input) - (1.*
      Abs(msq2(2,2))*((Abs(msu2(2,2))*((3.*Abs(msq2(2,2)))/Sqr(M3Input) - (3*Abs(
      msu2(2,2)))/Sqr(M3Input) - (2*Abs(msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(msq2(2,
      2))))/Sqr(M3Input) + (2.*Abs(msq2(2,2))*Log(Sqr(SCALE)/Abs(msu2(2,2))))/Sqr(
      M3Input) - (1.*Abs(msq2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(
      msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (Log((1.*Abs(msq2(2,
      2)))/Sqr(M3Input))*(-2 + (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2)
      )))/Quad(M3Input) + (1.*Abs(msq2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input)
      ))/Sqr((1.*Abs(msq2(2,2)))/Sqr(M3Input) - Abs(msu2(2,2))/Sqr(M3Input)) + (1.
      *(AtInput - MuInput/TanBeta)*Sqr(Abs(msq2(2,2)))*((2.*(AtInput - MuInput/
      TanBeta)*Cube(M3Input)*((-1.*Abs(msq2(2,2)))/Sqr(M3Input) + Abs(msu2(2,2))/
      Sqr(M3Input) + (1.*Abs(msq2(2,2))*Log((1.*Abs(msq2(2,2)))/Abs(msu2(2,2))))/
      Sqr(M3Input))*((ComplexLog(1 - (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(Abs(
      msu2(2,2)))*Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/Quad(M3Input) - (1.*
      Abs(msq2(2,2))*((Abs(msu2(2,2))*((3.*Abs(msq2(2,2)))/Sqr(M3Input) - (3*Abs(
      msu2(2,2)))/Sqr(M3Input) - (2*Abs(msu2(2,2))*Log((1.*Sqr(SCALE))/Abs(msq2(2,
      2))))/Sqr(M3Input) + (2.*Abs(msq2(2,2))*Log(Sqr(SCALE)/Abs(msu2(2,2))))/Sqr(
      M3Input) - (1.*Abs(msq2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-2 + Abs(
      msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) + (Log((1.*Abs(msq2(2,
      2)))/Sqr(M3Input))*(-2 + (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(Abs(msu2(2,2)
      )))/Quad(M3Input) + (1.*Abs(msq2(2,2))*ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input)
      ))/(Cube((1.*Abs(msq2(2,2)))/Sqr(M3Input) - Abs(msu2(2,2))/Sqr(M3Input))*Sqr
      (Abs(msq2(2,2)))) - ((Log((1.*Abs(msq2(2,2)))/Abs(msu2(2,2)))*((1.*Abs(msq2(
      2,2)))/Sqr(M3Input) + Abs(msu2(2,2))/Sqr(M3Input)) - (2.*Abs(msq2(2,2)))/Sqr
      (M3Input) + (2*Abs(msu2(2,2)))/Sqr(M3Input))*Sqr(M3Input)*((-1 + (1.*Abs(
      msq2(2,2)))/Abs(msu2(2,2)))*((-4.*Quad(M3Input)*((Abs(msu2(2,2))*ComplexLog(
      1 - (1.*Abs(msq2(2,2)))/Sqr(M3Input))*(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)
      ))/Sqr(M3Input) + (1.*Abs(msq2(2,2))*(ComplexLog(1 - Abs(msu2(2,2))/Sqr(
      M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)) - (2*Abs(msu2(2,2))*(2 + Log(
      Sqr(SCALE)/Sqr(M3Input))))/Sqr(M3Input)))/Sqr(M3Input)))/(Abs(msq2(2,2))*Abs
      (msu2(2,2))) - (4*(-((Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input))*(-1 +
      (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/Sqr(M3Input)) + (1.*Abs(msq2(2,2))*Log((
      1.*Abs(msq2(2,2)))/Sqr(M3Input))*(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(
      M3Input))*Sqr(AtInput - MuInput/TanBeta))/((-1 + (1.*Abs(msq2(2,2)))/Sqr(
      M3Input))*((1.*Abs(msq2(2,2)))/Sqr(M3Input) - Abs(msu2(2,2))/Sqr(M3Input))*(
      -1 + Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input)) - (1.*(AtInput - MuInput/
      TanBeta)*Cube(M3Input)*((1.*Abs(msq2(2,2))*(3 + 2*Log((1.*Sqr(SCALE))/Abs(
      msq2(2,2))) - Log((1.*Abs(msq2(2,2)))/Sqr(M3Input))*(-2 + (1.*Abs(msq2(2,2))
      )/Sqr(M3Input)) + (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/Sqr(M3Input) +
      ComplexLog(1 - (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(-1 + (1.*Abs(msq2(2,2))
      )/Sqr(M3Input))))/Sqr(Abs(msq2(2,2))) + ((AtInput - MuInput/TanBeta)*(-8*
      ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input)) - (7.*Abs(msq2(2,2)))/Sqr(
      M3Input) + (7*Abs(msu2(2,2)))/Sqr(M3Input) + (4*Abs(msu2(2,2))*ComplexLog(1
      - Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + (9.*Abs(msq2(2,2))*Log((1.*
      Abs(msq2(2,2)))/Sqr(SCALE)))/Sqr(M3Input) - (Abs(msu2(2,2))*Log((1.*Abs(msq2
      (2,2)))/Sqr(SCALE)))/Sqr(M3Input) + (1.*Abs(msq2(2,2))*Log(Abs(msu2(2,2))/
      Sqr(SCALE)))/Sqr(M3Input) - (9*Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(SCALE))
      )/Sqr(M3Input) + (4.*Abs(msq2(2,2))*Log(Sqr(SCALE)/Sqr(M3Input)))/Sqr(
      M3Input) - (4*Abs(msu2(2,2))*Log(Sqr(SCALE)/Sqr(M3Input)))/Sqr(M3Input) +
      Abs(msu2(2,2))/((1 - (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(M3Input)) + (1.*
      Abs(msq2(2,2)))/((-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input))*Sqr(M3Input)) + Abs
      (msu2(2,2))/((1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input)) + (1.*Abs(msq2(
      2,2)))/((-1 + Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input)) + (4*ComplexLog(1 -
      Abs(msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))/Abs(msu2(2,2)) + (Log((1.*Abs(
      msq2(2,2)))/Sqr(M3Input))*((-1.*Abs(msq2(2,2)))/Sqr(M3Input) + Abs(msu2(2,2)
      )/Sqr(M3Input)))/Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)) - (1.*Abs(msq2(2
      ,2))*Log(Abs(msu2(2,2))/Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + Abs(msu2(2,2))
      /Sqr(M3Input))) + (Abs(msu2(2,2))*Log(Abs(msu2(2,2))/Sqr(M3Input)))/(Sqr(
      M3Input)*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input))) - (4.*ComplexLog(1 - (1.*Abs
      (msq2(2,2)))/Sqr(M3Input))*Sqr(M3Input)*Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(
      M3Input)))/Abs(msq2(2,2))))/(M3Input*((1.*Abs(msq2(2,2)))/Sqr(M3Input) - Abs
      (msu2(2,2))/Sqr(M3Input))) + (3*(AtInput - MuInput/TanBeta)*Cube(M3Input)*((
      Abs(msu2(2,2))*(3 + 2*Log(Sqr(SCALE)/Abs(msu2(2,2))) - Log(Abs(msu2(2,2))/
      Sqr(M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)) + Abs(msu2(2,2))/Sqr(
      M3Input)))/Sqr(M3Input) + ComplexLog(1 - Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1
       + Abs(msu2(2,2))/Sqr(M3Input))))/Sqr(Abs(msu2(2,2)))) + (6.*(AtInput -
      MuInput/TanBeta)*Power7(M3Input)*((ComplexLog(1 - (1.*Abs(msq2(2,2)))/Sqr(
      M3Input))*Sqr(Abs(msu2(2,2)))*Sqr(-1 + (1.*Abs(msq2(2,2)))/Sqr(M3Input)))/
      Quad(M3Input) - (1.*Abs(msq2(2,2))*((Abs(msu2(2,2))*((3.*Abs(msq2(2,2)))/Sqr
      (M3Input) - (3*Abs(msu2(2,2)))/Sqr(M3Input) - (2*Abs(msu2(2,2))*Log((1.*Sqr(
      SCALE))/Abs(msq2(2,2))))/Sqr(M3Input) + (2.*Abs(msq2(2,2))*Log(Sqr(SCALE)/
      Abs(msu2(2,2))))/Sqr(M3Input) - (1.*Abs(msq2(2,2))*Log(Abs(msu2(2,2))/Sqr(
      M3Input))*(-2 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/Sqr(M3Input) +
      (Log((1.*Abs(msq2(2,2)))/Sqr(M3Input))*(-2 + (1.*Abs(msq2(2,2)))/Sqr(M3Input
      ))*Sqr(Abs(msu2(2,2))))/Quad(M3Input) + (1.*Abs(msq2(2,2))*ComplexLog(1 -
      Abs(msu2(2,2))/Sqr(M3Input))*Sqr(-1 + Abs(msu2(2,2))/Sqr(M3Input)))/Sqr(
      M3Input)))/Sqr(M3Input)))/(Abs(msq2(2,2))*Cube(Abs(msu2(2,2))))))/(Abs(msu2(
      2,2))*Quad(-1 + (1.*Abs(msq2(2,2)))/Abs(msu2(2,2))))))/Power5(M3Input)))/
      M3Input))/(Sqr(Abs(msq2(2,2)))*Sqr(Abs(msu2(2,2))))))/Quad(3.141592653589793
      )), 0) + IF(TwoLoopAtAt >= 1, (Power6(Yu(2,2))*(1 + Sqr(TanBeta))*((
      0.0009765625*FiniteLog(log(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))))*(-36
       - (36*Quad(MuInput))/(msq2(2,2)*msu2(2,2)) + (6*Quad(AtInput - MuInput/
      TanBeta))/(msq2(2,2)*msu2(2,2)) + (72*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)
      ) + (12*Quad(AtInput - MuInput/TanBeta)*Sqr(MuInput))/(msq2(2,2)*msu2(2,2)*
      Sqrt(msq2(2,2)*msu2(2,2))) - (36*Sqr(AtInput - MuInput/TanBeta))/Sqrt(msq2(2
      ,2)*msu2(2,2)) + (108*Quad(MuInput)*Sqr(AtInput - MuInput/TanBeta))/(msq2(2,
      2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) - (72*Sqr(MuInput)*Sqr(AtInput -
      MuInput/TanBeta))/(msq2(2,2)*msu2(2,2)) - (18*Quad(MuInput)*Quad(AtInput -
      MuInput/TanBeta))/(Sqr(msq2(2,2))*Sqr(msu2(2,2)))))/Quad(3.141592653589793)
      + (0.0009765625*(108 + (36*Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Quad(
      MuInput))/(msq2(2,2)*msu2(2,2)) + (9*Quad(AtInput - MuInput/TanBeta))/(msq2(
      2,2)*msu2(2,2)) - (108*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) - (30*Quad(
      AtInput - MuInput/TanBeta)*Sqr(MuInput))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)
      *msu2(2,2))) - (54*Sqr(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2)*msu2(2,2))
      - (108*Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Quad(MuInput)*Sqr(AtInput
       - MuInput/TanBeta))/(msq2(2,2)*msu2(2,2)*Sqrt(msq2(2,2)*msu2(2,2))) + (180*
      Sqr(MuInput)*Sqr(AtInput - MuInput/TanBeta))/(msq2(2,2)*msu2(2,2)) - (6*IF(
      Abs(-1 + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) < 0.00001, 0.5, (1 + (Log(
      Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2
      ,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))))/(1 - Sqr(MuInput)/Sqrt(
      msq2(2,2)*msu2(2,2))))*Sqr(MuInput)*Sqr(AtInput - MuInput/TanBeta)*(-6 + Sqr
      (AtInput - MuInput/TanBeta)/Sqrt(msq2(2,2)*msu2(2,2))))/(msq2(2,2)*msu2(2,2)
      ) + 3*Log(Sqr(SCALE)/Sqrt(msq2(2,2)*msu2(2,2)))*(-24*(-1 + Sqr(MuInput)/Sqrt
      (msq2(2,2)*msu2(2,2))) - (2*Quad(AtInput - MuInput/TanBeta)*(3 + (2*Sqr(
      MuInput))/Sqrt(msq2(2,2)*msu2(2,2))))/(msq2(2,2)*msu2(2,2)) + (12*(3 + (2*
      Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(AtInput - MuInput/TanBeta))/
      Sqrt(msq2(2,2)*msu2(2,2))) + (6*(2 + Log(Sqr(SCALE)/Sqrt(msq2(2,2)*msu2(2,2)
      )))*Sqr(AtInput - MuInput/TanBeta)*(30 + Quad(AtInput - MuInput/TanBeta)/(
      msq2(2,2)*msu2(2,2)) - (10*Sqr(AtInput - MuInput/TanBeta))/Sqrt(msq2(2,2)*
      msu2(2,2)))*Sqr(TanBeta))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 + Sqr(TanBeta))) + (
      18*Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Quad(MuInput)*Quad(AtInput -
      MuInput/TanBeta))/(Sqr(msq2(2,2))*Sqr(msu2(2,2))) + (36 - (
      4.4688152583787755*Cube(AtInput - MuInput/TanBeta)*((AtInput - MuInput/
      TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta))
      )/Sqrt(Sqrt(msq2(2,2)*msu2(2,2)))))/Power3(Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))
      + (26.812891550272653*(AtInput - MuInput/TanBeta)*((AtInput - MuInput/
      TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta))
      )/Sqrt(Sqrt(msq2(2,2)*msu2(2,2)))))/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (6*Sqr
      (AtInput - MuInput/TanBeta)*(-9 - 1.1172038145946939*Sqr((AtInput - MuInput/
      TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta))
      )/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))))/Sqrt(msq2(2,2)*msu2(2,2)) +
      6.70322288756816*Sqr((AtInput - MuInput/TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,
      2))) + (2*MuInput*Csc(2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2,2)*msu2(2,2)))) +
      (Quad(AtInput - MuInput/TanBeta)*(9 + 1.1172038145946939*Sqr((AtInput -
      MuInput/TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(
      TanBeta)))/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))))/(msq2(2,2)*msu2(2,2)) + 6*Log(
      Sqr(SCALE)/Sqrt(msq2(2,2)*msu2(2,2)))*((-4*Cube(AtInput - MuInput/TanBeta)*(
      (AtInput - MuInput/TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc
      (2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2,2)*msu2(2,2)))))/Power3(Sqrt(Sqrt(msq2
      (2,2)*msu2(2,2)))) + (24*(AtInput - MuInput/TanBeta)*((AtInput - MuInput/
      TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta))
      )/Sqrt(Sqrt(msq2(2,2)*msu2(2,2)))))/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + 6*(1 +
      Sqr((AtInput - MuInput/TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput
      *Csc(2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))) + (Quad(AtInput -
      MuInput/TanBeta)*(2 + Sqr((AtInput - MuInput/TanBeta)/Sqrt(Sqrt(msq2(2,2)*
      msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2,2)*msu2(2,
      2))))))/(msq2(2,2)*msu2(2,2)) - (6*Sqr(AtInput - MuInput/TanBeta)*(2 + Sqr((
      AtInput - MuInput/TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(
      2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))))/Sqrt(msq2(2,2)*msu2(2
      ,2))))/(1 + Sqr(TanBeta))))/Quad(3.141592653589793)))/Sqr(TanBeta), 0), 0) +
      IF(TwoLoopAbAs >= 1, WHICH(IsCloseRel(msu2(2,2),msq2(2,2),0.01) &&
      IsCloseRel(Sqrt(msu2(2,2)),M3Input,0.01), (0.00390625*Quad(Yd(2,2))*Sqr(g3)*
      (32*Log(msq2(2,2)) - 32*Log(Sqr(SCALE)) + 32*Log(msq2(2,2))*Log(Sqr(SCALE))
      + (5.333333333333333*Cube(AbInput - MuInput*TanBeta))/Power3(Sqrt(msq2(2,2))
      ) - (10.666666666666666*Cube(AbInput - MuInput*TanBeta)*Log(msq2(2,2)))/
      Power3(Sqrt(msq2(2,2))) + (10.666666666666666*Cube(AbInput - MuInput*TanBeta
      )*Log(Sqr(SCALE)))/Power3(Sqrt(msq2(2,2))) - (32*(AbInput - MuInput*TanBeta)
      )/Sqrt(msq2(2,2)) + (32*(AbInput - MuInput*TanBeta)*Log(msq2(2,2)))/Sqrt(
      msq2(2,2)) - (32*(AbInput - MuInput*TanBeta)*Log(Sqr(SCALE)))/Sqrt(msq2(2,2)
      ) + (16*Sqr(AbInput - MuInput*TanBeta))/msq2(2,2) - (32*Log(msq2(2,2))*Sqr(
      AbInput - MuInput*TanBeta))/msq2(2,2) + (32*Log(Sqr(SCALE))*Sqr(AbInput -
      MuInput*TanBeta))/msq2(2,2) - 16*Sqr(Log(Sqr(SCALE))) - 16*Sqr(Log(msq2(2,2)
      )) - (1.3333333333333333*Quad(AbInput - MuInput*TanBeta))/Sqr(msq2(2,2))))/(
      Quad(3.141592653589793)*Quad(1 + (0.0625*(1 + Sqr(TanBeta))*(0.25*Log(Sqr(
      MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1 +
      Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(
      TanBeta)))*Sqr(Yu(2,2)))/(Sqr(3.141592653589793)*Sqr(TanBeta)) + 0.5*((-
      0.03125*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2
      )/msd2(2,2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(3.141592653589793)) - (
      0.03125*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2
      )/msu2(2,2)))))/(Sqrt(Abs(msq2(2,2)*msu2(2,2)))*Sqr(3.141592653589793))) + (
      0.0625*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/Sqr(SCALE)) +
      (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)
      ) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2))/MuInput)))/
      Sqr(3.141592653589793) + (0.08333333333333333*Sqr(g3)*(1 + Log(Sqr(M3Input)/
      Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/
      M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(
      msd2(2,2))/M3Input))/M3Input))/Sqr(3.141592653589793) + (0.0625*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput
      - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/
      MuInput))/MuInput))/(Sqr(3.141592653589793)*Sqr(TanBeta)))), IsCloseRel(msq2
      (2,2),msu2(2,2),0.01), (-0.010416666666666666*Quad(Yd(2,2))*Sqr(g3)*(24*Cube
      (msq2(2,2))*PolyLog(2,1 - Sqr(M3Input)/msq2(2,2))*((-2*M3Input*(AbInput -
      MuInput*TanBeta))/Abs(M3Input) + Sqrt(Sqr(M3Input)))*Power3(Sqrt(Sqr(M3Input
      ))) + 6*Cube(msq2(2,2))*Sqr(Log(msq2(2,2)))*(3*Quad(M3Input) - 2*msq2(2,2)*
      Sqr(M3Input) - (4*M3Input*(AbInput - MuInput*TanBeta)*Power3(Sqrt(Sqr(
      M3Input))))/Abs(M3Input) + Sqr(msq2(2,2))) + Log(Sqr(M3Input))*Power3(Sqrt(
      Sqr(M3Input)))*(-12*Cube(msq2(2,2))*(Log(msq2(2,2)) - Log(Sqr(SCALE)))*((-2*
      M3Input*(AbInput - MuInput*TanBeta))/Abs(M3Input) + Sqrt(Sqr(M3Input))) + (4
      *M3Input*Cube(AbInput - MuInput*TanBeta)*msq2(2,2)*(-2*msq2(2,2) + Sqr(
      M3Input)))/Abs(M3Input) + Quad(AbInput - MuInput*TanBeta)*Sqrt(Sqr(M3Input))
      *(-3*msq2(2,2) + 2*Sqr(M3Input)) - 12*msq2(2,2)*Sqrt(Sqr(M3Input))*(-2*msq2(
      2,2) + Sqr(M3Input))*Sqr(AbInput - MuInput*TanBeta) - (24*M3Input*(AbInput -
      MuInput*TanBeta)*(msq2(2,2) + Sqr(M3Input))*Sqr(msq2(2,2)))/Abs(M3Input) + 6
      *Sqrt(Sqr(M3Input))*(msq2(2,2) + 2*Sqr(M3Input))*Sqr(msq2(2,2))) + Log(msq2(
      2,2))*msq2(2,2)*(12*msq2(2,2)*(-2*msq2(2,2) + Sqr(M3Input))*(-msq2(2,2) + 2*
      Sqr(M3Input))*Sqr(AbInput - MuInput*TanBeta) + Quad(AbInput - MuInput*
      TanBeta)*(-3*Quad(M3Input) + 6*msq2(2,2)*Sqr(M3Input) - 2*Sqr(msq2(2,2))) +
      (4*M3Input*Cube(AbInput - MuInput*TanBeta)*Sqrt(Sqr(M3Input))*Sqr(msq2(2,2))
      )/Abs(M3Input) + (48*M3Input*(AbInput - MuInput*TanBeta)*Power3(Sqrt(Sqr(
      M3Input)))*Sqr(msq2(2,2)))/Abs(M3Input) - 18*Sqr(msq2(2,2))*(2*Quad(M3Input)
      - 2*msq2(2,2)*Sqr(M3Input) + Sqr(msq2(2,2))) - 12*Log(Sqr(SCALE))*Sqr(msq2(2
      ,2))*(2*Quad(M3Input) - 2*msq2(2,2)*Sqr(M3Input) - (2*M3Input*(AbInput -
      MuInput*TanBeta)*Power3(Sqrt(Sqr(M3Input))))/Abs(M3Input) + Sqr(msq2(2,2))))
      - (-msq2(2,2) + Sqr(M3Input))*((4*M3Input*Cube(AbInput - MuInput*TanBeta)*
      msq2(2,2)*Sqrt(Sqr(M3Input))*(-2*msq2(2,2) + Sqr(M3Input)))/Abs(M3Input) + (
      24*M3Input*(AbInput - MuInput*TanBeta)*(msq2(2,2) - Sqr(M3Input))*Sqrt(Sqr(
      M3Input))*Sqr(msq2(2,2)))/Abs(M3Input) + Quad(AbInput - MuInput*TanBeta)*(2*
      Quad(M3Input) - 4*msq2(2,2)*Sqr(M3Input) + Sqr(msq2(2,2))) - 12*msq2(2,2)*
      Sqr(AbInput - MuInput*TanBeta)*(Quad(M3Input) - 3*msq2(2,2)*Sqr(M3Input) +
      Sqr(msq2(2,2))) + 3*Sqr(msq2(2,2))*(4*Quad(M3Input) - 9*msq2(2,2)*Sqr(
      M3Input) + 3*Sqr(msq2(2,2))) + 2*Log(Sqr(SCALE))*(3*Cube(msq2(2,2))*Log(Sqr(
      SCALE))*(msq2(2,2) - Sqr(M3Input)) + (2*M3Input*Cube(AbInput - MuInput*
      TanBeta)*msq2(2,2)*Sqrt(Sqr(M3Input))*(-msq2(2,2) + Sqr(M3Input)))/Abs(
      M3Input) - 6*msq2(2,2)*(-2*msq2(2,2) + Sqr(M3Input))*(-msq2(2,2) + Sqr(
      M3Input))*Sqr(AbInput - MuInput*TanBeta) + Quad(AbInput - MuInput*TanBeta)*
      Sqr(Sqr(M3Input) - msq2(2,2)) - (12*M3Input*(AbInput - MuInput*TanBeta)*
      Power3(Sqrt(Sqr(M3Input)))*Sqr(msq2(2,2)))/Abs(M3Input) + 3*Sqr(msq2(2,2))*(
      2*Quad(M3Input) - 3*msq2(2,2)*Sqr(M3Input) + 3*Sqr(msq2(2,2)))))))/(Cube(
      msq2(2,2))*Quad(3.141592653589793)*Quad(1 + (0.0625*(1 + Sqr(TanBeta))*(0.25
      *Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))
      )/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1
       + Sqr(TanBeta)))*Sqr(Yu(2,2)))/(Sqr(3.141592653589793)*Sqr(TanBeta)) + 0.5*
      ((-0.03125*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(
      2,2)/msd2(2,2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(3.141592653589793)) -
      (0.03125*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,
      2)/msu2(2,2)))))/(Sqrt(Abs(msq2(2,2)*msu2(2,2)))*Sqr(3.141592653589793))) +
      (0.0625*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/Sqr(SCALE)) +
      (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)
      ) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2))/MuInput)))/
      Sqr(3.141592653589793) + (0.08333333333333333*Sqr(g3)*(1 + Log(Sqr(M3Input)/
      Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/
      M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(
      msd2(2,2))/M3Input))/M3Input))/Sqr(3.141592653589793) + (0.0625*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput
      - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/
      MuInput))/MuInput))/(Sqr(3.141592653589793)*Sqr(TanBeta)))*Sqr(-msq2(2,2) +
      Sqr(M3Input))), True, (0.00390625*Quad(Yd(2,2))*Sqr(g3)*(-16*Sqr(Log(Sqr(
      SCALE))) + (AbInput - MuInput*TanBeta)*((64*M3Input*Log(0.985*msd2(2,2))*
      Sqrt(Sqr(M3Input)))/(Abs(M3Input)*(-0.985*msd2(2,2) + 1.02*msq2(2,2))) - (64
      *M3Input*Log(0.985*msd2(2,2))*Log(Sqr(SCALE))*Power3(Sqrt(Sqr(M3Input))))/(
      Abs(M3Input)*(0.985*msd2(2,2) - 1.02*msq2(2,2))*(-0.985*msd2(2,2) + Sqr(
      M3Input))) + (128*M3Input*PolyLog(2,1 - (1.015228426395939*Sqr(M3Input))/
      msd2(2,2))*Power3(Sqrt(Sqr(M3Input))))/(Abs(M3Input)*(0.985*msd2(2,2) - 1.02
      *msq2(2,2))*(-0.985*msd2(2,2) + Sqr(M3Input))) + (128*M3Input*PolyLog(2,1 -
      (0.9803921568627451*Sqr(M3Input))/msq2(2,2))*Power3(Sqrt(Sqr(M3Input))))/(
      Abs(M3Input)*(-0.985*msd2(2,2) + 1.02*msq2(2,2))*(-1.02*msq2(2,2) + Sqr(
      M3Input))) + Log(1.02*msq2(2,2))*((-64*M3Input*Sqrt(Sqr(M3Input)))/(Abs(
      M3Input)*(-0.985*msd2(2,2) + 1.02*msq2(2,2))) - (64*M3Input*Log(Sqr(SCALE))*
      Power3(Sqrt(Sqr(M3Input))))/(Abs(M3Input)*(-0.985*msd2(2,2) + 1.02*msq2(2,2)
      )*(-1.02*msq2(2,2) + Sqr(M3Input)))) + Log(Sqr(M3Input))*((-64*M3Input*Log(
      0.985*msd2(2,2))*Power3(Sqrt(Sqr(M3Input))))/(Abs(M3Input)*(0.985*msd2(2,2)
      - 1.02*msq2(2,2))*(-0.985*msd2(2,2) + Sqr(M3Input))) - (64*M3Input*Log(1.02*
      msq2(2,2))*Power3(Sqrt(Sqr(M3Input))))/(Abs(M3Input)*(-0.985*msd2(2,2) +
      1.02*msq2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))) + (64*M3Input*Log(Sqr(
      SCALE))*Power3(Sqrt(Sqr(M3Input))))/(Abs(M3Input)*(-0.985*msd2(2,2) + Sqr(
      M3Input))*(-1.02*msq2(2,2) + Sqr(M3Input)))) + (64*M3Input*Power3(Sqrt(Sqr(
      M3Input)))*Sqr(Log(0.985*msd2(2,2))))/(Abs(M3Input)*(0.985*msd2(2,2) - 1.02*
      msq2(2,2))*(-0.985*msd2(2,2) + Sqr(M3Input))) + (64*M3Input*Power3(Sqrt(Sqr(
      M3Input)))*Sqr(Log(1.02*msq2(2,2))))/(Abs(M3Input)*(-0.985*msd2(2,2) + 1.02*
      msq2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input)))) + (7.962575893301484*(2*(0.985
      *msd2(2,2) + 1.02*msq2(2,2))*Power6(M3Input) + 6.0282*msd2(2,2)*msq2(2,2)*(
      0.985*msd2(2,2) + 1.02*msq2(2,2))*Sqr(M3Input) - 3.02826627*Sqr(msd2(2,2))*
      Sqr(msq2(2,2)) - Quad(M3Input)*(9.0423*msd2(2,2)*msq2(2,2) + 1.94045*Sqr(
      msd2(2,2)) + 2.0808*Sqr(msq2(2,2)))))/(msd2(2,2)*msq2(2,2)*(-0.985*msd2(2,2)
      + Sqr(M3Input))*(-1.02*msq2(2,2) + Sqr(M3Input))) + Cube(AbInput - MuInput*
      TanBeta)*(Log(1.02*msq2(2,2))*((64*M3Input*Log(Sqr(SCALE))*(0.985*msd2(2,2)
      + 1.02*msq2(2,2))*Sqrt(Sqr(M3Input)))/(Abs(M3Input)*Cube(-0.985*msd2(2,2) +
      1.02*msq2(2,2))) + (64*M3Input*(0.985*msd2(2,2) + 3.06*msq2(2,2))*Sqrt(Sqr(
      M3Input)))/(Abs(M3Input)*Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2)))) - (64*
      M3Input*Log(0.985*msd2(2,2))*(2.955*msd2(2,2) + 1.02*msq2(2,2))*Sqrt(Sqr(
      M3Input)))/(Abs(M3Input)*Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2))) + (64*
      M3Input*PolyLog(2,1 - (1.015228426395939*Sqr(M3Input))/msd2(2,2))*Sqrt(Sqr(
      M3Input))*(-0.985*msd2(2,2) - 1.02*msq2(2,2) + 2*Sqr(M3Input)))/(Abs(M3Input
      )*Cube(0.985*msd2(2,2) - 1.02*msq2(2,2))) + (64*M3Input*PolyLog(2,1 - (
      0.9803921568627451*Sqr(M3Input))/msq2(2,2))*Sqrt(Sqr(M3Input))*(-0.985*msd2(
      2,2) - 1.02*msq2(2,2) + 2*Sqr(M3Input)))/(Abs(M3Input)*Cube(-0.985*msd2(2,2)
      + 1.02*msq2(2,2))) + Log(Sqr(M3Input))*((128*M3Input*Log(0.985*msd2(2,2))*
      Power3(Sqrt(Sqr(M3Input))))/(Abs(M3Input)*Cube(-0.985*msd2(2,2) + 1.02*msq2(
      2,2))) - (128*M3Input*Log(1.02*msq2(2,2))*Power3(Sqrt(Sqr(M3Input))))/(Abs(
      M3Input)*Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2)))) + (32*M3Input*Sqrt(Sqr(
      M3Input))*(-0.985*msd2(2,2) - 1.02*msq2(2,2) + 2*Sqr(M3Input))*Sqr(Log(0.985
      *msd2(2,2))))/(Abs(M3Input)*Cube(0.985*msd2(2,2) - 1.02*msq2(2,2))) + (32*
      M3Input*Sqrt(Sqr(M3Input))*(-0.985*msd2(2,2) - 1.02*msq2(2,2) + 2*Sqr(
      M3Input))*Sqr(Log(1.02*msq2(2,2))))/(Abs(M3Input)*Cube(-0.985*msd2(2,2) +
      1.02*msq2(2,2))) + Log(Sqr(SCALE))*((-64*M3Input*Log(0.985*msd2(2,2))*(0.985
      *msd2(2,2) + 1.02*msq2(2,2))*Sqrt(Sqr(M3Input)))/(Abs(M3Input)*Cube(-0.985*
      msd2(2,2) + 1.02*msq2(2,2))) - (128*M3Input*Sqrt(Sqr(M3Input)))/(Abs(M3Input
      )*Sqr(-0.985*msd2(2,2) + 1.02*msq2(2,2)))) - (256*M3Input*Sqrt(Sqr(M3Input))
      )/(Abs(M3Input)*Sqr(-0.985*msd2(2,2) + 1.02*msq2(2,2)))) + Sqr(AbInput -
      MuInput*TanBeta)*((-31.850303573205935*Sqr(M3Input))/(msd2(2,2)*msq2(2,2)) +
      (31.850303573205935*Log(Sqr(M3Input))*Quad(M3Input)*(-0.985*msd2(2,2) - 1.02
      *msq2(2,2) + Sqr(M3Input)))/(msd2(2,2)*msq2(2,2)*(-0.985*msd2(2,2) + Sqr(
      M3Input))*(-1.02*msq2(2,2) + Sqr(M3Input))) + (32*Log(0.985*msd2(2,2))*(-
      1.97*msd2(2,2) + 3*Sqr(M3Input)))/((0.985*msd2(2,2) - 1.02*msq2(2,2))*(-
      0.985*msd2(2,2) + Sqr(M3Input))) + Log(Sqr(SCALE))*((-64*Log(0.985*msd2(2,2)
      ))/(-0.985*msd2(2,2) + 1.02*msq2(2,2)) - (31.850303573205935*Sqr(M3Input))/(
      msd2(2,2)*msq2(2,2))) + Log(1.02*msq2(2,2))*((64*Log(Sqr(SCALE)))/(-0.985*
      msd2(2,2) + 1.02*msq2(2,2)) + (32*(-2.04*msq2(2,2) + 3*Sqr(M3Input)))/((-
      0.985*msd2(2,2) + 1.02*msq2(2,2))*(-1.02*msq2(2,2) + Sqr(M3Input))) + (32*
      Log(0.985*msd2(2,2))*(0.985*msd2(2,2) + 1.02*msq2(2,2)))/Sqr(-0.985*msd2(2,2
      ) + 1.02*msq2(2,2))) + (16*(-2.955*msd2(2,2) + 1.02*msq2(2,2))*Sqr(Log(0.985
      *msd2(2,2))))/Sqr(-0.985*msd2(2,2) + 1.02*msq2(2,2)) + (16*(0.985*msd2(2,2)
      - 3.06*msq2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Sqr(-0.985*msd2(2,2) + 1.02*msq2
      (2,2))) + Quad(AbInput - MuInput*TanBeta)*((16*PolyLog(2,1 - (
      1.015228426395939*Sqr(M3Input))/msd2(2,2))*(0.985*msd2(2,2) + 1.02*msq2(2,2)
      - 2*Sqr(M3Input)))/Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2)) + (16*PolyLog(2,1
       - (0.9803921568627451*Sqr(M3Input))/msq2(2,2))*(-0.985*msd2(2,2) - 1.02*
      msq2(2,2) + 2*Sqr(M3Input)))/Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2)) + (16*
      Log(0.985*msd2(2,2))*(6.895*msd2(2,2) + 3.06*msq2(2,2) + 2*Sqr(M3Input)))/
      Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2)) + (16*(1.97*msd2(2,2)*(0.985*msd2(2,
      2) + 1.02*msq2(2,2)) + (0.985*msd2(2,2) - 1.02*msq2(2,2))*Sqr(M3Input))*Sqr(
      Log(0.985*msd2(2,2))))/Quad(-0.985*msd2(2,2) + 1.02*msq2(2,2)) + (16*(2.04*
      msq2(2,2)*(0.985*msd2(2,2) + 1.02*msq2(2,2)) + (-0.985*msd2(2,2) + 1.02*msq2
      (2,2))*Sqr(M3Input))*Sqr(Log(1.02*msq2(2,2))))/Quad(-0.985*msd2(2,2) + 1.02*
      msq2(2,2)) + Log(Sqr(SCALE))*((32*Log(0.985*msd2(2,2))*(0.985*msd2(2,2) +
      1.02*msq2(2,2) + Sqr(M3Input)))/Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2)) + (
      15.925151786602967*(4.0188*msd2(2,2)*msq2(2,2) + (0.985*msd2(2,2) + 1.02*
      msq2(2,2))*Sqr(M3Input)))/(msd2(2,2)*msq2(2,2)*Sqr(-0.985*msd2(2,2) + 1.02*
      msq2(2,2)))) - (15.925151786602967*Log(Sqr(M3Input))*(0.985*msd2(2,2) + 1.02
      *msq2(2,2))*Sqr(M3Input))/(msd2(2,2)*msq2(2,2)*Sqr(-0.985*msd2(2,2) + 1.02*
      msq2(2,2))) + (15.925151786602967*(6.0282*msd2(2,2)*msq2(2,2) + (0.985*msd2(
      2,2) + 1.02*msq2(2,2))*Sqr(M3Input)))/(msd2(2,2)*msq2(2,2)*Sqr(-0.985*msd2(2
      ,2) + 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((-32*Log(Sqr(SCALE))*(0.985*
      msd2(2,2) + 1.02*msq2(2,2) + Sqr(M3Input)))/Cube(-0.985*msd2(2,2) + 1.02*
      msq2(2,2)) - (16*(2.955*msd2(2,2) + 7.140000000000001*msq2(2,2) + 2*Sqr(
      M3Input)))/Cube(-0.985*msd2(2,2) + 1.02*msq2(2,2)) - (32*Log(0.985*msd2(2,2)
      )*Sqr(0.985*msd2(2,2) + 1.02*msq2(2,2)))/Quad(-0.985*msd2(2,2) + 1.02*msq2(2
      ,2)))) + Log(Sqr(SCALE))*((15.925151786602967*((0.985*msd2(2,2) + 1.02*msq2(
      2,2))*Power6(M3Input) + 3.0141*msd2(2,2)*msq2(2,2)*(0.985*msd2(2,2) + 1.02*
      msq2(2,2))*Sqr(M3Input) - 3.02826627*Sqr(msd2(2,2))*Sqr(msq2(2,2)) - Quad(
      M3Input)*(3.0141*msd2(2,2)*msq2(2,2) + 0.970225*Sqr(msd2(2,2)) + 1.0404*Sqr(
      msq2(2,2)))))/(msd2(2,2)*msq2(2,2)*(-0.985*msd2(2,2) + Sqr(M3Input))*(-1.02*
      msq2(2,2) + Sqr(M3Input))) + (16*Log(0.985*msd2(2,2))*(2*Quad(M3Input) -
      1.97*msd2(2,2)*Sqr(M3Input) + 0.970225*Sqr(msd2(2,2))))/Sqr(-0.985*msd2(2,2)
      + Sqr(M3Input))) - (32*PolyLog(2,1 - (1.015228426395939*Sqr(M3Input))/msd2(2
      ,2))*Quad(M3Input))/Sqr(-0.985*msd2(2,2) + Sqr(M3Input)) + (24*Log(0.985*
      msd2(2,2))*(2*Quad(M3Input) - 1.97*msd2(2,2)*Sqr(M3Input) + 0.970225*Sqr(
      msd2(2,2))))/Sqr(-0.985*msd2(2,2) + Sqr(M3Input)) - (8*Sqr(Log(0.985*msd2(2,
      2)))*(3*Quad(M3Input) - 1.97*msd2(2,2)*Sqr(M3Input) + 0.970225*Sqr(msd2(2,2)
      )))/Sqr(-0.985*msd2(2,2) + Sqr(M3Input)) + Log(1.02*msq2(2,2))*((24*(2*Quad(
      M3Input) - 2.04*msq2(2,2)*Sqr(M3Input) + 1.0404*Sqr(msq2(2,2))))/Sqr(-1.02*
      msq2(2,2) + Sqr(M3Input)) + (16*Log(Sqr(SCALE))*(2*Quad(M3Input) - 2.04*msq2
      (2,2)*Sqr(M3Input) + 1.0404*Sqr(msq2(2,2))))/Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input))) + Log(Sqr(M3Input))*((16*Log(0.985*msd2(2,2))*Quad(M3Input))/Sqr(
      -0.985*msd2(2,2) + Sqr(M3Input)) + (16*Log(1.02*msq2(2,2))*Quad(M3Input))/
      Sqr(-1.02*msq2(2,2) + Sqr(M3Input)) - (16*Log(Sqr(SCALE))*Quad(M3Input)*(2*
      Quad(M3Input) - 2*(0.985*msd2(2,2) + 1.02*msq2(2,2))*Sqr(M3Input) + 0.970225
      *Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2))))/(Sqr(-0.985*msd2(2,2) + Sqr(
      M3Input))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input))) - (7.962575893301484*Quad(
      M3Input)*(2*(0.985*msd2(2,2) + 1.02*msq2(2,2))*Power6(M3Input) + Quad(
      M3Input)*(2.0094*msd2(2,2)*msq2(2,2) - 3.8809*Sqr(msd2(2,2)) - 4.1616*Sqr(
      msq2(2,2))) + 1.0047*msd2(2,2)*msq2(2,2)*(0.970225*Sqr(msd2(2,2)) + 1.0404*
      Sqr(msq2(2,2))) + 2*(0.985*msd2(2,2) + 1.02*msq2(2,2))*Sqr(M3Input)*Sqr(-
      0.985*msd2(2,2) + 1.02*msq2(2,2))))/(msd2(2,2)*msq2(2,2)*Sqr(-0.985*msd2(2,2
      ) + Sqr(M3Input))*Sqr(-1.02*msq2(2,2) + Sqr(M3Input)))) - (32*PolyLog(2,1 -
      (0.9803921568627451*Sqr(M3Input))/msq2(2,2))*Quad(M3Input))/Sqr(-1.02*msq2(2
      ,2) + Sqr(M3Input)) - (8*Sqr(Log(1.02*msq2(2,2)))*(3*Quad(M3Input) - 2.04*
      msq2(2,2)*Sqr(M3Input) + 1.0404*Sqr(msq2(2,2))))/Sqr(-1.02*msq2(2,2) + Sqr(
      M3Input))))/(Quad(3.141592653589793)*Quad(1 + (0.0625*(1 + Sqr(TanBeta))*(
      0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(
      SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/(Sqr(3.141592653589793)*Sqr(
      TanBeta)) + 0.5*((-0.03125*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5
      )(Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(
      3.141592653589793)) - (0.03125*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*
      TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/(Sqrt(Abs(msq2(2,2)*msu2(2,2)))*Sqr(
      3.141592653589793))) + (0.0625*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr
      (MuInput)/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(
      Sqrt(msq2(2,2))/MuInput)))/Sqr(3.141592653589793) + (0.08333333333333333*Sqr
      (g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) +
      TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(
      msq2(2,2))/M3Input,Sqrt(msd2(2,2))/M3Input))/M3Input))/Sqr(3.141592653589793
      ) + (0.0625*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/
      MuInput) + (0.5*(AtInput - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/
      MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/(Sqr(3.141592653589793)*Sqr(
      TanBeta))))), 0) + IF(TwoLoopAtAb >= 1, WHICH(IsCloseRel(Sqr(SCALE),msq2(2,2
      ),0.01) && IsCloseRel(Sqr(SCALE),msu2(2,2),0.01) && IsCloseRel(Sqr(SCALE),
      msd2(2,2),0.01) && IsCloseRel(SCALE,mAInput,0.01) && IsCloseRel(SCALE,Abs(
      MuInput),0.01), (0.00390625*((Power6(Yd(2,2))*(15 - Power6(AbInput - MuInput
      *TanBeta)/Power3(Sqrt(msq2(2,2)*msu2(2,2))) + (10.5*Quad(AbInput - MuInput*
      TanBeta))/(msq2(2,2)*msu2(2,2)) - (42*Sqr(AbInput - MuInput*TanBeta))/Sqrt(
      msq2(2,2)*msu2(2,2)) + (1 + Sqr(TanBeta))*(12 + (3*Quad(AbInput - MuInput*
      TanBeta))/(msq2(2,2)*msu2(2,2)) - (6*Sqr(AbInput - MuInput*TanBeta))/Sqrt(
      msq2(2,2)*msu2(2,2))) + Sqr(TanBeta)*(-10.049794796731925 + (
      2.6243712000000006*(AbInput + MuInput/TanBeta)*Cube(AbInput - MuInput*
      TanBeta))/(msq2(2,2)*msu2(2,2)) + (1.5025151999999977*(AbInput + MuInput/
      TanBeta)*(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + (
      0.37562879999999943*Sqr(AbInput + MuInput/TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)
      ) + Quad(AbInput - MuInput*TanBeta)*(-1.5/(msq2(2,2)*msu2(2,2)) - (
      0.06218560000000006*Sqr(AbInput + MuInput/TanBeta))/Power3(Sqrt(msq2(2,2)*
      msu2(2,2)))) + (8.063443199999998/Sqrt(msq2(2,2)*msu2(2,2)) - (
      0.06344319999999826*Sqr(AbInput + MuInput/TanBeta))/(msq2(2,2)*msu2(2,2)))*
      Sqr(AbInput - MuInput*TanBeta))))/Power6(1 + (0.0625*(1 + Sqr(TanBeta))*(
      0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(
      SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/(Sqr(3.141592653589793)*Sqr(
      TanBeta)) + 0.5*((-0.03125*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5
      )(Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(
      3.141592653589793)) - (0.03125*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*
      TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/(Sqrt(Abs(msq2(2,2)*msu2(2,2)))*Sqr(
      3.141592653589793))) + (0.0625*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr
      (MuInput)/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(
      Sqrt(msq2(2,2))/MuInput)))/Sqr(3.141592653589793) + (0.08333333333333333*Sqr
      (g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) +
      TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(
      msq2(2,2))/M3Input,Sqrt(msd2(2,2))/M3Input))/M3Input))/Sqr(3.141592653589793
      ) + (0.0625*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/
      MuInput) + (0.5*(AtInput - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/
      MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/(Sqr(3.141592653589793)*Sqr(
      TanBeta))) + (Quad(Yd(2,2))*(8.35153120435743 - (0.4585429333333333*(AbInput
       + MuInput/TanBeta)*Cube(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) +
      (6.751257599999999*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta))/
      Sqrt(msq2(2,2)*msu2(2,2)) + (AtInput + MuInput*TanBeta)*((-
      0.4585429333333333*Cube(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) +
      (6.751257599999999*(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2))) +
      (0.5*Quad(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) + ((-2*MuInput*
      Cube(AbInput - MuInput*TanBeta))/(Abs(MuInput)*Power3(Sqrt(Sqrt(msq2(2,2)*
      msu2(2,2))))) + (12*MuInput*(AtInput - MuInput/TanBeta))/(Abs(MuInput)*Sqrt(
      Sqrt(msq2(2,2)*msu2(2,2)))))/(Sqrt(1/(1 + Sqr(TanBeta)))*Sqrt(Sqr(TanBeta)/(
      1 + Sqr(TanBeta)))) - (11.248742400000001*Sqr(AbInput - MuInput*TanBeta))/
      Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AtInput - MuInput/TanBeta)*(12/Sqrt(msq2(2,2
      )*msu2(2,2)) + Quad(AbInput - MuInput*TanBeta)/Power3(Sqrt(msq2(2,2)*msu2(2,
      2))) - (9*Sqr(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2))) + (1 + Sqr(
      TanBeta))*(-24 - (3*Sqr(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)
      ) + Sqr(AtInput - MuInput/TanBeta)*(-3/Sqrt(msq2(2,2)*msu2(2,2)) + (1.5*Sqr(
      AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)))) + (AtInput - MuInput/
      TanBeta)*((-9.3756288*(AbInput + MuInput/TanBeta))/Sqrt(msq2(2,2)*msu2(2,2))
      - (9.3756288*(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + (2*(
      AbInput + MuInput/TanBeta)*Sqr(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2
      ,2)) + (AtInput + MuInput*TanBeta)*((4*(AbInput + MuInput/TanBeta)*(AbInput
      - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) - (0.819514311111111*(AbInput +
      MuInput/TanBeta)*Cube(AbInput - MuInput*TanBeta))/Power3(Sqrt(msq2(2,2)*msu2
      (2,2))) - 9.3756288/Sqrt(msq2(2,2)*msu2(2,2)) + (2*Sqr(AbInput - MuInput*
      TanBeta))/(msq2(2,2)*msu2(2,2)))) + Sqr(TanBeta)*(-14.063443199999998 - (
      9.3756288*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta))/Sqrt(msq2
      (2,2)*msu2(2,2)) - (4.6878144*Sqr(AbInput + MuInput/TanBeta))/Sqrt(msq2(2,2)
      *msu2(2,2)) + (1.6878143999999997/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AbInput +
      MuInput/TanBeta)/(msq2(2,2)*msu2(2,2)))*Sqr(AbInput - MuInput*TanBeta) + (
      AtInput - MuInput/TanBeta)*((-9.3756288*(AbInput + MuInput/TanBeta))/Sqrt(
      msq2(2,2)*msu2(2,2)) - (9.3756288*(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2
      )*msu2(2,2)) + (2*(AbInput + MuInput/TanBeta)*Sqr(AbInput - MuInput*TanBeta)
      )/(msq2(2,2)*msu2(2,2))) + Sqr(AtInput - MuInput/TanBeta)*((2*(AbInput +
      MuInput/TanBeta)*(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) +
      1.6878143999999997/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AbInput + MuInput/TanBeta
      )/(msq2(2,2)*msu2(2,2)) - (0.2707285333333336*Sqr(AbInput + MuInput/TanBeta)
      *Sqr(AbInput - MuInput*TanBeta))/Power3(Sqrt(msq2(2,2)*msu2(2,2))))) + ((1 +
      Sqr(TanBeta))*(12 + (0.5*Quad(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,
      2)) + (1.3378828010893575 + (AtInput + MuInput*TanBeta)*((-
      0.4585429333333333*Cube(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) +
      (6.751257599999999*(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2))) -
      (0.5*Quad(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) + (4.6878144*Sqr
      (AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + (1.6878143999999997
      /Sqrt(msq2(2,2)*msu2(2,2)) + (0.069514311111111*Quad(AbInput - MuInput*
      TanBeta))/Power3(Sqrt(msq2(2,2)*msu2(2,2))) - (0.6878143999999999*Sqr(
      AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)))*Sqr(AtInput + MuInput*
      TanBeta))/(1 + Sqr(TanBeta))))/Sqr(TanBeta))*Sqr(Yu(2,2)))/Quad(1 + (0.0625*
      (1 + Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(
      Sqr(mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(
      SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/(Sqr(
      3.141592653589793)*Sqr(TanBeta)) + 0.5*((-0.03125*Sqr(AbInput - MuInput*
      TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/(Sqrt(Abs(msd2
      (2,2)*msq2(2,2)))*Sqr(3.141592653589793)) - (0.03125*Sqr(AtInput - MuInput/
      TanBeta)*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/(Sqrt(Abs(msq2
      (2,2)*msu2(2,2)))*Sqr(3.141592653589793))) + (0.0625*(1 + Sqr(TanBeta))*Sqr(
      Yd(2,2))*(0.75*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput
      )/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))
      /MuInput) + TCF(6)(Sqrt(msq2(2,2))/MuInput)))/Sqr(3.141592653589793) + (
      0.08333333333333333*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)) + TCF(6)(Sqrt(
      msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((AbInput - MuInput*
      TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2))/M3Input))/M3Input))/
      Sqr(3.141592653589793) + (0.0625*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)
      (Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*TanBeta*TCF(9)(
      Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/(Sqr(
      3.141592653589793)*Sqr(TanBeta))) + (Quad(Yu(2,2))*(8.35153120435743 - (
      9.3756288*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta))/Sqrt(msq2
      (2,2)*msu2(2,2)) - (9.3756288*(AbInput - MuInput*TanBeta)*(AtInput + MuInput
      *TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + Cube(AtInput - MuInput/TanBeta)*((-
      0.4585429333333333*(AbInput + MuInput/TanBeta))/(msq2(2,2)*msu2(2,2)) + (
      AtInput + MuInput*TanBeta)*(-0.4585429333333333/(msq2(2,2)*msu2(2,2)) - (
      0.819514311111111*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta))/
      Power3(Sqrt(msq2(2,2)*msu2(2,2))))) + (AtInput - MuInput/TanBeta)*((
      6.751257599999999*(AbInput + MuInput/TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) - (
      9.3756288*(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + (AtInput
      + MuInput*TanBeta)*((4*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*
      TanBeta))/(msq2(2,2)*msu2(2,2)) + 6.751257599999999/Sqrt(msq2(2,2)*msu2(2,2)
      ))) + (1.3378828010893575 - (0.4585429333333333*(AbInput + MuInput/TanBeta)*
      Cube(AtInput - MuInput/TanBeta))/(msq2(2,2)*msu2(2,2)) + (6.751257599999999*
      (AtInput - MuInput/TanBeta)*(AbInput + MuInput/TanBeta))/Sqrt(msq2(2,2)*msu2
      (2,2)) + (1.6878143999999997*Sqr(AbInput + MuInput/TanBeta))/Sqrt(msq2(2,2)*
      msu2(2,2)) + Sqr(AtInput - MuInput/TanBeta)*(7.6878144/Sqrt(msq2(2,2)*msu2(2
      ,2)) - (0.6878143999999999*Sqr(AbInput + MuInput/TanBeta))/(msq2(2,2)*msu2(2
      ,2))) + Quad(AtInput - MuInput/TanBeta)*(-0.75/(msq2(2,2)*msu2(2,2)) + (
      0.069514311111111*Sqr(AbInput + MuInput/TanBeta))/Power3(Sqrt(msq2(2,2)*msu2
      (2,2)))))*Sqr(TanBeta) + ((-2*MuInput*Cube(AtInput - MuInput/TanBeta))/(Abs(
      MuInput)*Power3(Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))) + (12*MuInput*(AbInput -
      MuInput*TanBeta))/(Abs(MuInput)*Sqrt(Sqrt(msq2(2,2)*msu2(2,2)))) + (MuInput*
      (AbInput - MuInput*TanBeta)*Quad(AtInput - MuInput/TanBeta))/(Abs(MuInput)*
      Power5(Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))) - (12*MuInput*(AbInput - MuInput*
      TanBeta)*Sqr(AtInput - MuInput/TanBeta))/(Abs(MuInput)*Power3(Sqrt(Sqrt(msq2
      (2,2)*msu2(2,2))))))/(Sqrt(1/(1 + Sqr(TanBeta)))*Sqrt(Sqr(TanBeta)/(1 + Sqr(
      TanBeta)))) + (12 + (0.5*Quad(AtInput - MuInput/TanBeta))/(msq2(2,2)*msu2(2,
      2)))*(1 + Sqr(TanBeta)) + (12*Sqr(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)
      *msu2(2,2)) + Sqr(AtInput - MuInput/TanBeta)*((2*(AbInput + MuInput/TanBeta)
      *(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) + (2*(AbInput - MuInput*
      TanBeta)*(AtInput + MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)) +
      12.751257599999999/Sqrt(msq2(2,2)*msu2(2,2)) - (3*Sqr(AbInput - MuInput*
      TanBeta))/(msq2(2,2)*msu2(2,2))) + Quad(AtInput - MuInput/TanBeta)*(-1.5/(
      msq2(2,2)*msu2(2,2)) + (0.5*Sqr(AbInput - MuInput*TanBeta))/Power3(Sqrt(msq2
      (2,2)*msu2(2,2)))) + ((1 + Sqr(TanBeta))*(-24 - (3*Sqr(AbInput - MuInput*
      TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AtInput - MuInput/TanBeta)*(-3/
      Sqrt(msq2(2,2)*msu2(2,2)) + (1.5*Sqr(AbInput - MuInput*TanBeta))/(msq2(2,2)*
      msu2(2,2))) + (-14.063443199999998 - (9.3756288*(AbInput - MuInput*TanBeta)*
      (AtInput + MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + (1.6878143999999997
      *Sqr(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*msu2(2,2)) + (AtInput -
      MuInput/TanBeta)*((-9.3756288*(AbInput - MuInput*TanBeta))/Sqrt(msq2(2,2)*
      msu2(2,2)) + (AtInput + MuInput*TanBeta)*(-9.3756288/Sqrt(msq2(2,2)*msu2(2,2
      )) + (2*Sqr(AbInput - MuInput*TanBeta))/(msq2(2,2)*msu2(2,2)))) + (-
      4.6878144/Sqrt(msq2(2,2)*msu2(2,2)) + Sqr(AbInput - MuInput*TanBeta)/(msq2(2
      ,2)*msu2(2,2)))*Sqr(AtInput + MuInput*TanBeta) + Sqr(AtInput - MuInput/
      TanBeta)*((2*(AbInput - MuInput*TanBeta)*(AtInput + MuInput*TanBeta))/(msq2(
      2,2)*msu2(2,2)) + 1.6878143999999997/Sqrt(msq2(2,2)*msu2(2,2)) + (1/(msq2(2,
      2)*msu2(2,2)) - (0.2707285333333336*Sqr(AbInput - MuInput*TanBeta))/Power3(
      Sqrt(msq2(2,2)*msu2(2,2))))*Sqr(AtInput + MuInput*TanBeta)))/(1 + Sqr(
      TanBeta))))/Sqr(TanBeta))*Sqr(Yd(2,2)))/Sqr(1 + (0.006332573977646111*(1 +
      Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(
      SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/Sqr(TanBeta) + 0.5*
      ((-0.0031662869888230555*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*msq2(2,2))) - (
      0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))) +
      0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/
      Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1
      + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2
      ))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)
      ) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((
      AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2))/
      M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(
      0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*
      TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/
      Sqr(TanBeta))))/Quad(3.141592653589793), True, (0.00390625*((Power6(Yd(2,2))
      *(3*(5 + Log(1.02*msq2(2,2))*(-6 - 6*Log(Sqr(SCALE))) + Log(0.96*msd2(2,2))*
      (-4 - 4*Log(Sqr(SCALE))) + 10*Log(Sqr(SCALE)) + 5*Sqr(Log(Sqr(SCALE))) + 2*
      Sqr(Log(0.96*msd2(2,2))) + 3*Sqr(Log(1.02*msq2(2,2)))) + Power6(AbInput -
      MuInput*TanBeta)*(9*((-2*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Quad(0.96*msd2(2
      ,2) - 1.02*msq2(2,2)) - (1.9584*msd2(2,2)*msq2(2,2)*(0.96*msd2(2,2) + 1.02*
      msq2(2,2))*Sqr(Log(0.96*msd2(2,2))))/Power6(0.96*msd2(2,2) - 1.02*msq2(2,2))
      - (1.9584*msd2(2,2)*msq2(2,2)*(0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqr(Log(1.02
      *msq2(2,2))))/Power6(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (Log(1.02*msq2(2,2))
      *(5.8751999999999995*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) + 1.0404*
      Sqr(msq2(2,2))))/Power5(0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(0.96*msd2(2,2
      ))*((3.9168*Log(1.02*msq2(2,2))*msd2(2,2)*msq2(2,2)*(0.96*msd2(2,2) + 1.02*
      msq2(2,2)))/Power6(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (5.8751999999999995*
      msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2)))/Power5(
      0.96*msd2(2,2) - 1.02*msq2(2,2)))) + 3*((6*PolyLog(2,1 - (0.9411764705882352
      *msd2(2,2))/msq2(2,2)))/Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2)) + Log(0.96*
      msd2(2,2))*((-6.12*Log(Sqr(SCALE))*msq2(2,2))/Quad(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) + (0.9803921568627451*(0.96*msd2(2,2) - 5.1*msq2(2,2))*(0.96*msd2
      (2,2) + 2.04*msq2(2,2)))/(msq2(2,2)*Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) +
      (2*Log(1.02*msq2(2,2))*(5.76*msd2(2,2) + 5.1*msq2(2,2)))/Quad(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) - (2*(2.88*msd2(2,2) + 1.02*msq2(2,2))*Sqr(Log(0.96*msd2(
      2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2*(2.88*msd2(2,2) + 4.08*
      msq2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      (1.0212418300653596*Log(Sqr(SCALE))*(4.896*msd2(2,2)*msq2(2,2) - 0.9216*Sqr(
      msd2(2,2)) + 2.0808*Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      msd2(2,2)*msq2(2,2)) + (1.0212418300653596*(10.7712*msd2(2,2)*msq2(2,2) -
      0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*msd2(2,2)*msq2(2,2)) + Log(1.02*msq2(2,2))*((6.12*Log(Sqr(SCALE))
      *msq2(2,2))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.0416666666666667*(
      12.7296*msd2(2,2)*msq2(2,2) - 2.7648*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2)))
      )/(msd2(2,2)*Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)))))) + Quad(AbInput -
      MuInput*TanBeta)*(9*((2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2)
      ))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(0.96*msd2(2,2))*((-2*(0.96*
      msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (7.8336
      *Log(1.02*msq2(2,2))*msd2(2,2)*msq2(2,2))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,
      2))) + (3.9168*msd2(2,2)*msq2(2,2)*Sqr(Log(0.96*msd2(2,2))))/Quad(0.96*msd2(
      2,2) - 1.02*msq2(2,2)) + (3.9168*msd2(2,2)*msq2(2,2)*Sqr(Log(1.02*msq2(2,2))
      ))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) + 3*((-2*(0.96*msd2(2,2) + 3.06*
      msq2(2,2))*Sqr(Log(0.96*msd2(2,2))))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      (4*(2.88*msd2(2,2) + 1.02*msq2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Cube(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (6*PolyLog(2,1 - (0.9411764705882352*msd2(2,2)
      )/msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2.0424836601307192*(-
      26.438399999999998*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(
      msq2(2,2))))/(msd2(2,2)*msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      2.0424836601307192*Log(Sqr(SCALE))*(-14.687999999999999*msd2(2,2)*msq2(2,2)
      + 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2))))/(msd2(2,2)*msq2(2,2)*Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((-2*Log(Sqr(SCALE))
      *(6.72*msd2(2,2) + 5.1*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      4.166666666666667*(-9.792*msd2(2,2)*msq2(2,2) - 2.7648*Sqr(msd2(2,2)) +
      1.0404*Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msd2(2,2))) +
      Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2))*(-4.8*msd2(2,2) + 1.02*msq2(2,2)
      ))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(Sqr(SCALE))*(6.72*msd2(2,2
      ) + 5.1*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      0.9803921568627451*(43.0848*msd2(2,2)*msq2(2,2) - 1.8432*Sqr(msd2(2,2)) +
      6.2424*Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)))))
      + Sqr(AbInput - MuInput*TanBeta)*(9*((1.9584*msd2(2,2)*msq2(2,2)*Sqr(Log(
      0.96*msd2(2,2))))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (1.9584*msd2(2,2)*
      msq2(2,2)*Sqr(Log(1.02*msq2(2,2))))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      Log(0.96*msd2(2,2))*((-3.9168*Log(Sqr(SCALE))*msd2(2,2)*msq2(2,2))/Cube(0.96
      *msd2(2,2) - 1.02*msq2(2,2)) - (0.96*msd2(2,2) + 1.02*msq2(2,2))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((3.9168*Log(Sqr(SCALE))*
      msd2(2,2)*msq2(2,2))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (0.96*msd2(2,2)
      + 1.02*msq2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (2*Log(Sqr(SCALE))
      *(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) +
      3*(-2.0833333333333335/msd2(2,2) + Log(Sqr(SCALE))*(-2.0833333333333335/msd2
      (2,2) + 2/(0.96*msd2(2,2) - 1.02*msq2(2,2)) - 0.9803921568627451/msq2(2,2))
      + 4/(0.96*msd2(2,2) - 1.02*msq2(2,2)) - 0.9803921568627451/msq2(2,2) + (6*
      PolyLog(2,1 - (0.9411764705882352*msd2(2,2))/msq2(2,2)))/(-0.96*msd2(2,2) +
      1.02*msq2(2,2)) + (6*Sqr(Log(0.96*msd2(2,2))))/(0.96*msd2(2,2) - 1.02*msq2(2
      ,2)) - (2*(5.76*msd2(2,2) - 5.1*msq2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(1.02*msq2(2,2))*((2*Log(Sqr(SCALE))*(
      8.64*msd2(2,2) - 8.16*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      1.0416666666666667*(-16.6464*msd2(2,2)*msq2(2,2) + 17.5104*Sqr(msd2(2,2)) +
      2.0808*Sqr(msq2(2,2))))/(msd2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) +
      Log(0.96*msd2(2,2))*((-2*Log(Sqr(SCALE))*(8.64*msd2(2,2) - 8.16*msq2(2,2)))/
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(1.02*msq2(2,2))*(5.76*msd2(2,2)
      - 4.08*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      0.9803921568627451*(-20.563200000000002*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(
      msd2(2,2)) + 16.6464*Sqr(msq2(2,2))))/(msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2)))))) + (1 + Sqr(TanBeta))*(Sqr(AbInput - MuInput*TanBeta)*(9*(Log(
      1.02*msq2(2,2))*(4/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(
      0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*(4/(-0.96*msd2(2,2)
      + 1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) -
      (4*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (3.84*msd2(2,2)*Sqr(Log(0.96*msd2(2,2))))/Sqr(0.96*msd2
      (2,2) - 1.02*msq2(2,2)) + (4.08*msq2(2,2)*Sqr(Log(1.02*msq2(2,2))))/Sqr(0.96
      *msd2(2,2) - 1.02*msq2(2,2))) + 3*(4.12*Log(1.03*Sqr(MuInput))*Sqr(MuInput)*
      ((1.0416666666666667*(0.9803921568627451/msq2(2,2) + 1/(-0.96*msd2(2,2) +
      1.02*msq2(2,2))))/msd2(2,2) + 1/((0.96*msd2(2,2) - 1.03*Sqr(MuInput))*(-1.02
      *msq2(2,2) + 1.03*Sqr(MuInput)))) + (2.0424836601307192*Log(Sqr(SCALE))*(-
      1.92*msd2(2,2)*(1.02*msq2(2,2) + 1.03*Sqr(MuInput)) + 2.04*msq2(2,2)*(-1.02*
      msq2(2,2) + 2.06*Sqr(MuInput)) + 0.9216*Sqr(msd2(2,2))))/(msd2(2,2)*(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)) + Log(0.96*msd2(2,2))*(-
      1.9607843137254901/msq2(2,2) - (2*Log(1.02*msq2(2,2))*(1.92*msd2(2,2) + 3.06
      *msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(Sqr(SCALE))*(10.2*
      msq2(2,2) - 4*(0.96*msd2(2,2) + 1.03*Sqr(MuInput))))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) - (4*(0.96*msd2(2,2)*(-3.06*msq2(2,2) + 1.03*Sqr(MuInput)) -
      1.03*Sqr(MuInput)*(-2.04*msq2(2,2) + 1.03*Sqr(MuInput)) + 0.9216*Sqr(msd2(2,
      2))))/((0.96*msd2(2,2) - 1.03*Sqr(MuInput))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2
      ,2)))) + (4*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2))*(1.02
      *msq2(2,2) - 1.03*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*
      PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))*(-1.02*msq2(2,2)
      + 1.03*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(1.92*msd2(2
      ,2) - 1.02*msq2(2,2) + 1.03*Sqr(MuInput))*Sqr(Log(0.96*msd2(2,2))))/Sqr(0.96
      *msd2(2,2) - 1.02*msq2(2,2)) + ((8.16*msq2(2,2) - 2.06*Sqr(MuInput))*Sqr(Log
      (1.02*msq2(2,2))))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (-
      5.8751999999999995*msd2(2,2)*msq2(2,2) - 3.9552*msd2(2,2)*Sqr(MuInput) +
      8.4048*msq2(2,2)*Sqr(MuInput) + 1.8432*Sqr(msd2(2,2)) - 4.1616*Sqr(msq2(2,2)
      ))/(0.940032*msq2(2,2)*Sqr(msd2(2,2)) - 0.998784*msd2(2,2)*Sqr(msq2(2,2))) +
      Log(1.02*msq2(2,2))*((Log(Sqr(SCALE))*(-10.2*msq2(2,2) + 4*(0.96*msd2(2,2) +
      1.03*Sqr(MuInput))))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      1.0416666666666667*(-4.244832*Cube(msq2(2,2)) - 1.9776*msd2(2,2)*Sqr(MuInput
      )*(0.96*msd2(2,2) + 2.06*Sqr(MuInput)) + 1.9584*msd2(2,2)*msq2(2,2)*(2.88*
      msd2(2,2) + 5.15*Sqr(MuInput)) + 2.0808*(-4.8*msd2(2,2) + 2.06*Sqr(MuInput))
      *Sqr(msq2(2,2))))/(msd2(2,2)*(1.02*msq2(2,2) - 1.03*Sqr(MuInput))*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)))))) + Quad(AbInput - MuInput*TanBeta)*(3*((
      PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))*(3.9168*msd2(2,2)
      *msq2(2,2) + 6.365399999999999*Quad(MuInput) + 4.08*msq2(2,2)*(1.02*msq2(2,2
      ) - 3.09*Sqr(MuInput)) - 1.8432*Sqr(msd2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) + Log(1.03*Sqr(MuInput))*((6.365399999999999*Log(1.02*msq2(2,2))*
      Quad(MuInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2.103758169934641*
      Sqr(MuInput)*(-1.9584*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) - 2.0808*
      Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msd2(2,2)*msq2(2,2))
      ) + (Sqr(Log(1.02*msq2(2,2)))*(-7.8336*msd2(2,2)*msq2(2,2) -
      3.1826999999999996*Quad(MuInput) + 6.303599999999999*msq2(2,2)*Sqr(MuInput)
      + 0.9216*Sqr(msd2(2,2)) - 10.404*Sqr(msq2(2,2))))/Quad(0.96*msd2(2,2) - 1.02
      *msq2(2,2)) + (2*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2))*
      (-1.9584*msd2(2,2)*msq2(2,2) - 3.1826999999999996*Quad(MuInput) +
      6.303599999999999*msq2(2,2)*Sqr(MuInput) + 0.9216*Sqr(msd2(2,2)) - 2.0808*
      Sqr(msq2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Sqr(Log(0.96*msd2(
      2,2)))*(3.1826999999999996*Quad(MuInput) - 2.04*msq2(2,2)*(0.96*msd2(2,2) +
      3.09*Sqr(MuInput)) - 4.608*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2))))/Quad(
      0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.0212418300653596*(-0.884736*Cube(msd2(
      2,2)) + 0.9792*msd2(2,2)*msq2(2,2)*(27.54*msq2(2,2) - 16.48*Sqr(MuInput)) +
      1.8432*(-5.1*msq2(2,2) + 1.03*Sqr(MuInput))*Sqr(msd2(2,2)) + 2.0808*(1.02*
      msq2(2,2) - 2.06*Sqr(MuInput))*Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*msd2(2,2)*msq2(2,2)) + (1.0212418300653596*Log(Sqr(SCALE))*(-
      0.884736*Cube(msd2(2,2)) + 0.9792*msd2(2,2)*msq2(2,2)*(19.38*msq2(2,2) -
      10.3*Sqr(MuInput)) + 1.8432*(-4.08*msq2(2,2) + 1.03*Sqr(MuInput))*Sqr(msd2(2
      ,2)) + 2.0808*(1.02*msq2(2,2) - 2.06*Sqr(MuInput))*Sqr(msq2(2,2))))/(Cube(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*msd2(2,2)*msq2(2,2)) + Log(1.02*msq2(2,2))*
      ((-6*Log(Sqr(SCALE))*(-0.9792*msd2(2,2)*msq2(2,2) + 2.04*msq2(2,2)*(-1.02*
      msq2(2,2) + 1.03*Sqr(MuInput)) + 0.9216*Sqr(msd2(2,2))))/Quad(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) + (1.0416666666666667*(-6.1931519999999995*Cube(msd2(2,2))
      + 2.122416*Cube(msq2(2,2)) - 5.8751999999999995*msd2(2,2)*msq2(2,2)*(0.96*
      msd2(2,2) + 3.09*Sqr(MuInput)) + 28.964736*msd2(2,2)*Sqr(msq2(2,2))))/(msd2(
      2,2)*Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + Log(0.96*msd2(2,2))*((-
      6.365399999999999*Log(1.03*Sqr(MuInput))*Quad(MuInput))/Quad(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) + (6*Log(Sqr(SCALE))*(-0.9792*msd2(2,2)*msq2(2,2) + 2.04*
      msq2(2,2)*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput)) + 0.9216*Sqr(msd2(2,2))))/
      Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2))*(4.896*msd2(2
      ,2)*msq2(2,2) + 1.8432*Sqr(msd2(2,2)) + 4.1616*Sqr(msq2(2,2))))/Quad(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (0.9803921568627451*(0.884736*Cube(msd2(2,2))
      + 2.9375999999999998*msd2(2,2)*msq2(2,2)*(-7.140000000000001*msq2(2,2) +
      2.06*Sqr(MuInput)) + 13.160447999999999*msq2(2,2)*Sqr(msd2(2,2)) + 12.4848*(
      -1.02*msq2(2,2) + 1.03*Sqr(MuInput))*Sqr(msq2(2,2))))/(msq2(2,2)*Quad(0.96*
      msd2(2,2) - 1.02*msq2(2,2))))) + 9*(Log(1.02*msq2(2,2))*((-4*Log(Sqr(SCALE))
      *(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) -
      (4*(0.96*msd2(2,2) + 3.06*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)))
      - (3.84*msd2(2,2)*(0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqr(Log(0.96*msd2(2,2)))
      )/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (4.08*msq2(2,2)*(0.96*msd2(2,2) +
      1.02*msq2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,
      2)) - 8/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (8*Log(Sqr(SCALE)))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + Log(0.96*msd2(2,2))*((4*Log(Sqr(SCALE))*(0.96*
      msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*(
      2.88*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4
      *Log(1.02*msq2(2,2))*Sqr(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Quad(0.96*msd2(2,
      2) - 1.02*msq2(2,2))))) + 3*(-2.5 + 0.96*msd2(2,2)*(-0.9803921568627451/msq2
      (2,2) + 1/(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + 2.06*(0.9803921568627451/
      msq2(2,2) + 1/(1.02*msq2(2,2) - 1.03*Sqr(MuInput)))*Sqr(MuInput) + (
      1.0416666666666667*(-2.04*msq2(2,2) + 4.12*Sqr(MuInput)))/msd2(2,2) + Log(
      Sqr(SCALE))*(2 - (2.125*msq2(2,2))/msd2(2,2) - (0.9803921568627451*(0.96*
      msd2(2,2) - 2.06*Sqr(MuInput)))/msq2(2,2) + (4.291666666666667*Sqr(MuInput))
      /msd2(2,2) + (1.92*msd2(2,2))/(-0.96*msd2(2,2) + 1.03*Sqr(MuInput)) + (4.12*
      Sqr(MuInput))/(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + 4*Sqr(Log(Sqr(SCALE))
      ) + Log(0.96*msd2(2,2))*(2*Log(1.02*msq2(2,2)) + (0.9411764705882352*msd2(2,
      2))/msq2(2,2) + (2.1218*Log(1.03*Sqr(MuInput))*Quad(MuInput))/Sqr(0.96*msd2(
      2,2) - 1.03*Sqr(MuInput)) + (2.1218*Quad(MuInput) + 1.9776*msd2(2,2)*Sqr(
      MuInput) - 0.9216*Sqr(msd2(2,2)))/Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput)) -
      (2*Log(Sqr(SCALE))*(1.0609*Quad(MuInput) - 3.9552*msd2(2,2)*Sqr(MuInput) +
      1.8432*Sqr(msd2(2,2))))/Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (2*
      PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))*(-1.0609*Quad(
      MuInput) - 1.9776*msd2(2,2)*Sqr(MuInput) + 0.9216*Sqr(msd2(2,2))))/Sqr(0.96*
      msd2(2,2) - 1.03*Sqr(MuInput)) + (Sqr(Log(0.96*msd2(2,2)))*(-1.0609*Quad(
      MuInput) - 1.9776*msd2(2,2)*Sqr(MuInput) + 0.9216*Sqr(msd2(2,2))))/Sqr(0.96*
      msd2(2,2) - 1.03*Sqr(MuInput)) + PolyLog(2,1 - (1.0098039215686274*Sqr(
      MuInput))/msq2(2,2))*(2 - (8.4872*Quad(MuInput))/Sqr(-1.02*msq2(2,2) + 1.03*
      Sqr(MuInput))) + Sqr(Log(1.02*msq2(2,2)))*(1 - (4.2436*Quad(MuInput))/Sqr(-
      1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + Log(1.03*Sqr(MuInput))*(2.1218*Log(
      Sqr(SCALE))*Quad(MuInput)*(-(1/Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) - 2/
      Sqr(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + 1.03*Sqr(MuInput)*(-
      1.9607843137254901/msq2(2,2) + (1.0416666666666667*(-4.2436*Quad(MuInput) +
      2.9664*msd2(2,2)*Sqr(MuInput) - 1.8432*Sqr(msd2(2,2))))/(msd2(2,2)*Sqr(0.96*
      msd2(2,2) - 1.03*Sqr(MuInput))) - (6.18*Sqr(MuInput))/Sqr(-1.02*msq2(2,2) +
      1.03*Sqr(MuInput))) + (4.2436*Log(1.02*msq2(2,2))*Quad(MuInput))/Sqr(-1.02*
      msq2(2,2) + 1.03*Sqr(MuInput))) + Log(1.02*msq2(2,2))*((2.125*msq2(2,2))/
      msd2(2,2) + Log(Sqr(SCALE))*(-4 + (4.2436*Quad(MuInput))/Sqr(-1.02*msq2(2,2)
      + 1.03*Sqr(MuInput))) + (3.1826999999999996*Quad(MuInput) + 2.1012*msq2(2,2)
      *Sqr(MuInput) + 1.0404*Sqr(msq2(2,2)))/Sqr(-1.02*msq2(2,2) + 1.03*Sqr(
      MuInput))))) + Sqr(TanBeta)*(3*Sqr(AbInput + MuInput/TanBeta)*(-
      2.0833333333333335/msd2(2,2) + Log(Sqr(SCALE))*(-2.0833333333333335/msd2(2,2
      ) - 0.9803921568627451/msq2(2,2)) - 0.9803921568627451/msq2(2,2) + (
      0.5106209150326798*Log(0.96*msd2(2,2))*(0.9411919999999999*Power6(mAInput) +
      0.9603999999999999*(-4.8*msd2(2,2) + 1.02*msq2(2,2))*Quad(mAInput) + 4.9*Sqr
      (mAInput)*(-1.9584*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) - 1.0404*Sqr(
      msq2(2,2))) - (0.96*msd2(2,2) - 1.02*msq2(2,2))*(-5.8751999999999995*msd2(2,
      2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) + 3.1212*Sqr(msq2(2,2))) - (-2.88*msd2(
      2,2) + 3.06*msq2(2,2) + 0.98*Sqr(mAInput))*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2))))/(msd2(2,2)*msq2(2,2)*TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),0.96*msd2(2,2))) - (0.5425347222222222*Log(1.02*msq2(2,2))*(-
      1.8823839999999998*Power6(mAInput) + 0.9603999999999999*(6.72*msd2(2,2) +
      6.12*msq2(2,2))*Quad(mAInput) + (0.96*msd2(2,2) - 1.02*msq2(2,2))*(
      2.9375999999999998*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) - 2.0808*Sqr(
      msq2(2,2))) - 1.96*Sqr(mAInput)*(0.9792*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(
      msd2(2,2)) + 3.1212*Sqr(msq2(2,2))) + (-2.88*msd2(2,2) - 2.04*msq2(2,2) +
      1.96*Sqr(mAInput))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))
      /(Sqr(msd2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) +
      (0.5318967864923747*Log(0.98*Sqr(mAInput))*(-0.9411919999999999*(0.96*msd2(2
      ,2) + 2.04*msq2(2,2))*Power6(mAInput) + Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))
      *(-3.9168*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2)
      )) + 0.9603999999999999*Quad(mAInput)*(5.8751999999999995*msd2(2,2)*msq2(2,2
      ) + 4.608*Sqr(msd2(2,2)) + 6.2424*Sqr(msq2(2,2))) + 0.98*Sqr(mAInput)*(-
      4.42368*Cube(msd2(2,2)) - 6.367248*Cube(msq2(2,2)) + 7.520256*msq2(2,2)*Sqr(
      msd2(2,2)) + 2.996352*msd2(2,2)*Sqr(msq2(2,2))) + (0.96*msd2(2,2)*(-0.96*
      msd2(2,2) + 0.98*Sqr(mAInput)) + 2.04*msq2(2,2)*(1.92*msd2(2,2) + 0.98*Sqr(
      mAInput)) - 2.0808*Sqr(msq2(2,2)))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),
      0.96*msd2(2,2))))/(msq2(2,2)*Sqr(msd2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2))) + (0.5651403356481481*(-7.529535999999999*(0.96*
      msd2(2,2) + 1.02*msq2(2,2))*Power6(mAInput) + 1.8447363199999998*Power8(
      mAInput) - 1.96*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(mAInput)*(2.7648*Sqr(
      msd2(2,2)) - 4.1616*Sqr(msq2(2,2))) + Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*(
      -3.9168*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2)))
      + 0.9603999999999999*Quad(mAInput)*(7.8336*msd2(2,2)*msq2(2,2) + 6.4512*Sqr(
      msd2(2,2)) + 12.4848*Sqr(msq2(2,2))) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,
      2),0.96*msd2(2,2))*(11.750399999999999*msd2(2,2)*msq2(2,2) -
      5.7623999999999995*Quad(mAInput) + 11.76*(0.96*msd2(2,2) + 1.02*msq2(2,2))*
      Sqr(mAInput) - 4.608*Sqr(msd2(2,2)) - 6.2424*Sqr(msq2(2,2)) + 4*TDelta(0.98*
      Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))*TPhi(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))/(Cube(msd2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))) + 3*(1.5 + Log(1.02*msq2(2,2))*(-6 - 6*Log(Sqr(
      SCALE))) + Log(0.96*msd2(2,2))*(-4 - 4*Log(Sqr(SCALE))) + Sqr(
      3.141592653589793) - (2.041666666666667*Sqr(mAInput))/msd2(2,2) - (
      0.9607843137254901*Sqr(mAInput))/msq2(2,2) + Log(Sqr(SCALE))*(3 - 0.98*(
      2.0833333333333335/msd2(2,2) + 0.9803921568627451/msq2(2,2))*Sqr(mAInput)) +
      Log(0.98*Sqr(mAInput))*(7 - 6*Log(Sqr(SCALE)) + 0.98*(2.0833333333333335/
      msd2(2,2) + 0.9803921568627451/msq2(2,2))*Sqr(mAInput)) + 3*Sqr(Log(0.98*Sqr
      (mAInput))) + 8*Sqr(Log(Sqr(SCALE))) + 2*Sqr(Log(0.96*msd2(2,2))) + 3*Sqr(
      Log(1.02*msq2(2,2))) + (1.0850694444444444*(0.9603999999999999*Quad(mAInput)
      - 5.6448*msd2(2,2)*Sqr(mAInput) - TDelta(0.98*Sqr(mAInput),0.96*msd2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.96*msd2(2,2),0.96*msd2(2,2)))/Sqr(
      msd2(2,2)) + (0.9611687812379854*(1.9207999999999998*Quad(mAInput) - 10.9956
      *msq2(2,2)*Sqr(mAInput) - 2*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*
      msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/Sqr(msq2(
      2,2))) + 3*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta)*(Log(1.02
      *msq2(2,2))*(12/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (12*Log(Sqr(SCALE)))/(
      0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput))*((2*Log(0.96*msd2
      (2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2)))/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2)))/
      (0.96*msd2(2,2) - 1.02*msq2(2,2)) + 12/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) +
      (12*Log(Sqr(SCALE)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + (4*Sqr(Log(0.96*
      msd2(2,2))))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (6*Sqr(Log(1.02*msq2(2,2)))
      )/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) - (2.170138888888889*(-0.98*(-5.76*msd2
      (2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.98*Sqr(mAInput),0.96*msd2
      (2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.96*msd2(2,2),0.96*msd2(2,2))
      )/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))) + (2.170138888888889*(-
      2.9375999999999998*msd2(2,2)*msq2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(
      msd2(2,2)) + 1.0404*Sqr(msq2(2,2)) - TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2)
      ,0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))) + (0.9611687812379854*(-
      3.8415999999999997*Quad(mAInput) + 21.9912*msq2(2,2)*Sqr(mAInput) + 4*TDelta
      (0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(
      2,2)))) + 3*(AbInput + MuInput/TanBeta)*Cube(AbInput - MuInput*TanBeta)*(Log
      (1.02*msq2(2,2))*((-12*Log(Sqr(SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/
      Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (12*(0.96*msd2(2,2) + 3.06*msq2(2,2)
      ))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((12*Log(Sqr
      (SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(
      2,2)) + (12*(2.88*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) - (2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 3.06*msq2(2,2) - 1.96*
      Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput
      ))*((-2*Log(0.96*msd2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 5.88*Sqr(
      mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2))*(
      0.96*msd2(2,2) - 1.02*msq2(2,2) + 5.88*Sqr(mAInput)))/Cube(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + (4*(-0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput))
      *Sqr(Log(0.96*msd2(2,2))))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(2.88*
      msd2(2,2) + 5.1*msq2(2,2) - 3.92*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/
      Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - 48/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,
      2)) - (24*Log(Sqr(SCALE)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      2.170138888888889*(0.98*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(-5.76*msd2(2,2) +
      0.98*Sqr(mAInput))*Sqr(mAInput) + (-2.88*msd2(2,2) + 1.02*msq2(2,2))*TDelta(
      0.98*Sqr(mAInput),0.96*msd2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),
      0.96*msd2(2,2),0.96*msd2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      msd2(2,2))) + (2.170138888888889*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-
      2.9375999999999998*msd2(2,2)*msq2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(
      msd2(2,2)) + 1.0404*Sqr(msq2(2,2))) + (-2.88*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,
      2))*Sqr(msd2(2,2))) + (0.9611687812379854*(1.96*(0.96*msd2(2,2) - 1.02*msq2(
      2,2))*Sqr(mAInput)*(-11.22*msq2(2,2) + 1.96*Sqr(mAInput)) - 4*(0.96*msd2(2,2
      ) - 3.06*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))
      *TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(Cube(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Sqr(msq2(2,2)))) + Quad(AbInput - MuInput*TanBeta)*(3*(Log
      (0.96*msd2(2,2))*((-5.9976*msq2(2,2)*Sqr(mAInput))/Quad(0.96*msd2(2,2) -
      1.02*msq2(2,2)) - (5.9976*Log(Sqr(SCALE))*msq2(2,2)*Sqr(mAInput))/Quad(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((5.9976*msq2(2,2)*Sqr(
      mAInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (5.9976*Log(Sqr(SCALE))*
      msq2(2,2)*Sqr(mAInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.98*
      Sqr(mAInput))*((5.9976*Log(0.96*msd2(2,2))*msq2(2,2)*Sqr(mAInput))/Quad(0.96
      *msd2(2,2) - 1.02*msq2(2,2)) - (5.9976*Log(1.02*msq2(2,2))*msq2(2,2)*Sqr(
      mAInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.0008169934640525*Sqr(
      mAInput)*(-4.896*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) - 2.0808*Sqr(
      msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msd2(2,2)*msq2(2,2))) +
      (1.0008169934640525*Sqr(mAInput)*(4.896*msd2(2,2)*msq2(2,2) - 0.9216*Sqr(
      msd2(2,2)) + 2.0808*Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      msd2(2,2)*msq2(2,2)) + (1.0008169934640525*Log(Sqr(SCALE))*Sqr(mAInput)*(
      4.896*msd2(2,2)*msq2(2,2) - 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2))))/
      (Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msd2(2,2)*msq2(2,2))) + 3*Sqr(AbInput
       + MuInput/TanBeta)*((-2*(2.88*msd2(2,2) + 1.02*msq2(2,2) - 1.96*Sqr(mAInput
      ))*Sqr(Log(0.96*msd2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + ((-2.88
      *msd2(2,2) - 11.22*msq2(2,2) + 6.859999999999999*Sqr(mAInput))*Sqr(Log(1.02*
      msq2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.0212418300653596*Log
      (Sqr(SCALE))*(4.896*msd2(2,2)*msq2(2,2) - 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr
      (msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msd2(2,2)*msq2(2,2)) +
      (1.0212418300653596*(10.7712*msd2(2,2)*msq2(2,2) - 0.9216*Sqr(msd2(2,2)) +
      2.0808*Sqr(msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msd2(2,2)*
      msq2(2,2)) + Log(0.98*Sqr(mAInput))*((3*Log(0.96*msd2(2,2))*(0.96*msd2(2,2)
      - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))
      - (3*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput
      )))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (0.5318967864923747*(
      0.9411919999999999*(0.96*msd2(2,2) - 2.04*msq2(2,2))*Power6(mAInput) - Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*(3.9168*msd2(2,2)*msq2(2,2) + 2.7648*Sqr(
      msd2(2,2)) - 2.0808*Sqr(msq2(2,2))) + 0.9603999999999999*Quad(mAInput)*(
      1.9584*msd2(2,2)*msq2(2,2) - 2.7648*Sqr(msd2(2,2)) + 6.2424*Sqr(msq2(2,2)))
      + 0.98*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(mAInput)*(0.9792*msd2(2,2)*msq2
      (2,2) + 4.608*Sqr(msd2(2,2)) + 6.2424*Sqr(msq2(2,2))) + (3.9168*msd2(2,2)*
      msq2(2,2) - 0.9408*msd2(2,2)*Sqr(mAInput) + 1.9992*msq2(2,2)*Sqr(mAInput) +
      2.7648*Sqr(msd2(2,2)) - 2.0808*Sqr(msq2(2,2)))*TDelta(0.98*Sqr(mAInput),1.02
      *msq2(2,2),0.96*msd2(2,2))))/(msq2(2,2)*Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) +
      Log(1.02*msq2(2,2))*((6.12*Log(Sqr(SCALE))*msq2(2,2))/Quad(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + (0.5425347222222222*(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      (1.8823839999999998*Power6(mAInput) - 2.8811999999999998*(2.88*msd2(2,2) +
      2.04*msq2(2,2))*Quad(mAInput) - (0.96*msd2(2,2) - 1.02*msq2(2,2))*(0.9792*
      msd2(2,2)*msq2(2,2) + 10.137599999999999*Sqr(msd2(2,2)) - 2.0808*Sqr(msq2(2,
      2))) + 1.96*Sqr(mAInput)*(2.9375999999999998*msd2(2,2)*msq2(2,2) + 6.4512*
      Sqr(msd2(2,2)) + 3.1212*Sqr(msq2(2,2)))) + (4.42368*Cube(msd2(2,2)) +
      2.122416*Cube(msq2(2,2)) + 15.040512*msq2(2,2)*Sqr(msd2(2,2)) - 1.96*Sqr(
      mAInput)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + 0.998784*msd2(2,2)*Sqr(msq2(
      2,2)))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(Quad(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))) + Log(0.96*msd2(2,2))*((-6.12*Log(Sqr(SCALE))*
      msq2(2,2))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(1.02*msq2(2,2))*(
      8.64*msd2(2,2) + 13.26*msq2(2,2) - 10.78*Sqr(mAInput)))/Quad(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) + (0.5106209150326798*(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)
      )*(-0.9411919999999999*Power6(mAInput) + 0.9603999999999999*(2.88*msd2(2,2)
      + 7.140000000000001*msq2(2,2))*Quad(mAInput) + (0.96*msd2(2,2) - 1.02*msq2(2
      ,2))*(11.750399999999999*msd2(2,2)*msq2(2,2) + 2.7648*Sqr(msd2(2,2)) - 5.202
      *Sqr(msq2(2,2))) - 0.98*Sqr(mAInput)*(9.792*msd2(2,2)*msq2(2,2) + 4.608*Sqr(
      msd2(2,2)) + 11.4444*Sqr(msq2(2,2)))) - (0.884736*Cube(msd2(2,2)) +
      5.306039999999999*Cube(msq2(2,2)) + 14.10048*msq2(2,2)*Sqr(msd2(2,2)) - 0.98
      *Sqr(mAInput)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + 2.996352*msd2(2,2)*Sqr(
      msq2(2,2)))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(msd2(
      2,2)*msq2(2,2)*Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput
      ),1.02*msq2(2,2),0.96*msd2(2,2)))) + (1.0850694444444444*(0.98*(0.96*msd2(2,
      2) - 1.02*msq2(2,2))*(-5.76*msd2(2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput) + (-
      4.8*msd2(2,2) + 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.96*msd2(2,2),0.96
      *msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.96*msd2(2,2),0.96*msd2(2,2)))/(Quad(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))) + (0.5651403356481481*(Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*(-7.529535999999999*(0.96*msd2(2,2) + 1.02*
      msq2(2,2))*Power6(mAInput) + 1.8447363199999998*Power8(mAInput) - 1.96*(0.96
      *msd2(2,2) - 1.02*msq2(2,2))*Sqr(mAInput)*(2.7648*Sqr(msd2(2,2)) - 4.1616*
      Sqr(msq2(2,2))) - Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*(3.9168*msd2(2,2)*
      msq2(2,2) + 2.7648*Sqr(msd2(2,2)) - 2.0808*Sqr(msq2(2,2))) +
      0.9603999999999999*Quad(mAInput)*(7.8336*msd2(2,2)*msq2(2,2) +
      10.137599999999999*Sqr(msd2(2,2)) + 12.4848*Sqr(msq2(2,2)))) + 2*(-
      8.812800000000001*msd2(2,2)*msq2(2,2) + 16.5888*Sqr(msd2(2,2)) + 2.0808*Sqr(
      msq2(2,2)))*Sqr(TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) - (
      0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.9207999999999998*(7.68*msd2(2,2) - 3.06*
      msq2(2,2))*Quad(mAInput) - 3.92*Sqr(mAInput)*(4.896*msd2(2,2)*msq2(2,2) +
      8.2944*Sqr(msd2(2,2)) - 3.1212*Sqr(msq2(2,2))) + (0.96*msd2(2,2) - 1.02*msq2
      (2,2))*(-21.542399999999997*msd2(2,2)*msq2(2,2) + 24.8832*Sqr(msd2(2,2)) +
      6.2424*Sqr(msq2(2,2))))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,
      2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(Cube(msd2(2,2))
      *Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,
      2),0.96*msd2(2,2))) + (0.9611687812379854*(0.98*(-0.96*msd2(2,2) + 1.02*msq2
      (2,2))*Sqr(mAInput)*(-11.22*msq2(2,2) + 1.96*Sqr(mAInput)) + (1.92*msd2(2,2)
      - 9.18*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(Quad(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Sqr(msq2(2,2))))) + Sqr(AbInput - MuInput*TanBeta)*(3*Sqr(
      AbInput + MuInput/TanBeta)*((2*Sqr(Log(0.96*msd2(2,2))))/Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) + (3*Sqr(Log(1.02*msq2(2,2))))/Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) + (1.92*msd2(2,2) - 4.08*msq2(2,2))/(0.940032*msq2(2,2)*Sqr(msd2(
      2,2)) - 0.998784*msd2(2,2)*Sqr(msq2(2,2))) + (Log(Sqr(SCALE))*(1.92*msd2(2,2
      ) - 4.08*msq2(2,2)))/(0.940032*msq2(2,2)*Sqr(msd2(2,2)) - 0.998784*msd2(2,2)
      *Sqr(msq2(2,2))) + Log(0.98*Sqr(mAInput))*(-(Log(0.96*msd2(2,2))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))/Sqr(0.96*msd2(2,2) - 1.02
      *msq2(2,2)) - (1.0637935729847494*((0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*
      Sqr(mAInput))*(1.9592159999999998*msq2(2,2)*Quad(mAInput) + 0.98*Sqr(mAInput
      )*(-5.8751999999999995*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) - 4.1616*
      Sqr(msq2(2,2))) - (0.96*msd2(2,2) - 1.02*msq2(2,2))*(-3.9168*msd2(2,2)*msq2(
      2,2) + 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2)))) + (-3.9168*msd2(2,2)*
      msq2(2,2) + 2.04*msq2(2,2)*(1.02*msq2(2,2) - 0.98*Sqr(mAInput)) + 0.9216*Sqr
      (msd2(2,2)))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*Sqr(msd2(2,2))*TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + Log(0.96*msd2(2,2))*((-5*Log(
      1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(Sqr(SCALE)))/
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.0212418300653596*((0.96*msd2(2,2)
      - 1.02*msq2(2,2))*(0.9603999999999999*(0.96*msd2(2,2) + 4.08*msq2(2,2))*Quad
      (mAInput) - 1.9992*msq2(2,2)*(6.72*msd2(2,2) + 4.08*msq2(2,2))*Sqr(mAInput)
      - (0.96*msd2(2,2) - 1.02*msq2(2,2))*(-6.8544*msd2(2,2)*msq2(2,2) + 0.9216*
      Sqr(msd2(2,2)) + 4.1616*Sqr(msq2(2,2)))) - (0.9792*msd2(2,2)*msq2(2,2) +
      0.9216*Sqr(msd2(2,2)) - 4.1616*Sqr(msq2(2,2)))*TDelta(0.98*Sqr(mAInput),1.02
      *msq2(2,2),0.96*msd2(2,2))))/(msd2(2,2)*msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + Log(
      1.02*msq2(2,2))*((-2*Log(Sqr(SCALE)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      (1.0850694444444444*(-2*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(0.884736*Cube(
      msd2(2,2)) + 1.061208*Cube(msq2(2,2)) - 0.9411919999999999*Power6(mAInput) +
      0.9603999999999999*(3.84*msd2(2,2) + 3.06*msq2(2,2))*Quad(mAInput) -
      1.997568*msd2(2,2)*Sqr(msq2(2,2)) - 0.98*Sqr(mAInput)*(1.9584*msd2(2,2)*msq2
      (2,2) + 5.5296*Sqr(msd2(2,2)) + 3.1212*Sqr(msq2(2,2)))) - 2*(0.9792*msd2(2,2
      )*msq2(2,2) + 0.98*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(mAInput) - 0.9216*
      Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2)))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2
      ,2),0.96*msd2(2,2))))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + (
      1.0850694444444444*(0.98*(-5.76*msd2(2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput)
      - TDelta(0.98*Sqr(mAInput),0.96*msd2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.96*msd2(2,2),0.96*msd2(2,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (1.1302806712962963*((0.96*msd2(2,2) - 1.02*msq2(2,2))*
      (-7.529535999999999*(0.96*msd2(2,2) + 1.02*msq2(2,2))*Power6(mAInput) +
      1.8447363199999998*Power8(mAInput) - 1.96*(0.96*msd2(2,2) - 2.04*msq2(2,2))*
      (0.96*msd2(2,2) - 1.02*msq2(2,2))*(0.96*msd2(2,2) + 2.04*msq2(2,2))*Sqr(
      mAInput) + Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*(-3.9168*msd2(2,2)*msq2(2,2)
      + 0.9216*Sqr(msd2(2,2)) + 2.0808*Sqr(msq2(2,2))) + 0.9603999999999999*Quad(
      mAInput)*(7.8336*msd2(2,2)*msq2(2,2) + 10.137599999999999*Sqr(msd2(2,2)) +
      12.4848*Sqr(msq2(2,2)))) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2
      (2,2))*(0.9603999999999999*(-8.64*msd2(2,2) + 6.12*msq2(2,2))*Quad(mAInput)
      - 3*(2.88*msd2(2,2) - 2.04*msq2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      0.98*Sqr(mAInput)*(5.8751999999999995*msd2(2,2)*msq2(2,2) + 21.1968*Sqr(msd2
      (2,2)) - 12.4848*Sqr(msq2(2,2))) + (6.72*msd2(2,2) - 4.08*msq2(2,2))*TDelta(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))*TPhi(0.98*Sqr(mAInput),
      1.02*msq2(2,2),0.96*msd2(2,2)))/(Cube(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) + (
      0.9611687812379854*(0.98*Sqr(mAInput)*(-11.22*msq2(2,2) + 1.96*Sqr(mAInput))
      - 2*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2
      ))*Sqr(msq2(2,2)))) + 3*((2*Sqr(Log(0.96*msd2(2,2))))/(0.96*msd2(2,2) - 1.02
      *msq2(2,2)) + Log(0.98*Sqr(mAInput))*((-2.001633986928105*(0.96*msd2(2,2) -
      2.04*msq2(2,2))*Sqr(mAInput))/(msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      msq2(2,2)) - (Log(0.96*msd2(2,2))*(-0.96*msd2(2,2) + 1.02*msq2(2,2) + 0.98*
      Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(1.02*msq2(2,2))*(
      -0.96*msd2(2,2) + 1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((-2*(-2.88*msd2(2,2) + 1.02*msq2(2,2
      ) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2*Log(Sqr(
      SCALE))*(-2.88*msd2(2,2) + 2.04*msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((Log(1.02*msq2(2,2))*(
      0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + (2*(-3.84*msd2(2,2) + 2.04*msq2(2,2) + 0.98*Sqr(mAInput)))
      /Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(Sqr(SCALE))*(-2.88*msd2(2,2)
      + 2.04*msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))
      + ((-2.88*msd2(2,2) + 1.02*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(
      2,2))))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(Sqr(SCALE))*(1.9584*msd2
      (2,2)*msq2(2,2) + 1.8816*msd2(2,2)*Sqr(mAInput) - 3.9984*msq2(2,2)*Sqr(
      mAInput)))/(0.940032*msq2(2,2)*Sqr(msd2(2,2)) - 0.998784*msd2(2,2)*Sqr(msq2(
      2,2))) + (3.9168*msd2(2,2)*msq2(2,2) + 1.8816*msd2(2,2)*Sqr(mAInput) -
      3.9984*msq2(2,2)*Sqr(mAInput))/(0.940032*msq2(2,2)*Sqr(msd2(2,2)) - 0.998784
      *msd2(2,2)*Sqr(msq2(2,2))) + (1.0850694444444444*(0.98*(-5.76*msd2(2,2) +
      0.98*Sqr(mAInput))*Sqr(mAInput) - TDelta(0.98*Sqr(mAInput),0.96*msd2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.96*msd2(2,2),0.96*msd2(2,2)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))) + (1.0416666666666667*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + (0.9611687812379854*(0.98*(-0.96*msd2(2,2) + 1.02*msq2(2,
      2))*Sqr(mAInput)*(-11.22*msq2(2,2) + 1.96*Sqr(mAInput)) + (1.92*msd2(2,2) -
      3.06*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))*Sqr(msq2(2,2))))))))/Power6(1 + (0.0625*(1 + Sqr(TanBeta))*(
      0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(
      SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/(Sqr(3.141592653589793)*Sqr(
      TanBeta)) + 0.5*((-0.03125*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5
      )(Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(
      3.141592653589793)) - (0.03125*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*
      TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/(Sqrt(Abs(msq2(2,2)*msu2(2,2)))*Sqr(
      3.141592653589793))) + (0.0625*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr
      (MuInput)/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(
      Sqrt(msq2(2,2))/MuInput)))/Sqr(3.141592653589793) + (0.08333333333333333*Sqr
      (g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) +
      TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(
      msq2(2,2))/M3Input,Sqrt(msd2(2,2))/M3Input))/M3Input))/Sqr(3.141592653589793
      ) + (0.0625*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/
      MuInput) + (0.5*(AtInput - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/
      MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/(Sqr(3.141592653589793)*Sqr(
      TanBeta))) + (Quad(Yd(2,2))*Sqr(Yu(2,2))*((3*(AbInput - MuInput*TanBeta)*(
      Log(0.96*msd2(2,2))*((8.119113252073776*MuInput*Sqrt(Sqr(MuInput)))/(Abs(
      MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (8.119113252073776*MuInput*Log
      (Sqr(SCALE))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2
      ,2)))) + Log(1.02*msq2(2,2))*((8.119113252073776*MuInput*Sqrt(Sqr(MuInput)))
      /(Abs(MuInput)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + (8.119113252073776*
      MuInput*Log(Sqr(SCALE))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*(-0.96*msd2(2,2) +
      1.02*msq2(2,2)))) + (8.119113252073776*MuInput*PolyLog(2,1 - (
      1.0729166666666667*Sqr(MuInput))/msd2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(MuInput
      )*(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + (8.119113252073776*MuInput*PolyLog(2
      ,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(
      MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (4.059556626036888*MuInput*
      Sqrt(Sqr(MuInput))*Sqr(Log(0.96*msd2(2,2))))/(Abs(MuInput)*(-0.96*msd2(2,2)
      + 1.02*msq2(2,2))) + (4.059556626036888*MuInput*Sqrt(Sqr(MuInput))*Sqr(Log(
      1.02*msq2(2,2))))/(Abs(MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + 3*(
      AtInput - MuInput/TanBeta)*((8.281495517115252*MuInput*Log(1.02*msq2(2,2))*
      Log(Sqr(SCALE))*msq2(2,2)*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*(1.02*msq2(2,2)
      - 0.98*msu2(2,2))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + (
      8.119113252073776*MuInput*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/
      msq2(2,2))*Sqrt(Sqr(MuInput))*(1.02*msq2(2,2) + 1.03*Sqr(MuInput)))/(Abs(
      MuInput)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*(-1.02*msq2(2,2) + 1.03*Sqr(
      MuInput))) + (7.9567309870323*MuInput*Log(0.98*msu2(2,2))*Log(Sqr(SCALE))*
      msu2(2,2)*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*(-1.02*msq2(2,2) + 0.98*msu2(2,2
      ))*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput))) - (8.119113252073776*MuInput*
      PolyLog(2,1 - (1.0510204081632653*Sqr(MuInput))/msu2(2,2))*Sqrt(Sqr(MuInput)
      )*(0.98*msu2(2,2) + 1.03*Sqr(MuInput)))/(Abs(MuInput)*(-1.02*msq2(2,2) +
      0.98*msu2(2,2))*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput))) + Log(1.03*Sqr(
      MuInput))*((8.36268664963599*MuInput*Log(1.02*msq2(2,2))*Power3(Sqrt(Sqr(
      MuInput))))/(Abs(MuInput)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-1.02*msq2(2,2)
      + 1.03*Sqr(MuInput))) + (8.36268664963599*MuInput*Log(Sqr(SCALE))*Power3(
      Sqrt(Sqr(MuInput))))/(Abs(MuInput)*(0.98*msu2(2,2) - 1.03*Sqr(MuInput))*(-
      1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + (8.36268664963599*MuInput*Log(0.98*
      msu2(2,2))*Power3(Sqrt(Sqr(MuInput))))/(Abs(MuInput)*(-1.02*msq2(2,2) + 0.98
      *msu2(2,2))*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput)))) + (4.059556626036888*
      MuInput*Sqrt(Sqr(MuInput))*(1.02*msq2(2,2) + 1.03*Sqr(MuInput))*Sqr(Log(1.02
      *msq2(2,2))))/(Abs(MuInput)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*(-1.02*msq2(2
      ,2) + 1.03*Sqr(MuInput))) - (4.059556626036888*MuInput*Sqrt(Sqr(MuInput))*(
      0.98*msu2(2,2) + 1.03*Sqr(MuInput))*Sqr(Log(0.98*msu2(2,2))))/(Abs(MuInput)*
      (-1.02*msq2(2,2) + 0.98*msu2(2,2))*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput)))) +
      3*Cube(AbInput - MuInput*TanBeta)*(Log(1.02*msq2(2,2))*((8.119113252073776*
      MuInput*Log(Sqr(SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqrt(Sqr(MuInput))
      )/(Abs(MuInput)*Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (8.119113252073776*
      MuInput*(0.96*msd2(2,2) + 3.06*msq2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*
      Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + (8.119113252073776*MuInput*PolyLog
      (2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))*(0.96*msd2(2,2) + 1.02*
      msq2(2,2) - 2.06*Sqr(MuInput))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) - (8.119113252073776*MuInput*PolyLog(2,1 - (
      1.0098039215686274*Sqr(MuInput))/msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2)
      - 2.06*Sqr(MuInput))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + (16.72537329927198*MuInput*Log(1.02*msq2(2,2))*Log(1.03*
      Sqr(MuInput))*Power3(Sqrt(Sqr(MuInput))))/(Abs(MuInput)*Cube(-0.96*msd2(2,2)
      + 1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((-8.119113252073776*MuInput*Log(
      Sqr(SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(
      MuInput)*Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (8.119113252073776*MuInput
      *(2.88*msd2(2,2) + 1.02*msq2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(-
      0.96*msd2(2,2) + 1.02*msq2(2,2))) + (16.72537329927198*MuInput*Log(1.03*Sqr(
      MuInput))*Power3(Sqrt(Sqr(MuInput))))/(Abs(MuInput)*Cube(0.96*msd2(2,2) -
      1.02*msq2(2,2)))) + (4.059556626036888*MuInput*(0.96*msd2(2,2) + 1.02*msq2(2
      ,2) - 2.06*Sqr(MuInput))*Sqrt(Sqr(MuInput))*Sqr(Log(0.96*msd2(2,2))))/(Abs(
      MuInput)*Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (4.059556626036888*MuInput
      *(0.96*msd2(2,2) + 1.02*msq2(2,2) - 2.06*Sqr(MuInput))*Sqrt(Sqr(MuInput))*
      Sqr(Log(1.02*msq2(2,2))))/(Abs(MuInput)*Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)
      )) + (32.4764530082951*MuInput*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + (16.23822650414755*MuInput*Log(Sqr(SCALE))*
      Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))))/(
      Sqrt(1/(1 + Sqr(TanBeta)))*Sqrt(Sqr(TanBeta)/(1 + Sqr(TanBeta)))) + 3*Quad(
      AbInput - MuInput*TanBeta)*(Log(0.96*msd2(2,2))*((4*Log(1.02*msq2(2,2))*(
      0.96*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (4
      *Log(Sqr(SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) -
      1.02*msq2(2,2)) - (2*(6.72*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((4*Log(Sqr(SCALE))*(0.96*msd2(2,2)
      + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(2.88*msd2(2,2
      ) + 5.1*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (4*(0.96*msd2(2
      ,2) + 1.02*msq2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Cube(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) + 16/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (8*Log(Sqr(SCALE)))/
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (6*PolyLog(2,1 - (0.9411764705882352*
      msd2(2,2))/msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (1 + Sqr(
      TanBeta))*(3*(-4*PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))
      - 4*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2)) + 4.12*Log(
      1.03*Sqr(MuInput))*(1/(0.96*msd2(2,2) - 1.03*Sqr(MuInput)) + 1/(1.02*msq2(2,
      2) - 1.03*Sqr(MuInput)))*Sqr(MuInput) + Log(0.96*msd2(2,2))*(2*Log(1.02*msq2
      (2,2)) + 2*Log(Sqr(SCALE)) + (4.12*Sqr(MuInput))/(-0.96*msd2(2,2) + 1.03*Sqr
      (MuInput))) + Log(1.02*msq2(2,2))*(2*Log(Sqr(SCALE)) + (4.12*Sqr(MuInput))/(
      -1.02*msq2(2,2) + 1.03*Sqr(MuInput))) - 2*Sqr(Log(Sqr(SCALE))) - 2*Sqr(Log(
      0.96*msd2(2,2))) - 2*Sqr(Log(1.02*msq2(2,2)))) + 3*Sqr(AbInput - MuInput*
      TanBeta)*(4/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(Sqr(SCALE)))/(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (4.12*Log(1.03*Sqr(MuInput))*Sqr(MuInput))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + Log
      (1.02*msq2(2,2))*((2.04*Log(Sqr(SCALE))*msq2(2,2))/Sqr(0.96*msd2(2,2) - 1.02
      *msq2(2,2)) + (2*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((-2.04*Log(1.02*msq2(2,2))*msq2(2,2)
      )/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2.04*Log(Sqr(SCALE))*msq2(2,2))/
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (4*(-1.0506*msq2(2,2)*Sqr(MuInput) +
      0.9216*Sqr(msd2(2,2))))/((0.96*msd2(2,2) - 1.03*Sqr(MuInput))*Sqr(0.96*msd2(
      2,2) - 1.02*msq2(2,2)))) + (4.08*msq2(2,2)*PolyLog(2,1 - (1.0729166666666667
      *Sqr(MuInput))/msd2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (4.08*msq2
      (2,2)*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2)))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (2.04*msq2(2,2)*Sqr(Log(0.96*msd2(2,2))))/Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))) + Sqr(AtInput - MuInput/TanBeta)*(3*(4/(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4.12*Log(1.03*Sqr(MuInput))*Sqr(MuInput)
      )/((-1.02*msq2(2,2) + 0.98*msu2(2,2))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput)))
      + Log(Sqr(SCALE))*(2/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.96*Log(0.98*msu2
      (2,2))*msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2)
      )*(2/(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (1.96*Log(1.02*msq2(2,2))*msu2(2,2)
      )/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.96*Log(0.98*msu2(2,2))*msu2(2,2)
      )/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((-1.96*Log(
      Sqr(SCALE))*msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*(-1.0094*
      msu2(2,2)*Sqr(MuInput) + 1.0404*Sqr(msq2(2,2))))/((-1.02*msq2(2,2) + 1.03*
      Sqr(MuInput))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + (3.92*Log(0.98*msu2(2
      ,2))*msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (3.92*msu2(2,2)*
      PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2)))/Sqr(1.02*msq2(2,
      2) - 0.98*msu2(2,2)) - (3.92*msu2(2,2)*PolyLog(2,1 - (1.0510204081632653*Sqr
      (MuInput))/msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.96*msu2(2,2
      )*Sqr(Log(1.02*msq2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (1.96*
      msu2(2,2)*Sqr(Log(0.98*msu2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      3*Sqr(AbInput - MuInput*TanBeta)*(2/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02
      *msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((2.04*msq2(2,2))/((-
      1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      1.9992*Log(1.02*msq2(2,2))*msq2(2,2)*msu2(2,2))/(Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (1.9992*Log(0.98*msu2(2,2
      ))*msq2(2,2)*msu2(2,2))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(
      2,2) - 0.98*msu2(2,2)))) + Log(1.02*msq2(2,2))*((1.9992*Log(0.98*msu2(2,2))*
      msq2(2,2)*msu2(2,2))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2
      ) - 0.98*msu2(2,2))) + (2*(-0.9408*msd2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2
      ))))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,
      2)))) + (1.96*Log(0.98*msu2(2,2))*msu2(2,2))/((0.96*msd2(2,2) - 1.02*msq2(2,
      2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (1.9992*msq2(2,2)*msu2(2,2)*Sqr(
      Log(1.02*msq2(2,2))))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,
      2) - 0.98*msu2(2,2)))))) + Sqr(AtInput - MuInput/TanBeta)*(9*((1.9992*msq2(2
      ,2)*msu2(2,2)*Sqr(Log(1.02*msq2(2,2))))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)
      ) + Log(0.96*msd2(2,2))*((1.9992*Log(1.02*msq2(2,2))*msq2(2,2)*msu2(2,2))/
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (1.9992*Log(0.98*msu2(2,2))*msq2(2,2
      )*msu2(2,2))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (1.02*msq2(2,2) + 0.98*
      msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((-
      1.9992*Log(0.98*msu2(2,2))*msq2(2,2)*msu2(2,2))/Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) - (3.9984*Log(Sqr(SCALE))*msq2(2,2)*msu2(2,2))/Cube(1.02*msq2(2,2
      ) - 0.98*msu2(2,2)) - (1.02*msq2(2,2) + 0.98*msu2(2,2))/Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + Log(Sqr(SCALE))*((3.9984*Log(0.98*msu2(2,2))*msq2(2,2)*
      msu2(2,2))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*(1.02*msq2(2,2) + 0.98
      *msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + 3*((6*PolyLog(2,1 - (
      1.0408163265306123*msq2(2,2))/msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))
      + (3*Sqr(Log(0.98*msu2(2,2))))/(-1.02*msq2(2,2) + 0.98*msu2(2,2)) + (-5.1*
      msq2(2,2) + 0.98*msu2(2,2))/(-0.9996*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2
      ,2))) + Log(Sqr(SCALE))*((-3.06*msq2(2,2) + 0.98*msu2(2,2))/(-0.9996*msq2(2,
      2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2))) - (1.96*Log(0.98*msu2(2,2))*msu2(2,2))
      /Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((Log(0.98*msu2
      (2,2))*(6.12*msq2(2,2) - 3.92*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2
      )) + (3.06*msq2(2,2) + 0.98*msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))
      + (1.96*Log(Sqr(SCALE))*msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      0.9607843137254901*Log(0.98*msu2(2,2))*(-5.1*msq2(2,2) + 0.98*msu2(2,2))*
      msu2(2,2))/(msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + ((-3.06*msq2(2
      ,2) + 0.98*msu2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) + Sqr(AbInput - MuInput*TanBeta)*(9*((3.9984*msq2(2,2)*msu2(2,2)
      *Sqr(Log(1.02*msq2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.96*
      msd2(2,2) + 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((3.9984*Log(0.98*msu2(2,
      2))*msq2(2,2)*msu2(2,2))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2
      ,2) - 1.02*msq2(2,2))) + (2*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/((0.96*msd2(2
      ,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + Log(0.96*msd2
      (2,2))*((3.9984*Log(1.02*msq2(2,2))*msq2(2,2)*msu2(2,2))/(Cube(1.02*msq2(2,2
      ) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (3.9984*Log(0.98*
      msu2(2,2))*msq2(2,2)*msu2(2,2))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-
      0.96*msd2(2,2) + 1.02*msq2(2,2))) + (2*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/((
      -0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))))) +
      3*((1.9215686274509802*Log(0.98*msu2(2,2))*msu2(2,2))/((0.96*msd2(2,2) -
      1.02*msq2(2,2))*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + ((11.52*msd2(
      2,2) - 6*(1.02*msq2(2,2) + 0.98*msu2(2,2)))*PolyLog(2,1 - (
      1.0408163265306123*msq2(2,2))/msu2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))
      *Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + 2/(0.9792*msd2(2,2)*msq2(2,2) -
      1.0404*Sqr(msq2(2,2))) + (2*Log(Sqr(SCALE)))/(0.9792*msd2(2,2)*msq2(2,2) -
      1.0404*Sqr(msq2(2,2))) + (6*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 1.96*msu2(2,2
      ))*PolyLog(2,1 - (0.9411764705882352*msd2(2,2))/msq2(2,2)))/((0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (3*Sqr(Log(0.98*
      msu2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (Sqr(Log(1.02*msq2(2,2))
      )*(5.8751999999999995*msd2(2,2)*msq2(2,2) - 7.5264*msd2(2,2)*msu2(2,2) +
      3.9984*msq2(2,2)*msu2(2,2) + 2.7648*Sqr(msd2(2,2)) - 5.202*Sqr(msq2(2,2))))/
      (Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))
      - (6*PolyLog(2,1 - (0.9795918367346939*msd2(2,2))/msu2(2,2))*Sqr(0.96*msd2(2
      ,2) - 0.98*msu2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((-2*Log(Sqr(SCALE)))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (2*(0.96*msd2(2,2) - 0.98*msu2(2,2)))/((-1.02*
      msq2(2,2) + 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (Log(
      1.02*msq2(2,2))*(2.04*msq2(2,2)*(1.02*msq2(2,2) - 1.96*msu2(2,2)) + 7.5264*
      msd2(2,2)*msu2(2,2) - 5.5296*Sqr(msd2(2,2))))/(Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*Log(0.98*msu2(2,2))*(-
      3.7632*msd2(2,2)*msu2(2,2) + 2.7648*Sqr(msd2(2,2)) + 0.9603999999999999*Sqr(
      msu2(2,2))))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98
      *msu2(2,2)))) + Log(1.02*msq2(2,2))*(2/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))
      + (2*Log(Sqr(SCALE)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2*Log(0.98*
      msu2(2,2))*(5.8751999999999995*msd2(2,2)*msq2(2,2) - 3.7632*msd2(2,2)*msu2(2
      ,2) - 3.1212*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/(Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))))) + Quad(
      AbInput - MuInput*TanBeta)*(3*((0.9803921568627451*Log(Sqr(SCALE))*(0.96*
      msd2(2,2) + 5.1*msq2(2,2)))/(Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2))*msq2(2,2
      )) + (0.9803921568627451*(0.96*msd2(2,2) + 11.22*msq2(2,2)))/(Cube(-0.96*
      msd2(2,2) + 1.02*msq2(2,2))*msq2(2,2)) + (0.9607843137254901*Log(0.98*msu2(2
      ,2))*msu2(2,2))/(msq2(2,2)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(0.96*msd2(
      2,2) - 1.02*msq2(2,2))) + (6*PolyLog(2,1 - (0.9411764705882352*msd2(2,2))/
      msq2(2,2)))/((-1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))) + (6*PolyLog(2,1 - (0.9795918367346939*msd2(2,2))/msu2(2,2))*Sqr
      (0.96*msd2(2,2) - 0.98*msu2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Quad(
      0.96*msd2(2,2) - 1.02*msq2(2,2))) - (6*PolyLog(2,1 - (1.0408163265306123*
      msq2(2,2))/msu2(2,2))*Sqr(0.96*msd2(2,2) - 0.98*msu2(2,2)))/((1.02*msq2(2,2)
      - 0.98*msu2(2,2))*Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2
      ))*((-2*Log(Sqr(SCALE))*(1.92*msd2(2,2) + 1.02*msq2(2,2)))/Quad(0.96*msd2(2,
      2) - 1.02*msq2(2,2)) - (4*(0.96*msd2(2,2) + 2.04*msq2(2,2)))/Quad(0.96*msd2(
      2,2) - 1.02*msq2(2,2)) + 3/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(
      2,2) - 1.02*msq2(2,2))) - (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98*msu2
      (2,2))*(-2.9375999999999998*msd2(2,2)*msq2(2,2) + 1.8816*msd2(2,2)*msu2(2,2)
      + 0.9996*msq2(2,2)*msu2(2,2)))/(Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)))) + Log(0.96*msd2(2,2))*((2*Log(Sqr(SCALE))
      *(1.92*msd2(2,2) + 1.02*msq2(2,2)))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) -
      (2*(-6.8544*msd2(2,2)*msq2(2,2) + 4.704*msd2(2,2)*msu2(2,2) + 0.9996*msq2(2,
      2)*msu2(2,2) + 0.9216*Sqr(msd2(2,2))))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (Log(1.02*msq2(2,2))*(2.04*msq2(2,2
      )*(2.88*msd2(2,2) + 1.02*msq2(2,2)) - 3.92*(0.96*msd2(2,2) + 1.02*msq2(2,2))
      *msu2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98*msu2(2,2))*
      (-2.9375999999999998*msd2(2,2)*msq2(2,2) + 1.8816*msd2(2,2)*msu2(2,2) +
      0.9996*msq2(2,2)*msu2(2,2)))/(Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02
      *msq2(2,2) - 0.98*msu2(2,2)))) + ((-2.04*msq2(2,2)*(2.88*msd2(2,2) + 1.02*
      msq2(2,2)) + 3.92*(0.96*msd2(2,2) + 1.02*msq2(2,2))*msu2(2,2))*Sqr(Log(1.02*
      msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2)))) + 9*((1.9992*msq2(2,2)*(0.96*msd2(2,2) + 1.02*msq2(2,2))*
      msu2(2,2)*Sqr(Log(1.02*msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (3.9984*Log(0.98*msu2(2,2))*msq2(2,
      2)*msu2(2,2))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((1.9992*Log(1.02*msq2(2,2))*msq2(2,2
      )*(0.96*msd2(2,2) + 1.02*msq2(2,2))*msu2(2,2))/(Cube(-0.96*msd2(2,2) + 1.02*
      msq2(2,2))*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (1.9992*Log(0.98*msu2(2,
      2))*msq2(2,2)*(0.96*msd2(2,2) + 1.02*msq2(2,2))*msu2(2,2))/(Cube(0.96*msd2(2
      ,2) - 1.02*msq2(2,2))*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + ((0.96*msd2(2
      ,2) + 1.02*msq2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/(Cube(0.96*msd2(2,2
      ) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) - (2*(1.02*msq2(2
      ,2) + 0.98*msu2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((1.9992*Log(0.98*msu2(2,2))*
      msq2(2,2)*(0.96*msd2(2,2) + 1.02*msq2(2,2))*msu2(2,2))/(Cube(-0.96*msd2(2,2)
      + 1.02*msq2(2,2))*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (3.9984*(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*msu2(2,2) - 1.0404*(0.96*msd2(2,2) +
      1.02*msq2(2,2))*Sqr(msq2(2,2)) + 0.9603999999999999*(0.96*msd2(2,2) + 1.02*
      msq2(2,2))*Sqr(msu2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*Cube(1.02*
      msq2(2,2) - 0.98*msu2(2,2))))))) + 3*(AbInput + MuInput/TanBeta)*(AbInput -
      MuInput*TanBeta)*(Log(1.02*msq2(2,2))*(4/(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      (4*Log(Sqr(SCALE)))/(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.98*Sqr(
      mAInput))*((2*Log(0.96*msd2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*
      Log(1.02*msq2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + Log(0.96*msd2(2,2
      ))*((2*Log(1.02*msq2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + 4/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(-0.96*msd2(2,2) + 1.02*
      msq2(2,2))) + (2*Sqr(Log(1.02*msq2(2,2))))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)
      ) + (2.170138888888889*(-2.9375999999999998*msd2(2,2)*msq2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*
      msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2)) -
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*
      Sqr(msd2(2,2))) + (1.9223375624759709*(-0.98*(-5.1*msq2(2,2) + 0.98*Sqr(
      mAInput))*Sqr(mAInput) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2
      ,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2
      ) - 1.02*msq2(2,2))*Sqr(msq2(2,2)))) + 3*(AbInput + MuInput/TanBeta)*Cube(
      AbInput - MuInput*TanBeta)*(Log(1.02*msq2(2,2))*((-4*Log(Sqr(SCALE))*(0.96*
      msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (4*(
      0.96*msd2(2,2) + 3.06*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) +
      Log(0.96*msd2(2,2))*((4*Log(Sqr(SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/
      Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*(2.88*msd2(2,2) + 1.02*msq2(2,2))
      )/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2*Log(1.02*msq2(2,2))*(0.96*msd2(
      2,2) + 3.06*msq2(2,2) - 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(
      2,2))) + Log(0.98*Sqr(mAInput))*((-2*Log(0.96*msd2(2,2))*(0.96*msd2(2,2) -
      1.02*msq2(2,2) + 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      (2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 1.96*Sqr(mAInput))
      )/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (2*(0.96*msd2(2,2) + 3.06*msq2(2,
      2) - 1.96*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/Cube(0.96*msd2(2,2) - 1.02
      *msq2(2,2)) - 16/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (8*Log(Sqr(SCALE)))/
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2.170138888888889*((0.96*msd2(2,2) -
      1.02*msq2(2,2))*(-2.9375999999999998*msd2(2,2)*msq2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*
      msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2))) + (-
      2.88*msd2(2,2) + 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(
      Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))) + (1.9223375624759709*
      (0.98*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))
      *Sqr(mAInput) - (0.96*msd2(2,2) - 3.06*msq2(2,2))*TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*
      msq2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2)))) + 3*(-1
      + Log(0.98*Sqr(mAInput))*(2*Log(0.96*msd2(2,2)) - 2*Log(0.98*msu2(2,2)) - 8*
      Log(Sqr(SCALE))) + Log(0.96*msd2(2,2))*(-4 + 2*Log(0.98*msu2(2,2)) - 4*Log(
      Sqr(SCALE))) + Log(1.02*msq2(2,2))*(-2 - 2*Log(Sqr(SCALE))) + 6*Log(Sqr(
      SCALE)) + 1.3333333333333333*Sqr(3.141592653589793) + 4*Sqr(Log(0.98*Sqr(
      mAInput))) + 7*Sqr(Log(Sqr(SCALE))) + Sqr(Log(1.02*msq2(2,2))) - (
      1.9223375624759709*(-0.98*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput)
      + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/Sqr(msq2(2,2)) + (2.170138888888889
      *(-2.8224*msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*
      msd2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*
      Sqr(msd2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput
      ),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96
      *msd2(2,2)))/Sqr(msd2(2,2))) + 3*Sqr(AbInput - MuInput*TanBeta)*(8/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(-0.96*msd2(2,2) + 1.02*
      msq2(2,2)) + (12*PolyLog(2,1 - (0.9411764705882352*msd2(2,2))/msq2(2,2)))/(
      0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(0.98*Sqr(mAInput))*((2*Log(0.96*msd2(
      2,2))*(0.96*msd2(2,2) - 2.04*msq2(2,2) + 0.98*msu2(2,2)))/Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) - (2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) - 2.04*msq2(2,2)
      + 0.98*msu2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2
      ))*((-3.84*Log(Sqr(SCALE))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) -
      (2*(0.96*msd2(2,2) + 3.06*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      (2*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput))
      )/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((3.84*Log(Sqr
      (SCALE))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(2.88*msd2(2,2
      ) + 1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2*Log(0.98*msu2
      (2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(
      2,2) - 1.02*msq2(2,2)) - (2*Log(1.02*msq2(2,2))*(2.88*msd2(2,2) - 3.06*msq2(
      2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (2*(2.88*
      msd2(2,2) - 3.06*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2.0833333333333335*TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2
      ),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      1.9223375624759709*(0.98*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-5.1*msq2(2,2)
      + 0.98*Sqr(mAInput))*Sqr(mAInput) + (0.96*msd2(2,2) - 2.04*msq2(2,2))*TDelta
      (0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      msq2(2,2))) + (2.170138888888889*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.8224
      *msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*
      Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2
      ,2)) + 0.9603999999999999*Sqr(msu2(2,2))) + (-1.92*msd2(2,2) + 1.02*msq2(2,2
      ))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (1.9607843137254901*TDelta(0.98*Sqr(mAInput),0.98*msu2(
      2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/
      (msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + (AtInput + MuInput*
      TanBeta)*(3*(AbInput - MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(0.96
      *msd2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2)))/(-
      0.96*msd2(2,2) + 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*(4/(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2
      *Log(0.98*msu2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + Log(0.96*msd2(2,
      2))*((2*Log(0.98*msu2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + 4/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(-0.96*msd2(2,2) + 1.02*
      msq2(2,2))) + (2.170138888888889*(-2.8224*msd2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))) - (1.9223375624759709*(-
      2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(
      2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(
      msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2)))) + 3*Cube(
      AbInput - MuInput*TanBeta)*(Log(0.96*msd2(2,2))*((4*Log(Sqr(SCALE))*(0.96*
      msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*(
      2.88*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2
      *Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) + 1.96*msu2(2,2) -
      1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2
      ,2))*((-4*Log(Sqr(SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(
      2,2) - 1.02*msq2(2,2)) - (4*(0.96*msd2(2,2) + 3.06*msq2(2,2)))/Cube(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*
      msq2(2,2) + 1.96*msu2(2,2) - 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*
      msq2(2,2))) + Log(0.98*Sqr(mAInput))*((-2*Log(0.96*msd2(2,2))*(0.96*msd2(2,2
      ) + 1.02*msq2(2,2) - 1.96*msu2(2,2) + 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2
      ) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2
      ) - 1.96*msu2(2,2) + 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2
      ))) - 16/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (8*Log(Sqr(SCALE)))/Sqr(0.96
      *msd2(2,2) - 1.02*msq2(2,2)) + (2.170138888888889*((0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(-2.8224*msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) +
      1.8432*Sqr(msd2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) + (-2.88*msd2(2,2
      ) + 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))
      *TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Cube(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Sqr(msd2(2,2))) + (1.9223375624759709*((0.96*msd2(2,2) -
      1.02*msq2(2,2))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(
      mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(
      mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) - (
      0.96*msd2(2,2) - 3.06*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(
      Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2))))) + Sqr(TanBeta)*(3*
      Sqr(AbInput + MuInput/TanBeta)*((-2*Log(0.96*msd2(2,2))*(0.96*msd2(2,2) -
      1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2)
      ,0.96*msd2(2,2)) + (1.0416666666666667*Log(0.98*Sqr(mAInput))*(Sqr(0.98*Sqr(
      mAInput) + 0.96*msd2(2,2) - 1.02*msq2(2,2)) - TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2))))/(msd2(2,2)*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,
      2),0.96*msd2(2,2))) + (1.0416666666666667*Log(1.02*msq2(2,2))*(-
      0.9603999999999999*Quad(mAInput) + 1.9992*msq2(2,2)*Sqr(mAInput) + 0.9216*
      Sqr(msd2(2,2)) - 1.0404*Sqr(msq2(2,2)) + TDelta(0.98*Sqr(mAInput),1.02*msq2(
      2,2),0.96*msd2(2,2))))/(msd2(2,2)*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),
      0.96*msd2(2,2))) + (1.0850694444444444*(-1.02*msq2(2,2) + 0.98*Sqr(mAInput)
      + ((0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput))*(1.02*(0.96*msd2(2,
      2) - 1.02*msq2(2,2))*msq2(2,2) - 0.9603999999999999*Quad(mAInput) + 0.98*(
      0.96*msd2(2,2) + 2.04*msq2(2,2))*Sqr(mAInput)))/TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*
      msd2(2,2)))/Sqr(msd2(2,2))) + 3*(AbInput + MuInput/TanBeta)*(AbInput -
      MuInput*TanBeta)*((2*Log(0.96*msd2(2,2))*Log(1.02*msq2(2,2)))/(0.96*msd2(2,2
      ) - 1.02*msq2(2,2)) + Log(0.98*Sqr(mAInput))*((2*Log(1.02*msq2(2,2)))/(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(0.96*msd2(2,2)))/(-0.96*msd2(2,2) +
      1.02*msq2(2,2))) + (2*Sqr(Log(1.02*msq2(2,2))))/(-0.96*msd2(2,2) + 1.02*msq2
      (2,2)) - (2.0833333333333335*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(
      mAInput))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*
      (0.96*msd2(2,2) - 1.02*msq2(2,2))) + (1.9215686274509802*Sqr(mAInput)*TPhi(
      0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*
      msq2(2,2))*msq2(2,2))) + 3*((4 + Log(0.96*msd2(2,2)) - Log(0.98*msu2(2,2)))*
      Log(0.98*Sqr(mAInput)) + Log(1.02*msq2(2,2))*(-2 - 2*Log(Sqr(SCALE))) + Log(
      0.96*msd2(2,2))*(-2 + Log(0.98*msu2(2,2)) - 2*Log(Sqr(SCALE))) + 2*Sqr(Log(
      Sqr(SCALE))) + Sqr(Log(1.02*msq2(2,2))) - (0.9607843137254901*Sqr(mAInput)*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/msq2(2,2) + (
      1.0850694444444444*(-2.8224*msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(
      mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(
      mAInput) + 1.8432*Sqr(msd2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)) -
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/Sqr(msd2(2,2))) + Sqr(AbInput -
      MuInput*TanBeta)*(3*Sqr(AbInput + MuInput/TanBeta)*(Sqr(Log(1.02*msq2(2,2)))
      /Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(0.96*msd2(2,2))*(-(Log(1.02*msq2
      (2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (2*(0.96*msd2(2,2) - 1.02*
      msq2(2,2) + 0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*TDelta(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + Log(0.98*Sqr(mAInput))*
      (Log(0.96*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - Log(1.02*msq2(2,
      2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.0416666666666667*(Sqr(0.98*Sqr
      (mAInput) + 0.96*msd2(2,2) - 1.02*msq2(2,2)) - TDelta(0.98*Sqr(mAInput),1.02
      *msq2(2,2),0.96*msd2(2,2))))/(msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + (
      1.0416666666666667*Log(1.02*msq2(2,2))*(-0.9603999999999999*Quad(mAInput) +
      1.9992*msq2(2,2)*Sqr(mAInput) + 0.9216*Sqr(msd2(2,2)) - 1.0404*Sqr(msq2(2,2)
      ) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(msd2(2,2)*(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),
      0.96*msd2(2,2))) + (1.0850694444444444*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(
      0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput))*(1.02*(0.96*msd2(2,2) -
      1.02*msq2(2,2))*msq2(2,2) - 0.9603999999999999*Quad(mAInput) + 0.98*(0.96*
      msd2(2,2) + 2.04*msq2(2,2))*Sqr(mAInput)) + (0.98*(1.92*msd2(2,2) - 1.02*
      msq2(2,2))*Sqr(mAInput) + Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))*TDelta(0.98*
      Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,
      2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) - (
      0.9607843137254901*Sqr(mAInput)*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*
      msq2(2,2)))/(msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + 3*(4/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2)) + (2*Log(Sqr(SCALE)))/(-0.96*msd2(2,2) + 1.02*
      msq2(2,2)) + Log(0.96*msd2(2,2))*((2.04*Log(Sqr(SCALE))*msq2(2,2))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (2*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Sqr(0.96
      *msd2(2,2) - 1.02*msq2(2,2)) - (Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*
      msu2(2,2) - 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(
      1.02*msq2(2,2))*((-4.08*msq2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (
      2.04*Log(Sqr(SCALE))*msq2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(
      0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput))*(-((Log(0.96*msd2
      (2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(
      2,2) - 1.02*msq2(2,2))) + (Log(1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2
      ,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      1.0850694444444444*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.8224*msd2(2,2)*
      msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput)
      - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (-1.92*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (0.9803921568627451*TDelta(0.98*Sqr(mAInput),0.98*msu2(
      2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/
      (msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))))) + (AtInput - MuInput/
      TanBeta)*(3*(AbInput + MuInput/TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(0.98
      *msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(1.02*msq2(2,2)))/(-
      1.02*msq2(2,2) + 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,
      2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2)))/(-1.02*msq2
      (2,2) + 0.98*msu2(2,2))) + (2.170138888888889*(-2.9375999999999998*msd2(2,2)
      *msq2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput
      ) - 1.9992*msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(
      2,2)) - TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*
      Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,
      2))*Sqr(msd2(2,2))) - (2.170138888888889*(-2.8224*msd2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msd2(2,2)))) + 3*(AbInput - MuInput*
      TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(0.96*msd2(2,2)))/(-0.96*msd2(2,2) +
      1.02*msq2(2,2)) + 2*Log(1.02*msq2(2,2))*(1/(0.96*msd2(2,2) - 1.02*msq2(2,2))
      + 1/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (2*Log(0.98*msu2(2,2)))/(1.02*msq2
      (2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2))*(0.96
      *msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2*Log(0.98*msu2(2,2))*(0.96
      *msd2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + (2*Log(1.02*msq2(2,2))*Log(
      0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*(-
      2.04*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/((0.96*msd2(2,
      2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*TDelta(0.98*Sqr
      (mAInput),1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,
      2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.9792*msd2(2,2)*
      msq2(2,2) + 0.9216*Sqr(msd2(2,2)))) - (1.9607843137254901*TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2
      ),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))) - (2*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*
      msd2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/((1.02*
      msq2(2,2) - 0.98*msu2(2,2))*(-0.9792*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2
      ,2)))) + (1.9607843137254901*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2)))) +
      3*(AbInput + MuInput/TanBeta)*Sqr(AbInput - MuInput*TanBeta)*(Log(0.98*Sqr(
      mAInput))*((2*Log(0.96*msd2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2
      *Log(1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.96*msd2(
      2,2))*((2*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(
      mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2
      (2,2))) + (2*Log(1.02*msq2(2,2))*(-2.04*msq2(2,2) + 0.98*Sqr(mAInput)))/((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) - (2
      *Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + ((4.08*msq2(2,2) - 1.96*Sqr(mAInput))*Sqr(Log(1.02*msq2(2
      ,2))))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2
      ))) + (2.170138888888889*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-
      2.9375999999999998*msd2(2,2)*msq2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(
      msd2(2,2)) + 1.0404*Sqr(msq2(2,2))) + (-1.92*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (1.9607843137254901*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (2.170138888888889*((0.96
      *msd2(2,2) - 1.02*msq2(2,2))*(-2.8224*msd2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (-1.92*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (1.9607843137254901*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))))) + Sqr(AtInput - MuInput/
      TanBeta)*(3*Sqr(AbInput + MuInput/TanBeta)*(Log(0.96*msd2(2,2))*(-(Log(1.02*
      msq2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*msu2(2,2))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2*(0.96*msd2(2,2) - 1.02*msq2(2,2) +
      0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + Log(0.98*Sqr(mAInput))*(Log(1.02
      *msq2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - Log(0.98*msu2(2,2))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.0416666666666667*(Sqr(0.98*Sqr(mAInput
      ) + 0.96*msd2(2,2) - 1.02*msq2(2,2)) - TDelta(0.98*Sqr(mAInput),1.02*msq2(2,
      2),0.96*msd2(2,2))))/(msd2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + (1.0416666666666667*Log
      (1.02*msq2(2,2))*(-0.9603999999999999*Quad(mAInput) + 1.9992*msq2(2,2)*Sqr(
      mAInput) + 0.9216*Sqr(msd2(2,2)) - 1.0404*Sqr(msq2(2,2)) + TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(msd2(2,2)*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) + (
      1.0850694444444444*(-((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) -
      1.02*msq2(2,2) + 0.98*Sqr(mAInput))*(1.02*msq2(2,2)*(-0.96*msd2(2,2) + 1.02*
      msq2(2,2)) + 0.9603999999999999*Quad(mAInput) - 0.98*(0.96*msd2(2,2) + 2.04*
      msq2(2,2))*Sqr(mAInput))) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*
      msd2(2,2))*(2.9375999999999998*msd2(2,2)*msq2(2,2) + 1.02*msq2(2,2)*(-2.04*
      msq2(2,2) + 0.98*msu2(2,2)) - 0.9603999999999999*Quad(mAInput) + 2.94*(0.96*
      msd2(2,2) + 1.02*msq2(2,2))*Sqr(mAInput) - 0.9603999999999999*msu2(2,2)*Sqr(
      mAInput) - 1.8432*Sqr(msd2(2,2)) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),
      0.96*msd2(2,2))))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(
      Sqr(msd2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput)
      ,1.02*msq2(2,2),0.96*msd2(2,2))) + (1.0850694444444444*(-2.8224*msd2(2,2)*
      msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput)
      - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Sqr
      (msd2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + 3*(4/(-1.02*msq2(2,2) +
      0.98*msu2(2,2)) + Log(Sqr(SCALE))*(2/(-1.02*msq2(2,2) + 0.98*msu2(2,2)) - (
      1.96*Log(0.98*msu2(2,2))*msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      Log(1.02*msq2(2,2))*((2*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2
      ) - 0.98*msu2(2,2)) + (1.96*Log(Sqr(SCALE))*msu2(2,2))/Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2)) + (Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*Sqr(
      mAInput))*(-((Log(1.02*msq2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2(2,2) + 0.98*
      Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (Log(0.98*msu2(2,2))*
      (-1.02*msq2(2,2) + 0.98*msu2(2,2) + 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) - (3.92*Log(0.98*msu2(2,2))*msu2(2,2))/Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2)) - ((1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput))*Sqr
      (Log(1.02*msq2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (
      0.9803921568627451*(0.98*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput) +
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (0.9803921568627451*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,
      2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(
      msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + 3*(AbInput + MuInput/
      TanBeta)*(AbInput - MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(1.02*
      msq2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2*Log(0.98*msu2(2,2)))/
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((-2*Log(1.02*
      msq2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*Log(
      0.98*msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) - (2
      *Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput))*
      Sqr(Log(1.02*msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))) - (2.0833333333333335*((1.02*msq2(2,2) - 0.98*msu2(2,
      2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)) + TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,
      2),0.96*msd2(2,2)))/(msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2))) + (1.9607843137254901*(0.98*(1.02*msq2(2,2) -
      0.98*msu2(2,2))*Sqr(mAInput) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*
      msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))
      + (2.0833333333333335*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)
      )*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (
      1.9607843137254901*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) -
      1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + 3*Sqr(
      AbInput + MuInput/TanBeta)*Sqr(AbInput - MuInput*TanBeta)*(-(((1.02*msq2(2,2
      ) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/(Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + Log(
      0.96*msd2(2,2))*((Log(1.02*msq2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) - (Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) - (2*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)
      ))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + (
      1.0416666666666667*Log(0.98*Sqr(mAInput))*(Sqr(0.98*Sqr(mAInput) + 0.96*msd2
      (2,2) - 1.02*msq2(2,2)) - TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(
      2,2))))/(msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) + Log(
      1.02*msq2(2,2))*((Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + (-0.9603999999999999*Quad(mAInput) + 1.9992*msq2(2,2)*
      Sqr(mAInput) + 0.9216*Sqr(msd2(2,2)) - 1.0404*Sqr(msq2(2,2)) + TDelta(0.98*
      Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,
      2))*(-0.9792*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)))*TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + (1.0850694444444444*((0.96*msd2(
      2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) -
      1.02*msq2(2,2) + 0.98*Sqr(mAInput))*(1.02*(0.96*msd2(2,2) - 1.02*msq2(2,2))*
      msq2(2,2) - 0.9603999999999999*Quad(mAInput) + 0.98*(0.96*msd2(2,2) + 2.04*
      msq2(2,2))*Sqr(mAInput)) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2
      (2,2))*(0.9603999999999999*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*Quad(mAInput)
      + 0.98*Sqr(mAInput)*(0.96*msd2(2,2)*(1.02*msq2(2,2) - 1.96*msu2(2,2)) + 1.02
      *msq2(2,2)*(-3.06*msq2(2,2) + 0.98*msu2(2,2)) + 2.7648*Sqr(msd2(2,2))) - (
      1.92*msd2(2,2) - 2.04*msq2(2,2) + 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) + (1.92*msd2(2,2) - 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02
      *msq2(2,2),0.96*msd2(2,2))))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2
      (2,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)
      )) - (0.9803921568627451*(0.98*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput
      ) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (1.0850694444444444*
      ((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.8224*msd2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (-1.92*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      0.9803921568627451*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))))) + (
      AtInput - MuInput/TanBeta)*(3*(AbInput + MuInput/TanBeta)*(Log(0.98*Sqr(
      mAInput))*((2*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*
      Log(1.02*msq2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + Log(0.96*msd2(2,2
      ))*((2*Log(1.02*msq2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*
      msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (2.170138888888889*(-
      2.9375999999999998*msd2(2,2)*msq2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(
      msd2(2,2)) + 1.0404*Sqr(msq2(2,2)) - TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2)
      ,0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msd2(2,2))) - (2.170138888888889*(-
      2.8224*msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(
      2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(
      msd2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*
      msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msd2(2,2)))) + 3*(AbInput
       - MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(0.96*msd2(2,2)))/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2)) + 2*Log(1.02*msq2(2,2))*(1/(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + 1/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (2*Log(0.98*msu2(2
      ,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((2*Log(1.02*
      msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2*Log(0.98
      *msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + (2*Log(
      1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 0.98*
      Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) + (2*(-2.04*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2
      ))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      (2*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      (-0.9792*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)))) - (1.9607843137254901
      *TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*
      msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2*TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*
      msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.9792*msd2(2,2)*msq2(2,2)
      + 0.9216*Sqr(msd2(2,2)))) + (1.9607843137254901*TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*(1.02*msq2(2,2) -
      0.98*msu2(2,2)))) + 3*(AbInput + MuInput/TanBeta)*Sqr(AbInput - MuInput*
      TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(0.96*msd2(2,2)))/Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) - (2*Log(1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2
      ,2))) + Log(0.96*msd2(2,2))*((2*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*
      msu2(2,2) - 0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + (2*Log(1.02*msq2(2,2))*(-2.04*msq2(2,2) +
      0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)))) - (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(1.02*msq2(2,
      2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))
      *Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + ((4.08*msq2(2,2) - 1.96*Sqr(mAInput
      ))*Sqr(Log(1.02*msq2(2,2))))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + (2.170138888888889*((0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(-2.9375999999999998*msd2(2,2)*msq2(2,2) + 0.9603999999999999*
      Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*msq2(2,2)*Sqr(mAInput
      ) + 1.8432*Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2))) + (-1.92*msd2(2,2) + 1.02
      *msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*
      msu2(2,2))*Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      1.9607843137254901*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*(1.02*msq2
      (2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (
      2.170138888888889*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.8224*msd2(2,2)*msu2
      (2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (-1.92*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (1.9607843137254901*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + (AtInput + MuInput*
      TanBeta)*(3*(Log(0.98*Sqr(mAInput))*((2*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2)
      - 0.98*msu2(2,2)) + (2*Log(1.02*msq2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2
      ))) + (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*
      msu2(2,2)) + (2*Sqr(Log(1.02*msq2(2,2))))/(1.02*msq2(2,2) - 0.98*msu2(2,2))
      - (1.9223375624759709*(-0.98*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(
      mAInput) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(
      0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((1.02*msq2(2,2) - 0.98*
      msu2(2,2))*Sqr(msq2(2,2))) - (1.9223375624759709*(-2.9988*msq2(2,2)*msu2(2,2
      ) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msq2(2,2)))) + 3*(AbInput + MuInput/
      TanBeta)*(AbInput - MuInput*TanBeta)*(Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(
      2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2)))
      + (2*Log(0.98*msu2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-1.02*msq2(2,2
      ) + 0.98*msu2(2,2)))) + (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*Sqr(Log(
      1.02*msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 0.98
      *msu2(2,2))) + (2.170138888888889*(-2.9375999999999998*msd2(2,2)*msq2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*
      msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2)) -
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*
      (1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msd2(2,2))) + (1.9223375624759709*(-
      0.98*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,
      2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) -
      0.98*msu2(2,2))*Sqr(msq2(2,2))) - (2.170138888888889*(-2.8224*msd2(2,2)*msu2
      (2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msd2(
      2,2))) + (1.9223375624759709*(-2.9988*msq2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msq2(
      2,2)))) + 3*(AbInput + MuInput/TanBeta)*Cube(AbInput - MuInput*TanBeta)*(((4
      *Log(0.96*msd2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*Log(1.02*
      msq2(2,2)))/Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2)))*Log(0.98*Sqr(mAInput)) +
      Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 3.06*msq2(2,2)
      - 1.96*Sqr(mAInput)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*(-1.02*msq2(2,2
      ) + 0.98*msu2(2,2))) + (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,
      2) + 1.96*msu2(2,2) - 1.96*Sqr(mAInput)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2
      ,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2)))) - (2*Log(1.02*msq2(2,2))*Log(0.98*
      msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) + 1.96*msu2(2,2) - 1.96*Sqr(
      mAInput)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) + (2*(0.96*msd2(2,2) + 3.06*msq2(2,2) - 1.96*Sqr(mAInput))*Sqr(
      Log(1.02*msq2(2,2))))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + (2.170138888888889*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(
      -2.9375999999999998*msd2(2,2)*msq2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9992*msq2(2,2)*Sqr(mAInput) + 1.8432*Sqr(
      msd2(2,2)) + 1.0404*Sqr(msq2(2,2))) + (-2.88*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,
      2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msd2(2,2))) + (1.9223375624759709*
      (0.98*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))
      *Sqr(mAInput) - (0.96*msd2(2,2) - 3.06*msq2(2,2))*TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*
      msq2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*Sqr(msq2(2,2))) - (2.170138888888889*((0.96*msd2(2,2) - 1.02*msq2
      (2,2))*(-2.8224*msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) +
      1.8432*Sqr(msd2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) + (-2.88*msd2(2,2
      ) + 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))
      *TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Cube(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msd2(2,2))) + (
      1.9223375624759709*(-((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.9988*msq2(2,2)*
      msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput)
      - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)))) + (0.96*msd2(2,2) - 3.06*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,
      2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msq2(2,2)))) + 3*Sqr(AbInput -
      MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((-2*Log(0.96*msd2(2,2)))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2))*(0.96*msd2(2,
      2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))
      *Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (2*Log(0.98*msu2(2,2))*(0.96*msd2(2
      ,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2)
      )*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + (2*Log(1.02*msq2(2,2))*Log(0.98*
      msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((1.02*
      msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (2*(-
      0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))
      ))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))
      + (2.0833333333333335*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)
      )*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(1.02*
      msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      1.9223375624759709*(0.98*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-5.1*msq2(2,2)
      + 0.98*Sqr(mAInput))*Sqr(mAInput) + (0.96*msd2(2,2) - 2.04*msq2(2,2))*TDelta
      (0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2))) - (2.0833333333333335*TDelta(
      0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),0.98
      *msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr
      (0.96*msd2(2,2) - 1.02*msq2(2,2))) + (1.9223375624759709*((0.96*msd2(2,2) -
      1.02*msq2(2,2))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(
      mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(
      mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) - (
      0.96*msd2(2,2) - 2.04*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      msq2(2,2)))))) + ((1 + Sqr(TanBeta))*(3*Quad(AbInput - MuInput*TanBeta)*((
      0.9607843137254901*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 5.1*msq2(2,2))*msu2
      (2,2))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)) + (
      0.9803921568627451*Log(Sqr(SCALE))*(1.02*msq2(2,2)*(1.02*msq2(2,2) + 4.9*
      msu2(2,2) - 10.3*Sqr(MuInput)) + 0.96*msd2(2,2)*(5.1*msq2(2,2) + 0.98*msu2(2
      ,2) - 2.06*Sqr(MuInput))))/(Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2))*msq2(2,2)
      ) + (0.9803921568627451*(1.02*msq2(2,2)*(3.06*msq2(2,2) + 4.9*msu2(2,2) -
      16.48*Sqr(MuInput)) + 0.96*msd2(2,2)*(9.18*msq2(2,2) + 0.98*msu2(2,2) - 2.06
      *Sqr(MuInput))))/(Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2))*msq2(2,2)) - (2*
      PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))*(0.96*msd2(2,2) +
      2.04*msq2(2,2) - 3.09*Sqr(MuInput))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput)))/
      Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*PolyLog(2,1 - (1.0098039215686274
      *Sqr(MuInput))/msq2(2,2))*(0.96*msd2(2,2) + 2.04*msq2(2,2) - 3.09*Sqr(
      MuInput))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput)))/Quad(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) + Log(1.03*Sqr(MuInput))*((-6.365399999999999*Log(1.02*msq2(2,2))
      *Quad(MuInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2.019607843137255*(
      0.96*msd2(2,2) + 2.04*msq2(2,2))*Sqr(MuInput))/(Cube(-0.96*msd2(2,2) + 1.02*
      msq2(2,2))*msq2(2,2))) + ((0.96*msd2(2,2) + 2.04*msq2(2,2) - 3.09*Sqr(
      MuInput))*(-0.96*msd2(2,2) + 1.03*Sqr(MuInput))*Sqr(Log(0.96*msd2(2,2))))/
      Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + ((0.96*msd2(2,2) + 2.04*msq2(2,2) -
      3.09*Sqr(MuInput))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))*Sqr(Log(1.02*msq2(2,
      2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(1.02*msq2(2,2))*((1.96*Log
      (0.98*msu2(2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2))*msu2(2,2))/Quad(0.96*msd2
      (2,2) - 1.02*msq2(2,2)) - (1.02*msq2(2,2)*(1.02*msq2(2,2) + 1.96*msu2(2,2) -
      10.3*Sqr(MuInput)) + 0.96*msd2(2,2)*(10.2*msq2(2,2) + 3.92*msu2(2,2) - 8.24*
      Sqr(MuInput)) + 0.9216*Sqr(msd2(2,2)))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))
      - (2*Log(Sqr(SCALE))*(1.02*msq2(2,2)*(0.98*msu2(2,2) - 2.06*Sqr(MuInput)) +
      1.92*msd2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(MuInput)) +
      0.9216*Sqr(msd2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.96*
      msd2(2,2))*((-1.96*Log(0.98*msu2(2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2))*
      msu2(2,2))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (6.365399999999999*Log(
      1.03*Sqr(MuInput))*Quad(MuInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2
      *Log(Sqr(SCALE))*(1.02*msq2(2,2)*(0.98*msu2(2,2) - 2.06*Sqr(MuInput)) + 1.92
      *msd2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(MuInput)) + 0.9216*
      Sqr(msd2(2,2))))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(0.96*msd2(2,2)*
      (3.06*msq2(2,2) + 1.96*msu2(2,2) - 7.21*Sqr(MuInput)) + 1.02*msq2(2,2)*(0.98
      *msu2(2,2) - 2.06*Sqr(MuInput)) + 2.7648*Sqr(msd2(2,2))))/Quad(0.96*msd2(2,2
      ) - 1.02*msq2(2,2)))) + 3*Sqr(AbInput - MuInput*TanBeta)*((
      1.9215686274509802*Log(0.98*msu2(2,2))*msu2(2,2))/(msq2(2,2)*(-0.96*msd2(2,2
      ) + 1.02*msq2(2,2))) + (1.9607843137254901*Log(Sqr(SCALE))*(2.04*msq2(2,2) +
      0.98*msu2(2,2) - 2.06*Sqr(MuInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2
      (2,2)) + (1.9607843137254901*(3.06*msq2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(
      MuInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)) + Log(0.96*msd2(2,2
      ))*((1.96*Log(0.98*msu2(2,2))*msu2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)
      ) - (2*Log(Sqr(SCALE))*(1.92*msd2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(MuInput))
      )/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2*(2.88*msd2(2,2) + 0.98*msu2(2,2)
      - 2.06*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(
      2,2))*((-1.96*Log(0.98*msu2(2,2))*msu2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(
      2,2)) + (2*Log(Sqr(SCALE))*(1.92*msd2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(
      MuInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(0.96*msd2(2,2) + 2.04*
      msq2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))) + (4*PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))*
      (0.96*msd2(2,2) - 1.03*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) -
      (4*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2))*(0.96*msd2(2,2
      ) - 1.03*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(0.96*msd2
      (2,2) - 1.03*Sqr(MuInput))*Sqr(Log(0.96*msd2(2,2))))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) - (2*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))*Sqr(Log(1.02*msq2(
      2,2))))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4.12*Log(1.03*Sqr(MuInput))*
      Sqr(MuInput))/(0.9792*msd2(2,2)*msq2(2,2) - 1.0404*Sqr(msq2(2,2)))) + 3*(-
      0.5 + Log(1.02*msq2(2,2))*(-1 - 2*Log(Sqr(SCALE))) - (0.9607843137254901*
      msu2(2,2))/msq2(2,2) + 2*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/
      msq2(2,2)) + 1.03*(1.9607843137254901/msq2(2,2) + 1/(0.98*msu2(2,2) - 1.03*
      Sqr(MuInput)))*Sqr(MuInput) + 2*Sqr(Log(Sqr(SCALE))) + Sqr(Log(1.02*msq2(2,2
      ))) + Log(Sqr(SCALE))*((0.9803921568627451*(-3.0282*msu2(2,2)*Sqr(MuInput) +
      2.06*Sqr(MuInput)*(1.02*msq2(2,2) + 1.03*Sqr(MuInput)) + 0.9603999999999999*
      Sqr(msu2(2,2))))/(msq2(2,2)*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput))) + Log(
      0.98*msu2(2,2))*(-2 + (2.1218*Quad(MuInput))/Sqr(-0.98*msu2(2,2) + 1.03*Sqr(
      MuInput)))) + PolyLog(2,1 - (1.0510204081632653*Sqr(MuInput))/msu2(2,2))*(2
      - (4.2436*Quad(MuInput))/Sqr(-0.98*msu2(2,2) + 1.03*Sqr(MuInput))) + Sqr(Log
      (0.98*msu2(2,2)))*(1 - (2.1218*Quad(MuInput))/Sqr(-0.98*msu2(2,2) + 1.03*Sqr
      (MuInput))) + Log(1.03*Sqr(MuInput))*((2.1218*Log(0.98*msu2(2,2))*Quad(
      MuInput))/Sqr(-0.98*msu2(2,2) + 1.03*Sqr(MuInput)) - (2.1218*Log(Sqr(SCALE))
      *Quad(MuInput))/Sqr(-0.98*msu2(2,2) + 1.03*Sqr(MuInput)) - (
      1.0098039215686274*Sqr(MuInput)*(1.96*(1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2
      (2,2) + 2.1218*Quad(MuInput) + 1.03*(1.02*msq2(2,2) - 3.92*msu2(2,2))*Sqr(
      MuInput)))/(msq2(2,2)*Sqr(-0.98*msu2(2,2) + 1.03*Sqr(MuInput)))) + (
      0.9607843137254901*Log(0.98*msu2(2,2))*msu2(2,2)*(0.98*(1.02*msq2(2,2) +
      0.98*msu2(2,2))*msu2(2,2) + 1.0609*Quad(MuInput) + 2.06*(1.02*msq2(2,2) -
      0.98*msu2(2,2))*Sqr(MuInput)))/(msq2(2,2)*Sqr(-0.98*msu2(2,2) + 1.03*Sqr(
      MuInput)))) + (3*Quad(AbInput - MuInput*TanBeta)*((0.9607843137254901*(0.96*
      msd2(2,2) + 5.1*msq2(2,2))*Sqr(mAInput))/(Cube(-0.96*msd2(2,2) + 1.02*msq2(2
      ,2))*msq2(2,2)) + (0.9607843137254901*Log(Sqr(SCALE))*(0.96*msd2(2,2) + 5.1*
      msq2(2,2))*Sqr(mAInput))/(Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2))*msq2(2,2))
      + Log(0.98*Sqr(mAInput))*((0.9607843137254901*(0.96*msd2(2,2) + 5.1*msq2(2,2
      ))*Sqr(mAInput))/(Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)) - (1.96*
      Log(0.96*msd2(2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2))*Sqr(mAInput))/Quad(
      0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.96*Log(1.02*msq2(2,2))*(1.92*msd2(2,2)
      + 1.02*msq2(2,2))*Sqr(mAInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log
      (1.02*msq2(2,2))*((-1.96*(1.92*msd2(2,2) + 1.02*msq2(2,2))*Sqr(mAInput))/
      Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (1.96*Log(Sqr(SCALE))*(1.92*msd2(2,2
      ) + 1.02*msq2(2,2))*Sqr(mAInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) +
      Log(0.96*msd2(2,2))*((1.96*(1.92*msd2(2,2) + 1.02*msq2(2,2))*Sqr(mAInput))/
      Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (1.96*Log(Sqr(SCALE))*(1.92*msd2(2,2
      ) + 1.02*msq2(2,2))*Sqr(mAInput))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)))) +
      3*(0.5 + Log(0.96*msd2(2,2))*(-2 + Log(0.98*msu2(2,2)) - 2*Log(Sqr(SCALE)))
      + 0.3333333333333333*Sqr(3.141592653589793) + Log(0.98*Sqr(mAInput))*(Log(
      0.96*msd2(2,2)) - Log(0.98*msu2(2,2)) - 2*Log(Sqr(SCALE)) + (
      0.9803921568627451*(1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/msq2(2,2)) - (
      0.9607843137254901*Sqr(mAInput))/msq2(2,2) + Log(Sqr(SCALE))*(1 - (
      0.9607843137254901*Sqr(mAInput))/msq2(2,2)) + Sqr(Log(0.98*Sqr(mAInput))) +
      2*Sqr(Log(Sqr(SCALE))) + (1.0850694444444444*(-2.8224*msd2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/Sqr(
      msd2(2,2))) + 3*Sqr(AbInput - MuInput*TanBeta)*((1.9607843137254901*(-2.04*
      msq2(2,2) + 0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)
      ) + (1.9607843137254901*Log(Sqr(SCALE))*(-1.02*msq2(2,2) + 0.98*Sqr(mAInput)
      ))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)) + Log(0.96*msd2(2,2))*((2*(
      0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) - (Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2*Log(Sqr(SCALE)
      )*(-1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)
      )) + Log(1.02*msq2(2,2))*((Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2
      ,2) - 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(-2.04*
      msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*
      Log(Sqr(SCALE))*(-1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput))*((1.9215686274509802*Sqr(mAInput))
      /(msq2(2,2)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + (Log(0.96*msd2(2,2))*(-
      1.02*msq2(2,2) + 0.98*msu2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) - (Log(1.02*msq2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2(2,2) +
      0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      1.0850694444444444*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.8224*msd2(2,2)*
      msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput)
      - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (-1.92*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (0.9803921568627451*TDelta(0.98*Sqr(mAInput),0.98*msu2(
      2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/
      (msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + (AtInput + MuInput*
      TanBeta)*(3*(AbInput - MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(0.96
      *msd2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2)))/(-
      0.96*msd2(2,2) + 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*(4/(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2
      *Log(0.98*msu2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + Log(0.96*msd2(2,
      2))*((2*Log(0.98*msu2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + 4/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(-0.96*msd2(2,2) + 1.02*
      msq2(2,2))) + (2.170138888888889*(-2.8224*msd2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msd2(2,2))) - (1.9223375624759709*(-
      2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(
      2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(
      msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2)))) + 3*Cube(
      AbInput - MuInput*TanBeta)*(Log(0.96*msd2(2,2))*((4*Log(Sqr(SCALE))*(0.96*
      msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*(
      2.88*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (2
      *Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) + 1.96*msu2(2,2) -
      1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2
      ,2))*((-4*Log(Sqr(SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2)))/Cube(0.96*msd2(
      2,2) - 1.02*msq2(2,2)) - (4*(0.96*msd2(2,2) + 3.06*msq2(2,2)))/Cube(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*
      msq2(2,2) + 1.96*msu2(2,2) - 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*
      msq2(2,2))) + Log(0.98*Sqr(mAInput))*((-2*Log(0.96*msd2(2,2))*(0.96*msd2(2,2
      ) + 1.02*msq2(2,2) - 1.96*msu2(2,2) + 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2
      ) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2
      ) - 1.96*msu2(2,2) + 1.96*Sqr(mAInput)))/Cube(0.96*msd2(2,2) - 1.02*msq2(2,2
      ))) - 16/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (8*Log(Sqr(SCALE)))/Sqr(0.96
      *msd2(2,2) - 1.02*msq2(2,2)) + (2.170138888888889*((0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(-2.8224*msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.8224*msd2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) +
      1.8432*Sqr(msd2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) + (-2.88*msd2(2,2
      ) + 1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))
      *TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Cube(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Sqr(msd2(2,2))) + (1.9223375624759709*(-((0.96*msd2(2,2) -
      1.02*msq2(2,2))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(
      mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(
      mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)))) + (
      0.96*msd2(2,2) - 3.06*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(
      Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqr(msq2(2,2))))) + Sqr(AtInput +
      MuInput*TanBeta)*(3*Quad(AbInput - MuInput*TanBeta)*((0.9803921568627451*Log
      (Sqr(SCALE))*(0.96*msd2(2,2) + 5.1*msq2(2,2)))/(Cube(-0.96*msd2(2,2) + 1.02*
      msq2(2,2))*msq2(2,2)) + (0.9803921568627451*(0.96*msd2(2,2) + 11.22*msq2(2,2
      )))/(Cube(-0.96*msd2(2,2) + 1.02*msq2(2,2))*msq2(2,2)) + Log(0.96*msd2(2,2))
      *((2*Log(Sqr(SCALE))*(1.92*msd2(2,2) + 1.02*msq2(2,2)))/Quad(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) + (2*(4.8*msd2(2,2) + 1.02*msq2(2,2)))/Quad(0.96*msd2(2,2)
      - 1.02*msq2(2,2)) - (Log(0.98*msu2(2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2) +
      2.94*msu2(2,2) - 2.94*Sqr(mAInput)))/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)))
      + Log(1.02*msq2(2,2))*((-2*Log(Sqr(SCALE))*(1.92*msd2(2,2) + 1.02*msq2(2,2))
      )/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(0.98*msu2(2,2))*(1.92*msd2(2,
      2) + 1.02*msq2(2,2) + 2.94*msu2(2,2) - 2.94*Sqr(mAInput)))/Quad(0.96*msd2(2,
      2) - 1.02*msq2(2,2)) - (0.9803921568627451*(1.9584*msd2(2,2)*msq2(2,2) +
      0.9216*Sqr(msd2(2,2)) + 9.3636*Sqr(msq2(2,2)) - (Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) -
      2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) +
      2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))/(msq2(2,2)*Quad(0.96*msd2(2,2) -
      1.02*msq2(2,2)))) + (0.4805843906189927*Log(0.98*msu2(2,2))*((-1.02*msq2(2,2
      ) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*(-2.9988*msq2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (2.04*msq2(2,2) + 0.98*msu2(2,2) - 0.98
      *Sqr(mAInput))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))/(
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2))*TDelta(0.98*Sqr(mAInput)
      ,0.98*msu2(2,2),1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput))*(-((Log(0.96*msd2(
      2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2) - 2.94*msu2(2,2) + 2.94*Sqr(mAInput))
      )/Quad(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (Log(1.02*msq2(2,2))*(1.92*msd2(2
      ,2) + 1.02*msq2(2,2) - 2.94*msu2(2,2) + 2.94*Sqr(mAInput)))/Quad(0.96*msd2(2
      ,2) - 1.02*msq2(2,2)) + (0.4805843906189927*(-((1.02*msq2(2,2) - 0.98*msu2(2
      ,2) + 0.98*Sqr(mAInput))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*
      Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)
      *Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))
      + (2.04*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,
      2))*Sqr(msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))
      + (1.0850694444444444*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.8224*msd2(2,2)*
      msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*Sqr(mAInput)
      - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (-3.84*msd2(2,2) + 1.02*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(Quad(0.96*msd2(2,2) - 1.02*msq2(2,
      2))*Sqr(msd2(2,2))) + (0.4711611672735223*(-((0.96*msd2(2,2) - 1.02*msq2(2,2
      ))*(-5.9976*(0.96*msd2(2,2) - 3.06*msq2(2,2))*msq2(2,2)*msu2(2,2) +
      0.9603999999999999*(2.88*msd2(2,2) - 7.140000000000001*msq2(2,2))*Quad(
      mAInput) + 1.96*(-2.88*msd2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2)) + 1.02*
      msq2(2,2)*(9.18*msq2(2,2) + 6.859999999999999*msu2(2,2)))*Sqr(mAInput) +
      2.0808*(0.96*msd2(2,2) - 5.1*msq2(2,2))*Sqr(msq2(2,2)) + 0.9603999999999999*
      (2.88*msd2(2,2) - 7.140000000000001*msq2(2,2))*Sqr(msu2(2,2)))) + ((0.98*(-
      1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2) + 0.9603999999999999*Quad(mAInput
      ) - 0.98*(1.02*msq2(2,2) + 1.96*msu2(2,2))*Sqr(mAInput))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(
      mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(
      mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)) + 2*(-3.9168*msd2(2,
      2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)) + 6.2424*Sqr(msq2(2,2)))*TDelta(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2)))/(Cube(msq2(2,2))*Quad(0.96*msd2(2,2) - 1.02*msq2(
      2,2)))) + 3*(-0.9803921568627451/msq2(2,2) - (0.9803921568627451*Log(Sqr(
      SCALE)))/msq2(2,2) + (0.9803921568627451*Log(1.02*msq2(2,2))*(-2.9988*msq2(2
      ,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(
      mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2))
      + 0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2
      ),1.02*msq2(2,2))))/(msq2(2,2)*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2))) + (0.4805843906189927*Log(0.98*msu2(2,2))*((-1.02*msq2(2,2) -
      0.98*msu2(2,2) + 0.98*Sqr(mAInput))*(-2.9988*msq2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + (2.04*msq2(2,2) + 0.98*msu2(2,2) - 0.98
      *Sqr(mAInput))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))/(
      Sqr(msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (
      0.4805843906189927*Log(0.98*Sqr(mAInput))*(-((1.02*msq2(2,2) - 0.98*msu2(2,2
      ) + 0.98*Sqr(mAInput))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*
      Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)
      *Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))
      + (2.04*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))/(Sqr(msq2(2,2))*TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (0.4711611672735223*((0.98*(-1.02
      *msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2) + 0.9603999999999999*Quad(mAInput) -
      0.98*(1.02*msq2(2,2) + 1.96*msu2(2,2))*Sqr(mAInput))*(-2.9988*msq2(2,2)*msu2
      (2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2)
      ,1.02*msq2(2,2))*(5.9976*msq2(2,2)*msu2(2,2) - 2.8811999999999998*Quad(
      mAInput) + 5.88*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(mAInput) - 2.0808*Sqr(
      msq2(2,2)) - 2.8811999999999998*Sqr(msu2(2,2)) + 2*TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2))))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))/(Cube(msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))) + 3*Sqr(AbInput - MuInput*TanBeta)*(1.9607843137254901/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)) + Log(0.96*msd2(2,2))*(-2/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + Log(0.98*msu2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) - (2*Log(Sqr(SCALE)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (2*
      Log(Sqr(SCALE)))/(0.9792*msd2(2,2)*msq2(2,2) - 1.0404*Sqr(msq2(2,2))) + (
      0.9611687812379854*Log(0.98*msu2(2,2))*(-((-1.02*msq2(2,2) - 0.98*msu2(2,2)
      + 0.98*Sqr(mAInput))*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(
      mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(
      mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)))) + (-
      2.04*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TDelta(0.98*Sqr(mAInput
      ),0.98*msu2(2,2),1.02*msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + Log(
      1.02*msq2(2,2))*(-(Log(0.98*msu2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))
      + (2*Log(Sqr(SCALE)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      1.9607843137254901*(-((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.9988*msq2(2,2)*
      msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput)
      - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)))) + 0.96*msd2(2,2)*TDelta(0.98*Sqr(mAInput
      ),0.98*msu2(2,2),1.02*msq2(2,2))))/(msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2
      (2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + Log(0.98*
      Sqr(mAInput))*(Log(0.96*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) -
      Log(1.02*msq2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      0.9611687812379854*((1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*(-
      2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(
      2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(
      msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) + (-2.04*msq2(2,2) + 0.98*
      msu2(2,2) - 0.98*Sqr(mAInput))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2))*TDelta(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + (1.0850694444444444*(-2.8224
      *msd2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.8224*msd2(2,2)*
      Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.8432*Sqr(msd2(2
      ,2)) + 0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*
      msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2
      ,2)))/(Sqr(msd2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      0.9423223345470446*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(0.98*(1.02*msq2(2,2)
      - 0.98*msu2(2,2))*msu2(2,2) - 0.9603999999999999*Quad(mAInput) + 0.98*(1.02*
      msq2(2,2) + 1.96*msu2(2,2))*Sqr(mAInput))*(-2.9988*msq2(2,2)*msu2(2,2) +
      0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2)
      ,1.02*msq2(2,2))*(2.9988*msq2(2,2)*(-1.92*msd2(2,2) + 3.06*msq2(2,2))*msu2(2
      ,2) + 0.9603999999999999*(2.88*msd2(2,2) - 4.08*msq2(2,2))*Quad(mAInput) +
      0.98*(-5.76*msd2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2)) + 1.02*msq2(2,2)*(
      9.18*msq2(2,2) + 7.84*msu2(2,2)))*Sqr(mAInput) + 2.0808*(0.96*msd2(2,2) -
      2.04*msq2(2,2))*Sqr(msq2(2,2)) + 0.9603999999999999*(2.88*msd2(2,2) - 4.08*
      msq2(2,2))*Sqr(msu2(2,2)) + (-1.92*msd2(2,2) + 3.06*msq2(2,2))*TDelta(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))*TPhi(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2)))/(Cube(msq2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2
      ,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))))/(1 + Sqr(
      TanBeta))))/Sqr(TanBeta)))/Quad(1 + (0.0625*(1 + Sqr(TanBeta))*(0.25*Log(Sqr
      (MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1 +
      Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(
      TanBeta)))*Sqr(Yu(2,2)))/(Sqr(3.141592653589793)*Sqr(TanBeta)) + 0.5*((-
      0.03125*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2
      )/msd2(2,2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(3.141592653589793)) - (
      0.03125*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2
      )/msu2(2,2)))))/(Sqrt(Abs(msq2(2,2)*msu2(2,2)))*Sqr(3.141592653589793))) + (
      0.0625*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/Sqr(SCALE)) +
      (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)
      ) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2))/MuInput)))/
      Sqr(3.141592653589793) + (0.08333333333333333*Sqr(g3)*(1 + Log(Sqr(M3Input)/
      Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/
      M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(
      msd2(2,2))/M3Input))/M3Input))/Sqr(3.141592653589793) + (0.0625*(1 + Sqr(
      TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput
      - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/
      MuInput))/MuInput))/(Sqr(3.141592653589793)*Sqr(TanBeta))) + (Quad(Yu(2,2))*
      Sqr(Yd(2,2))*(3*Sqr(AbInput - MuInput*TanBeta)*((0.9803921568627451*(0.96*
      msd2(2,2) - 5.1*msq2(2,2)))/(msq2(2,2)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))) +
      (0.9803921568627451*Log(Sqr(SCALE))*(0.96*msd2(2,2) - 3.06*msq2(2,2)))/(msq2
      (2,2)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + (6*PolyLog(2,1 - (
      0.9411764705882352*msd2(2,2))/msq2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2))
      + Log(0.96*msd2(2,2))*((1.92*Log(1.02*msq2(2,2))*msd2(2,2))/Sqr(0.96*msd2(2,
      2) - 1.02*msq2(2,2)) - (1.92*Log(Sqr(SCALE))*msd2(2,2))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + (0.9411764705882352*msd2(2,2)*(0.96*msd2(2,2) - 5.1*msq2(2
      ,2)))/(msq2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) + Log(1.02*msq2(2,2)
      )*((1.92*Log(Sqr(SCALE))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      0.96*msd2(2,2) + 3.06*msq2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (
      1.92*msd2(2,2)*Sqr(Log(1.02*msq2(2,2))))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)
      )) + (3*(AbInput - MuInput*TanBeta)*Sqr(AtInput - MuInput/TanBeta)*(Log(0.96
      *msd2(2,2))*((7.794348721990825*MuInput*Log(0.98*msu2(2,2))*msd2(2,2)*Sqrt(
      Sqr(MuInput)))/(Abs(MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,
      2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (
      7.794348721990825*MuInput*Log(1.02*msq2(2,2))*msd2(2,2)*Sqrt(Sqr(MuInput)))/
      (Abs(MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2
      (2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput)))) + (8.281495517115252*MuInput*
      Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*msq2(2,2)*Sqrt(Sqr(MuInput)))/(Abs(
      MuInput)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2(2,2
      ))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + Log(1.03*Sqr(MuInput))*((
      8.36268664963599*MuInput*Log(1.02*msq2(2,2))*Power3(Sqrt(Sqr(MuInput))))/(
      Abs(MuInput)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(
      MuInput))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + (8.36268664963599*MuInput
      *Log(0.98*msu2(2,2))*Power3(Sqrt(Sqr(MuInput))))/(Abs(MuInput)*(1.02*msq2(2,
      2) - 0.98*msu2(2,2))*(-0.96*msd2(2,2) + 1.03*Sqr(MuInput))*(-1.02*msq2(2,2)
      + 1.03*Sqr(MuInput)))) + (8.281495517115252*MuInput*msq2(2,2)*Sqrt(Sqr(
      MuInput))*Sqr(Log(1.02*msq2(2,2))))/(Abs(MuInput)*(-0.96*msd2(2,2) + 1.02*
      msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-1.02*msq2(2,2) + 1.03*Sqr(
      MuInput)))) + 3*(AbInput - MuInput*TanBeta)*((8.119113252073776*MuInput*
      PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput))/msd2(2,2))*Sqrt(Sqr(MuInput)
      )*(0.96*msd2(2,2) + 1.03*Sqr(MuInput)))/(Abs(MuInput)*(0.96*msd2(2,2) - 1.02
      *msq2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (4.140747758557626*
      MuInput*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*msq2(2,2)*Sqrt(Sqr(MuInput))
      )/(Abs(MuInput)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 1.03*
      Sqr(MuInput))) - (8.119113252073776*MuInput*PolyLog(2,1 - (
      1.0098039215686274*Sqr(MuInput))/msq2(2,2))*Sqrt(Sqr(MuInput))*(1.02*msq2(2,
      2) + 1.03*Sqr(MuInput)))/(Abs(MuInput)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-
      1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + Log(0.96*msd2(2,2))*((
      3.8971743609954124*MuInput*Log(1.02*msq2(2,2))*msd2(2,2)*Sqrt(Sqr(MuInput)))
      /(Abs(MuInput)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr
      (MuInput))) + (3.8971743609954124*MuInput*Log(0.98*msu2(2,2))*msd2(2,2)*Sqrt
      (Sqr(MuInput)))/(Abs(MuInput)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(0.96*msd2(
      2,2) - 1.03*Sqr(MuInput))) + (8.36268664963599*MuInput*Log(1.03*Sqr(MuInput)
      )*Power3(Sqrt(Sqr(MuInput))))/(Abs(MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2)
      )*(-0.96*msd2(2,2) + 1.03*Sqr(MuInput)))) + Log(1.03*Sqr(MuInput))*((
      4.181343324817995*MuInput*Log(0.98*msu2(2,2))*Power3(Sqrt(Sqr(MuInput))))/(
      Abs(MuInput)*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))*(-1.02*msq2(2,2) + 1.03*
      Sqr(MuInput))) - (4.181343324817995*MuInput*Log(1.02*msq2(2,2))*(0.96*msd2(2
      ,2) + 1.02*msq2(2,2) - 2.06*Sqr(MuInput))*Power3(Sqrt(Sqr(MuInput))))/(Abs(
      MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(
      MuInput))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput)))) + (4.059556626036888*
      MuInput*Sqrt(Sqr(MuInput))*(0.96*msd2(2,2) + 1.03*Sqr(MuInput))*Sqr(Log(0.96
      *msd2(2,2))))/(Abs(MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(0.96*msd2(2,2
      ) - 1.03*Sqr(MuInput))) + (4.181343324817995*MuInput*Power3(Sqrt(Sqr(MuInput
      )))*Sqr(Log(1.02*msq2(2,2))))/(Abs(MuInput)*(0.96*msd2(2,2) - 1.02*msq2(2,2)
      )*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput)))) + 3*(AtInput - MuInput/TanBeta)*(
      Log(1.02*msq2(2,2))*((8.119113252073776*MuInput*Sqrt(Sqr(MuInput)))/(Abs(
      MuInput)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (8.119113252073776*MuInput*Log
      (Sqr(SCALE))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*(1.02*msq2(2,2) - 0.98*msu2(2
      ,2)))) + (8.119113252073776*MuInput*Log(0.98*msu2(2,2))*Sqrt(Sqr(MuInput)))/
      (Abs(MuInput)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (8.119113252073776*
      MuInput*Log(0.98*msu2(2,2))*Log(Sqr(SCALE))*Sqrt(Sqr(MuInput)))/(Abs(MuInput
      )*(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (8.119113252073776*MuInput*PolyLog(2
      ,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(
      MuInput)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (8.119113252073776*MuInput*
      PolyLog(2,1 - (1.0510204081632653*Sqr(MuInput))/msu2(2,2))*Sqrt(Sqr(MuInput)
      ))/(Abs(MuInput)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (4.059556626036888*
      MuInput*Sqrt(Sqr(MuInput))*Sqr(Log(1.02*msq2(2,2))))/(Abs(MuInput)*(-1.02*
      msq2(2,2) + 0.98*msu2(2,2))) + (4.059556626036888*MuInput*Sqrt(Sqr(MuInput))
      *Sqr(Log(0.98*msu2(2,2))))/(Abs(MuInput)*(1.02*msq2(2,2) - 0.98*msu2(2,2))))
      + 3*(AbInput - MuInput*TanBeta)*Quad(AtInput - MuInput/TanBeta)*((-
      4.140747758557626*MuInput*msq2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqrt(
      Sqr(MuInput))*Sqr(Log(1.02*msq2(2,2))))/(Abs(MuInput)*Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2))*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 1.03*
      Sqr(MuInput))) + Log(0.96*msd2(2,2))*((3.8971743609954124*MuInput*Log(1.02*
      msq2(2,2))*msd2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqrt(Sqr(MuInput)))/(
      Abs(MuInput)*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) - (3.8971743609954124*
      MuInput*Log(0.98*msu2(2,2))*msd2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqrt
      (Sqr(MuInput)))/(Abs(MuInput)*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) - (
      7.794348721990825*MuInput*msd2(2,2)*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))*Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2)))) + Log(1.02*msq2(2,2))*((4.140747758557626*
      MuInput*Log(0.98*msu2(2,2))*msq2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqrt
      (Sqr(MuInput)))/(Abs(MuInput)*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.96*
      msd2(2,2) + 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + (
      8.281495517115252*MuInput*msq2(2,2)*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*(-0.96
      *msd2(2,2) + 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))*Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2)))) + Log(1.03*Sqr(MuInput))*((-4.181343324817995*
      MuInput*Log(1.02*msq2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Power3(Sqrt(
      Sqr(MuInput))))/(Abs(MuInput)*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*
      msd2(2,2) - 1.03*Sqr(MuInput))*(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) + (
      4.181343324817995*MuInput*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,
      2))*Power3(Sqrt(Sqr(MuInput))))/(Abs(MuInput)*Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))*(-1.02*msq2(2,2) + 1.03*Sqr(
      MuInput))) + (8.36268664963599*MuInput*Power3(Sqrt(Sqr(MuInput))))/(Abs(
      MuInput)*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))*(-1.02*msq2(2,2) + 1.03*Sqr(
      MuInput))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))))) + 3*Cube(AtInput - MuInput
      /TanBeta)*(Log(1.02*msq2(2,2))*((-8.119113252073776*MuInput*Log(Sqr(SCALE))*
      (1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) - (8.119113252073776*MuInput*(3.06*msq2(2,
      2) + 0.98*msu2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2)))) + (8.119113252073776*MuInput*Log(0.98*msu2(2,2))*(1.02*
      msq2(2,2) + 2.94*msu2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(1.02*msq2
      (2,2) - 0.98*msu2(2,2))) + (8.119113252073776*MuInput*PolyLog(2,1 - (
      1.0098039215686274*Sqr(MuInput))/msq2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)
      - 2.06*Sqr(MuInput))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2))) - (8.119113252073776*MuInput*PolyLog(2,1 - (
      1.0510204081632653*Sqr(MuInput))/msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)
      - 2.06*Sqr(MuInput))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + Log(1.03*Sqr(MuInput))*((16.72537329927198*MuInput*Log(
      1.02*msq2(2,2))*Power3(Sqrt(Sqr(MuInput))))/(Abs(MuInput)*Cube(1.02*msq2(2,2
      ) - 0.98*msu2(2,2))) + (16.72537329927198*MuInput*Log(0.98*msu2(2,2))*Power3
      (Sqrt(Sqr(MuInput))))/(Abs(MuInput)*Cube(-1.02*msq2(2,2) + 0.98*msu2(2,2))))
      + (4.059556626036888*MuInput*(1.02*msq2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(
      MuInput))*Sqrt(Sqr(MuInput))*Sqr(Log(1.02*msq2(2,2))))/(Abs(MuInput)*Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) - (4.059556626036888*MuInput*(1.02*msq2(2,
      2) + 0.98*msu2(2,2) - 2.06*Sqr(MuInput))*Sqrt(Sqr(MuInput))*Sqr(Log(0.98*
      msu2(2,2))))/(Abs(MuInput)*Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(Sqr(
      SCALE))*((8.119113252073776*MuInput*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) +
      0.98*msu2(2,2))*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*Cube(1.02*msq2(2,2) - 0.98
      *msu2(2,2))) + (16.23822650414755*MuInput*Sqrt(Sqr(MuInput)))/(Abs(MuInput)*
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + (32.4764530082951*MuInput*Sqrt(Sqr(
      MuInput)))/(Abs(MuInput)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))))/(Sqrt(1/(1
      + Sqr(TanBeta)))*Sqrt(Sqr(TanBeta)/(1 + Sqr(TanBeta)))) + Quad(AtInput -
      MuInput/TanBeta)*(3*((6*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 2.94*msu2(2,2)
      ))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(1.02*msq2(2,2))*((2*Log(0.98*
      msu2(2,2))*(-5.1*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) - (8*Log(Sqr(SCALE))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02
      *msq2(2,2) - 0.98*msu2(2,2)) - (2*(7.140000000000001*msq2(2,2) + 4.9*msu2(2,
      2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + ((7.140000000000001*msq2(2,2)
      + 0.98*msu2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Cube(1.02*msq2(2,2) - 0.98*msu2(
      2,2)) + Log(0.98*Sqr(mAInput))*((4*Log(1.02*msq2(2,2))*(1.02*msq2(2,2) +
      0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (4*Log(0.98*msu2(2,
      2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))
      - 8/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(Sqr(SCALE))*((8*Log(0.98*
      msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + 16/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + 24/Sqr(1.02*msq2(2,2
      ) - 0.98*msu2(2,2)) + (6*PolyLog(2,1 - (1.0408163265306123*msq2(2,2))/msu2(2
      ,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (3*Sqr(Log(0.98*msu2(2,2))))/
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + 3*Sqr(AbInput - MuInput*TanBeta)*((
      0.9803921568627451*(11.22*msq2(2,2) + 0.98*msu2(2,2)))/(Cube(1.02*msq2(2,2)
      - 0.98*msu2(2,2))*msq2(2,2)) + Log(Sqr(SCALE))*((0.9803921568627451*(5.1*
      msq2(2,2) + 0.98*msu2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2
      )) + (2*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 1.96*msu2(2,2)))/Quad(1.02*
      msq2(2,2) - 0.98*msu2(2,2))) + (2*Log(0.98*msu2(2,2))*(0.98*(-
      7.140000000000001*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2) + 0.96*msd2(2,2)*(
      1.02*msq2(2,2) + 4.9*msu2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + (3*(0.96*msd2(2,2) + 1.02*msq2(2,2) -
      1.96*msu2(2,2))*Sqr(Log(0.98*msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,
      2)) + (Sqr(Log(1.02*msq2(2,2)))*(-3*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 1.96*
      msu2(2,2)) + (2*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2)*(-5.1*msq2
      (2,2) + 0.98*msu2(2,2)) + 4.1616*Sqr(msq2(2,2))))/Sqr(0.96*msd2(2,2) - 1.02*
      msq2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (6*PolyLog(2,1 - (
      0.9411764705882352*msd2(2,2))/msq2(2,2))*Sqr(0.96*msd2(2,2) - 0.98*msu2(2,2)
      ))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)))
      + (6*PolyLog(2,1 - (0.9795918367346939*msd2(2,2))/msu2(2,2))*Sqr(0.96*msd2(2
      ,2) - 0.98*msu2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Quad(1.02*msq2(2,2
      ) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((1.92*Log(1.02*msq2(2,2))*msd2(2
      ,2)*(0.96*msd2(2,2) - 0.98*msu2(2,2))*(2.88*msd2(2,2) - 2.04*msq2(2,2) -
      0.98*msu2(2,2)))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2))) - (1.92*Log(0.98*msu2(2,2))*msd2(2,2)*(0.96*msd2(2,2) -
      0.98*msu2(2,2))*(2.88*msd2(2,2) - 2.04*msq2(2,2) - 0.98*msu2(2,2)))/(Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      0.9411764705882352*msd2(2,2))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) - (6*PolyLog(2,1 - (
      1.0408163265306123*msq2(2,2))/msu2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))
      *Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((-2*Log(Sqr(
      SCALE))*(1.02*msq2(2,2) + 1.96*msu2(2,2)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2
      ,2)) + (2*Log(0.98*msu2(2,2))*(4.896*msd2(2,2)*msq2(2,2) - 0.9408*msd2(2,2)*
      msu2(2,2) - 4.1616*Sqr(msq2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (-3.84*msd2(2,2)*(2.04*msq2(2,2) +
      0.98*msu2(2,2)) + 9.996*msq2(2,2)*msu2(2,2) + 5.202*Sqr(msq2(2,2)) -
      2.8811999999999998*Sqr(msu2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2)))))) + (1 + Sqr(TanBeta))*(3*Sqr(AtInput -
      MuInput/TanBeta)*((1.9607843137254901*(0.96*msd2(2,2) + 3.06*msq2(2,2) -
      2.06*Sqr(MuInput)))/(msq2(2,2)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + Log(Sqr
      (SCALE))*((1.9607843137254901*(0.96*msd2(2,2) + 2.04*msq2(2,2) - 2.06*Sqr(
      MuInput)))/(msq2(2,2)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))) - (2*Log(0.98*msu2
      (2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(MuInput)
      ))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((2*Log(Sqr(
      SCALE))*(0.96*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) - 2.06*Sqr(MuInput
      )))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*(0.96*msd2(2,2) + 3.06*msq2(2,
      2) - 2.06*Sqr(MuInput)) + ((-1.02*msq2(2,2) + 0.98*msu2(2,2))*(0.96*msd2(2,2
      ) + 1.03*Sqr(MuInput)))/(0.96*msd2(2,2) - 1.03*Sqr(MuInput)))/Sqr(1.02*msq2(
      2,2) - 0.98*msu2(2,2))) + (Log(0.98*msu2(2,2))*(-2*(0.96*msd2(2,2) + 3.06*
      msq2(2,2) - 2.06*Sqr(MuInput)) + ((1.02*msq2(2,2) - 0.98*msu2(2,2))*(4.8*
      msd2(2,2) - 3.09*Sqr(MuInput)))/(0.96*msd2(2,2) - 1.03*Sqr(MuInput))))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*PolyLog(2,1 - (1.0098039215686274*Sqr(
      MuInput))/msq2(2,2))*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput)))/Sqr(1.02*msq2(2,
      2) - 0.98*msu2(2,2)) - (4*PolyLog(2,1 - (1.0510204081632653*Sqr(MuInput))/
      msu2(2,2))*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput)))/Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + (2*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput))*Sqr(Log(1.02*msq2(2,2)
      )))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2*(-0.98*msu2(2,2) + 1.03*Sqr(
      MuInput))*Sqr(Log(0.98*msu2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      Log(0.96*msd2(2,2))*((1.92*msd2(2,2))/(-0.9996*msq2(2,2)*msu2(2,2) + 1.0404*
      Sqr(msq2(2,2))) + (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + (0.96*msd2(2,2)*(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 2.06*Sqr(MuInput)))/Sqr(
      0.96*msd2(2,2) - 1.03*Sqr(MuInput))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      (1.92*Log(1.02*msq2(2,2))*msd2(2,2)*(-1 + ((-1.02*msq2(2,2) + 0.98*msu2(2,2)
      )*(0.96*msd2(2,2) - 2.06*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.03*Sqr(
      MuInput))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.03*Sqr(MuInput))*(
      (-4.12*Sqr(MuInput))/(-0.9996*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2))) +
      (2.1218*Log(0.98*msu2(2,2))*Quad(MuInput))/((1.02*msq2(2,2) - 0.98*msu2(2,2)
      )*Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (2.1218*Log(1.02*msq2(2,2))*
      Quad(MuInput))/((-1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.03
      *Sqr(MuInput))))) + 3*Quad(AtInput - MuInput/TanBeta)*((-2*PolyLog(2,1 - (
      1.0510204081632653*Sqr(MuInput))/msu2(2,2))*(2.04*msq2(2,2) + 0.98*msu2(2,2)
      - 3.09*Sqr(MuInput))*(0.98*msu2(2,2) - 1.03*Sqr(MuInput)))/Quad(1.02*msq2(2,
      2) - 0.98*msu2(2,2)) + (2*PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/
      msq2(2,2))*(-0.98*msu2(2,2) + 1.03*Sqr(MuInput))*(-2.04*msq2(2,2) - 0.98*
      msu2(2,2) + 3.09*Sqr(MuInput)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + ((-
      0.98*msu2(2,2) + 1.03*Sqr(MuInput))*(-2.04*msq2(2,2) - 0.98*msu2(2,2) + 3.09
      *Sqr(MuInput))*Sqr(Log(1.02*msq2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2
      )) + ((0.98*msu2(2,2) - 1.03*Sqr(MuInput))*(-2.04*msq2(2,2) - 0.98*msu2(2,2)
      + 3.09*Sqr(MuInput))*Sqr(Log(0.98*msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + (0.9803921568627451*(0.96*msd2(2,2)*(4.08*msq2(2,2)*(1.02*msq2(
      2,2) + 1.96*msu2(2,2)) - 3.09*(7.140000000000001*msq2(2,2) + 0.98*msu2(2,2))
      *Sqr(MuInput)) + 2.06*Sqr(MuInput)*(-3.06*msq2(2,2)*(1.02*msq2(2,2) + 0.98*
      msu2(2,2)) + 1.03*(8.16*msq2(2,2) + 0.98*msu2(2,2))*Sqr(MuInput)) + 0.9216*(
      5.1*msq2(2,2) + 0.98*msu2(2,2))*Sqr(msd2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98
      *msu2(2,2))*msq2(2,2)*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (0.5*Log(0.98*
      msu2(2,2))*(0.96*msd2(2,2)*((1.02*msq2(2,2) + 0.98*msu2(2,2))*(1.02*msq2(2,2
      ) + 10.78*msu2(2,2)) - 12.36*(1.02*msq2(2,2) + 2.94*msu2(2,2))*Sqr(MuInput))
      + 1.03*Sqr(MuInput)*(-3*(1.02*msq2(2,2) + 0.98*msu2(2,2))*(1.02*msq2(2,2) +
      2.94*msu2(2,2)) + 4.12*(2.04*msq2(2,2) + 6.859999999999999*msu2(2,2))*Sqr(
      MuInput)) + 3.6864*(1.02*msq2(2,2) + 1.96*msu2(2,2))*Sqr(msd2(2,2))))/(Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + Log
      (Sqr(SCALE))*((0.9803921568627451*(3.06*msq2(2,2)*(1.02*msq2(2,2) + 0.98*
      msu2(2,2)) + 0.96*msd2(2,2)*(5.1*msq2(2,2) + 0.98*msu2(2,2)) - 2.06*(5.1*
      msq2(2,2) + 0.98*msu2(2,2))*Sqr(MuInput)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(
      2,2))*msq2(2,2)) + (Log(0.98*msu2(2,2))*(1.02*msq2(2,2)*(1.92*msd2(2,2) +
      1.02*msq2(2,2) - 4.12*Sqr(MuInput)) + 3.92*msu2(2,2)*(0.96*msd2(2,2) + 1.02*
      msq2(2,2) - 2.06*Sqr(MuInput)) + 0.9603999999999999*Sqr(msu2(2,2))))/Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*(-((Log(Sqr(SCALE))*
      (1.02*msq2(2,2)*(1.92*msd2(2,2) + 1.02*msq2(2,2) - 4.12*Sqr(MuInput)) + 3.92
      *msu2(2,2)*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 2.06*Sqr(MuInput)) +
      0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      (0.5*(-3.92*msu2(2,2)*(1.92*msd2(2,2) + 5.1*msq2(2,2) - 4.12*Sqr(MuInput))*(
      0.96*msd2(2,2) - 1.03*Sqr(MuInput)) + 1.02*msq2(2,2)*(5.15*(1.02*msq2(2,2) -
      4.12*Sqr(MuInput))*Sqr(MuInput) + 2.88*msd2(2,2)*(-1.02*msq2(2,2) + 8.24*Sqr
      (MuInput)) - 3.6864*Sqr(msd2(2,2))) - 0.9603999999999999*(0.96*msd2(2,2) +
      1.03*Sqr(MuInput))*Sqr(msu2(2,2))))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*(
      0.96*msd2(2,2) - 1.03*Sqr(MuInput)))) + Log(1.03*Sqr(MuInput))*((2.06*Sqr(
      MuInput)*(-3 + (0.9803921568627451*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(1.0506
      *msq2(2,2)*Sqr(MuInput) + Sqr(-1.03*Sqr(MuInput) + 0.96*msd2(2,2))))/(msq2(2
      ,2)*Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput)))))/Cube(-1.02*msq2(2,2) + 0.98*
      msu2(2,2)) + (1.0609*Log(1.02*msq2(2,2))*Quad(MuInput)*((1.02*msq2(2,2) -
      0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)) - 6*Sqr(-1.03*Sqr(MuInput)
      + 0.96*msd2(2,2))))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2
      ) - 1.03*Sqr(MuInput))) + (1.0609*Log(0.98*msu2(2,2))*Quad(MuInput)*(6*Sqr(-
      1.03*Sqr(MuInput) + 0.96*msd2(2,2)) - 1.0404*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput)))) + Log(0.96*msd2(2,2))*((0.96*msd2(
      2,2)*(6 - (0.9803921568627451*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-1.92*msd2(
      2,2)*(1.02*msq2(2,2) + 1.03*Sqr(MuInput)) + 1.03*Sqr(MuInput)*(4.08*msq2(2,2
      ) + 1.03*Sqr(MuInput)) + 0.9216*Sqr(msd2(2,2))))/(msq2(2,2)*Sqr(0.96*msd2(2,
      2) - 1.03*Sqr(MuInput)))))/Cube(-1.02*msq2(2,2) + 0.98*msu2(2,2)) + (0.96*
      Log(1.02*msq2(2,2))*msd2(2,2)*(3.92*msu2(2,2)*Sqr(-1.03*Sqr(MuInput) + 0.96*
      msd2(2,2)) + 1.02*msq2(2,2)*(1.02*msq2(2,2)*(0.96*msd2(2,2) - 2.06*Sqr(
      MuInput)) + 2*Sqr(-1.03*Sqr(MuInput) + 0.96*msd2(2,2))) - 0.9603999999999999
      *(0.96*msd2(2,2) - 2.06*Sqr(MuInput))*Sqr(msu2(2,2))))/(Quad(1.02*msq2(2,2)
      - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (0.96*Log(0.98*
      msu2(2,2))*msd2(2,2)*(1.02*msq2(2,2)*(-1.02*msq2(2,2)*(0.96*msd2(2,2) - 2.06
      *Sqr(MuInput)) - 2*Sqr(-1.03*Sqr(MuInput) + 0.96*msd2(2,2))) - 3.92*msu2(2,2
      )*Sqr(-1.03*Sqr(MuInput) + 0.96*msd2(2,2)) + 0.9603999999999999*(0.96*msd2(2
      ,2) - 2.06*Sqr(MuInput))*Sqr(msu2(2,2))))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2
      ,2))*Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput))))) + 3*(-1.5 + 2*PolyLog(2,1 -
      (1.0098039215686274*Sqr(MuInput))/msq2(2,2)) + Log(Sqr(SCALE))*(Log(0.98*
      msu2(2,2)) - (0.9803921568627451*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 2.06*Sqr
      (MuInput)))/msq2(2,2)) - (0.9803921568627451*(0.96*msd2(2,2) - 2.06*Sqr(
      MuInput)))/msq2(2,2) + (0.96*msd2(2,2))/(0.96*msd2(2,2) - 1.03*Sqr(MuInput))
      + Log(0.98*msu2(2,2))*(1.5 + (0.96*msd2(2,2))/(-0.96*msd2(2,2) + 1.03*Sqr(
      MuInput))) + Log(1.02*msq2(2,2))*(0.5 - Log(Sqr(SCALE)) + (0.96*msd2(2,2))/(
      -0.96*msd2(2,2) + 1.03*Sqr(MuInput))) + Sqr(Log(1.02*msq2(2,2))) + Log(0.96*
      msd2(2,2))*(0.96*msd2(2,2)*(0.9803921568627451/msq2(2,2) + (0.96*msd2(2,2) +
      2.06*Sqr(MuInput))/Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (2.1218*Log(
      1.03*Sqr(MuInput))*Quad(MuInput))/Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput)) -
      (0.96*Log(1.02*msq2(2,2))*msd2(2,2)*(0.96*msd2(2,2) - 2.06*Sqr(MuInput)))/
      Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput)) - (0.96*Log(0.98*msu2(2,2))*msd2(2,2
      )*(0.96*msd2(2,2) - 2.06*Sqr(MuInput)))/Sqr(0.96*msd2(2,2) - 1.03*Sqr(
      MuInput))) + Log(1.03*Sqr(MuInput))*((-2.019607843137255*Sqr(MuInput))/msq2(
      2,2) - (1.0609*Log(1.02*msq2(2,2))*Quad(MuInput))/Sqr(0.96*msd2(2,2) - 1.03*
      Sqr(MuInput)) - (1.0609*Log(0.98*msu2(2,2))*Quad(MuInput))/Sqr(0.96*msd2(2,2
      ) - 1.03*Sqr(MuInput)) - (1.03*Sqr(MuInput)*(1.92*msd2(2,2) + 1.03*Sqr(
      MuInput)))/Sqr(0.96*msd2(2,2) - 1.03*Sqr(MuInput))) + (2*PolyLog(2,1 - (
      1.0729166666666667*Sqr(MuInput))/msd2(2,2))*(-1.0609*Quad(MuInput) - 1.9776*
      msd2(2,2)*Sqr(MuInput) + 0.9216*Sqr(msd2(2,2))))/Sqr(0.96*msd2(2,2) - 1.03*
      Sqr(MuInput)) + (Sqr(Log(0.96*msd2(2,2)))*(-1.0609*Quad(MuInput) - 1.9776*
      msd2(2,2)*Sqr(MuInput) + 0.9216*Sqr(msd2(2,2))))/Sqr(0.96*msd2(2,2) - 1.03*
      Sqr(MuInput)))) + 3*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta)*
      ((2*Log(0.96*msd2(2,2))*Log(1.02*msq2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2
      )) + Log(0.98*Sqr(mAInput))*((2*Log(1.02*msq2(2,2)))/(0.96*msd2(2,2) - 1.02*
      msq2(2,2)) + (2*Log(0.96*msd2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2))) + (
      2*Sqr(Log(1.02*msq2(2,2))))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) - (
      2.0833333333333335*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(0.96*msd2
      (2,2) - 1.02*msq2(2,2))) + (1.9215686274509802*Sqr(mAInput)*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*
      msq2(2,2))) + 3*(-1 + 2*Log(0.96*msd2(2,2))*Log(0.98*msu2(2,2)) + (-2*Log(
      0.96*msd2(2,2)) - 4*Log(1.02*msq2(2,2)) - 2*Log(0.98*msu2(2,2)))*Log(0.98*
      Sqr(mAInput)) - 2*Log(Sqr(SCALE)) + Log(1.02*msq2(2,2))*(2 + 2*Log(Sqr(SCALE
      ))) + 1.3333333333333333*Sqr(3.141592653589793) + 4*Sqr(Log(0.98*Sqr(mAInput
      ))) - Sqr(Log(Sqr(SCALE))) + Sqr(Log(1.02*msq2(2,2))) - (1.9215686274509802*
      Sqr(mAInput)*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/msq2(2,2
      ) - (2.0833333333333335*(0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)
      )*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/msd2(2,2)) + Sqr(
      TanBeta)*(3*Sqr(AbInput + MuInput/TanBeta)*(-0.9803921568627451/msq2(2,2) -
      (0.9803921568627451*Log(Sqr(SCALE)))/msq2(2,2) - (Log(1.02*msq2(2,2))*(0.96*
      msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),0.96*msd2(2,2)) + (0.49019607843137253*Log(0.98*Sqr(mAInput))
      *(0.9603999999999999*Quad(mAInput) - Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(msq2(2,2)*TDelta(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) + (0.49019607843137253*Log
      (0.96*msd2(2,2))*(-0.9603999999999999*Quad(mAInput) + 1.9992*msq2(2,2)*Sqr(
      mAInput) + 0.9216*Sqr(msd2(2,2)) - 1.0404*Sqr(msq2(2,2)) + TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(msq2(2,2)*TDelta(0.98*Sqr(mAInput
      ),1.02*msq2(2,2),0.96*msd2(2,2))) + (0.5208333333333334*(-Sqr(0.98*Sqr(
      mAInput) + 0.96*msd2(2,2) - 1.02*msq2(2,2)) + TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2
      ,2)))/(msd2(2,2)*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) +
      3*(0.5 - 1.5*Log(0.98*msu2(2,2)) + Log(0.96*msd2(2,2))*Log(0.98*msu2(2,2)) +
      Log(1.02*msq2(2,2))*(0.5 + Log(Sqr(SCALE))) + 0.3333333333333333*Sqr(
      3.141592653589793) + Log(0.98*Sqr(mAInput))*(-Log(0.96*msd2(2,2)) - Log(1.02
      *msq2(2,2)) + (0.9803921568627451*(1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/msq2
      (2,2)) - (0.9607843137254901*Sqr(mAInput))/msq2(2,2) + Log(Sqr(SCALE))*(-Log
      (0.98*msu2(2,2)) - (0.9607843137254901*Sqr(mAInput))/msq2(2,2)) + Sqr(Log(
      0.98*Sqr(mAInput))) - (1.0416666666666667*(0.96*msd2(2,2) - 0.98*msu2(2,2) +
      0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/
      msd2(2,2)) + 3*(AtInput - MuInput/TanBeta)*(AbInput + MuInput/TanBeta)*(Log(
      0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      (2*Log(0.98*msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + Log(0.98*Sqr(
      mAInput))*((2*Log(1.02*msq2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*
      Log(0.98*msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + Log(1.02*msq2(2,2
      ))*(4/(-1.02*msq2(2,2) + 0.98*msu2(2,2)) + (4*Log(Sqr(SCALE)))/(-1.02*msq2(2
      ,2) + 0.98*msu2(2,2))) + (4*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2
      (2,2)) + (4*Log(0.98*msu2(2,2))*Log(Sqr(SCALE)))/(1.02*msq2(2,2) - 0.98*msu2
      (2,2)) + (2.0833333333333335*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(
      mAInput))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*
      (-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (2.0833333333333335*(0.96*msd2(2,2) -
      0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))/(msd2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + 3*(AbInput
       + MuInput/TanBeta)*Cube(AtInput - MuInput/TanBeta)*(Log(1.02*msq2(2,2))*((4
      *Log(Sqr(SCALE))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2)) + (4*(3.06*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) - (4*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 2.94*msu2(2,2)
      ))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(0.96*msd2(2,2))*((-2*Log(1.02
      *msq2(2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) - 1.96*Sqr(
      mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2))*(
      1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) - 1.96*Sqr(mAInput)))/Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*Sqr(mAInput))*((-2*Log(1.02*
      msq2(2,2))*(-1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) + 1.96*Sqr(
      mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2))*(-
      1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) + 1.96*Sqr(mAInput)))/Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(Sqr(SCALE))*((-4*Log(0.98*msu2(2,2))
      *(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) -
      8/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - 16/Sqr(1.02*msq2(2,2) - 0.98*msu2(
      2,2)) - (2.0833333333333335*((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,
      2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)) + 2*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2
      ,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msd2(2,2)) + (
      1.0416666666666667*(2*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*(0.96*msd2(2,2) -
      0.98*msu2(2,2) + 0.98*Sqr(mAInput)) + 4*TDelta(0.98*Sqr(mAInput),0.98*msu2(2
      ,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/
      (Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msd2(2,2))) + Quad(AtInput - MuInput/
      TanBeta)*(3*((0.9803921568627451*(1.02*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2
      (2,2)) + 0.98*(5.1*msq2(2,2) + 0.98*msu2(2,2))*Sqr(mAInput)))/(Cube(1.02*
      msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2)) + Log(Sqr(SCALE))*((
      0.9803921568627451*(2.04*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2)) + 0.98*
      (5.1*msq2(2,2) + 0.98*msu2(2,2))*Sqr(mAInput)))/(Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*msq2(2,2)) + (Log(0.98*msu2(2,2))*(1.9992*msq2(2,2)*Sqr(mAInput)
      + 3.8415999999999997*msu2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      Log(0.98*Sqr(mAInput))*((0.9803921568627451*(2.04*msq2(2,2)*(-1.02*msq2(2,2)
      + 0.98*msu2(2,2)) - 0.98*(5.1*msq2(2,2) + 0.98*msu2(2,2))*Sqr(mAInput)))/(
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2)) + (Log(1.02*msq2(2,2))*(
      1.9992*msq2(2,2)*Sqr(mAInput) + 3.8415999999999997*msu2(2,2)*Sqr(mAInput) +
      1.0404*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2
      ,2) - 0.98*msu2(2,2)) + (Log(0.98*msu2(2,2))*(-3.8415999999999997*msu2(2,2)*
      Sqr(mAInput) - 1.02*msq2(2,2)*(1.02*msq2(2,2) + 1.96*Sqr(mAInput)) +
      0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      Log(1.02*msq2(2,2))*((Log(Sqr(SCALE))*(-3.8415999999999997*msu2(2,2)*Sqr(
      mAInput) - 1.02*msq2(2,2)*(1.02*msq2(2,2) + 1.96*Sqr(mAInput)) +
      0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      (0.5*(-7.683199999999999*msu2(2,2)*Sqr(mAInput) - 1.02*msq2(2,2)*(1.02*msq2(
      2,2) + 3.92*Sqr(mAInput)) + 0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*
      msq2(2,2) - 0.98*msu2(2,2))) + (0.5*Log(0.98*msu2(2,2))*(3.9984*msq2(2,2)*
      Sqr(mAInput) + 7.683199999999999*msu2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,
      2)) - 0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,
      2))) + 3*Sqr(AbInput + MuInput/TanBeta)*((0.9803921568627451*(11.22*msq2(2,2
      ) + 0.98*msu2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2)) + Log
      (Sqr(SCALE))*((0.9803921568627451*(5.1*msq2(2,2) + 0.98*msu2(2,2)))/(Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2)) + (2*Log(0.98*msu2(2,2))*(1.02*
      msq2(2,2) + 1.96*msu2(2,2)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*Log
      (0.98*msu2(2,2))*(1.02*msq2(2,2) + 4.9*msu2(2,2)))/Quad(1.02*msq2(2,2) -
      0.98*msu2(2,2)) + Log(1.02*msq2(2,2))*((-4*(2.04*msq2(2,2) + 0.98*msu2(2,2))
      )/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2*Log(Sqr(SCALE))*(1.02*msq2(2,2)
      + 1.96*msu2(2,2)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (0.96*msd2(2,2) -
      1.02*msq2(2,2) + 0.98*Sqr(mAInput))/(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + Log(0.98*Sqr(
      mAInput))*((Log(1.02*msq2(2,2))*(-2.88*msd2(2,2) + 1.02*msq2(2,2) + 1.96*
      msu2(2,2) + 2.94*Sqr(mAInput)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (Log
      (0.98*msu2(2,2))*(-2.88*msd2(2,2) + 1.02*msq2(2,2) + 1.96*msu2(2,2) + 2.94*
      Sqr(mAInput)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (0.49019607843137253*
      (0.9603999999999999*Quad(mAInput) - Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) +
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(msq2(2,2)*Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),
      0.96*msd2(2,2)))) + Log(0.96*msd2(2,2))*((Log(1.02*msq2(2,2))*(2.88*msd2(2,2
      ) + 1.02*msq2(2,2) + 1.96*msu2(2,2) - 2.94*Sqr(mAInput)))/Quad(1.02*msq2(2,2
      ) - 0.98*msu2(2,2)) - (Log(0.98*msu2(2,2))*(2.88*msd2(2,2) + 1.02*msq2(2,2)
      + 1.96*msu2(2,2) - 2.94*Sqr(mAInput)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))
      + (0.49019607843137253*(-0.9603999999999999*Quad(mAInput) + 1.9992*msq2(2,2)
      *Sqr(mAInput) + 0.9216*Sqr(msd2(2,2)) - 1.0404*Sqr(msq2(2,2)) + TDelta(0.98*
      Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))/(msq2(2,2)*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))))
      + (0.5208333333333334*(-(Sqr(0.98*Sqr(mAInput) + 0.96*msd2(2,2) - 1.02*msq2(
      2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + 6*Sqr(TDelta(0.98*Sqr(mAInput)
      ,1.02*msq2(2,2),0.96*msd2(2,2))) + (1.02*msq2(2,2) - 0.98*msu2(2,2))*(3.84*
      msd2(2,2) - 3.06*msq2(2,2) - 0.98*msu2(2,2) + 3.92*Sqr(mAInput))*TDelta(0.98
      *Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) + (
      1.0416666666666667*((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98
      *msu2(2,2) + 0.98*Sqr(mAInput)) - 3*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(
      msd2(2,2)*Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))))) + Sqr(AtInput - MuInput/
      TanBeta)*(3*Sqr(AbInput + MuInput/TanBeta)*(-2/(-0.9996*msq2(2,2)*msu2(2,2)
      + 1.0404*Sqr(msq2(2,2))) + Log(Sqr(SCALE))*(-2/(-0.9996*msq2(2,2)*msu2(2,2)
      + 1.0404*Sqr(msq2(2,2))) - (2*Log(0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98
      *msu2(2,2))) - (2*Log(0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))
      + Log(1.02*msq2(2,2))*((2*Log(Sqr(SCALE)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,
      2)) + (2*((-1.02*msq2(2,2) + 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2
      ) + 0.98*Sqr(mAInput)) + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2
      ,2))))/(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2)))) + Log(0.98*Sqr(mAInput))*(-(Log(1.02*msq2(2,2))/
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*msu2(2,2))/Sqr(1.02*msq2(2,
      2) - 0.98*msu2(2,2)) + (0.9803921568627451*(0.9603999999999999*Quad(mAInput)
      - Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + TDelta(0.98*Sqr(mAInput),1.02*msq2(
      2,2),0.96*msd2(2,2))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + Log(0.96*msd2(2,2))*(-(
      Log(1.02*msq2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*msu2(2,
      2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (0.9803921568627451*(-
      0.9603999999999999*Quad(mAInput) + 1.9992*msq2(2,2)*Sqr(mAInput) + 0.9216*
      Sqr(msd2(2,2)) - 1.0404*Sqr(msq2(2,2)) + TDelta(0.98*Sqr(mAInput),1.02*msq2(
      2,2),0.96*msd2(2,2))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))) + (1.0416666666666667*(-(
      (1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.98*Sqr(mAInput) + 0.96*msd2(2,2) -
      1.02*msq2(2,2))) + (0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))) - (
      1.0416666666666667*(0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2)))) + 3*((0.9803921568627451*(4.08*msq2(2,2) -
      1.96*Sqr(mAInput)))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96
      *msd2(2,2))*((Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*
      Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (Log(0.98*msu2(2,2))*(
      0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + Log(0.98*Sqr(mAInput))*((1.96*Sqr(mAInput))/(-0.9996*msq2
      (2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2))) - (Log(1.02*msq2(2,2))*(0.96*msd2(2
      ,2) + 1.02*msq2(2,2) - 1.96*msu2(2,2) + 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,
      2) - 0.98*msu2(2,2)) + (Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2)
      - 1.96*msu2(2,2) + 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))
      + Log(Sqr(SCALE))*((1.9607843137254901*(1.02*msq2(2,2) - 0.98*Sqr(mAInput)))
      /(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2*Log(0.98*msu2(2,2))*(-
      0.98*msu2(2,2) + 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      Log(1.02*msq2(2,2))*(-((3.06*msq2(2,2) + 0.98*msu2(2,2) - 1.96*Sqr(mAInput))
      /Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*Log(Sqr(SCALE))*(-0.98*msu2(2,2)
      + 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (Log(0.98*msu2
      (2,2))*(1.02*msq2(2,2) + 2.94*msu2(2,2) - 1.96*Sqr(mAInput)))/Sqr(1.02*msq2(
      2,2) - 0.98*msu2(2,2)) + (1.0416666666666667*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,
      2)))/(msd2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (1.0416666666666667*
      ((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*
      Sqr(mAInput)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2)))))) + 3*(AbInput - MuInput*TanBeta)*(AtInput +
      MuInput*TanBeta)*((2*Log(0.96*msd2(2,2))*Log(0.98*msu2(2,2)))/(0.96*msd2(2,2
      ) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2)))/(-0.96*
      msd2(2,2) + 1.02*msq2(2,2)) + Log(0.98*Sqr(mAInput))*((2*Log(1.02*msq2(2,2))
      )/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(0.96*msd2(2,2)))/(-0.96*msd2(2,
      2) + 1.02*msq2(2,2))) - (2.0833333333333335*(0.96*msd2(2,2) - 0.98*msu2(2,2)
      + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/
      (msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (1.9607843137254901*(1.02*
      msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2))) +
      Sqr(AtInput - MuInput/TanBeta)*(3*Sqr(AbInput - MuInput*TanBeta)*((6*(-1.92*
      msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2))*PolyLog(2,1 - (
      1.0408163265306123*msq2(2,2))/msu2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))
      *Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - 2/(-0.9996*msq2(2,2)*msu2(2,2) +
      1.0404*Sqr(msq2(2,2))) + Log(Sqr(SCALE))*(-2/(-0.9996*msq2(2,2)*msu2(2,2) +
      1.0404*Sqr(msq2(2,2))) - (2*Log(0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) + Log(0.96*msd2(2,2))*((1.8823529411764703*msd2(2,2))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      3.84*Log(1.02*msq2(2,2))*msd2(2,2)*(0.96*msd2(2,2) - 0.98*msu2(2,2)))/(Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      3.84*Log(0.98*msu2(2,2))*msd2(2,2)*(-0.96*msd2(2,2) + 0.98*msu2(2,2)))/(Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) +
      Log(1.02*msq2(2,2))*(2/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(Sqr(
      SCALE)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (Log(0.98*msu2(2,2))*(7.8336
      *msd2(2,2)*msq2(2,2) - 3.7632*msd2(2,2)*msu2(2,2) - 4.1616*Sqr(msq2(2,2))))/
      (Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))))
      - (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98*msu2(2,2)))/((0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (6*(0.96*msd2(2,2)
      + 1.02*msq2(2,2) - 1.96*msu2(2,2))*PolyLog(2,1 - (0.9411764705882352*msd2(2,
      2))/msq2(2,2)))/((-0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (3*Sqr(Log(0.98*msu2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + (Sqr(Log(1.02*msq2(2,2)))*(-1.9584*msd2(2,2)*msq2(2,2) + 3.7632
      *msd2(2,2)*msu2(2,2) - 2.7648*Sqr(msd2(2,2)) + 1.0404*Sqr(msq2(2,2))))/(Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (6*
      PolyLog(2,1 - (0.9795918367346939*msd2(2,2))/msu2(2,2))*Sqr(0.96*msd2(2,2) -
      0.98*msu2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2)))) + 3*(8/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (12*PolyLog(2,1
       - (1.0408163265306123*msq2(2,2))/msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,
      2)) + (6*Sqr(Log(0.98*msu2(2,2))))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(
      0.98*Sqr(mAInput))*((-2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 2.04*msq2(2,2)
      - 2.94*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2
      ,2))*(0.96*msd2(2,2) + 2.04*msq2(2,2) - 2.94*msu2(2,2)))/Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + Log(Sqr(SCALE))*(4/(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      (Log(0.98*msu2(2,2))*(-8.16*msq2(2,2) + 11.76*msu2(2,2)))/Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2))*(0.96*msd2(
      2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2
      ,2)) + (2*Log(0.98*msu2(2,2))*(-0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(
      mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((2*(
      1.02*msq2(2,2) - 4.9*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*
      Log(Sqr(SCALE))*(2.04*msq2(2,2) - 2.94*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98
      *msu2(2,2)) - (2*Log(0.98*msu2(2,2))*(3.06*msq2(2,2) - 2.94*msu2(2,2) + 0.98
      *Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (Log(0.98*msu2(2,2))
      *(-6.12*msq2(2,2) + 13.719999999999999*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98
      *msu2(2,2)) + (1.96*Sqr(mAInput)*Sqr(Log(1.02*msq2(2,2))))/Sqr(1.02*msq2(2,2
      ) - 0.98*msu2(2,2)) + (2.0833333333333335*TDelta(0.98*Sqr(mAInput),1.02*msq2
      (2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))
      /(msd2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (1.9607843137254901*(
      0.98*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput) + TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,
      2),1.02*msq2(2,2)))/(msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      2.0833333333333335*((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98
      *msu2(2,2) + 0.98*Sqr(mAInput)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(
      msd2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (1.9607843137254901*TDelta
      (0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,
      2)))) + 3*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta)*(Log(0.98*
      Sqr(mAInput))*((2*Log(1.02*msq2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))
      - (2*Log(0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*
      msd2(2,2))*((-2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98*
      Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) + (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98
      *Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98
      *msu2(2,2)))) - (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) +
      0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,
      2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2.0833333333333335*((1.02*msq2(
      2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput))
      + TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(0.96*msd2(2,2) - 1.02*
      msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (1.9607843137254901*(0.98
      *(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput) + TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*
      msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + (2.0833333333333335*TDelta(0.98*Sqr(mAInput),0.98*msu2(
      2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/
      (msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(
      2,2))) - (1.9607843137254901*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))))
      + 3*(AbInput - MuInput*TanBeta)*(AtInput + MuInput*TanBeta)*(Log(0.98*Sqr(
      mAInput))*((-2*Log(1.02*msq2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (
      2*Log(0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2
      (2,2))*((2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(
      mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2
      (2,2))) + (2*Log(0.98*msu2(2,2))*(-0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*
      Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2)))) + (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(-2.04*msq2(2,2) +
      0.98*Sqr(mAInput)))/((-0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2*(-2.04*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*
      msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) + (2.0833333333333335*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),
      0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2
      (2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))
      ) - (1.9607843137254901*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,
      2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2)
      - 1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      2.0833333333333335*((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98
      *msu2(2,2) + 0.98*Sqr(mAInput)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(
      msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2
      ,2))) + (0.9803921568627451*(-2*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(1.02*msq2
      (2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)) + 2*TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))))) + (AtInput - MuInput/TanBeta)*(3*(AbInput + MuInput/
      TanBeta)*(Log(0.96*msd2(2,2))*((2*Log(1.02*msq2(2,2)))/(1.02*msq2(2,2) -
      0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))
      ) + Log(0.98*Sqr(mAInput))*((2*Log(1.02*msq2(2,2)))/(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + (2*Log(0.98*msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) +
      Log(1.02*msq2(2,2))*(4/(-1.02*msq2(2,2) + 0.98*msu2(2,2)) + (4*Log(Sqr(SCALE
      )))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (4*Log(0.98*msu2(2,2)))/(1.02*msq2
      (2,2) - 0.98*msu2(2,2)) + (4*Log(0.98*msu2(2,2))*Log(Sqr(SCALE)))/(1.02*msq2
      (2,2) - 0.98*msu2(2,2)) + (2.0833333333333335*(0.96*msd2(2,2) - 1.02*msq2(2,
      2) + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)
      ))/(msd2(2,2)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (2.0833333333333335*(
      0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),
      0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))
      )) + 3*(AbInput - MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((2*Log(0.96*msd2
      (2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) + 2*Log(1.02*msq2(2,2))*(1/(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + 1/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (2*Log
      (0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*(
      (2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput))
      )/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2
      *Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/
      ((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + (2*
      Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2) -
      0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2*(-2.04*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*
      msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(
      2,2))) + (2*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*
      msu2(2,2))*(-0.9792*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)))) - (
      1.9607843137254901*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) -
      1.02*msq2(2,2))*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2*TDelta(
      0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),0.98
      *msu2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.9792*msd2
      (2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)))) + (1.9607843137254901*TDelta(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2
      (2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*(1.02*
      msq2(2,2) - 0.98*msu2(2,2)))) + (AtInput + MuInput*TanBeta)*(3*(Log(0.98*Sqr
      (mAInput))*((2*Log(1.02*msq2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*
      Log(0.98*msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + Log(1.02*msq2(2,2
      ))*(4/(-1.02*msq2(2,2) + 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2)))/(-1.02*
      msq2(2,2) + 0.98*msu2(2,2)) + (4*Log(Sqr(SCALE)))/(-1.02*msq2(2,2) + 0.98*
      msu2(2,2))) + (4*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4
      *Log(0.98*msu2(2,2))*Log(Sqr(SCALE)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2
      *Sqr(Log(1.02*msq2(2,2))))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (1.96*Sqr(
      mAInput)*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(-0.9996*
      msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2))) + (1.9607843137254901*(1.02*
      msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2)))) +
      3*(AbInput + MuInput/TanBeta)*(AbInput - MuInput*TanBeta)*(Log(0.96*msd2(2,2
      ))*((2*Log(1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,
      2) - 0.98*msu2(2,2))) + (2*Log(0.98*msu2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2
      (2,2))*(-1.02*msq2(2,2) + 0.98*msu2(2,2)))) + (2*Log(1.02*msq2(2,2))*Log(
      0.98*msu2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) + (2*Sqr(Log(1.02*msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2)
      )*(-1.02*msq2(2,2) + 0.98*msu2(2,2))) - (2*(0.96*msd2(2,2) - 1.02*msq2(2,2)
      + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/
      ((1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.9792*msd2(2,2)*msq2(2,2) + 0.9216*Sqr
      (msd2(2,2)))) + (1.9215686274509802*Sqr(mAInput)*TPhi(0.98*Sqr(mAInput),1.02
      *msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2.0833333333333335*(0.96*msd2(2,2) -
      0.98*msu2(2,2) + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))/(msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2
      ) - 0.98*msu2(2,2))) + (1.9607843137254901*(1.02*msq2(2,2) - 0.98*msu2(2,2)
      + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/
      (msq2(2,2)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,
      2)))))) + Cube(AtInput - MuInput/TanBeta)*(3*(AbInput + MuInput/TanBeta)*(
      Log(1.02*msq2(2,2))*((4*Log(Sqr(SCALE))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*(3.06*msq2(2,2) + 0.98*msu2(2,2))
      )/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (4*Log(0.98*msu2(2,2))*(1.02*msq2
      (2,2) + 2.94*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(0.96*
      msd2(2,2))*((-2*Log(1.02*msq2(2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*
      msu2(2,2) - 1.96*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*
      Log(0.98*msu2(2,2))*(1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) - 1.96
      *Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*Sqr(
      mAInput))*((-2*Log(1.02*msq2(2,2))*(-1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*
      msu2(2,2) + 1.96*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*
      Log(0.98*msu2(2,2))*(-1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) +
      1.96*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(Sqr(SCALE))
      *((-4*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(
      2,2) - 0.98*msu2(2,2)) - 8/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - 16/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2.0833333333333335*((1.02*msq2(2,2) -
      0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)) + 2*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,
      2))*msd2(2,2)) + (1.0416666666666667*(2*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*(
      0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)) + 4*TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,
      2),0.96*msd2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msd2(2,2))) + (
      AtInput + MuInput*TanBeta)*(3*(AbInput + MuInput/TanBeta)*(AbInput - MuInput
      *TanBeta)*(((4*Log(1.02*msq2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      (4*Log(0.98*msu2(2,2)))/Cube(-1.02*msq2(2,2) + 0.98*msu2(2,2)))*Log(0.98*Sqr
      (mAInput)) + Log(0.96*msd2(2,2))*((-2*Log(1.02*msq2(2,2))*(1.92*msd2(2,2) +
      1.02*msq2(2,2) + 0.98*msu2(2,2) - 1.96*Sqr(mAInput)))/(Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (2*Log(0.98*msu2(2,2))*
      (1.92*msd2(2,2) + 1.02*msq2(2,2) + 0.98*msu2(2,2) - 1.96*Sqr(mAInput)))/(
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2)))) -
      (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(3.06*msq2(2,2) + 0.98*msu2(2,2)
      - 1.96*Sqr(mAInput)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (2*(3.06*msq2(2,2) + 0.98*msu2(2,2) - 1.96*Sqr(mAInput)
      )*Sqr(Log(1.02*msq2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) - (2.0833333333333335*((1.02*msq2(2,2) - 0.98*
      msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)) + 2*TDelta(
      0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),
      1.02*msq2(2,2),0.96*msd2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msd2(
      2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (1.9607843137254901*(0.98*(1.02*
      msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput) + 2*TDelta(0.98*Sqr(mAInput),1.02*
      msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2
      ,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2
      ))*msq2(2,2)) + (1.0416666666666667*(2*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*(
      0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)) + 4*TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,
      2),0.96*msd2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msd2(2,2)*(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + (1.9607843137254901*((1.02*msq2(2,2) - 0.98*
      msu2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)) - 2*TDelta(
      0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96
      *msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2))) + 3*((-4*Log(0.98*msu2(2,2))*(1.02*
      msq2(2,2) + 2.94*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(
      1.02*msq2(2,2))*((4*Log(Sqr(SCALE))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*(3.06*msq2(2,2) + 0.98*msu2(2,2)))/
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2))*(3.06*msq2(2,
      2) + 0.98*msu2(2,2) - 1.96*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,
      2))) + Log(0.98*Sqr(mAInput))*((2*Log(1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98
      *msu2(2,2) - 1.96*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*
      Log(0.98*msu2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2(2,2) + 1.96*Sqr(mAInput)))/
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2*(3.06*msq2(2,2) + 0.98*msu2(2,2)
      - 1.96*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + Log(Sqr(SCALE))*((-4*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98
      *msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) - 8/Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) - 16/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (
      1.9607843137254901*(0.98*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput) + 2*
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,
      2))*msq2(2,2)) + (0.9803921568627451*(-2*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(
      1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)) + 4*TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,
      2),1.02*msq2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2))))) + (
      (1 + Sqr(TanBeta))*(3*(2*Log(0.98*msu2(2,2))*Log(Sqr(SCALE)) - 4*PolyLog(2,1
       - (1.0098039215686274*Sqr(MuInput))/msq2(2,2)) - 4*PolyLog(2,1 - (
      1.0510204081632653*Sqr(MuInput))/msu2(2,2)) + 4.12*Log(1.03*Sqr(MuInput))*(1
      /(1.02*msq2(2,2) - 1.03*Sqr(MuInput)) + 1/(0.98*msu2(2,2) - 1.03*Sqr(MuInput
      )))*Sqr(MuInput) + (4.12*Log(0.98*msu2(2,2))*Sqr(MuInput))/(-0.98*msu2(2,2)
      + 1.03*Sqr(MuInput)) + Log(1.02*msq2(2,2))*(2*Log(0.98*msu2(2,2)) + 2*Log(
      Sqr(SCALE)) + (4.12*Sqr(MuInput))/(-1.02*msq2(2,2) + 1.03*Sqr(MuInput))) - 2
      *Sqr(Log(Sqr(SCALE))) - 2*Sqr(Log(1.02*msq2(2,2))) - 2*Sqr(Log(0.98*msu2(2,2
      )))) + 3*Sqr(AbInput - MuInput*TanBeta)*(4/(-0.96*msd2(2,2) + 1.02*msq2(2,2)
      ) + (2*Log(0.98*msu2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) + (2*Log(Sqr(
      SCALE)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) + (4.12*Log(1.03*Sqr(MuInput))*
      Sqr(MuInput))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-1.02*msq2(2,2) + 1.03*Sqr
      (MuInput))) + Log(0.96*msd2(2,2))*((3.84*msd2(2,2))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + (1.92*Log(0.98*msu2(2,2))*msd2(2,2))/Sqr(0.96*msd2(2,2) -
      1.02*msq2(2,2)) + (1.92*Log(Sqr(SCALE))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02
      *msq2(2,2))) - (3.84*msd2(2,2)*PolyLog(2,1 - (1.0729166666666667*Sqr(MuInput
      ))/msd2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (3.84*msd2(2,2)*
      PolyLog(2,1 - (1.0098039215686274*Sqr(MuInput))/msq2(2,2)))/Sqr(0.96*msd2(2,
      2) - 1.02*msq2(2,2)) - (1.92*msd2(2,2)*Sqr(Log(0.96*msd2(2,2))))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (1.92*msd2(2,2)*Sqr(Log(1.02*msq2(2,2))))/Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(1.02*msq2(2,2))*((-1.92*Log(0.98*msu2
      (2,2))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (1.92*Log(Sqr(SCALE
      ))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (4*(0.9888*msd2(2,2)*
      Sqr(MuInput) - 1.0404*Sqr(msq2(2,2))))/((1.02*msq2(2,2) - 1.03*Sqr(MuInput))
      *Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))))) + Sqr(AtInput - MuInput/TanBeta)*(3
      *Sqr(AbInput - MuInput*TanBeta)*(2/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*
      msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((1.92*msd2(2,2))/((-1.02
      *msq2(2,2) + 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (1.9584
      *Log(1.02*msq2(2,2))*msd2(2,2)*msq2(2,2))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,
      2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (1.9584*Log(0.98*msu2(2,2))*msd2
      (2,2)*msq2(2,2))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2)))) + Log(1.02*msq2(2,2))*((1.9584*Log(0.98*msu2(2,2))*msd2(2,
      2)*msq2(2,2))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2*(-0.9408*msd2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2))))
      /(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))
      ) + (2.04*Log(0.98*msu2(2,2))*msq2(2,2))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (1.9584*msd2(2,2)*msq2(2,2)*Sqr(Log(
      1.02*msq2(2,2))))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2)))) + 3*(4/(-1.02*msq2(2,2) + 0.98*msu2(2,2)) + (4.12*Log(1.03
      *Sqr(MuInput))*Sqr(MuInput))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.98*msu2(
      2,2) + 1.03*Sqr(MuInput))) + Log(Sqr(SCALE))*(2/(-1.02*msq2(2,2) + 0.98*msu2
      (2,2)) - (2.04*Log(0.98*msu2(2,2))*msq2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2
      (2,2))) + Log(1.02*msq2(2,2))*((-2.04*Log(0.98*msu2(2,2))*msq2(2,2))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2.04*Log(Sqr(SCALE))*msq2(2,2))/Sqr(1.02
      *msq2(2,2) - 0.98*msu2(2,2)) + (2*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) - (4.08*msq2(2,2)*PolyLog(2,1 - (
      1.0098039215686274*Sqr(MuInput))/msq2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(
      2,2)) + (4.08*msq2(2,2)*PolyLog(2,1 - (1.0510204081632653*Sqr(MuInput))/msu2
      (2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2.04*msq2(2,2)*Sqr(Log(0.98*
      msu2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (Log(0.98*msu2(2,2))*(-
      4.2024*msq2(2,2)*Sqr(MuInput) + 3.8415999999999997*Sqr(msu2(2,2))))/((-0.98*
      msu2(2,2) + 1.03*Sqr(MuInput))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))))) + (3*
      Sqr(AbInput - MuInput*TanBeta)*(4/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log
      (Sqr(SCALE)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(1.02*msq2(2,2))*((1.92
      *Log(Sqr(SCALE))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*(0.96*
      msd2(2,2) + 1.02*msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(
      0.96*msd2(2,2))*((-3.84*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) - (
      1.92*Log(Sqr(SCALE))*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (Log(
      1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput))*((Log(0.96*msd2(2
      ,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,
      2) - 1.02*msq2(2,2)) - (Log(1.02*msq2(2,2))*(0.96*msd2(2,2) - 1.02*msq2(2,2)
      + 0.98*Sqr(mAInput)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - ((0.96*msd2(2,
      2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2)) + (1.0416666666666667*TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*
      msd2(2,2)))/(msd2(2,2)*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      0.9611687812379854*(0.98*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-5.1*msq2(2,2)
      + 0.98*Sqr(mAInput))*Sqr(mAInput) + (0.96*msd2(2,2) - 2.04*msq2(2,2))*TDelta
      (0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      msq2(2,2)))) + 3*(-2*Log(0.98*msu2(2,2)) + Log(0.96*msd2(2,2))*Log(0.98*msu2
      (2,2)) + (4 - Log(0.96*msd2(2,2)) + Log(0.98*msu2(2,2)))*Log(0.98*Sqr(
      mAInput)) + Log(1.02*msq2(2,2))*(-2 - 2*Log(Sqr(SCALE))) - 2*Log(0.98*msu2(2
      ,2))*Log(Sqr(SCALE)) + 2*Sqr(Log(Sqr(SCALE))) + Sqr(Log(1.02*msq2(2,2))) + (
      0.9611687812379854*(0.9603999999999999*Quad(mAInput) - 4.998*msq2(2,2)*Sqr(
      mAInput) - TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(
      0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/Sqr(msq2(2,2)) - (
      1.0416666666666667*(0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/msd2(2,2)) + 3*(
      AbInput - MuInput*TanBeta)*(AtInput + MuInput*TanBeta)*((2*Log(0.96*msd2(2,2
      ))*Log(0.98*msu2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02*msq2
      (2,2))*Log(0.98*msu2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) + Log(0.98*
      Sqr(mAInput))*((2*Log(1.02*msq2(2,2)))/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (
      2*Log(0.96*msd2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2))) - (
      2.0833333333333335*(0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(0.96*msd2
      (2,2) - 1.02*msq2(2,2))) + (1.9607843137254901*(1.02*msq2(2,2) - 0.98*msu2(2
      ,2) + 0.98*Sqr(mAInput))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2
      )))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2))) + Sqr(AtInput + MuInput*
      TanBeta)*(3*((-2*Log(1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98
      *Sqr(mAInput)))/TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)) + (
      0.9803921568627451*Log(0.98*Sqr(mAInput))*(Sqr(0.98*Sqr(mAInput) + 1.02*msq2
      (2,2) - 0.98*msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(
      2,2))))/(msq2(2,2)*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))
      + (0.9803921568627451*Log(0.98*msu2(2,2))*(-0.9603999999999999*Quad(mAInput)
      + 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2)) + TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2))))/(msq2(2,2)*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2))) + (0.9611687812379854*(-0.98*msu2(2,2) + 0.98*Sqr(mAInput) + ((
      1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*(0.98*(1.02*msq2(2,2) -
      0.98*msu2(2,2))*msu2(2,2) - 0.9603999999999999*Quad(mAInput) + 0.98*(1.02*
      msq2(2,2) + 1.96*msu2(2,2))*Sqr(mAInput)))/TDelta(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2
      ,2)))/Sqr(msq2(2,2))) + 3*Sqr(AbInput - MuInput*TanBeta)*((Log(0.96*msd2(2,2
      ))*Log(0.98*msu2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + Log(1.02*msq2
      (2,2))*(-(Log(0.98*msu2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (2*(
      1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)))/((0.96*msd2(2,2) -
      1.02*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) +
      (0.9803921568627451*Log(0.98*msu2(2,2))*(-0.9603999999999999*Quad(mAInput) +
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2)) + TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2))))/(msq2(2,2)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*TDelta(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput))*(-(Log
      (0.96*msd2(2,2))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(1.02*msq2(2,2))
      /Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (0.9803921568627451*(-Sqr(0.98*Sqr(
      mAInput) + 1.02*msq2(2,2) - 0.98*msu2(2,2)) + TDelta(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) - (
      1.0416666666666667*(0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + (0.9611687812379854*((0.96*msd2(2,2) - 1.02*
      msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*(0.98*(-
      1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2) + 0.9603999999999999*Quad(mAInput
      ) - 0.98*(1.02*msq2(2,2) + 1.96*msu2(2,2))*Sqr(mAInput)) + (0.9408*msd2(2,2)
      *msu2(2,2) - 1.9992*msq2(2,2)*msu2(2,2) - 0.9408*msd2(2,2)*Sqr(mAInput) +
      1.9992*msq2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)))*TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,
      2),1.02*msq2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))) + (AtInput -
      MuInput/TanBeta)*(3*(AbInput - MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((2*
      Log(0.96*msd2(2,2)))/(-0.96*msd2(2,2) + 1.02*msq2(2,2)) + 2*Log(1.02*msq2(2,
      2))*(1/(0.96*msd2(2,2) - 1.02*msq2(2,2)) + 1/(-1.02*msq2(2,2) + 0.98*msu2(2,
      2))) + (2*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96
      *msd2(2,2))*((2*Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*
      Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))) - (2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98
      *Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*
      msu2(2,2)))) + (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) +
      0.98*msu2(2,2) - 0.98*Sqr(mAInput)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*(-2.04*msq2(2,2) + 0.98*Sqr(mAInput))
      *Sqr(Log(1.02*msq2(2,2))))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2
      ) - 0.98*msu2(2,2))) + (2*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(
      2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2
      ) - 0.98*msu2(2,2))*(-0.9792*msd2(2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)))) -
      (1.9607843137254901*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) -
      1.02*msq2(2,2))*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (2*TDelta(
      0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),0.98
      *msu2(2,2),0.96*msd2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.9792*msd2
      (2,2)*msq2(2,2) + 0.9216*Sqr(msd2(2,2)))) + (1.9607843137254901*TDelta(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*TPhi(0.98*Sqr(mAInput),0.98*msu2
      (2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*(1.02*
      msq2(2,2) - 0.98*msu2(2,2)))) + (AtInput + MuInput*TanBeta)*(3*(Log(0.98*Sqr
      (mAInput))*((2*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*
      Log(1.02*msq2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2))) + (2*Log(1.02*msq2(
      2,2))*Log(0.98*msu2(2,2)))/(-1.02*msq2(2,2) + 0.98*msu2(2,2)) + (2*Sqr(Log(
      1.02*msq2(2,2))))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (1.9223375624759709*(-
      0.98*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,
      2),1.02*msq2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msq2(2,2))) - (
      1.9223375624759709*(-2.9988*msq2(2,2)*msu2(2,2) + 0.9603999999999999*Quad(
      mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput) - 1.9207999999999998*msu2(2,2)*Sqr(
      mAInput) + 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)) -
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(msq2(2,2)))) + 3*Sqr(AbInput - MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*
      ((-2*Log(0.96*msd2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)) + (2*Log(1.02
      *msq2(2,2)))/Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + Log(0.96*msd2(2,2))*((2
      *Log(1.02*msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/
      ((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) - (
      2*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2) - 0.98*Sqr(mAInput)))
      /((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2)))) +
      (2*Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 0.98*msu2(2,2)
      - 0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (2*(-0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput
      ))*Sqr(Log(1.02*msq2(2,2))))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))) + (2.0833333333333335*TDelta(0.98*Sqr(mAInput),
      1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*
      msd2(2,2)))/(msd2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2)
      - 1.02*msq2(2,2))) + (1.9223375624759709*(0.98*(-0.96*msd2(2,2) + 1.02*msq2(
      2,2))*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput) + (0.96*msd2(2,2) -
      2.04*msq2(2,2))*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((1.02*msq2(2,2) -
      0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2))) - (
      2.0833333333333335*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2))*
      TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(1.02*msq2
      (2,2) - 0.98*msu2(2,2))*Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))) + (
      1.9223375624759709*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(-2.9988*msq2(2,2)*
      msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 2.9988*msq2(2,2)*Sqr(mAInput)
      - 1.9207999999999998*msu2(2,2)*Sqr(mAInput) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) - (0.96*msd2(2,2) - 2.04*msq2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2)))))) + Sqr(AtInput -
      MuInput/TanBeta)*(3*(4/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(Sqr(SCALE))*(
      2/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2.04*Log(0.98*msu2(2,2))*msq2(2,2))/
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((-4.08*msq2(2,2
      ))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2.04*Log(Sqr(SCALE))*msq2(2,2))/
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((Log(1.02*msq2(
      2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2
      ,2) - 0.98*msu2(2,2)) - (Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2
      ) - 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*Sqr
      (mAInput))*((Log(1.02*msq2(2,2))*(-0.96*msd2(2,2) + 1.02*msq2(2,2) + 0.98*
      Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (Log(0.98*msu2(2,2))*(
      -0.96*msd2(2,2) + 1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))
      /Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.0416666666666667*TDelta(0.98*Sqr(
      mAInput),1.02*msq2(2,2),0.96*msd2(2,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2
      ),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      1.0416666666666667*((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98
      *msu2(2,2) + 0.98*Sqr(mAInput)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(
      msd2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + 3*(AbInput - MuInput*
      TanBeta)*(AtInput + MuInput*TanBeta)*(Log(0.98*Sqr(mAInput))*((-2*Log(1.02*
      msq2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2)))/
      Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.96*msd2(2,2))*((2*Log(1.02*
      msq2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/((0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*Log(
      0.98*msu2(2,2))*(-0.96*msd2(2,2) - 1.02*msq2(2,2) + 0.98*Sqr(mAInput)))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + (2
      *Log(1.02*msq2(2,2))*Log(0.98*msu2(2,2))*(-2.04*msq2(2,2) + 0.98*Sqr(mAInput
      )))/((-0.96*msd2(2,2) + 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))
      ) + (2*(-2.04*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/((
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      2.0833333333333335*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*(0.96*msd2
      (2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (
      1.9607843137254901*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2))*
      TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/((0.96*msd2(2,2) -
      1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      2.0833333333333335*((1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98
      *msu2(2,2) + 0.98*Sqr(mAInput)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      0.96*msd2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(
      msd2(2,2)*(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2
      ,2))) + (0.9803921568627451*(-2*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(1.02*msq2
      (2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)) + 2*TDelta(0.98*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))/((0.96*msd2(2,2) - 1.02*msq2(2,2))*msq2(2,2)*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2)))) + Sqr(AtInput + MuInput*TanBeta)*(3*(Sqr(Log(1.02*msq2(2
      ,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(1.02*msq2(2,2))*(-(Log(0.98
      *msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (2*(1.02*msq2(2,2) -
      0.98*msu2(2,2) + 0.98*Sqr(mAInput)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + (
      0.9803921568627451*Log(0.98*msu2(2,2))*(0.9603999999999999*Quad(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) - 1.0404*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + Log(0.98*Sqr(mAInput))*(-(Log
      (1.02*msq2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(0.98*msu2(2,2))
      /Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (0.9803921568627451*(-Sqr(0.98*Sqr(
      mAInput) + 1.02*msq2(2,2) - 0.98*msu2(2,2)) + TDelta(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + (
      0.9611687812379854*(0.98*(-5.1*msq2(2,2) + 0.98*Sqr(mAInput))*Sqr(mAInput) -
      TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(Sqr(msq2(2,2))*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + (0.9611687812379854*((1.02*msq2(2,2) - 0.98*msu2(2,2))*
      (1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*(0.98*(-1.02*msq2(2,2)
      + 0.98*msu2(2,2))*msu2(2,2) + 0.9603999999999999*Quad(mAInput) - 0.98*(1.02*
      msq2(2,2) + 1.96*msu2(2,2))*Sqr(mAInput)) + TDelta(0.98*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2))*(-0.9603999999999999*Quad(mAInput) + 0.98*(2.04*
      msq2(2,2) + 2.94*msu2(2,2))*Sqr(mAInput) - 2*Sqr(1.02*msq2(2,2) - 0.98*msu2(
      2,2)) + TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))*TPhi(0.98*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(Sqr(msq2(2,2))*Sqr(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)
      ))) + 3*Sqr(AbInput - MuInput*TanBeta)*(Log(0.96*msd2(2,2))*((Log(1.02*msq2(
      2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/(Sqr(0.96*msd2(
      2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) - (Log(0.98*
      msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput)))/(Sqr(0.96*
      msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))) - ((0.96*
      msd2(2,2) + 1.02*msq2(2,2) - 0.98*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/(
      Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      Log(1.02*msq2(2,2))*((Log(0.98*msu2(2,2))*(0.96*msd2(2,2) + 1.02*msq2(2,2) -
      0.98*Sqr(mAInput)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + (2*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput)
      ))/((-0.96*msd2(2,2) + 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + (
      0.9803921568627451*Log(0.98*msu2(2,2))*(0.9603999999999999*Quad(mAInput) -
      1.9207999999999998*msu2(2,2)*Sqr(mAInput) - 1.0404*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2))))/(msq2(2,2)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)
      )) + (0.9803921568627451*Log(0.98*Sqr(mAInput))*(-Sqr(0.98*Sqr(mAInput) +
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2))))/(msq2(2,2)*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)
      )) + (1.0416666666666667*TDelta(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2
      ,2))*TPhi(0.98*Sqr(mAInput),1.02*msq2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(
      0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      0.9611687812379854*(0.98*(-0.96*msd2(2,2) + 1.02*msq2(2,2))*(-5.1*msq2(2,2)
      + 0.98*Sqr(mAInput))*Sqr(mAInput) + (0.96*msd2(2,2) - 2.04*msq2(2,2))*TDelta
      (0.98*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.98*Sqr(mAInput),
      1.02*msq2(2,2),1.02*msq2(2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(
      msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (1.0416666666666667*((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.96*msd2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr
      (mAInput)) - TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))*TPhi(
      0.98*Sqr(mAInput),0.98*msu2(2,2),0.96*msd2(2,2)))/(msd2(2,2)*Sqr(0.96*msd2(2
      ,2) - 1.02*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      0.9611687812379854*((0.96*msd2(2,2) - 1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98
      *msu2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.98*Sqr(mAInput))*(0.98*(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*msu2(2,2) - 0.9603999999999999*Quad(mAInput
      ) + 0.98*(1.02*msq2(2,2) + 1.96*msu2(2,2))*Sqr(mAInput)) + TDelta(0.98*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*(0.9603999999999999*(0.96*msd2(2,2)
      - 1.02*msq2(2,2))*Quad(mAInput) + 0.98*Sqr(mAInput)*(-1.9584*msd2(2,2)*msq2(
      2,2) - 2.8224*msd2(2,2)*msu2(2,2) + 3.9984*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(
      msq2(2,2))) + (1.92*msd2(2,2) - 3.06*msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) - (0.96*msd2(2,2) - 2.04*msq2(2,2))*TDelta(0.98*Sqr(mAInput),0.98
      *msu2(2,2),1.02*msq2(2,2))))*TPhi(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2
      (2,2)))/(Sqr(0.96*msd2(2,2) - 1.02*msq2(2,2))*Sqr(msq2(2,2))*Sqr(1.02*msq2(2
      ,2) - 0.98*msu2(2,2))*TDelta(0.98*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)
      ))))))/(1 + Sqr(TanBeta))))/Sqr(TanBeta)))/Sqr(1 + (0.006332573977646111*(1
      + Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(
      SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/Sqr(TanBeta) + 0.5*
      ((-0.0031662869888230555*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*msq2(2,2))) - (
      0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))) +
      0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/
      Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1
      + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2
      ))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)
      ) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((
      AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2))/
      M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(
      0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*
      TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/
      Sqr(TanBeta))))/Quad(3.141592653589793)), 0) + IF(TwoLoopAtAs >= 1, WHICH(
      IsCloseRel(Sqr(SCALE),msq2(2,2),0.01) && IsCloseRel(Sqr(SCALE),msu2(2,2),
      0.01) && IsCloseRel(SCALE,M3Input,0.01), (0.010416666666666666*Quad(Yu(2,2))
      *Sqr(g3)*((-12*(AtInput - MuInput/TanBeta))/SCALE + (14*Cube(AtInput -
      MuInput/TanBeta))/Cube(SCALE) - Power5(AtInput - MuInput/TanBeta)/Power5(
      SCALE) + (0.5*Quad(AtInput - MuInput/TanBeta))/Quad(SCALE) - (6*Sqr(AtInput
      - MuInput/TanBeta))/Sqr(SCALE)))/Quad(3.141592653589793), IsCloseRel(Sqr(
      M3Input),msq2(2,2),0.01) && IsCloseRel(Sqr(M3Input),msu2(2,2),0.01), (-
      0.005208333333333333*Quad(Yu(2,2))*Sqr(g3)*(((AtInput - MuInput/TanBeta)*(24
       + (12*(AtInput - MuInput/TanBeta))/M3Input - Cube(AtInput - MuInput/TanBeta
      )/Cube(M3Input) + (2*Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) - (28*
      Sqr(AtInput - MuInput/TanBeta))/Sqr(M3Input)))/M3Input - (2*(AtInput -
      MuInput/TanBeta)*Log(Sqr(M3Input)/Sqr(SCALE))*(24 - (24*(AtInput - MuInput/
      TanBeta))/M3Input + Cube(AtInput - MuInput/TanBeta)/Cube(M3Input) - (4*Sqr(
      AtInput - MuInput/TanBeta))/Sqr(M3Input)))/M3Input + 36*Sqr(Log(Sqr(M3Input)
      /Sqr(SCALE)))))/Quad(3.141592653589793), IsCloseRel(Sqr(M3Input),msq2(2,2),
      0.01), (-0.015625*Quad(Yu(2,2))*Sqr(g3)*Sqr(M3Input)*(4 - (32*(AtInput -
      MuInput/TanBeta)*msu2(2,2))/Cube(M3Input) + (24*(AtInput - MuInput/TanBeta)*
      Log(msu2(2,2)/Sqr(M3Input))*msu2(2,2))/Cube(M3Input) - (64*Cube(AtInput -
      MuInput/TanBeta)*msu2(2,2))/Power5(M3Input) - (32*Cube(AtInput - MuInput/
      TanBeta)*Log(msu2(2,2)/Sqr(M3Input))*msu2(2,2))/Power5(M3Input) + (16*Cube(
      msu2(2,2))*Power5(AtInput - MuInput/TanBeta))/Power11(M3Input) - (24*Cube(
      msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input))*Power5(AtInput - MuInput/TanBeta))/
      Power11(M3Input) - (21*Power5(msu2(2,2)))/Power10(M3Input) + (17*Log(msu2(2,
      2)/Sqr(M3Input))*Power5(msu2(2,2)))/Power10(M3Input) - (32*(AtInput -
      MuInput/TanBeta)*Power5(msu2(2,2)))/Power11(M3Input) + (24*(AtInput -
      MuInput/TanBeta)*Log(msu2(2,2)/Sqr(M3Input))*Power5(msu2(2,2)))/Power11(
      M3Input) - (82*Cube(msu2(2,2)))/Power6(M3Input) + (60*Cube(msu2(2,2))*Log(
      msu2(2,2)/Sqr(M3Input)))/Power6(M3Input) + (3*Power6(msu2(2,2)))/Power12(
      M3Input) - (3*Log(msu2(2,2)/Sqr(M3Input))*Power6(msu2(2,2)))/Power12(M3Input
      ) - (192*(AtInput - MuInput/TanBeta)*Cube(msu2(2,2)))/Power7(M3Input) + (144
      *(AtInput - MuInput/TanBeta)*Cube(msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input)))/
      Power7(M3Input) + (16*msu2(2,2)*Power5(AtInput - MuInput/TanBeta))/Power7(
      M3Input) + (8*Log(msu2(2,2)/Sqr(M3Input))*msu2(2,2)*Power5(AtInput - MuInput
      /TanBeta))/Power7(M3Input) - (192*Cube(AtInput - MuInput/TanBeta)*Cube(msu2(
      2,2)))/Power9(M3Input) + (32*Cube(AtInput - MuInput/TanBeta)*Cube(msu2(2,2))
      *Log(msu2(2,2)/Sqr(M3Input)))/Power9(M3Input) + (66*Cube(msu2(2,2))*Quad(
      AtInput - MuInput/TanBeta))/Power10(M3Input) - (47*Cube(msu2(2,2))*Log(msu2(
      2,2)/Sqr(M3Input))*Quad(AtInput - MuInput/TanBeta))/Power10(M3Input) + (14*
      msu2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input) + (19*Log(msu2(2,
      2)/Sqr(M3Input))*msu2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input)
      + (4*Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) + (64*Cube(AtInput -
      MuInput/TanBeta)*Quad(msu2(2,2)))/Power11(M3Input) - (32*Cube(AtInput -
      MuInput/TanBeta)*Log(msu2(2,2)/Sqr(M3Input))*Quad(msu2(2,2)))/Power11(
      M3Input) + (58*Quad(msu2(2,2)))/Power8(M3Input) - (44*Log(msu2(2,2)/Sqr(
      M3Input))*Quad(msu2(2,2)))/Power8(M3Input) + (128*(AtInput - MuInput/TanBeta
      )*Quad(msu2(2,2)))/Power9(M3Input) - (96*(AtInput - MuInput/TanBeta)*Log(
      msu2(2,2)/Sqr(M3Input))*Quad(msu2(2,2)))/Power9(M3Input) - (22*Quad(AtInput
      - MuInput/TanBeta)*Quad(msu2(2,2)))/Power12(M3Input) + (29*Log(msu2(2,2)/Sqr
      (M3Input))*Quad(AtInput - MuInput/TanBeta)*Quad(msu2(2,2)))/Power12(M3Input)
      - (25*msu2(2,2))/Sqr(M3Input) + (11*Log(msu2(2,2)/Sqr(M3Input))*msu2(2,2))/
      Sqr(M3Input) + (4*Cube(-1 + msu2(2,2)/Sqr(M3Input))*msu2(2,2)*PolyLog(2,((-1
       + msu2(2,2)/Sqr(M3Input))*Sqr(M3Input))/msu2(2,2))*(2 + (8*(AtInput -
      MuInput/TanBeta))/M3Input + (4*Cube(AtInput - MuInput/TanBeta))/Cube(M3Input
      ) + Quad(AtInput - MuInput/TanBeta)/Quad(M3Input)))/Sqr(M3Input) - (10*Log(
      msu2(2,2)/Sqr(M3Input))*Power5(msu2(2,2))*Sqr(AtInput - MuInput/TanBeta))/
      Power12(M3Input) + (32*Cube(msu2(2,2))*Sqr(AtInput - MuInput/TanBeta))/
      Power8(M3Input) - (96*Cube(msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input))*Sqr(
      AtInput - MuInput/TanBeta))/Power8(M3Input) + (32*msu2(2,2)*Sqr(AtInput -
      MuInput/TanBeta))/Quad(M3Input) - (22*Log(msu2(2,2)/Sqr(M3Input))*msu2(2,2)*
      Sqr(AtInput - MuInput/TanBeta))/Quad(M3Input) - (8*Quad(msu2(2,2))*Sqr(
      AtInput - MuInput/TanBeta))/Power10(M3Input) + (52*Log(msu2(2,2)/Sqr(M3Input
      ))*Quad(msu2(2,2))*Sqr(AtInput - MuInput/TanBeta))/Power10(M3Input) - (8*Sqr
      (AtInput - MuInput/TanBeta))/Sqr(M3Input) - (16*(AtInput - MuInput/TanBeta)*
      msu2(2,2)*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Cube(M3Input) - (8*Cube(AtInput
      - MuInput/TanBeta)*msu2(2,2)*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power5(
      M3Input) + (8*Cube(msu2(2,2))*Power5(AtInput - MuInput/TanBeta)*Sqr(Log(msu2
      (2,2)/Sqr(M3Input))))/Power11(M3Input) - (20*Power5(msu2(2,2))*Sqr(Log(msu2(
      2,2)/Sqr(M3Input))))/Power10(M3Input) - (8*(AtInput - MuInput/TanBeta)*
      Power5(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power11(M3Input) - (46*
      Cube(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power6(M3Input) + (4*
      Power6(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power12(M3Input) - (72*(
      AtInput - MuInput/TanBeta)*Cube(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))
      /Power7(M3Input) + (8*Cube(AtInput - MuInput/TanBeta)*Cube(msu2(2,2))*Sqr(
      Log(msu2(2,2)/Sqr(M3Input))))/Power9(M3Input) - (2*Cube(msu2(2,2))*Quad(
      AtInput - MuInput/TanBeta)*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power10(M3Input
      ) - (4*msu2(2,2)*Quad(AtInput - MuInput/TanBeta)*Sqr(Log(msu2(2,2)/Sqr(
      M3Input))))/Power6(M3Input) - (8*Cube(AtInput - MuInput/TanBeta)*Quad(msu2(2
      ,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power11(M3Input) + (42*Quad(msu2(2,2)
      )*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power8(M3Input) + (40*(AtInput - MuInput
      /TanBeta)*Quad(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power9(M3Input)
      - (10*Quad(AtInput - MuInput/TanBeta)*Quad(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(
      M3Input))))/Power12(M3Input) - (6*msu2(2,2)*Sqr(Log(msu2(2,2)/Sqr(M3Input)))
      )/Sqr(M3Input) + (16*Power5(msu2(2,2))*Sqr(AtInput - MuInput/TanBeta)*Sqr(
      Log(msu2(2,2)/Sqr(M3Input))))/Power12(M3Input) + (68*Cube(msu2(2,2))*Sqr(
      AtInput - MuInput/TanBeta)*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power8(M3Input)
      + (4*msu2(2,2)*Sqr(AtInput - MuInput/TanBeta)*Sqr(Log(msu2(2,2)/Sqr(M3Input)
      )))/Quad(M3Input) - (56*Quad(msu2(2,2))*Sqr(AtInput - MuInput/TanBeta)*Sqr(
      Log(msu2(2,2)/Sqr(M3Input))))/Power10(M3Input) + (12*msu2(2,2)*Power5(-1 +
      msu2(2,2)/Sqr(M3Input))*Sqr(Log(Sqr(M3Input)/Sqr(SCALE))))/Sqr(M3Input) + (
      128*(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power5(M3Input) - (96*(
      AtInput - MuInput/TanBeta)*Log(msu2(2,2)/Sqr(M3Input))*Sqr(msu2(2,2)))/
      Power5(M3Input) + (192*Cube(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/
      Power7(M3Input) + (32*Cube(AtInput - MuInput/TanBeta)*Log(msu2(2,2)/Sqr(
      M3Input))*Sqr(msu2(2,2)))/Power7(M3Input) - (32*Power5(AtInput - MuInput/
      TanBeta)*Sqr(msu2(2,2)))/Power9(M3Input) + (16*Log(msu2(2,2)/Sqr(M3Input))*
      Power5(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power9(M3Input) + (63*Sqr(
      msu2(2,2)))/Quad(M3Input) - (41*Log(msu2(2,2)/Sqr(M3Input))*Sqr(msu2(2,2)))/
      Quad(M3Input) - (62*Quad(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power8(
      M3Input) - (Log(msu2(2,2)/Sqr(M3Input))*Quad(AtInput - MuInput/TanBeta)*Sqr(
      msu2(2,2)))/Power8(M3Input) - (48*Sqr(AtInput - MuInput/TanBeta)*Sqr(msu2(2,
      2)))/Power6(M3Input) + (76*Log(msu2(2,2)/Sqr(M3Input))*Sqr(AtInput - MuInput
      /TanBeta)*Sqr(msu2(2,2)))/Power6(M3Input) + (56*(AtInput - MuInput/TanBeta)*
      Sqr(Log(msu2(2,2)/Sqr(M3Input)))*Sqr(msu2(2,2)))/Power5(M3Input) + (8*Cube(
      AtInput - MuInput/TanBeta)*Sqr(Log(msu2(2,2)/Sqr(M3Input)))*Sqr(msu2(2,2)))/
      Power7(M3Input) + (8*Power5(AtInput - MuInput/TanBeta)*Sqr(Log(msu2(2,2)/Sqr
      (M3Input)))*Sqr(msu2(2,2)))/Power9(M3Input) + (26*Sqr(Log(msu2(2,2)/Sqr(
      M3Input)))*Sqr(msu2(2,2)))/Quad(M3Input) + (20*Quad(AtInput - MuInput/
      TanBeta)*Sqr(Log(msu2(2,2)/Sqr(M3Input)))*Sqr(msu2(2,2)))/Power8(M3Input) -
      (32*Sqr(AtInput - MuInput/TanBeta)*Sqr(Log(msu2(2,2)/Sqr(M3Input)))*Sqr(msu2
      (2,2)))/Power6(M3Input) + 4*Log(Sqr(M3Input)/Sqr(SCALE))*(-(Quad(msu2(2,2))/
      Power8(M3Input)) + (msu2(2,2)*(4 + (8*Cube(AtInput - MuInput/TanBeta))/Cube(
      M3Input) - (6*Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) - (6*Sqr(
      AtInput - MuInput/TanBeta))/Sqr(M3Input)))/Sqr(M3Input) - (2*Cube(msu2(2,2))
      *(-2 + Sqr(AtInput - MuInput/TanBeta)/Sqr(M3Input)))/Power6(M3Input) + ((-6
      - (8*Cube(AtInput - MuInput/TanBeta))/Cube(M3Input) + (7*Quad(AtInput -
      MuInput/TanBeta))/Quad(M3Input) + (6*Sqr(AtInput - MuInput/TanBeta))/Sqr(
      M3Input))*Sqr(msu2(2,2)))/Quad(M3Input) + (Log(msu2(2,2)/Sqr(M3Input))*msu2(
      2,2)*(-3 - (4*(AtInput - MuInput/TanBeta))/M3Input + (4*Cube(AtInput -
      MuInput/TanBeta))/Cube(M3Input) + (3*Cube(msu2(2,2)))/Power6(M3Input) - (5*
      Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) + (6*Sqr(AtInput - MuInput/
      TanBeta))/Sqr(M3Input) + (msu2(2,2)*(9 + (8*(AtInput - MuInput/TanBeta))/
      M3Input + (4*Cube(AtInput - MuInput/TanBeta))/Cube(M3Input) - (3*Quad(
      AtInput - MuInput/TanBeta))/Quad(M3Input) - (12*Sqr(AtInput - MuInput/
      TanBeta))/Sqr(M3Input)))/Sqr(M3Input) + ((-9 - (4*(AtInput - MuInput/TanBeta
      ))/M3Input + (6*Sqr(AtInput - MuInput/TanBeta))/Sqr(M3Input))*Sqr(msu2(2,2))
      )/Quad(M3Input)))/Sqr(M3Input) - Sqr(-1 + Sqr(AtInput - MuInput/TanBeta)/Sqr
      (M3Input)))*Sqr(-1 + msu2(2,2)/Sqr(M3Input))))/(msu2(2,2)*Power5(-1 + msu2(2
      ,2)/Sqr(M3Input))*Quad(3.141592653589793)), IsCloseRel(Sqr(M3Input),msu2(2,2
      ),0.01), (-0.015625*Quad(Yu(2,2))*Sqr(g3)*Sqr(M3Input)*(4 - (32*(AtInput -
      MuInput/TanBeta)*msq2(2,2))/Cube(M3Input) + (24*(AtInput - MuInput/TanBeta)*
      Log(msq2(2,2)/Sqr(M3Input))*msq2(2,2))/Cube(M3Input) - (64*Cube(AtInput -
      MuInput/TanBeta)*msq2(2,2))/Power5(M3Input) - (32*Cube(AtInput - MuInput/
      TanBeta)*Log(msq2(2,2)/Sqr(M3Input))*msq2(2,2))/Power5(M3Input) + (16*Cube(
      msq2(2,2))*Power5(AtInput - MuInput/TanBeta))/Power11(M3Input) - (24*Cube(
      msq2(2,2))*Log(msq2(2,2)/Sqr(M3Input))*Power5(AtInput - MuInput/TanBeta))/
      Power11(M3Input) - (21*Power5(msq2(2,2)))/Power10(M3Input) + (17*Log(msq2(2,
      2)/Sqr(M3Input))*Power5(msq2(2,2)))/Power10(M3Input) - (32*(AtInput -
      MuInput/TanBeta)*Power5(msq2(2,2)))/Power11(M3Input) + (24*(AtInput -
      MuInput/TanBeta)*Log(msq2(2,2)/Sqr(M3Input))*Power5(msq2(2,2)))/Power11(
      M3Input) - (82*Cube(msq2(2,2)))/Power6(M3Input) + (60*Cube(msq2(2,2))*Log(
      msq2(2,2)/Sqr(M3Input)))/Power6(M3Input) + (3*Power6(msq2(2,2)))/Power12(
      M3Input) - (3*Log(msq2(2,2)/Sqr(M3Input))*Power6(msq2(2,2)))/Power12(M3Input
      ) - (192*(AtInput - MuInput/TanBeta)*Cube(msq2(2,2)))/Power7(M3Input) + (144
      *(AtInput - MuInput/TanBeta)*Cube(msq2(2,2))*Log(msq2(2,2)/Sqr(M3Input)))/
      Power7(M3Input) + (16*msq2(2,2)*Power5(AtInput - MuInput/TanBeta))/Power7(
      M3Input) + (8*Log(msq2(2,2)/Sqr(M3Input))*msq2(2,2)*Power5(AtInput - MuInput
      /TanBeta))/Power7(M3Input) - (192*Cube(AtInput - MuInput/TanBeta)*Cube(msq2(
      2,2)))/Power9(M3Input) + (32*Cube(AtInput - MuInput/TanBeta)*Cube(msq2(2,2))
      *Log(msq2(2,2)/Sqr(M3Input)))/Power9(M3Input) + (66*Cube(msq2(2,2))*Quad(
      AtInput - MuInput/TanBeta))/Power10(M3Input) - (47*Cube(msq2(2,2))*Log(msq2(
      2,2)/Sqr(M3Input))*Quad(AtInput - MuInput/TanBeta))/Power10(M3Input) + (14*
      msq2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input) + (19*Log(msq2(2,
      2)/Sqr(M3Input))*msq2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input)
      + (4*Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) + (64*Cube(AtInput -
      MuInput/TanBeta)*Quad(msq2(2,2)))/Power11(M3Input) - (32*Cube(AtInput -
      MuInput/TanBeta)*Log(msq2(2,2)/Sqr(M3Input))*Quad(msq2(2,2)))/Power11(
      M3Input) + (58*Quad(msq2(2,2)))/Power8(M3Input) - (44*Log(msq2(2,2)/Sqr(
      M3Input))*Quad(msq2(2,2)))/Power8(M3Input) + (128*(AtInput - MuInput/TanBeta
      )*Quad(msq2(2,2)))/Power9(M3Input) - (96*(AtInput - MuInput/TanBeta)*Log(
      msq2(2,2)/Sqr(M3Input))*Quad(msq2(2,2)))/Power9(M3Input) - (22*Quad(AtInput
      - MuInput/TanBeta)*Quad(msq2(2,2)))/Power12(M3Input) + (29*Log(msq2(2,2)/Sqr
      (M3Input))*Quad(AtInput - MuInput/TanBeta)*Quad(msq2(2,2)))/Power12(M3Input)
      - (25*msq2(2,2))/Sqr(M3Input) + (11*Log(msq2(2,2)/Sqr(M3Input))*msq2(2,2))/
      Sqr(M3Input) + (4*Cube(-1 + msq2(2,2)/Sqr(M3Input))*msq2(2,2)*PolyLog(2,((-1
       + msq2(2,2)/Sqr(M3Input))*Sqr(M3Input))/msq2(2,2))*(2 + (8*(AtInput -
      MuInput/TanBeta))/M3Input + (4*Cube(AtInput - MuInput/TanBeta))/Cube(M3Input
      ) + Quad(AtInput - MuInput/TanBeta)/Quad(M3Input)))/Sqr(M3Input) - (10*Log(
      msq2(2,2)/Sqr(M3Input))*Power5(msq2(2,2))*Sqr(AtInput - MuInput/TanBeta))/
      Power12(M3Input) + (32*Cube(msq2(2,2))*Sqr(AtInput - MuInput/TanBeta))/
      Power8(M3Input) - (96*Cube(msq2(2,2))*Log(msq2(2,2)/Sqr(M3Input))*Sqr(
      AtInput - MuInput/TanBeta))/Power8(M3Input) + (32*msq2(2,2)*Sqr(AtInput -
      MuInput/TanBeta))/Quad(M3Input) - (22*Log(msq2(2,2)/Sqr(M3Input))*msq2(2,2)*
      Sqr(AtInput - MuInput/TanBeta))/Quad(M3Input) - (8*Quad(msq2(2,2))*Sqr(
      AtInput - MuInput/TanBeta))/Power10(M3Input) + (52*Log(msq2(2,2)/Sqr(M3Input
      ))*Quad(msq2(2,2))*Sqr(AtInput - MuInput/TanBeta))/Power10(M3Input) - (8*Sqr
      (AtInput - MuInput/TanBeta))/Sqr(M3Input) - (16*(AtInput - MuInput/TanBeta)*
      msq2(2,2)*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Cube(M3Input) - (8*Cube(AtInput
      - MuInput/TanBeta)*msq2(2,2)*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power5(
      M3Input) + (8*Cube(msq2(2,2))*Power5(AtInput - MuInput/TanBeta)*Sqr(Log(msq2
      (2,2)/Sqr(M3Input))))/Power11(M3Input) - (20*Power5(msq2(2,2))*Sqr(Log(msq2(
      2,2)/Sqr(M3Input))))/Power10(M3Input) - (8*(AtInput - MuInput/TanBeta)*
      Power5(msq2(2,2))*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power11(M3Input) - (46*
      Cube(msq2(2,2))*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power6(M3Input) + (4*
      Power6(msq2(2,2))*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power12(M3Input) - (72*(
      AtInput - MuInput/TanBeta)*Cube(msq2(2,2))*Sqr(Log(msq2(2,2)/Sqr(M3Input))))
      /Power7(M3Input) + (8*Cube(AtInput - MuInput/TanBeta)*Cube(msq2(2,2))*Sqr(
      Log(msq2(2,2)/Sqr(M3Input))))/Power9(M3Input) - (2*Cube(msq2(2,2))*Quad(
      AtInput - MuInput/TanBeta)*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power10(M3Input
      ) - (4*msq2(2,2)*Quad(AtInput - MuInput/TanBeta)*Sqr(Log(msq2(2,2)/Sqr(
      M3Input))))/Power6(M3Input) - (8*Cube(AtInput - MuInput/TanBeta)*Quad(msq2(2
      ,2))*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power11(M3Input) + (42*Quad(msq2(2,2)
      )*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power8(M3Input) + (40*(AtInput - MuInput
      /TanBeta)*Quad(msq2(2,2))*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power9(M3Input)
      - (10*Quad(AtInput - MuInput/TanBeta)*Quad(msq2(2,2))*Sqr(Log(msq2(2,2)/Sqr(
      M3Input))))/Power12(M3Input) - (6*msq2(2,2)*Sqr(Log(msq2(2,2)/Sqr(M3Input)))
      )/Sqr(M3Input) + (16*Power5(msq2(2,2))*Sqr(AtInput - MuInput/TanBeta)*Sqr(
      Log(msq2(2,2)/Sqr(M3Input))))/Power12(M3Input) + (68*Cube(msq2(2,2))*Sqr(
      AtInput - MuInput/TanBeta)*Sqr(Log(msq2(2,2)/Sqr(M3Input))))/Power8(M3Input)
      + (4*msq2(2,2)*Sqr(AtInput - MuInput/TanBeta)*Sqr(Log(msq2(2,2)/Sqr(M3Input)
      )))/Quad(M3Input) - (56*Quad(msq2(2,2))*Sqr(AtInput - MuInput/TanBeta)*Sqr(
      Log(msq2(2,2)/Sqr(M3Input))))/Power10(M3Input) + (12*msq2(2,2)*Power5(-1 +
      msq2(2,2)/Sqr(M3Input))*Sqr(Log(Sqr(M3Input)/Sqr(SCALE))))/Sqr(M3Input) + (
      128*(AtInput - MuInput/TanBeta)*Sqr(msq2(2,2)))/Power5(M3Input) - (96*(
      AtInput - MuInput/TanBeta)*Log(msq2(2,2)/Sqr(M3Input))*Sqr(msq2(2,2)))/
      Power5(M3Input) + (192*Cube(AtInput - MuInput/TanBeta)*Sqr(msq2(2,2)))/
      Power7(M3Input) + (32*Cube(AtInput - MuInput/TanBeta)*Log(msq2(2,2)/Sqr(
      M3Input))*Sqr(msq2(2,2)))/Power7(M3Input) - (32*Power5(AtInput - MuInput/
      TanBeta)*Sqr(msq2(2,2)))/Power9(M3Input) + (16*Log(msq2(2,2)/Sqr(M3Input))*
      Power5(AtInput - MuInput/TanBeta)*Sqr(msq2(2,2)))/Power9(M3Input) + (63*Sqr(
      msq2(2,2)))/Quad(M3Input) - (41*Log(msq2(2,2)/Sqr(M3Input))*Sqr(msq2(2,2)))/
      Quad(M3Input) - (62*Quad(AtInput - MuInput/TanBeta)*Sqr(msq2(2,2)))/Power8(
      M3Input) - (Log(msq2(2,2)/Sqr(M3Input))*Quad(AtInput - MuInput/TanBeta)*Sqr(
      msq2(2,2)))/Power8(M3Input) - (48*Sqr(AtInput - MuInput/TanBeta)*Sqr(msq2(2,
      2)))/Power6(M3Input) + (76*Log(msq2(2,2)/Sqr(M3Input))*Sqr(AtInput - MuInput
      /TanBeta)*Sqr(msq2(2,2)))/Power6(M3Input) + (56*(AtInput - MuInput/TanBeta)*
      Sqr(Log(msq2(2,2)/Sqr(M3Input)))*Sqr(msq2(2,2)))/Power5(M3Input) + (8*Cube(
      AtInput - MuInput/TanBeta)*Sqr(Log(msq2(2,2)/Sqr(M3Input)))*Sqr(msq2(2,2)))/
      Power7(M3Input) + (8*Power5(AtInput - MuInput/TanBeta)*Sqr(Log(msq2(2,2)/Sqr
      (M3Input)))*Sqr(msq2(2,2)))/Power9(M3Input) + (26*Sqr(Log(msq2(2,2)/Sqr(
      M3Input)))*Sqr(msq2(2,2)))/Quad(M3Input) + (20*Quad(AtInput - MuInput/
      TanBeta)*Sqr(Log(msq2(2,2)/Sqr(M3Input)))*Sqr(msq2(2,2)))/Power8(M3Input) -
      (32*Sqr(AtInput - MuInput/TanBeta)*Sqr(Log(msq2(2,2)/Sqr(M3Input)))*Sqr(msq2
      (2,2)))/Power6(M3Input) + 4*Log(Sqr(M3Input)/Sqr(SCALE))*(-(Quad(msq2(2,2))/
      Power8(M3Input)) + (msq2(2,2)*(4 + (8*Cube(AtInput - MuInput/TanBeta))/Cube(
      M3Input) - (6*Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) - (6*Sqr(
      AtInput - MuInput/TanBeta))/Sqr(M3Input)))/Sqr(M3Input) - (2*Cube(msq2(2,2))
      *(-2 + Sqr(AtInput - MuInput/TanBeta)/Sqr(M3Input)))/Power6(M3Input) + ((-6
      - (8*Cube(AtInput - MuInput/TanBeta))/Cube(M3Input) + (7*Quad(AtInput -
      MuInput/TanBeta))/Quad(M3Input) + (6*Sqr(AtInput - MuInput/TanBeta))/Sqr(
      M3Input))*Sqr(msq2(2,2)))/Quad(M3Input) + (Log(msq2(2,2)/Sqr(M3Input))*msq2(
      2,2)*(-3 - (4*(AtInput - MuInput/TanBeta))/M3Input + (4*Cube(AtInput -
      MuInput/TanBeta))/Cube(M3Input) + (3*Cube(msq2(2,2)))/Power6(M3Input) - (5*
      Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) + (6*Sqr(AtInput - MuInput/
      TanBeta))/Sqr(M3Input) + (msq2(2,2)*(9 + (8*(AtInput - MuInput/TanBeta))/
      M3Input + (4*Cube(AtInput - MuInput/TanBeta))/Cube(M3Input) - (3*Quad(
      AtInput - MuInput/TanBeta))/Quad(M3Input) - (12*Sqr(AtInput - MuInput/
      TanBeta))/Sqr(M3Input)))/Sqr(M3Input) + ((-9 - (4*(AtInput - MuInput/TanBeta
      ))/M3Input + (6*Sqr(AtInput - MuInput/TanBeta))/Sqr(M3Input))*Sqr(msq2(2,2))
      )/Quad(M3Input)))/Sqr(M3Input) - Sqr(-1 + Sqr(AtInput - MuInput/TanBeta)/Sqr
      (M3Input)))*Sqr(-1 + msq2(2,2)/Sqr(M3Input))))/(msq2(2,2)*Power5(-1 + msq2(2
      ,2)/Sqr(M3Input))*Quad(3.141592653589793)), !IsClose(Sqr(M3Input),0) &&
      IsCloseRel(msu2(2,2)/Sqr(M3Input),msq2(2,2)/Sqr(M3Input),0.01), (
      0.005208333333333333*Power6(M3Input)*Quad(Yu(2,2))*Sqr(g3)*((32*Cube(AtInput
       - MuInput/TanBeta)*Log(msu2(2,2)/Sqr(M3Input)))/Cube(M3Input) + (32*Cube(
      AtInput - MuInput/TanBeta)*Log(1 - ((-1 + msu2(2,2)/Sqr(M3Input))*Sqr(
      M3Input))/msu2(2,2)))/Cube(M3Input) + (8*Cube(AtInput - MuInput/TanBeta)*
      msu2(2,2))/Power5(M3Input) - (72*Cube(AtInput - MuInput/TanBeta)*Log(msu2(2,
      2)/Sqr(M3Input))*msu2(2,2))/Power5(M3Input) - (8*Cube(AtInput - MuInput/
      TanBeta)*Log(Sqr(M3Input))*msu2(2,2))/Power5(M3Input) - (72*Cube(AtInput -
      MuInput/TanBeta)*Log(1 - ((-1 + msu2(2,2)/Sqr(M3Input))*Sqr(M3Input))/msu2(2
      ,2))*msu2(2,2))/Power5(M3Input) + (8*Cube(AtInput - MuInput/TanBeta)*Log(Sqr
      (SCALE))*msu2(2,2))/Power5(M3Input) - (18*Power5(msu2(2,2)))/Power10(M3Input
      ) + (24*Log(msu2(2,2)/Sqr(M3Input))*Power5(msu2(2,2)))/Power10(M3Input) + (
      24*Log(Sqr(M3Input))*Power5(msu2(2,2)))/Power10(M3Input) - (72*Log(msu2(2,2)
      /Sqr(M3Input))*Log(Sqr(M3Input))*Power5(msu2(2,2)))/Power10(M3Input) - (24*
      Log(Sqr(SCALE))*Power5(msu2(2,2)))/Power10(M3Input) + (72*Log(msu2(2,2)/Sqr(
      M3Input))*Log(Sqr(SCALE))*Power5(msu2(2,2)))/Power10(M3Input) + (72*Log(Sqr(
      M3Input))*Log(Sqr(SCALE))*Power5(msu2(2,2)))/Power10(M3Input) - (78*Cube(
      msu2(2,2)))/Power6(M3Input) + (84*Cube(msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input)
      ))/Power6(M3Input) + (72*Cube(msu2(2,2))*Log(Sqr(M3Input)))/Power6(M3Input)
      - (72*Cube(msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input))*Log(Sqr(M3Input)))/Power6(
      M3Input) - (72*Cube(msu2(2,2))*Log(Sqr(SCALE)))/Power6(M3Input) + (72*Cube(
      msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input))*Log(Sqr(SCALE)))/Power6(M3Input) + (
      72*Cube(msu2(2,2))*Log(Sqr(M3Input))*Log(Sqr(SCALE)))/Power6(M3Input) - (48*
      Cube(msu2(2,2))*PolyLog(2,((-1 + msu2(2,2)/Sqr(M3Input))*Sqr(M3Input))/msu2(
      2,2)))/Power6(M3Input) + (96*(AtInput - MuInput/TanBeta)*Cube(msu2(2,2)))/
      Power7(M3Input) - (144*(AtInput - MuInput/TanBeta)*Cube(msu2(2,2))*Log(msu2(
      2,2)/Sqr(M3Input)))/Power7(M3Input) - (96*(AtInput - MuInput/TanBeta)*Cube(
      msu2(2,2))*Log(Sqr(M3Input)))/Power7(M3Input) + (96*(AtInput - MuInput/
      TanBeta)*Cube(msu2(2,2))*Log(Sqr(SCALE)))/Power7(M3Input) + (96*(AtInput -
      MuInput/TanBeta)*Cube(msu2(2,2))*PolyLog(2,((-1 + msu2(2,2)/Sqr(M3Input))*
      Sqr(M3Input))/msu2(2,2)))/Power7(M3Input) + (4*msu2(2,2)*Power5(AtInput -
      MuInput/TanBeta))/Power7(M3Input) + (4*Log(msu2(2,2)/Sqr(M3Input))*msu2(2,2)
      *Power5(AtInput - MuInput/TanBeta))/Power7(M3Input) + (64*Cube(AtInput -
      MuInput/TanBeta)*Cube(msu2(2,2)))/Power9(M3Input) - (8*Cube(AtInput -
      MuInput/TanBeta)*Cube(msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input)))/Power9(M3Input
      ) - (8*Cube(AtInput - MuInput/TanBeta)*Cube(msu2(2,2))*Log(Sqr(M3Input)))/
      Power9(M3Input) + (8*Cube(AtInput - MuInput/TanBeta)*Cube(msu2(2,2))*Log(Sqr
      (SCALE)))/Power9(M3Input) - (Cube(msu2(2,2))*Quad(AtInput - MuInput/TanBeta)
      )/Power10(M3Input) + (6*Cube(msu2(2,2))*Log(msu2(2,2)/Sqr(M3Input))*Quad(
      AtInput - MuInput/TanBeta))/Power10(M3Input) + (6*Cube(msu2(2,2))*Log(Sqr(
      M3Input))*Quad(AtInput - MuInput/TanBeta))/Power10(M3Input) - (6*Cube(msu2(2
      ,2))*Log(Sqr(SCALE))*Quad(AtInput - MuInput/TanBeta))/Power10(M3Input) - (13
      *msu2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input) - (12*Log(msu2(2
      ,2)/Sqr(M3Input))*msu2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input)
      + (14*Log(Sqr(M3Input))*msu2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(
      M3Input) - (18*Log(1 - ((-1 + msu2(2,2)/Sqr(M3Input))*Sqr(M3Input))/msu2(2,2
      ))*msu2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input) - (14*Log(Sqr(
      SCALE))*msu2(2,2)*Quad(AtInput - MuInput/TanBeta))/Power6(M3Input) + (4*Quad
      (AtInput - MuInput/TanBeta))/Quad(M3Input) + (8*Log(msu2(2,2)/Sqr(M3Input))*
      Quad(AtInput - MuInput/TanBeta))/Quad(M3Input) - (4*Log(Sqr(M3Input))*Quad(
      AtInput - MuInput/TanBeta))/Quad(M3Input) + (8*Log(1 - ((-1 + msu2(2,2)/Sqr(
      M3Input))*Sqr(M3Input))/msu2(2,2))*Quad(AtInput - MuInput/TanBeta))/Quad(
      M3Input) + (4*Log(Sqr(SCALE))*Quad(AtInput - MuInput/TanBeta))/Quad(M3Input)
      + (72*Quad(msu2(2,2)))/Power8(M3Input) - (72*Log(msu2(2,2)/Sqr(M3Input))*
      Quad(msu2(2,2)))/Power8(M3Input) - (72*Log(Sqr(M3Input))*Quad(msu2(2,2)))/
      Power8(M3Input) + (144*Log(msu2(2,2)/Sqr(M3Input))*Log(Sqr(M3Input))*Quad(
      msu2(2,2)))/Power8(M3Input) + (72*Log(Sqr(SCALE))*Quad(msu2(2,2)))/Power8(
      M3Input) - (144*Log(msu2(2,2)/Sqr(M3Input))*Log(Sqr(SCALE))*Quad(msu2(2,2)))
      /Power8(M3Input) - (144*Log(Sqr(M3Input))*Log(Sqr(SCALE))*Quad(msu2(2,2)))/
      Power8(M3Input) - (48*(AtInput - MuInput/TanBeta)*Quad(msu2(2,2)))/Power9(
      M3Input) + (48*(AtInput - MuInput/TanBeta)*Log(msu2(2,2)/Sqr(M3Input))*Quad(
      msu2(2,2)))/Power9(M3Input) + (48*(AtInput - MuInput/TanBeta)*Log(Sqr(
      M3Input))*Quad(msu2(2,2)))/Power9(M3Input) - (48*(AtInput - MuInput/TanBeta)
      *Log(Sqr(SCALE))*Quad(msu2(2,2)))/Power9(M3Input) - (96*Cube(msu2(2,2))*Sqr(
      AtInput - MuInput/TanBeta))/Power8(M3Input) + (168*Cube(msu2(2,2))*Log(msu2(
      2,2)/Sqr(M3Input))*Sqr(AtInput - MuInput/TanBeta))/Power8(M3Input) + (168*
      Cube(msu2(2,2))*Log(Sqr(M3Input))*Sqr(AtInput - MuInput/TanBeta))/Power8(
      M3Input) - (168*Cube(msu2(2,2))*Log(Sqr(SCALE))*Sqr(AtInput - MuInput/
      TanBeta))/Power8(M3Input) - (24*msu2(2,2)*Sqr(AtInput - MuInput/TanBeta))/
      Quad(M3Input) + (24*Log(Sqr(M3Input))*msu2(2,2)*Sqr(AtInput - MuInput/
      TanBeta))/Quad(M3Input) - (24*Log(Sqr(SCALE))*msu2(2,2)*Sqr(AtInput -
      MuInput/TanBeta))/Quad(M3Input) + (12*Quad(msu2(2,2))*Sqr(AtInput - MuInput/
      TanBeta))/Power10(M3Input) - (72*Log(msu2(2,2)/Sqr(M3Input))*Quad(msu2(2,2))
      *Sqr(AtInput - MuInput/TanBeta))/Power10(M3Input) - (72*Log(Sqr(M3Input))*
      Quad(msu2(2,2))*Sqr(AtInput - MuInput/TanBeta))/Power10(M3Input) + (72*Log(
      Sqr(SCALE))*Quad(msu2(2,2))*Sqr(AtInput - MuInput/TanBeta))/Power10(M3Input)
      - (36*Power5(msu2(2,2))*Sqr(Log(Sqr(M3Input))))/Power10(M3Input) - (36*Cube(
      msu2(2,2))*Sqr(Log(Sqr(M3Input))))/Power6(M3Input) + (72*Quad(msu2(2,2))*Sqr
      (Log(Sqr(M3Input))))/Power8(M3Input) - (36*Power5(msu2(2,2))*Sqr(Log(Sqr(
      SCALE))))/Power10(M3Input) - (36*Cube(msu2(2,2))*Sqr(Log(Sqr(SCALE))))/
      Power6(M3Input) + (72*Quad(msu2(2,2))*Sqr(Log(Sqr(SCALE))))/Power8(M3Input)
      - (36*Power5(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power10(M3Input) -
      (36*Cube(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power6(M3Input) + (72*
      Quad(msu2(2,2))*Sqr(Log(msu2(2,2)/Sqr(M3Input))))/Power8(M3Input) - (48*(
      AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power5(M3Input) + (96*(AtInput -
      MuInput/TanBeta)*Log(msu2(2,2)/Sqr(M3Input))*Sqr(msu2(2,2)))/Power5(M3Input)
      + (48*(AtInput - MuInput/TanBeta)*Log(Sqr(M3Input))*Sqr(msu2(2,2)))/Power5(
      M3Input) + (96*(AtInput - MuInput/TanBeta)*Log(1 - ((-1 + msu2(2,2)/Sqr(
      M3Input))*Sqr(M3Input))/msu2(2,2))*Sqr(msu2(2,2)))/Power5(M3Input) - (48*(
      AtInput - MuInput/TanBeta)*Log(Sqr(SCALE))*Sqr(msu2(2,2)))/Power5(M3Input) -
      (72*Cube(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power7(M3Input) + (16*
      Cube(AtInput - MuInput/TanBeta)*Log(Sqr(M3Input))*Sqr(msu2(2,2)))/Power7(
      M3Input) + (48*Cube(AtInput - MuInput/TanBeta)*Log(1 - ((-1 + msu2(2,2)/Sqr(
      M3Input))*Sqr(M3Input))/msu2(2,2))*Sqr(msu2(2,2)))/Power7(M3Input) - (16*
      Cube(AtInput - MuInput/TanBeta)*Log(Sqr(SCALE))*Sqr(msu2(2,2)))/Power7(
      M3Input) - (4*Power5(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power9(
      M3Input) + (24*Sqr(msu2(2,2)))/Quad(M3Input) - (24*Log(Sqr(M3Input))*Sqr(
      msu2(2,2)))/Quad(M3Input) + (24*Log(Sqr(SCALE))*Sqr(msu2(2,2)))/Quad(M3Input
      ) + (10*Quad(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power8(M3Input) - (4
      *Log(msu2(2,2)/Sqr(M3Input))*Quad(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))
      /Power8(M3Input) - (16*Log(Sqr(M3Input))*Quad(AtInput - MuInput/TanBeta)*Sqr
      (msu2(2,2)))/Power8(M3Input) + (12*Log(1 - ((-1 + msu2(2,2)/Sqr(M3Input))*
      Sqr(M3Input))/msu2(2,2))*Quad(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/
      Power8(M3Input) + (16*Log(Sqr(SCALE))*Quad(AtInput - MuInput/TanBeta)*Sqr(
      msu2(2,2)))/Power8(M3Input) + (108*Sqr(AtInput - MuInput/TanBeta)*Sqr(msu2(2
      ,2)))/Power6(M3Input) - (48*Log(msu2(2,2)/Sqr(M3Input))*Sqr(AtInput -
      MuInput/TanBeta)*Sqr(msu2(2,2)))/Power6(M3Input) - (120*Log(Sqr(M3Input))*
      Sqr(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power6(M3Input) + (120*Log(
      Sqr(SCALE))*Sqr(AtInput - MuInput/TanBeta)*Sqr(msu2(2,2)))/Power6(M3Input)))
      /(Cube(msu2(2,2))*Quad(3.141592653589793)*Sqr(-1 + msu2(2,2)/Sqr(M3Input))),
      True, (0.015625*Quad(Yu(2,2))*Sqr(g3)*(Log(Sqr(M3Input)/Sqr(SCALE))*(8 - 12*
      Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input)) - 12*Log((1.01007550314386*
      msu2(2,2))/Sqr(M3Input)) - (4.028144723618089*Sqr(M3Input))/msq2(2,2) - (
      3.9601*Sqr(M3Input))/msu2(2,2)) + ((AtInput - MuInput/TanBeta)*(Log(Sqr(
      M3Input)/Sqr(SCALE))*((16*Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input)))/
      ((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)) - (16*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input)))/((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))) + Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(-16/((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)) - (8*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input)))/((-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,2
      ))/Sqr(M3Input)))) + (16*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input)))/((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)) + (32*(PolyLog(2,(0.990025*(-1 + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))*Sqr(M3Input))/msu2(2,2))*(-1 + (0.9930129810249694*msq2(2,2))/
      Sqr(M3Input)) - PolyLog(2,(1.0070361809045223*(-1 + (0.9930129810249694*msq2
      (2,2))/Sqr(M3Input))*Sqr(M3Input))/msq2(2,2))*(-1 + (1.01007550314386*msu2(2
      ,2))/Sqr(M3Input))))/((-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))) + (8*(-2 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(Log((0.9930129810249694*msq2
      (2,2))/Sqr(M3Input))))/((-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(
      (0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))) - (8*(-2 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(Log
      ((1.01007550314386*msu2(2,2))/Sqr(M3Input))))/(((0.9930129810249694*msq2(2,2
      ))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))))/M3Input + (Power5(AtInput -
      MuInput/TanBeta)*(Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((8*Log((
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*((0.9930129810249694*msq2(2,2))/
      Sqr(M3Input) + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*((-
      1.01007550314386*msu2(2,2))/Sqr(M3Input) + (0.9930129810249694*msq2(2,2)*(-1
       + (2.02015100628772*msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/(Quad((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))*(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))) + (15.88820769639951*msq2(2,2))/(
      Cube((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,
      2))/Sqr(M3Input))*(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(
      M3Input))) + (16.16120805030176*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input
      ))*msu2(2,2))/(Cube((-0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,2))
      /Sqr(M3Input))*Sqr(M3Input)) - (7.944103848199755*msq2(2,2)*((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))*Sqr(Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input))))/(Quad((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))*(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(M3Input
      )) - (8.08060402515088*msu2(2,2)*((0.9930129810249694*msq2(2,2))/Sqr(M3Input
      ) + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(Log((1.01007550314386*
      msu2(2,2))/Sqr(M3Input))))/(Quad((0.9930129810249694*msq2(2,2))/Sqr(M3Input)
      - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,
      2))/Sqr(M3Input))*Sqr(M3Input))))/Power5(M3Input) - 12*Sqr(Log(Sqr(M3Input)/
      Sqr(SCALE))) + (Sqr(Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input)))*(-6 +
      (7.944103848199755*msq2(2,2))/Sqr(M3Input) - (3.9442991219363845*Sqr(msq2(2,
      2)))/Quad(M3Input)))/Sqr(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input)) +
      (Cube(AtInput - MuInput/TanBeta)*((-16*(4*((0.9930129810249694*msq2(2,2))/
      Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input)) + PolyLog(2,(
      1.0070361809045223*(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(
      M3Input))/msq2(2,2))*(-2 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)) - PolyLog(2,(0.990025*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))/msu2(2,2))*(-2 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))))/Cube((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)) - (16*Log((1.01007550314386*msu2(2
      ,2))/Sqr(M3Input))*((0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      3.03022650943158*msu2(2,2))/Sqr(M3Input)))/Cube((0.9930129810249694*msq2(2,2
      ))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input)) + (8*Sqr(Log((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input)))*(-2 + (1.01007550314386*msu2(2,
      2))/Sqr(M3Input) - (2.979038943074908*msq2(2,2)*(-1 + (1.01007550314386*msu2
      (2,2))/Sqr(M3Input)))/Sqr(M3Input) + (0.9860747804840961*Sqr(msq2(2,2)))/
      Quad(M3Input)))/(Cube((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (0.9930129810249694*msq2(2,2
      ))/Sqr(M3Input))) + (8*Sqr(Log((1.01007550314386*msu2(2,2))/Sqr(M3Input)))*(
      -2 - (3.0090542593115406*msq2(2,2)*msu2(2,2))/Quad(M3Input) + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (3.03022650943158*msu2(2,2))/
      Sqr(M3Input) + (1.0202525220513219*Sqr(msu2(2,2)))/Quad(M3Input)))/(Cube((-
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))) + Log(Sqr(
      M3Input)/Sqr(SCALE))*((-16*Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input))*
      ((0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)))/Cube((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)) + (16*Log((1.01007550314386*msu2(2
      ,2))/Sqr(M3Input))*((0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Cube((0.9930129810249694*msq2(2,2
      ))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input)) + 32/Sqr((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))) + Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((16*((
      2.979038943074908*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/Sqr
      (M3Input)))/Cube((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)) + (16*Log((1.01007550314386*msu2(2
      ,2))/Sqr(M3Input))*((-2.0060361728743605*msq2(2,2)*msu2(2,2))/Quad(M3Input)
      + (0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))
      /Sqr(M3Input)))/((-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr((0.9930129810249694*msq2(2,2))
      /Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))))))/Cube(M3Input)
      + Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((2*(-7 + (
      6.06045301886316*msu2(2,2))/Sqr(M3Input) + (0.9930129810249694*msq2(2,2)*(7
      - (5.050377515719299*msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + (
      0.9860747804840961*(-3 + (2.02015100628772*msu2(2,2))/Sqr(M3Input))*Sqr(msq2
      (2,2)))/Quad(M3Input)))/((-1 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*
      Sqr(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))) - (2*Log((
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*((1.01007550314386*msu2(2,2)*(-2 +
      (1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + (
      0.9930129810249694*msq2(2,2)*(-2 + (8.08060402515088*msu2(2,2))/Sqr(M3Input)
      - (4.0810100882052875*Sqr(msu2(2,2)))/Quad(M3Input)))/Sqr(M3Input) + (
      0.9860747804840961*Sqr(msq2(2,2))*(1 - (4.04030201257544*msu2(2,2))/Sqr(
      M3Input) + (2.0405050441026438*Sqr(msu2(2,2)))/Quad(M3Input)))/Quad(M3Input)
      ))/(Sqr(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))) + (Sqr(AtInput - MuInput/TanBeta
      )*((-7.975927959999997*Quad(M3Input))/(msq2(2,2)*msu2(2,2)) + Log(Sqr(
      M3Input)/Sqr(SCALE))*((7.975927959999997*Quad(M3Input))/(msq2(2,2)*msu2(2,2)
      ) - (24*Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input)))/((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)) + (24*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input)))/((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input))) - (4*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input))*(7 - (
      4.04030201257544*msu2(2,2))/Sqr(M3Input) + (2.979038943074908*msq2(2,2)*(-2
      + (1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/((-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((0.9930129810249694*msq2(2,2))/
      Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))) - (4*Sqr(Log((0.9930129810249694*
      msq2(2,2))/Sqr(M3Input)))*((3.9167402291282185*Cube(msq2(2,2)))/Power6(
      M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input) + (0.9930129810249694*
      msq2(2,2)*(3 + (4.04030201257544*msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) - (
      1.9721495609681923*(4 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(msq2(
      2,2)))/Quad(M3Input)))/(Sqr(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input)
      )*Sqr((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2
      ,2))/Sqr(M3Input))) + Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((4*(
      7 - (6.06045301886316*msu2(2,2))/Sqr(M3Input) + (0.9930129810249694*msq2(2,2
      )*(-4 + (3.03022650943158*msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/((-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((0.9930129810249694*msq2(2,2))/
      Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))) + (4*Log((1.01007550314386*msu2(2
      ,2))/Sqr(M3Input))*(2*((0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)) + (0.9930129810249694*msq2(2,2)*(-
      2 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((0.9930129810249694*msq2(2
      ,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input)))/(Sqr(M3Input
      )*Sqr(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))) + (1.01007550314386
      *msu2(2,2)*(-2 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*((-
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + (1.01007550314386*msu2(2,2))/Sqr(
      M3Input)))))/Sqr((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))) + (4*Sqr(Log((1.01007550314386*
      msu2(2,2))/Sqr(M3Input)))*((1.01007550314386*msu2(2,2)*(-3 + (
      8.08060402515088*msu2(2,2))/Sqr(M3Input) - (4.0810100882052875*Sqr(msu2(2,2)
      ))/Quad(M3Input)))/Sqr(M3Input) + (0.9930129810249694*msq2(2,2)*(1 - (
      4.04030201257544*msu2(2,2))/Sqr(M3Input) + (2.0405050441026438*Sqr(msu2(2,2)
      ))/Quad(M3Input)))/Sqr(M3Input)))/(Sqr((0.9930129810249694*msq2(2,2))/Sqr(
      M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))))/Sqr(M3Input) + (Sqr(Log((
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))*(-6 + (8.08060402515088*msu2(2,2)
      )/Sqr(M3Input) - (4.0810100882052875*Sqr(msu2(2,2)))/Quad(M3Input)))/Sqr(-1
      + (1.01007550314386*msu2(2,2))/Sqr(M3Input)) + (2*Log((1.01007550314386*msu2
      (2,2))/Sqr(M3Input))*(-7 + (7.07052852200702*msu2(2,2))/Sqr(M3Input) - (
      3.0607575661539657*Sqr(msu2(2,2)))/Quad(M3Input) + (0.9930129810249694*msq2(
      2,2)*(6 - (5.050377515719299*msu2(2,2))/Sqr(M3Input) + (2.0405050441026438*
      Sqr(msu2(2,2)))/Quad(M3Input)))/Sqr(M3Input)))/((-1 + (0.9930129810249694*
      msq2(2,2))/Sqr(M3Input))*Sqr(-1 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))
      ) - (1.9939819899999993*Quad(M3Input)*((-1 + (0.9930129810249694*msq2(2,2))/
      Sqr(M3Input))*((4.012072345748721*msq2(2,2)*msu2(2,2)*PolyLog(2,(0.990025*(-
      1 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))/msu2(2,2))*(-1
      + (0.9930129810249694*msq2(2,2))/Sqr(M3Input)))/Quad(M3Input) + (-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*((2.02015100628772*msu2(2,2)*(-1 +
      (1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input) + (
      0.9930129810249694*msq2(2,2)*(-2 + (9.09067952829474*msu2(2,2))/Sqr(M3Input)
      - (6.121515132307931*Sqr(msu2(2,2)))/Quad(M3Input)))/Sqr(M3Input) + (
      0.9860747804840961*Sqr(msq2(2,2))*(2 - (6.06045301886316*msu2(2,2))/Sqr(
      M3Input) + (3.0607575661539657*Sqr(msu2(2,2)))/Quad(M3Input)))/Quad(M3Input)
      )) + (4.012072345748721*msq2(2,2)*msu2(2,2)*PolyLog(2,(1.0070361809045223*(-
      1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(M3Input))/msq2(2,2))*
      Sqr(-1 + (1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Quad(M3Input)))/(msq2(2
      ,2)*msu2(2,2)*Sqr(-1 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(-1 +
      (1.01007550314386*msu2(2,2))/Sqr(M3Input))) + (Quad(AtInput - MuInput/
      TanBeta)*((-3.9879639799999986*Quad(M3Input)*((0.9791850572820546*Cube(msq2(
      2,2)))/Power6(M3Input) - (1.030532079544781*Cube(msu2(2,2)))/Power6(M3Input)
      - (5.116658661753119*Cube(msu2(2,2))*msq2(2,2))/Power8(M3Input) + (
      4.945254197025603*Cube(msq2(2,2))*msu2(2,2))/Power8(M3Input) + (
      1.0030180864371803*msq2(2,2)*msu2(2,2)*PolyLog(2,(1.0070361809045223*(-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(M3Input))/msq2(2,2))*(-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,2
      ))/Sqr(M3Input))*(-2 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Quad(M3Input) - (
      1.0030180864371803*msq2(2,2)*msu2(2,2)*PolyLog(2,(0.990025*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(M3Input))/msu2(2,2))*(-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,2
      ))/Sqr(M3Input))*(-2 + (0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Quad(M3Input) + (
      5.080908470594695*Cube(msu2(2,2))*Sqr(msq2(2,2)))/Power10(M3Input) - (
      5.976059880209667*msu2(2,2)*Sqr(msq2(2,2)))/Power6(M3Input) - (
      0.9860747804840961*Sqr(msq2(2,2)))/Quad(M3Input) - (4.995080121234921*Cube(
      msq2(2,2))*Sqr(msu2(2,2)))/Power10(M3Input) + (6.078743989922559*msq2(2,2)*
      Sqr(msu2(2,2)))/Power6(M3Input) + (1.0202525220513219*Sqr(msu2(2,2)))/Quad(
      M3Input)))/(Cube((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*msq2(2,2)*msu2(2,2)*(-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(-1 + (1.01007550314386*msu2(2,2
      ))/Sqr(M3Input))) + (2*Sqr(Log((0.9930129810249694*msq2(2,2))/Sqr(M3Input)))
      *((4.861717363533792*Quad(msq2(2,2)))/Power8(M3Input) + (3.9167402291282185*
      Cube(msq2(2,2))*(-2 + (1.01007550314386*msu2(2,2))/Sqr(M3Input)))/Power6(
      M3Input) - (2.02015100628772*msu2(2,2))/Sqr(M3Input) - (0.9960099800349447*
      msu2(2,2)*(10 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(msq2(2,2)))/
      Power6(M3Input) + (1.9860259620499388*msq2(2,2)*(1 + (4.04030201257544*msu2(
      2,2))/Sqr(M3Input) + (1.0202525220513219*Sqr(msu2(2,2)))/Quad(M3Input)))/Sqr
      (M3Input)))/(Quad((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(-1 + (0.9930129810249694*msq2(
      2,2))/Sqr(M3Input))) + Log(Sqr(M3Input)/Sqr(SCALE))*((4*Log((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*(2 + (2.979038943074908*msq2(2,2
      ))/Sqr(M3Input) + (3.03022650943158*msu2(2,2))/Sqr(M3Input)))/Cube((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)) - (4*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input))*(2 + (
      2.979038943074908*msq2(2,2))/Sqr(M3Input) + (3.03022650943158*msu2(2,2))/Sqr
      (M3Input)))/Cube((0.9930129810249694*msq2(2,2))/Sqr(M3Input) - (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)) - (3.9879639799999986*Quad(M3Input
      )*((6.018108518623082*msq2(2,2)*msu2(2,2))/Quad(M3Input) + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)))/(msq2(2,2)*msu2(2,2)*Sqr((0.9930129810249694*msq2(2,2))/Sqr(
      M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input)))) + Log((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((2*(4 + (0.9791850572820546*
      Cube(msq2(2,2))*(16 - (15.1511325471579*msu2(2,2))/Sqr(M3Input)))/Power6(
      M3Input) + (3.03022650943158*msu2(2,2))/Sqr(M3Input) - (6.121515132307931*
      Sqr(msu2(2,2)))/Quad(M3Input) + (0.9860747804840961*Sqr(msq2(2,2))*(-29 + (
      31.31234059745966*msu2(2,2))/Sqr(M3Input) - (3.0607575661539657*Sqr(msu2(2,2
      )))/Quad(M3Input)))/Quad(M3Input) + (0.9930129810249694*msq2(2,2)*(7 - (
      15.1511325471579*msu2(2,2))/Sqr(M3Input) + (7.141767654359253*Sqr(msu2(2,2))
      )/Quad(M3Input)))/Sqr(M3Input)))/(Cube((0.9930129810249694*msq2(2,2))/Sqr(
      M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(-1 + (0.9930129810249694*msq2(
      2,2))/Sqr(M3Input))) + (2*Log((1.01007550314386*msu2(2,2))/Sqr(M3Input))*((-
      8.024144691497442*msq2(2,2)*msu2(2,2))/Quad(M3Input) - (3.9442991219363845*
      Sqr(msq2(2,2)))/Quad(M3Input) - (0.9930129810249694*msq2(2,2)*(-2 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*((0.9930129810249694*msq2(2,2))/
      Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*((
      0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))/
      Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + (0.9930129810249694*msq2(2,2))/Sqr(
      M3Input))) - (4.0810100882052875*Sqr(msu2(2,2)))/Quad(M3Input) - (
      1.01007550314386*msu2(2,2)*(-2 + (1.01007550314386*msu2(2,2))/Sqr(M3Input))*
      ((-0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (1.01007550314386*msu2(2,2))
      /Sqr(M3Input))*((0.9930129810249694*msq2(2,2))/Sqr(M3Input) + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))/(Sqr(M3Input)*Sqr(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))))/Quad((0.9930129810249694*msq2(2
      ,2))/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))) + (2*Log((
      1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-4 - (16.488513272716496*Cube(
      msu2(2,2)))/Power6(M3Input) - (7.07052852200702*msu2(2,2))/Sqr(M3Input) + (
      29.587323139488333*Sqr(msu2(2,2)))/Quad(M3Input) + (0.9930129810249694*msq2(
      2,2)*(-3 + (15.457981193171715*Cube(msu2(2,2)))/Power6(M3Input) + (
      15.1511325471579*msu2(2,2))/Sqr(M3Input) - (31.627828183590978*Sqr(msu2(2,2)
      ))/Quad(M3Input)))/Sqr(M3Input) + (0.9860747804840961*Sqr(msq2(2,2))*(6 - (
      7.07052852200702*msu2(2,2))/Sqr(M3Input) + (3.0607575661539657*Sqr(msu2(2,2)
      ))/Quad(M3Input)))/Quad(M3Input)))/(Cube((0.9930129810249694*msq2(2,2))/Sqr(
      M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*(-1 + (
      0.9930129810249694*msq2(2,2))/Sqr(M3Input))*Sqr(-1 + (1.01007550314386*msu2(
      2,2))/Sqr(M3Input))) - (2*Sqr(Log((1.01007550314386*msu2(2,2))/Sqr(M3Input))
      )*((8.244256636358248*Cube(msu2(2,2)))/Power6(M3Input) - (5.204576043760415*
      Quad(msu2(2,2)))/Power8(M3Input) - (2.02015100628772*msu2(2,2))/Sqr(M3Input)
      + (0.9960099800349447*msu2(2,2)*(-2 + (1.01007550314386*msu2(2,2))/Sqr(
      M3Input))*Sqr(msq2(2,2)))/Power6(M3Input) - (1.9860259620499388*msq2(2,2)*(-
      1 + (2.02015100628772*msu2(2,2))/Sqr(M3Input))*Sqr(-1 + (1.01007550314386*
      msu2(2,2))/Sqr(M3Input)))/Sqr(M3Input)))/(Quad((0.9930129810249694*msq2(2,2)
      )/Sqr(M3Input) - (1.01007550314386*msu2(2,2))/Sqr(M3Input))*Sqr(-1 + (
      1.01007550314386*msu2(2,2))/Sqr(M3Input)))))/Quad(M3Input)))/Quad(
      3.141592653589793)), 0) + IF(TwoLoopAtAt >= 1, WHICH(IsCloseRel(msu2(2,2),
      msq2(2,2),0.01) && IsCloseRel(Sqr(mAInput),msu2(2,2),0.01), (0.01171875*
      Power6(Yu(2,2))*(1 + Sqr(TanBeta))*(0.5 - 8*IF(Abs(-1 + Sqr(MuInput)/Sqrt(
      msq2(2,2)*msu2(2,2))) < 0.00001, -1, (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2
      ,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2
      ,2)*msu2(2,2))))) + 4*IF(Abs(-1 + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) <
      0.00001, -2.25, Re(((-1 + (2*Quad(MuInput))/(msq2(2,2)*msu2(2,2)) + (2*Sqr(
      MuInput))/Sqrt(msq2(2,2)*msu2(2,2)))*(Log(Abs(1 - Sqr(MuInput)/Sqrt(msq2(2,2
      )*msu2(2,2))))*Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) + PolyLog(2,Sqr(
      MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) - 0.16666666666666666*Sqr(
      3.141592653589793) - (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(
      MuInput))/Sqrt(msq2(2,2)*msu2(2,2))))/Sqr(Abs(1 - Sqr(MuInput)/Sqrt(msq2(2,2
      )*msu2(2,2)))))) - 4*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE)) + (6*Sqr(
      MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) - (2*IF(Abs(-1 + Sqr(MuInput)/Sqrt(msq2(
      2,2)*msu2(2,2))) < 0.00001, -1, (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))
      *Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*
      msu2(2,2)))))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) + (3*IF(Abs(-1 + Sqr(
      MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) < 0.00001, 0.5, (1 + (Log(Sqr(MuInput)/
      Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr
      (MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))))/(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2
      (2,2))))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) - (6*Log(Sqrt(msq2(2,2)*
      msu2(2,2))/Sqr(SCALE))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) -
      8.34993159891064/(1 + Sqr(TanBeta)) + (13*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(
      SCALE)))/(1 + Sqr(TanBeta)) + (Power6(AtInput - MuInput/TanBeta)*(-0.5 + 0.5
      *Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE)) + (0.5 - 0.5*Log(Sqrt(msq2(2,2)*
      msu2(2,2))/Sqr(SCALE)))/(1 + Sqr(TanBeta))))/Power3(Sqrt(msq2(2,2)*msu2(2,2)
      )) + (Cube(AtInput - MuInput/TanBeta)*((AtInput - MuInput/TanBeta)/Sqrt(Sqrt
      (msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2
      ,2)*msu2(2,2))))*(0.8747904000000002/(1 + Sqr(TanBeta)) - (2*Log(Sqrt(msq2(2
      ,2)*msu2(2,2))/Sqr(SCALE)))/(1 + Sqr(TanBeta))))/Power3(Sqrt(Sqrt(msq2(2,2)*
      msu2(2,2)))) + ((AtInput - MuInput/TanBeta)*((AtInput - MuInput/TanBeta)/
      Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta)))/Sqrt(
      Sqrt(msq2(2,2)*msu2(2,2))))*(0.5008383999999992/(1 + Sqr(TanBeta)) + (12*Log
      (Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE)))/(1 + Sqr(TanBeta))))/Sqrt(Sqrt(msq2(
      2,2)*msu2(2,2))) + 3*Sqr(Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE))) - (3*Sqr
      (Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + (
      0.1252095999999998/(1 + Sqr(TanBeta)) + (3*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr
      (SCALE)))/(1 + Sqr(TanBeta)))*Sqr((AtInput - MuInput/TanBeta)/Sqrt(Sqrt(msq2
      (2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2,2)*
      msu2(2,2)))) + (Sqr(AtInput - MuInput/TanBeta)*(-7 + 4*IF(Abs(-1 + Sqr(
      MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) < 0.00001, -1, (Log(Sqr(MuInput)/Sqrt(
      msq2(2,2)*msu2(2,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr(
      MuInput)/Sqrt(msq2(2,2)*msu2(2,2))))) - 4*IF(Abs(-1 + Sqr(MuInput)/Sqrt(msq2
      (2,2)*msu2(2,2))) < 0.00001, 0.5, (1 + (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2
      (2,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2
      (2,2)*msu2(2,2)))))/(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))) + 27*Log(
      Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE)) - (6*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2
      (2,2)) - (6*IF(Abs(-1 + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) < 0.00001, -
      1, (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2
      )*msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))))*Sqr(MuInput))/
      Sqrt(msq2(2,2)*msu2(2,2)) - (6*IF(Abs(-1 + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(
      2,2))) < 0.00001, 0.5, (1 + (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr
      (MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(
      2,2)))))/(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))))*Sqr(MuInput))/Sqrt(
      msq2(2,2)*msu2(2,2)) + (6*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE))*Sqr(
      MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) + 19.6878144/(1 + Sqr(TanBeta)) - (24*
      Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE)))/(1 + Sqr(TanBeta)) + (-
      0.021147733333332752/(1 + Sqr(TanBeta)) - (3*Log(Sqrt(msq2(2,2)*msu2(2,2))/
      Sqr(SCALE)))/(1 + Sqr(TanBeta)))*Sqr((AtInput - MuInput/TanBeta)/Sqrt(Sqrt(
      msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta)))/Sqrt(Sqrt(msq2(2,
      2)*msu2(2,2))))))/Sqrt(msq2(2,2)*msu2(2,2)) + (Quad(AtInput - MuInput/
      TanBeta)*(5.5 - 0.5*IF(Abs(-1 + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) <
      0.00001, -1, (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(MuInput))/(
      Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))))) +
      0.5*IF(Abs(-1 + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) < 0.00001, 0.5, (1 +
      (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*
      msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))))/(1 - Sqr(MuInput)/
      Sqrt(msq2(2,2)*msu2(2,2)))) - 6.5*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE))
      + Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)) + (IF(Abs(-1 + Sqr(MuInput)/Sqrt(
      msq2(2,2)*msu2(2,2))) < 0.00001, -1, (Log(Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2
      ,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))*(1 - Sqr(MuInput)/Sqrt(msq2(2
      ,2)*msu2(2,2)))))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) + (0.5*IF(Abs(-1 +
      Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2))) < 0.00001, 0.5, (1 + (Log(Sqr(
      MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))*Sqr(MuInput))/(Sqrt(msq2(2,2)*msu2(2,2))
      *(1 - Sqr(MuInput)/Sqrt(msq2(2,2)*msu2(2,2)))))/(1 - Sqr(MuInput)/Sqrt(msq2(
      2,2)*msu2(2,2))))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) - (Log(Sqrt(msq2(2
      ,2)*msu2(2,2))/Sqr(SCALE))*Sqr(MuInput))/Sqrt(msq2(2,2)*msu2(2,2)) - 6.25/(1
       + Sqr(TanBeta)) + (6*Log(Sqrt(msq2(2,2)*msu2(2,2))/Sqr(SCALE)))/(1 + Sqr(
      TanBeta)) + (-0.020728533333333354/(1 + Sqr(TanBeta)) + (0.5*Log(Sqrt(msq2(2
      ,2)*msu2(2,2))/Sqr(SCALE)))/(1 + Sqr(TanBeta)))*Sqr((AtInput - MuInput/
      TanBeta)/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))) + (2*MuInput*Csc(2*ArcTan(TanBeta))
      )/Sqrt(Sqrt(msq2(2,2)*msu2(2,2))))))/(msq2(2,2)*msu2(2,2))))/(Quad(
      3.141592653589793)*Sqr(TanBeta)), True, (0.00390625*Power6(Yu(2,2))*(3*(5 -
      4*Log(0.98*msu2(2,2)) + Log(1.02*msq2(2,2))*(-6 - 6*Log(Sqr(SCALE))) + (10 -
      4*Log(0.98*msu2(2,2)))*Log(Sqr(SCALE)) + 5*Sqr(Log(Sqr(SCALE))) + 3*Sqr(Log(
      1.02*msq2(2,2))) + 2*Sqr(Log(0.98*msu2(2,2)))) + 3*Power6(AtInput - MuInput/
      TanBeta)*((-2*PolyLog(2,1 - (1.0408163265306123*msq2(2,2))/msu2(2,2)))/Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*PolyLog(2,(0.9803921568627451*(1.02*
      msq2(2,2) - 0.98*msu2(2,2)))/msq2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2
      )) + ((-9.18*msq2(2,2) - 4.9*msu2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Quad(1.02*
      msq2(2,2) - 0.98*msu2(2,2)) + ((-3.06*msq2(2,2) - 4.9*msu2(2,2))*Sqr(Log(
      0.98*msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(1.02*msq2(2,2)
      )*((6.12*Log(Sqr(SCALE))*msq2(2,2))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      (2*Log(0.98*msu2(2,2))*(6.12*msq2(2,2) + 4.9*msu2(2,2)))/Quad(1.02*msq2(2,2)
      - 0.98*msu2(2,2)) + (1.0204081632653061*(12.994800000000001*msq2(2,2)*msu2(2
      ,2) + 2.0808*Sqr(msq2(2,2)) - 2.8811999999999998*Sqr(msu2(2,2))))/(msu2(2,2)
      *Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)))) + Log(Sqr(SCALE))*((-6.12*Log(0.98*
      msu2(2,2))*msq2(2,2))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (
      1.0004001600640255*(-4.998*msq2(2,2)*msu2(2,2) - 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      msq2(2,2)*msu2(2,2))) + (0.9803921568627451*Log(0.98*msu2(2,2))*(-2.9988*
      msq2(2,2)*msu2(2,2) - 10.404*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,
      2))))/(msq2(2,2)*Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      1.0004001600640255*(-10.9956*msq2(2,2)*msu2(2,2) - 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      msq2(2,2)*msu2(2,2))) + 3*Quad(AtInput - MuInput/TanBeta)*(((-
      7.140000000000001*msq2(2,2) - 8.82*msu2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Cube
      (1.02*msq2(2,2) - 0.98*msu2(2,2)) + ((3.06*msq2(2,2) + 4.9*msu2(2,2))*Sqr(
      Log(0.98*msu2(2,2))))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (6*PolyLog(2,1
       - (1.0408163265306123*msq2(2,2))/msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2
      (2,2)) + Log(1.02*msq2(2,2))*((4*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*
      msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(Sqr(SCALE))*(5.1*
      msq2(2,2) + 6.859999999999999*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,
      2)) - (4.081632653061225*(-9.996*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2))
      - 2.8811999999999998*Sqr(msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))
      *msu2(2,2))) + Log(Sqr(SCALE))*((-2*Log(0.98*msu2(2,2))*(5.1*msq2(2,2) +
      6.859999999999999*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (
      2.000800320128051*(-14.994*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))))/(msq2(2,2)*msu2(2,2)*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2)))) + (2.000800320128051*(-26.9892*msq2(2,2)*msu2(2,2) +
      2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/(msq2(2,2)*msu2(
      2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (0.9803921568627451*Log(0.98*
      msu2(2,2))*(-43.9824*msq2(2,2)*msu2(2,2) - 6.2424*Sqr(msq2(2,2)) +
      1.9207999999999998*Sqr(msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      msq2(2,2))) + 3*Sqr(AtInput - MuInput/TanBeta)*((-6*PolyLog(2,1 - (
      1.0408163265306123*msq2(2,2))/msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2))
      - (9*Sqr(Log(0.98*msu2(2,2))))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + ((
      7.140000000000001*msq2(2,2) - 8.82*msu2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(Sqr(SCALE))*((2*Log(0.98*msu2(2,2))*(
      8.16*msq2(2,2) - 8.82*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (
      1.0004001600640255*(-0.9996*msq2(2,2)*msu2(2,2) - 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2
      ,2))*msu2(2,2))) + (1.0004001600640255*(-2.9988*msq2(2,2)*msu2(2,2) - 2.0808
      *Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/(msq2(2,2)*(1.02*msq2(
      2,2) - 0.98*msu2(2,2))*msu2(2,2)) + (0.9803921568627451*Log(0.98*msu2(2,2))*
      (-20.9916*msq2(2,2)*msu2(2,2) + 16.6464*Sqr(msq2(2,2)) + 0.9603999999999999*
      Sqr(msu2(2,2))))/(msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02
      *msq2(2,2))*((2.04*Log(0.98*msu2(2,2))*msq2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) - (2*Log(Sqr(SCALE))*(8.16*msq2(2,2) - 8.82*msu2(2,2)))/Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2)) + (1.0204081632653061*(-16.9932*msq2(2,2)*msu2(2
      ,2) + 2.0808*Sqr(msq2(2,2)) + 18.2476*Sqr(msu2(2,2))))/(msu2(2,2)*Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2))))) + ((1 + Sqr(TanBeta))*(Sqr(AtInput - MuInput/
      TanBeta)*(9*((4*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*
      Log(0.98*msu2(2,2))*Log(Sqr(SCALE)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log
      (1.02*msq2(2,2))*(-4/(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (4*Log(Sqr(SCALE)))
      /(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (4*Log(0.98*msu2(2,2))*(1.02*msq2(2,2)
      + 0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (4.08*msq2(2,2)*
      Sqr(Log(1.02*msq2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (3.92*msu2(
      2,2)*Sqr(Log(0.98*msu2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + 3*((-
      4*PolyLog(2,1 - (0.9702970297029703*msu2(2,2))/Sqr(MuInput))*(-1.02*msq2(2,2
      ) + 1.01*Sqr(MuInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*PolyLog(2,
      (0.9900990099009901*(-1.02*msq2(2,2) + 1.01*Sqr(MuInput)))/Sqr(MuInput))*(-
      1.02*msq2(2,2) + 1.01*Sqr(MuInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      Log(Sqr(SCALE))*((2*Log(0.98*msu2(2,2))*(2.04*msq2(2,2) + 0.98*msu2(2,2) -
      2.02*Sqr(MuInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (
      1.0004001600640255*(3.9984*msq2(2,2)*msu2(2,2) - 8.2416*msq2(2,2)*Sqr(
      MuInput) + 3.9592*msu2(2,2)*Sqr(MuInput) + 4.1616*Sqr(msq2(2,2)) -
      1.9207999999999998*Sqr(msu2(2,2))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2
      ,2))*msu2(2,2))) + (1.0004001600640255*(5.9976*msq2(2,2)*msu2(2,2) - 8.2416*
      msq2(2,2)*Sqr(MuInput) + 3.9592*msu2(2,2)*Sqr(MuInput) + 4.1616*Sqr(msq2(2,2
      )) - 1.9207999999999998*Sqr(msu2(2,2))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*msu2(2,2)) + (0.9803921568627451*Log(0.98*msu2(2,2))*(-4.20362808
      *msq2(2,2)*Power6(MuInput) + 1.0201*Quad(MuInput)*(16.9932*msq2(2,2)*msu2(2,
      2) + 1.0404*Sqr(msq2(2,2)) - 1.9207999999999998*Sqr(msu2(2,2))) + 0.9996*
      msq2(2,2)*msu2(2,2)*(2.9988*msq2(2,2)*msu2(2,2) + 7.2828*Sqr(msq2(2,2)) -
      1.9207999999999998*Sqr(msu2(2,2))) - 1.01*Sqr(MuInput)*(1.061208*Cube(msq2(2
      ,2)) - 1.8823839999999998*Cube(msu2(2,2)) + 16.313472*msu2(2,2)*Sqr(msq2(2,2
      )) + 4.89804*msq2(2,2)*Sqr(msu2(2,2)))))/(msq2(2,2)*(-1.02*msq2(2,2) + 1.01*
      Sqr(MuInput))*(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2.04*msq2(2,2)*(1.02*msq2(2,2)*(1.02*msq2(2,2) + 1.96*
      msu2(2,2)) + 3.0603*Quad(MuInput) - 2.02*(1.02*msq2(2,2) + 1.96*msu2(2,2))*
      Sqr(MuInput))*Sqr(Log(1.02*msq2(2,2))))/(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)
      )*Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))) + Log(1.02*msq2(2,2))*((-2*Log(
      Sqr(SCALE))*(2.04*msq2(2,2) + 0.98*msu2(2,2) - 2.02*Sqr(MuInput)))/Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2)) + (1.0204081632653061*(4.03877992*msu2(2,2)*
      Power6(MuInput) + 0.9996*msq2(2,2)*msu2(2,2)*(-6.997199999999999*msq2(2,2)*
      msu2(2,2) - 4.1616*Sqr(msq2(2,2)) + 2.8811999999999998*Sqr(msu2(2,2))) -
      1.0201*Quad(MuInput)*(0.9996*msq2(2,2)*msu2(2,2) + 4.1616*Sqr(msq2(2,2)) +
      10.5644*Sqr(msu2(2,2))) + 1.01*Sqr(MuInput)*(4.244832*Cube(msq2(2,2)) +
      4.705959999999999*Cube(msu2(2,2)) + 9.176328*msu2(2,2)*Sqr(msq2(2,2)) +
      1.9592159999999998*msq2(2,2)*Sqr(msu2(2,2)))))/(msu2(2,2)*(-1.02*msq2(2,2) +
      1.01*Sqr(MuInput))*(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + Log(0.98*msu2(2,2))*((-6.12*msq2(2,2))/Sqr(1.02*msq2(2,
      2) - 0.98*msu2(2,2)) - (3.92*msu2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))
      + (4.08*msq2(2,2)*(1.02*msq2(2,2) - 2.02*Sqr(MuInput)))/((1.02*msq2(2,2) -
      0.98*msu2(2,2))*Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))) + (1.96*msu2(2,2)*
      (0.98*msu2(2,2) - 2.02*Sqr(MuInput)))/((-1.02*msq2(2,2) + 0.98*msu2(2,2))*
      Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))))) + Log(1.01*Sqr(MuInput))*((
      4.041616646658663*Sqr(MuInput)*(1.0201*(2.04*msq2(2,2) - 0.98*msu2(2,2))*
      Quad(MuInput) + 1.019592*msu2(2,2)*Sqr(msq2(2,2)) + 1.01*Sqr(MuInput)*(-
      0.9996*msq2(2,2)*msu2(2,2) - 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(
      msu2(2,2)))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*msu2(2,2)*(-1.02*
      msq2(2,2) + 1.01*Sqr(MuInput))*(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))) + (2*
      Log(0.98*msu2(2,2))*(2.1020201002*Power10(MuInput) - 1.04060401*(3.06*msq2(2
      ,2) + 6.859999999999999*msu2(2,2))*Power8(MuInput) + 2.05957584*msu2(2,2)*(
      2.04*msq2(2,2) + 2.94*msu2(2,2))*Sqr(MuInput)*Sqr(msq2(2,2)) -
      2.0383683263999997*Cube(msq2(2,2))*Sqr(msu2(2,2)) + 2.060602*Power6(MuInput)
      *(4.998*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)) + 2.8811999999999998*Sqr
      (msu2(2,2))) - 1.0201*Quad(MuInput)*(1.061208*Cube(msq2(2,2)) +
      1.8823839999999998*Cube(msu2(2,2)) + 13.254696000000001*msu2(2,2)*Sqr(msq2(2
      ,2)) + 3.9184319999999997*msq2(2,2)*Sqr(msu2(2,2)))))/(Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))*Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*msu2(2,2)
      + 1.01*Sqr(MuInput))) + (2*Log(1.02*msq2(2,2))*(-2.1020201002*Power10(
      MuInput) + 1.04060401*(3.06*msq2(2,2) + 6.859999999999999*msu2(2,2))*Power8(
      MuInput) - 2.05957584*msu2(2,2)*(2.04*msq2(2,2) + 2.94*msu2(2,2))*Sqr(
      MuInput)*Sqr(msq2(2,2)) + 2.0383683263999997*Cube(msq2(2,2))*Sqr(msu2(2,2))
      - 2.060602*Power6(MuInput)*(4.998*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)
      ) + 2.8811999999999998*Sqr(msu2(2,2))) + 1.0201*Quad(MuInput)*(1.061208*Cube
      (msq2(2,2)) + 1.8823839999999998*Cube(msu2(2,2)) + 13.254696000000001*msu2(2
      ,2)*Sqr(msq2(2,2)) + 3.9184319999999997*msq2(2,2)*Sqr(msu2(2,2)))))/(Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*
      Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput)))) + (1.96*msu2(2,2)*(0.98*(1.02*
      msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2) + 2.0402*Quad(MuInput) - 2.02*(1.02*
      msq2(2,2) + 0.98*msu2(2,2))*Sqr(MuInput))*Sqr(Log(0.98*msu2(2,2))))/(Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput)))))
      + Quad(AtInput - MuInput/TanBeta)*(9*((-4*Log(0.98*msu2(2,2))*(1.02*msq2(2,2
      ) + 2.94*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (4.08*msq2(2,2)
      *(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(Log(1.02*msq2(2,2))))/Quad(1.02*msq2(
      2,2) - 0.98*msu2(2,2)) - (3.92*(1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2)*
      Sqr(Log(0.98*msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(Sqr(
      SCALE))*((-4*Log(0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) - 8/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) -
      8/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(1.02*msq2(2,2))*((4*Log(Sqr(
      SCALE))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2
      ,2)) + (4*(3.06*msq2(2,2) + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2
      (2,2)) + (4*Log(0.98*msu2(2,2))*Sqr(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2)))) + 3*((2*PolyLog(2,(0.9900990099009901*(-
      1.02*msq2(2,2) + 1.01*Sqr(MuInput)))/Sqr(MuInput))*(1.9992*msq2(2,2)*msu2(2,
      2) + 3.0603*Quad(MuInput) - 6.1812000000000005*msq2(2,2)*Sqr(MuInput) +
      2.0808*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2
      ,2) - 0.98*msu2(2,2)) + (2*PolyLog(2,1 - (0.9702970297029703*msu2(2,2))/Sqr(
      MuInput))*(-1.9992*msq2(2,2)*msu2(2,2) - 3.0603*Quad(MuInput) +
      6.1812000000000005*msq2(2,2)*Sqr(MuInput) - 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) +
      Log(Sqr(SCALE))*((3*Log(0.98*msu2(2,2))*(-1.9992*msq2(2,2)*msu2(2,2) +
      4.1208*msq2(2,2)*Sqr(MuInput) - 3.1212*Sqr(msq2(2,2)) + 0.9603999999999999*
      Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.0004001600640255
      *(-2.122416*Cube(msq2(2,2)) + 0.9411919999999999*Cube(msu2(2,2)) -
      13.254696000000001*msu2(2,2)*Sqr(msq2(2,2)) + 2.02*Sqr(MuInput)*(4.998*msq2(
      2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2)))
      + 1.9592159999999998*msq2(2,2)*Sqr(msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*msq2(2,2)*msu2(2,2))) + (1.0004001600640255*(2.060602*Power6(
      MuInput)*(7.9968*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2))) + 0.9996*msq2(2,2)*msu2(2,2)*(-2.122416*
      Cube(msq2(2,2)) + 0.9411919999999999*Cube(msu2(2,2)) - 24.470208*msu2(2,2)*
      Sqr(msq2(2,2)) + 6.857256*msq2(2,2)*Sqr(msu2(2,2))) - 1.0201*Quad(MuInput)*(
      6.367248*Cube(msq2(2,2)) - 2.8235759999999996*Cube(msu2(2,2)) +
      38.744496000000005*msu2(2,2)*Sqr(msq2(2,2)) + 12.734903999999998*msq2(2,2)*
      Sqr(msu2(2,2))) + 1.01*Sqr(MuInput)*(-5.7600950399999995*Cube(msu2(2,2))*
      msq2(2,2) + 29.119547519999998*Cube(msq2(2,2))*msu2(2,2) + 2.16486432*Quad(
      msq2(2,2)) - 0.9223681599999999*Quad(msu2(2,2)) + 30.975204959999996*Sqr(
      msq2(2,2))*Sqr(msu2(2,2)))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2
      )*msu2(2,2)*(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*(-0.98*msu2(2,2) + 1.01*
      Sqr(MuInput))) - (2.04*msq2(2,2)*Sqr(Log(1.02*msq2(2,2)))*(1.0201*(4.08*msq2
      (2,2) + 2.94*msu2(2,2))*Quad(MuInput) + 1.02*msq2(2,2)*(2.9988*msq2(2,2)*
      msu2(2,2) + 3.1212*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) -
      2.02*Sqr(MuInput)*(2.9988*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)))))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))) + Log(1.02*msq2(2,2))*((-3*Log(Sqr
      (SCALE))*(-1.9992*msq2(2,2)*msu2(2,2) + 4.1208*msq2(2,2)*Sqr(MuInput) -
      3.1212*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2
      ,2) - 0.98*msu2(2,2)) + (0.5102040816326531*(-37.446759662256*msq2(2,2)*msu2
      (2,2)*Power8(MuInput) - 1.019592*msu2(2,2)*Sqr(msq2(2,2))*(4.244832*Cube(
      msq2(2,2)) - 17.882648*Cube(msu2(2,2)) + 47.920824*msu2(2,2)*Sqr(msq2(2,2))
      + 3.9184319999999997*msq2(2,2)*Sqr(msu2(2,2))) + 1.030301*Power6(MuInput)*(
      4.244832*Cube(msq2(2,2)) - 4.705959999999999*Cube(msu2(2,2)) + 123.370632*
      msu2(2,2)*Sqr(msq2(2,2)) + 23.510592*msq2(2,2)*Sqr(msu2(2,2))) + 1.0201*Quad
      (MuInput)*(40.32066528*Cube(msu2(2,2))*msq2(2,2) - 131.03796384*Cube(msq2(2,
      2))*msu2(2,2) - 8.65945728*Quad(msq2(2,2)) + 6.4565771199999995*Quad(msu2(2,
      2)) - 130.89522096*Sqr(msq2(2,2))*Sqr(msu2(2,2))) + 1.0302*msq2(2,2)*Sqr(
      MuInput)*(-8.64014256*Cube(msu2(2,2))*msq2(2,2) + 55.119143519999994*Cube(
      msq2(2,2))*msu2(2,2) + 4.32972864*Quad(msq2(2,2)) - 31.360517439999995*Quad(
      msu2(2,2)) + 129.89602079999997*Sqr(msq2(2,2))*Sqr(msu2(2,2)))))/(msu2(2,2)*
      Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))*
      Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))) + Log(0.98*msu2(2,2))*((10.2*msq2(
      2,2))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (5.88*msu2(2,2))/Cube(-1.02*
      msq2(2,2) + 0.98*msu2(2,2)) + (21.9912*msq2(2,2)*msu2(2,2))/Quad(1.02*msq2(2
      ,2) - 0.98*msu2(2,2)) - 2/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2.04*msq2(
      2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2))*(-1.02*msq2(2,2) + 2.02*Sqr(MuInput))
      )/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(-1.02*msq2(2,2) + 1.01*Sqr(
      MuInput))) + (0.98*(1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2)*(-0.98*msu2(2
      ,2) + 2.02*Sqr(MuInput)))/(Cube(-1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(-0.98*
      msu2(2,2) + 1.01*Sqr(MuInput))))) + Log(1.01*Sqr(MuInput))*((Log(1.02*msq2(2
      ,2))*(12.864363013223999*msq2(2,2)*Power10(MuInput) - 1.04060401*Power8(
      MuInput)*(27.9888*msq2(2,2)*msu2(2,2) + 26.009999999999998*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + 1.9984003199999998*Sqr(msq2(2,2))*(-
      1.9992*msq2(2,2)*msu2(2,2) - 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(
      msu2(2,2)))*Sqr(msu2(2,2)) + 6.3054421199999995*msq2(2,2)*Power6(MuInput)*(
      9.996*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) + 2.8811999999999998*Sqr(
      msu2(2,2))) + 4.038384*msq2(2,2)*msu2(2,2)*Sqr(MuInput)*(2.122416*Cube(msq2(
      2,2)) - 0.9411919999999999*Cube(msu2(2,2)) + 7.137144*msu2(2,2)*Sqr(msq2(2,2
      )) + 0.9796079999999999*msq2(2,2)*Sqr(msu2(2,2))) - 1.040502*msq2(2,2)*Quad(
      MuInput)*(3.183624*Cube(msq2(2,2)) - 3.7647679999999997*Cube(msu2(2,2)) +
      44.862048*msu2(2,2)*Sqr(msq2(2,2)) + 40.163928*msq2(2,2)*Sqr(msu2(2,2)))))/(
      Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput
      ))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))) + (Log(0.98*msu2(2,2))*(-
      12.864363013223999*msq2(2,2)*Power10(MuInput) + 1.04060401*Power8(MuInput)*(
      27.9888*msq2(2,2)*msu2(2,2) + 26.009999999999998*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))) + 1.9984003199999998*Sqr(msq2(2,2))*(
      1.9992*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(
      msu2(2,2)))*Sqr(msu2(2,2)) - 6.3054421199999995*msq2(2,2)*Power6(MuInput)*(
      9.996*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) + 2.8811999999999998*Sqr(
      msu2(2,2))) - 4.038384*msq2(2,2)*msu2(2,2)*Sqr(MuInput)*(2.122416*Cube(msq2(
      2,2)) - 0.9411919999999999*Cube(msu2(2,2)) + 7.137144*msu2(2,2)*Sqr(msq2(2,2
      )) + 0.9796079999999999*msq2(2,2)*Sqr(msu2(2,2))) + 1.040502*msq2(2,2)*Quad(
      MuInput)*(3.183624*Cube(msq2(2,2)) - 3.7647679999999997*Cube(msu2(2,2)) +
      44.862048*msu2(2,2)*Sqr(msq2(2,2)) + 40.163928*msq2(2,2)*Sqr(msu2(2,2)))))/(
      Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput
      ))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))) - (2.0208083233293315*Sqr(
      MuInput)*(-3.1511510352*Cube(msq2(2,2))*msu2(2,2)*(1.02*msq2(2,2) + 2.94*
      msu2(2,2))*Sqr(MuInput) + 1.04060401*Power8(MuInput)*(1.9992*msq2(2,2)*msu2(
      2,2) + 2.0808*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2))) +
      0.9992001599999999*Sqr(msq2(2,2))*(1.9992*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(
      msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2)))*Sqr(msu2(2,2)) - 1.030301*
      Power6(MuInput)*(4.244832*Cube(msq2(2,2)) - 1.8823839999999998*Cube(msu2(2,2
      )) + 5.0979600000000005*msu2(2,2)*Sqr(msq2(2,2)) + 4.89804*msq2(2,2)*Sqr(
      msu2(2,2))) + 1.0201*Quad(MuInput)*(1.92003168*Cube(msu2(2,2))*msq2(2,2) +
      8.319870719999999*Cube(msq2(2,2))*msu2(2,2) + 2.16486432*Quad(msq2(2,2)) -
      0.9223681599999999*Quad(msu2(2,2)) + 6.994401119999999*Sqr(msq2(2,2))*Sqr(
      msu2(2,2)))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2)*msu2(2,2)*Sqr
      (-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput
      )))) - (0.98*(1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2)*(0.98*msu2(2,2)*(
      1.02*msq2(2,2) + 2.94*msu2(2,2)) + 4.0804*Quad(MuInput) - 2.02*(1.02*msq2(2,
      2) + 2.94*msu2(2,2))*Sqr(MuInput))*Sqr(Log(0.98*msu2(2,2))))/(Quad(1.02*msq2
      (2,2) - 0.98*msu2(2,2))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))) + (
      0.49019607843137253*Log(0.98*msu2(2,2))*(12.736993082400001*msq2(2,2)*(2.04*
      msq2(2,2) + 0.98*msu2(2,2))*Power8(MuInput) + 0.9796079999999999*msq2(2,2)*
      Sqr(msu2(2,2))*(26.530199999999997*Cube(msq2(2,2)) - 1.8823839999999998*Cube
      (msu2(2,2)) + 34.666128*msu2(2,2)*Sqr(msq2(2,2)) - 20.571768*msq2(2,2)*Sqr(
      msu2(2,2))) - 1.030301*Power6(MuInput)*(41.387111999999995*Cube(msq2(2,2)) -
      1.8823839999999998*Cube(msu2(2,2)) + 103.998384*msu2(2,2)*Sqr(msq2(2,2)) +
      4.89804*msq2(2,2)*Sqr(msu2(2,2))) + 0.9898*msu2(2,2)*Sqr(MuInput)*(
      20.16033264*Cube(msu2(2,2))*msq2(2,2) - 117.51817391999998*Cube(msq2(2,2))*
      msu2(2,2) - 51.95674368*Quad(msq2(2,2)) + 1.8447363199999998*Quad(msu2(2,2))
      - 5.995200959999999*Sqr(msq2(2,2))*Sqr(msu2(2,2))) + 1.0201*Quad(MuInput)*(-
      21.12034848*Cube(msu2(2,2))*msq2(2,2) + 135.1978992*Cube(msq2(2,2))*msu2(2,2
      ) + 20.56621104*Quad(msq2(2,2)) - 3.6894726399999995*Quad(msu2(2,2)) +
      92.92561487999998*Sqr(msq2(2,2))*Sqr(msu2(2,2)))))/(msq2(2,2)*Quad(1.02*msq2
      (2,2) - 0.98*msu2(2,2))*(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*msu2
      (2,2) + 1.01*Sqr(MuInput))))) + 3*(Log(Sqr(SCALE))*(Log(0.98*msu2(2,2)) - (
      1.0004001600640255*(2.04*msq2(2,2) + 0.98*msu2(2,2))*(1.02*msq2(2,2) + 0.98*
      msu2(2,2) - 2.02*Sqr(MuInput)))/(msq2(2,2)*msu2(2,2))) - 2*Sqr(Log(Sqr(SCALE
      ))) + (0.5002000800320128*(4.121204*(2.04*msq2(2,2) + 0.98*msu2(2,2))*Power6
      (MuInput) - 0.9996*msq2(2,2)*msu2(2,2)*(2.9988*msq2(2,2)*msu2(2,2) + 4.1616*
      Sqr(msq2(2,2)) + 1.9207999999999998*Sqr(msu2(2,2))) - 3.0603*Quad(MuInput)*(
      6.9972*msq2(2,2)*msu2(2,2) + 4.1616*Sqr(msq2(2,2)) + 1.9207999999999998*Sqr(
      msu2(2,2))) + 1.01*Sqr(MuInput)*(4.244832*Cube(msq2(2,2)) +
      1.8823839999999998*Cube(msu2(2,2)) + 17.333064*msu2(2,2)*Sqr(msq2(2,2)) +
      12.734903999999998*msq2(2,2)*Sqr(msu2(2,2)))))/(msq2(2,2)*msu2(2,2)*(-1.02*
      msq2(2,2) + 1.01*Sqr(MuInput))*(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))) + (
      4.0804*PolyLog(2,(0.9900990099009901*(-1.02*msq2(2,2) + 1.01*Sqr(MuInput)))/
      Sqr(MuInput))*Quad(MuInput))/Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput)) + (
      2.04*msq2(2,2)*(-1.02*msq2(2,2) + 2.02*Sqr(MuInput))*Sqr(Log(1.02*msq2(2,2))
      ))/Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput)) + (2*PolyLog(2,1 - (
      1.00990099009901*msq2(2,2))/Sqr(MuInput))*(1.0201*Quad(MuInput) + 2.0604*
      msq2(2,2)*Sqr(MuInput) - 1.0404*Sqr(msq2(2,2))))/Sqr(-1.02*msq2(2,2) + 1.01*
      Sqr(MuInput)) + Log(1.02*msq2(2,2))*(3*Log(Sqr(SCALE)) + (0.5102040816326531
      *(1.030301*(4.08*msq2(2,2) + 14.7*msu2(2,2))*Power6(MuInput) - 1.019592*msu2
      (2,2)*(4.08*msq2(2,2) + 4.9*msu2(2,2))*Sqr(msq2(2,2)) + 1.0302*msq2(2,2)*Sqr
      (MuInput)*(14.994*msq2(2,2)*msu2(2,2) + 4.1616*Sqr(msq2(2,2)) +
      5.7623999999999995*Sqr(msu2(2,2))) - 1.0201*Quad(MuInput)*(13.9944*msq2(2,2)
      *msu2(2,2) + 8.3232*Sqr(msq2(2,2)) + 12.485199999999999*Sqr(msu2(2,2)))))/(
      msu2(2,2)*(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))*Sqr(-1.02*msq2(2,2) + 1.01*
      Sqr(MuInput))) + (Log(0.98*msu2(2,2))*(-2.01938996*msu2(2,2)*Power6(MuInput)
      + 2.08120802*Power8(MuInput) + 0.999698*(-4.08*msq2(2,2) + 0.98*msu2(2,2))*
      msu2(2,2)*Quad(MuInput) + 2.019192*msq2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2
      ))*msu2(2,2)*Sqr(MuInput) - 0.9992001599999999*Sqr(msq2(2,2))*Sqr(msu2(2,2))
      ))/(Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(
      MuInput)))) + Log(1.01*Sqr(MuInput))*((Log(0.98*msu2(2,2))*(2.10181404*msq2(
      2,2)*Power6(MuInput) - 3.12181203*Power8(MuInput) - 1.040502*msq2(2,2)*(1.02
      *msq2(2,2) - 7.84*msu2(2,2))*Quad(MuInput) - 4.038384*msq2(2,2)*(1.02*msq2(2
      ,2) + 0.98*msu2(2,2))*msu2(2,2)*Sqr(MuInput) + 1.9984003199999998*Sqr(msq2(2
      ,2))*Sqr(msu2(2,2))))/(Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*
      msu2(2,2) + 1.01*Sqr(MuInput))) + (Log(1.02*msq2(2,2))*(-2.060602*(1.02*msq2
      (2,2) - 3.92*msu2(2,2))*Power6(MuInput) - 5.20302005*Power8(MuInput) -
      4.038384*msq2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2)*Sqr(MuInput)
      + 1.0201*Quad(MuInput)*(7.9968*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2)) -
      3.8415999999999997*Sqr(msu2(2,2))) + 1.9984003199999998*Sqr(msq2(2,2))*Sqr(
      msu2(2,2))))/(Sqr(-1.02*msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*msu2(2,2) +
      1.01*Sqr(MuInput))) + (1.0104041616646657*Sqr(MuInput)*(-2.08120802*(2.04*
      msq2(2,2) + 0.98*msu2(2,2))*Power8(MuInput) + 1.009596*msq2(2,2)*msu2(2,2)*
      Sqr(MuInput)*(7.9968*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) -
      1.9207999999999998*Sqr(msu2(2,2))) - 1.9984003199999998*(1.02*msq2(2,2) +
      0.98*msu2(2,2))*Sqr(msq2(2,2))*Sqr(msu2(2,2)) + 1.030301*Power6(MuInput)*(
      0.9996*msq2(2,2)*msu2(2,2) + 8.3232*Sqr(msq2(2,2)) + 3.8415999999999997*Sqr(
      msu2(2,2))) - 2.0402*Quad(MuInput)*(2.122416*Cube(msq2(2,2)) +
      0.9411919999999999*Cube(msu2(2,2)) + 4.078368*msu2(2,2)*Sqr(msq2(2,2)) -
      0.9796079999999999*msq2(2,2)*Sqr(msu2(2,2)))))/(msq2(2,2)*msu2(2,2)*Sqr(-
      1.02*msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))
      )) + (0.98*msu2(2,2)*(-0.98*msu2(2,2) + 2.02*Sqr(MuInput))*Sqr(Log(0.98*msu2
      (2,2))))/Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput)) + (2*PolyLog(2,1 - (
      0.9702970297029703*msu2(2,2))/Sqr(MuInput))*(1.0201*Quad(MuInput) + 1.9796*
      msu2(2,2)*Sqr(MuInput) - 0.9603999999999999*Sqr(msu2(2,2))))/Sqr(-0.98*msu2(
      2,2) + 1.01*Sqr(MuInput)) + (0.49019607843137253*Log(0.98*msu2(2,2))*(
      1.030301*(13.26*msq2(2,2) + 1.96*msu2(2,2))*Power6(MuInput) -
      0.9796079999999999*msq2(2,2)*(1.02*msq2(2,2) + 1.96*msu2(2,2))*Sqr(msu2(2,2)
      ) + 0.9898*msu2(2,2)*Sqr(MuInput)*(8.996400000000001*msq2(2,2)*msu2(2,2) +
      4.1616*Sqr(msq2(2,2)) + 1.9207999999999998*Sqr(msu2(2,2))) - 1.0201*Quad(
      MuInput)*(13.9944*msq2(2,2)*msu2(2,2) + 9.3636*Sqr(msq2(2,2)) +
      3.8415999999999997*Sqr(msu2(2,2)))))/(msq2(2,2)*(-1.02*msq2(2,2) + 1.01*Sqr(
      MuInput))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput))) + (Sqr(Log(1.01*Sqr(
      MuInput)))*(-4.03877992*msu2(2,2)*Power6(MuInput) + 4.16241604*Power8(
      MuInput) + 1.999396*(-4.08*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2)*Quad(
      MuInput) + 4.038384*msq2(2,2)*(1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2)*
      Sqr(MuInput) - 1.9984003199999998*Sqr(msq2(2,2))*Sqr(msu2(2,2))))/(Sqr(-1.02
      *msq2(2,2) + 1.01*Sqr(MuInput))*Sqr(-0.98*msu2(2,2) + 1.01*Sqr(MuInput)))) +
      (3*Sqr(AtInput + MuInput*TanBeta)*(-0.9803921568627451/msq2(2,2) + Log(Sqr(
      SCALE))*(-0.9803921568627451/msq2(2,2) - 2.0408163265306123/msu2(2,2)) -
      2.0408163265306123/msu2(2,2) + (0.5206164098292378*Log(1.02*msq2(2,2))*((-
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + 0.9702989999999999*Power6(mAInput) -
      2.9402999999999997*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Quad(mAInput) +
      2.9699999999999998*Sqr(mAInput)*(1.0404*Sqr(msq2(2,2)) - 0.9603999999999999*
      Sqr(msu2(2,2))))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)) +
      TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*(-0.98*msu2(2,2)*(
      0.9801*Quad(mAInput) - 1.9404*msu2(2,2)*Sqr(mAInput) - 1.0404*Sqr(msq2(2,2))
      + 0.9603999999999999*Sqr(msu2(2,2))) + (1.02*msq2(2,2) + 1.96*msu2(2,2) -
      0.99*Sqr(mAInput))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))))
      )/(Sqr(msu2(2,2))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*
      TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (
      0.5002000800320128*Log(0.98*msu2(2,2))*((Cube(1.02*msq2(2,2) - 0.98*msu2(2,2
      )) + 0.9702989999999999*Power6(mAInput) - 0.9801*(1.02*msq2(2,2) + 4.9*msu2(
      2,2))*Quad(mAInput) - 0.99*Sqr(mAInput)*(3.9984*msq2(2,2)*msu2(2,2) + 1.0404
      *Sqr(msq2(2,2)) - 4.802*Sqr(msu2(2,2))))*TDelta(0.99*Sqr(mAInput),0.98*msu2(
      2,2),1.02*msq2(2,2)) - TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2
      ))*(1.9992*msq2(2,2)*msu2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.99*Sqr(
      mAInput)) + (1.02*msq2(2,2) - 2.94*msu2(2,2) + 0.99*Sqr(mAInput))*TDelta(
      0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))))/(msq2(2,2)*msu2(2,2)*
      TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (0.5104082449306253*Log(0.99*Sqr(
      mAInput))*((-0.9702989999999999*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Power6(
      mAInput) + Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) - 0.99*(3.06*msq2(2,2) +
      4.9*msu2(2,2))*Sqr(mAInput)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + 0.9801*
      Quad(mAInput)*(3.9984*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) + 4.802*
      Sqr(msu2(2,2))))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)) +
      TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*(-0.9996*msq2(2,2)*
      msu2(2,2)*(-0.9801*Quad(mAInput) + Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      2.9988*msq2(2,2)*msu2(2,2) + 0.99*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(
      mAInput) - 1.0404*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2)))*TDelta
      (0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))))/(msq2(2,2)*Sqr(msu2(2,2
      ))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (0.5312412345196305*(-
      3.8811959999999996*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Power6(mAInput) +
      0.96059601*Power8(mAInput) + Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) - 3.96*(
      1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(mAInput)*Sqr(1.02*msq2(2,2) - 0.98*msu2
      (2,2)) + 1.9602*Quad(mAInput)*(1.9992*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(
      2,2)) + 0.9603999999999999*Sqr(msu2(2,2))) + 2*Sqr(TDelta(0.99*Sqr(mAInput),
      1.02*msq2(2,2),0.98*msu2(2,2))) - 3*(0.9801*Quad(mAInput) - 1.98*(1.02*msq2(
      2,2) + 0.98*msu2(2,2))*Sqr(mAInput) + Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)))*
      TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2)))*TPhi(0.99*Sqr(
      mAInput),1.02*msq2(2,2),0.98*msu2(2,2)))/(Cube(msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),1.02*msq2(2,2),0.98*msu2(2,2))) + (0.49019607843137253*(-Sqr(0.99*
      Sqr(mAInput) + 1.02*msq2(2,2) - 0.98*msu2(2,2)) + TDelta(0.99*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*
      msq2(2,2)))/(msq2(2,2)*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2
      )))) + 3*(1.5 - 2.5*Log(0.98*msu2(2,2)) + Log(1.02*msq2(2,2))*(-4.5 - 3*Log(
      Sqr(SCALE))) + Sqr(3.141592653589793) - (0.9705882352941176*Sqr(mAInput))/
      msq2(2,2) - (2.0204081632653064*Sqr(mAInput))/msu2(2,2) + Log(0.99*Sqr(
      mAInput))*(7 - 3*Log(1.02*msq2(2,2)) - 3*Log(0.98*msu2(2,2)) + 0.99*(
      0.9803921568627451/msq2(2,2) + 2.0408163265306123/msu2(2,2))*Sqr(mAInput)) +
      Log(Sqr(SCALE))*(-Log(0.98*msu2(2,2)) - (0.9903961584633852*(2.04*msq2(2,2)
      + 0.98*msu2(2,2))*Sqr(mAInput))/(msq2(2,2)*msu2(2,2))) + 3*Sqr(Log(0.99*Sqr(
      mAInput))) + 2*Sqr(Log(Sqr(SCALE))) + 3*Sqr(Log(1.02*msq2(2,2))) + 2*Sqr(Log
      (0.98*msu2(2,2))) + (0.9611687812379854*(0.9801*Quad(mAInput) - 7.0686*msq2(
      2,2)*Sqr(mAInput) - TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))
      *TPhi(0.99*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/Sqr(msq2(2,2)) + (
      1.0412328196584757*(0.9801*Quad(mAInput) - 5.8212*msu2(2,2)*Sqr(mAInput) -
      TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))*TPhi(0.99*Sqr(
      mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))/Sqr(msu2(2,2))) + 3*(AtInput -
      MuInput/TanBeta)*(AtInput + MuInput*TanBeta)*(Log(0.99*Sqr(mAInput))*((2*Log
      (1.02*msq2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2*Log(0.98*msu2(2,2))
      )/(1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*(-12/(1.02*msq2(2
      ,2) - 0.98*msu2(2,2)) - (2*Log(0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(
      2,2)) - (12*Log(Sqr(SCALE)))/(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (12*Log(
      0.98*msu2(2,2)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (12*Log(0.98*msu2(2,2))
      *Log(Sqr(SCALE)))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (6*Sqr(Log(1.02*msq2(2
      ,2))))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (4*Sqr(Log(0.98*msu2(2,2))))/(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) - (1.9223375624759709*(-0.99*(-
      7.140000000000001*msq2(2,2) + 0.99*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.99*
      Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(mAInput),1.02*
      msq2(2,2),1.02*msq2(2,2)))/((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msq2(2,2))
      ) + (1.9607843137254901*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.99*Sqr(mAInput)
      )*TPhi(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*(1.02*
      msq2(2,2) - 0.98*msu2(2,2))) + (2.0824656393169514*(-0.99*(-5.88*msu2(2,2) +
      0.99*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),
      0.98*msu2(2,2)))*TPhi(0.99*Sqr(mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))/((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msu2(2,2)))) + 3*(AtInput + MuInput*
      TanBeta)*Cube(AtInput - MuInput/TanBeta)*((-12*Log(0.98*msu2(2,2))*(1.02*
      msq2(2,2) + 2.94*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(
      1.02*msq2(2,2))*((12*Log(Sqr(SCALE))*(1.02*msq2(2,2) + 0.98*msu2(2,2)))/Cube
      (1.02*msq2(2,2) - 0.98*msu2(2,2)) + (12*(3.06*msq2(2,2) + 0.98*msu2(2,2)))/
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (2*Log(0.98*msu2(2,2))*(3.06*msq2(2,
      2) + 0.98*msu2(2,2) - 1.98*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,
      2))) + Log(0.99*Sqr(mAInput))*((2*Log(1.02*msq2(2,2))*(1.02*msq2(2,2) - 0.98
      *msu2(2,2) - 5.9399999999999995*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + (2*Log(0.98*msu2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2(2,2) +
      5.9399999999999995*Sqr(mAInput)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (
      2*(-5.1*msq2(2,2) - 2.94*msu2(2,2) + 3.96*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,
      2))))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4*(1.02*msq2(2,2) + 0.98*msu2
      (2,2) - 0.99*Sqr(mAInput))*Sqr(Log(0.98*msu2(2,2))))/Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2)) + Log(Sqr(SCALE))*((-12*Log(0.98*msu2(2,2))*(1.02*msq2(2,2)
      + 0.98*msu2(2,2)))/Cube(1.02*msq2(2,2) - 0.98*msu2(2,2)) - 24/Sqr(1.02*msq2(
      2,2) - 0.98*msu2(2,2))) - 48/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (
      1.9223375624759709*(0.99*(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-
      7.140000000000001*msq2(2,2) + 0.99*Sqr(mAInput))*Sqr(mAInput) + (-5.1*msq2(2
      ,2) + 0.98*msu2(2,2))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)
      ))*TPhi(0.99*Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))/(Cube(1.02*msq2(2,
      2) - 0.98*msu2(2,2))*Sqr(msq2(2,2))) - (1.9607843137254901*((1.02*msq2(2,2)
      - 0.98*msu2(2,2))*(1.02*msq2(2,2) - 0.98*msu2(2,2) + 0.99*Sqr(mAInput)) - 2*
      TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,
      2))*msq2(2,2)) + (2.0824656393169514*(0.99*(1.02*msq2(2,2) - 0.98*msu2(2,2))
      *(-5.88*msu2(2,2) + 0.99*Sqr(mAInput))*Sqr(mAInput) - (1.02*msq2(2,2) - 2.94
      *msu2(2,2))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))*TPhi(
      0.99*Sqr(mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))/(Cube(1.02*msq2(2,2) -
      0.98*msu2(2,2))*Sqr(msu2(2,2)))) + Sqr(AtInput - MuInput/TanBeta)*(3*Sqr(
      AtInput + MuInput*TanBeta)*((3*Sqr(Log(1.02*msq2(2,2))))/Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2)) + (2*Sqr(Log(0.98*msu2(2,2))))/Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + (4.08*msq2(2,2) - 1.96*msu2(2,2))/(1.019592*msu2(2,2)*Sqr(msq2(
      2,2)) - 0.9796079999999999*msq2(2,2)*Sqr(msu2(2,2))) + Log(Sqr(SCALE))*((2*
      Log(0.98*msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (4.08*msq2(2,2)
      - 1.96*msu2(2,2))/(1.019592*msu2(2,2)*Sqr(msq2(2,2)) - 0.9796079999999999*
      msq2(2,2)*Sqr(msu2(2,2)))) + (1.0004001600640255*Log(0.98*msu2(2,2))*(-((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*(0.9801*(2.04*msq2(2,2) + 0.98*msu2(2,2))*
      Quad(mAInput) - 4.0392*msq2(2,2)*(1.02*msq2(2,2) + 1.96*msu2(2,2))*Sqr(
      mAInput) + (2.04*msq2(2,2) - 0.98*msu2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(
      2,2)))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + TDelta(
      0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*(1.9992*msq2(2,2)*(1.02*
      msq2(2,2) - 0.98*msu2(2,2))*msu2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2) +
      0.99*Sqr(mAInput)) + (0.9996*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2)))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))))/(msq2(2,2)*msu2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))
      *TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + Log(1.02*msq2(2,2))*((-5*Log(0.98
      *msu2(2,2)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2*Log(Sqr(SCALE)))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.0412328196584757*((1.02*msq2(2,2) -
      0.98*msu2(2,2))*(-0.9702989999999999*Power6(mAInput) + 0.9801*(3.06*msq2(2,2
      ) + 3.92*msu2(2,2))*Quad(mAInput) + 1.02*msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98
      *msu2(2,2)) - 0.99*Sqr(mAInput)*(1.9992*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(
      msq2(2,2)) + 6.722799999999999*Sqr(msu2(2,2))))*TDelta(0.99*Sqr(mAInput),
      0.98*msu2(2,2),1.02*msq2(2,2)) + TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),
      0.98*msu2(2,2))*(-0.98*(1.02*msq2(2,2) - 0.98*msu2(2,2))*msu2(2,2)*(-0.9801*
      Quad(mAInput) + 1.9404*msu2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2))) + (-1.9992*msq2(2,2)*msu2(2,2) + 0.99*(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput) - 1.0404*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))))/(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msu2(2,2))*
      TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + Log(0.99*Sqr(mAInput))*(Log(1.02
      *msq2(2,2))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - Log(0.98*msu2(2,2))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.0208164898612506*((0.98970498*msq2(2,2
      )*Power6(mAInput) - Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + 1.0098*msq2(2,2)
      *Sqr(mAInput)*(-1.9992*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2))) + 0.9801*Quad(mAInput)*(-1.9992*msq2(2,2)
      *msu2(2,2) - 3.1212*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))*
      TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)) + TDelta(0.99*Sqr(
      mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*(0.9996*msq2(2,2)*msu2(2,2)*(-0.9801
      *Quad(mAInput) + Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (-2.9988*msq2(2,2)*
      msu2(2,2) - 1.0098*msq2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2)))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)))))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msu2(2,2
      ))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + (0.9611687812379854*(0.99*(-
      7.140000000000001*msq2(2,2) + 0.99*Sqr(mAInput))*Sqr(mAInput) - TDelta(0.99*
      Sqr(mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(mAInput),1.02*
      msq2(2,2),1.02*msq2(2,2)))/(Sqr(msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,
      2))) + (1.062482469039261*(-((1.02*msq2(2,2) - 0.98*msu2(2,2))*(-
      3.8811959999999996*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Power6(mAInput) +
      0.96059601*Power8(mAInput) + Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) -
      4.119984*(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(mAInput)*Sqr(msq2(2,2)) +
      0.9801*Quad(mAInput)*(3.9984*msq2(2,2)*msu2(2,2) + 6.2424*Sqr(msq2(2,2)) +
      5.7623999999999995*Sqr(msu2(2,2))))) - 2*(1.02*msq2(2,2) - 1.96*msu2(2,2))*
      Sqr(TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))) + (0.9801*(3.06
      *msq2(2,2) - 4.9*msu2(2,2))*Quad(mAInput) + (3.06*msq2(2,2) - 4.9*msu2(2,2))
      *Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + 0.99*Sqr(mAInput)*(3.9984*msq2(2,2)*
      msu2(2,2) - 6.2424*Sqr(msq2(2,2)) + 13.445599999999999*Sqr(msu2(2,2))))*
      TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2)))*TPhi(0.99*Sqr(
      mAInput),1.02*msq2(2,2),0.98*msu2(2,2)))/(Cube(msu2(2,2))*Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))) +
      (0.9803921568627451*((1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(0.99*Sqr(mAInput)
      + 1.02*msq2(2,2) - 0.98*msu2(2,2)) + 0.99*Sqr(mAInput)*TDelta(0.99*Sqr(
      mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(mAInput),0.98*msu2(2,
      2),1.02*msq2(2,2)))/(msq2(2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(
      0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (1.0412328196584757*(
      0.99*(-5.88*msu2(2,2) + 0.99*Sqr(mAInput))*Sqr(mAInput) - TDelta(0.99*Sqr(
      mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))*TPhi(0.99*Sqr(mAInput),0.98*msu2(2,
      2),0.98*msu2(2,2)))/(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msu2(2,2)))) +
      3*((-2*Sqr(Log(0.98*msu2(2,2))))/(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(
      1.02*msq2(2,2))*((1.02*msq2(2,2) + 2.94*msu2(2,2) - 1.98*Sqr(mAInput))/Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2)) + (Log(0.98*msu2(2,2))*(1.02*msq2(2,2) +
      0.98*msu2(2,2) - 0.99*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) -
      (2*Log(Sqr(SCALE))*(-1.02*msq2(2,2) + 0.99*Sqr(mAInput)))/Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + Log(0.99*Sqr(mAInput))*((1.9807923169267705*(-2.04*msq2
      (2,2) + 0.98*msu2(2,2))*Sqr(mAInput))/(msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2
      (2,2))*msu2(2,2)) + (Log(1.02*msq2(2,2))*(-5.1*msq2(2,2) + 4.9*msu2(2,2) +
      0.99*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (Log(0.98*msu2(2,
      2))*(-5.1*msq2(2,2) + 4.9*msu2(2,2) + 0.99*Sqr(mAInput)))/Sqr(1.02*msq2(2,2)
      - 0.98*msu2(2,2))) + (Log(0.98*msu2(2,2))*(1.02*msq2(2,2) - 4.9*msu2(2,2) +
      1.98*Sqr(mAInput)))/Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + ((1.02*msq2(2,2)
      - 2.94*msu2(2,2) + 0.99*Sqr(mAInput))*Sqr(Log(1.02*msq2(2,2))))/Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2)) + (-3.9984*msq2(2,2)*msu2(2,2) + 4.0392*msq2(2,2
      )*Sqr(mAInput) - 1.9404*msu2(2,2)*Sqr(mAInput))/(1.019592*msu2(2,2)*Sqr(msq2
      (2,2)) - 0.9796079999999999*msq2(2,2)*Sqr(msu2(2,2))) + Log(Sqr(SCALE))*((2*
      Log(0.98*msu2(2,2))*(-1.02*msq2(2,2) + 0.99*Sqr(mAInput)))/Sqr(1.02*msq2(2,2
      ) - 0.98*msu2(2,2)) + (-1.9992*msq2(2,2)*msu2(2,2) + 4.0392*msq2(2,2)*Sqr(
      mAInput) - 1.9404*msu2(2,2)*Sqr(mAInput))/(1.019592*msu2(2,2)*Sqr(msq2(2,2))
      - 0.9796079999999999*msq2(2,2)*Sqr(msu2(2,2)))) + (0.9611687812379854*(0.99*
      (1.02*msq2(2,2) - 0.98*msu2(2,2))*(-7.140000000000001*msq2(2,2) + 0.99*Sqr(
      mAInput))*Sqr(mAInput) + (-2.04*msq2(2,2) + 0.98*msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(mAInput),1.02*msq2(2,
      2),1.02*msq2(2,2)))/(Sqr(msq2(2,2))*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      (0.9803921568627451*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))*
      TPhi(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*Sqr(1.02*
      msq2(2,2) - 0.98*msu2(2,2))) + (1.0412328196584757*(-0.99*(-5.88*msu2(2,2) +
      0.99*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),
      0.98*msu2(2,2)))*TPhi(0.99*Sqr(mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))/((
      1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msu2(2,2))))) + Quad(AtInput - MuInput/
      TanBeta)*(3*(Log(Sqr(SCALE))*((1.0004001600640255*(5.9976*msq2(2,2)*(1.02*
      msq2(2,2) - 0.98*msu2(2,2))*msu2(2,2) + 0.99*Sqr(mAInput)*(-4.998*msq2(2,2)*
      msu2(2,2) - 2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2)))))/(
      Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*msq2(2,2)*msu2(2,2)) + (3*Log(0.98*
      msu2(2,2))*(-2.0196*msq2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) +
      (1.0004001600640255*(2.9988*msq2(2,2)*(1.02*msq2(2,2) - 0.98*msu2(2,2))*msu2
      (2,2) + 0.99*Sqr(mAInput)*(-4.998*msq2(2,2)*msu2(2,2) - 2.0808*Sqr(msq2(2,2)
      ) + 0.9603999999999999*Sqr(msu2(2,2)))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,
      2))*msq2(2,2)*msu2(2,2)) + Log(0.99*Sqr(mAInput))*((1.0004001600640255*(
      5.9976*msq2(2,2)*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*msu2(2,2) + 0.99*Sqr(
      mAInput)*(4.998*msq2(2,2)*msu2(2,2) + 2.0808*Sqr(msq2(2,2)) -
      0.9603999999999999*Sqr(msu2(2,2)))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      msq2(2,2)*msu2(2,2)) + (3*Log(1.02*msq2(2,2))*(-2.0196*msq2(2,2)*Sqr(mAInput
      ) + 1.0404*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2))))/Quad(1.02*
      msq2(2,2) - 0.98*msu2(2,2)) + (3*Log(0.98*msu2(2,2))*(2.0196*msq2(2,2)*Sqr(
      mAInput) - 1.0404*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2))) + Log(1.02*msq2(2,2))*((3*Log(Sqr(SCALE))*
      (2.0196*msq2(2,2)*Sqr(mAInput) - 1.0404*Sqr(msq2(2,2)) + 0.9603999999999999*
      Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (1.5*(4.0392*msq2(2
      ,2)*Sqr(mAInput) - 1.0404*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))
      ))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) + (1.5*Log(0.98*msu2(2,2))*(-
      4.0392*msq2(2,2)*Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) - 0.9603999999999999*
      Sqr(msu2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))) + 3*Sqr(AtInput +
      MuInput*TanBeta)*(((-11.22*msq2(2,2) - 2.94*msu2(2,2) + 6.93*Sqr(mAInput))*
      Sqr(Log(1.02*msq2(2,2))))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) - (2*(1.02*
      msq2(2,2) + 2.94*msu2(2,2) - 1.98*Sqr(mAInput))*Sqr(Log(0.98*msu2(2,2))))/
      Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Log(Sqr(SCALE))*((-6.12*Log(0.98*
      msu2(2,2))*msq2(2,2))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (
      1.0004001600640255*(-4.998*msq2(2,2)*msu2(2,2) - 2.0808*Sqr(msq2(2,2)) +
      0.9603999999999999*Sqr(msu2(2,2))))/(Cube(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      msq2(2,2)*msu2(2,2))) + (1.0004001600640255*(-10.9956*msq2(2,2)*msu2(2,2) -
      2.0808*Sqr(msq2(2,2)) + 0.9603999999999999*Sqr(msu2(2,2))))/(Cube(1.02*msq2(
      2,2) - 0.98*msu2(2,2))*msq2(2,2)*msu2(2,2)) - (0.5002000800320128*Log(0.98*
      msu2(2,2))*(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-3.183624*Cube(msq2(2,2))
      - 2.8235759999999996*Cube(msu2(2,2)) + 0.9702989999999999*Power6(mAInput) -
      0.9801*(5.1*msq2(2,2) + 2.94*msu2(2,2))*Quad(mAInput) + 11.215512*msu2(2,2)*
      Sqr(msq2(2,2)) - 4.89804*msq2(2,2)*Sqr(msu2(2,2)) + 0.99*Sqr(mAInput)*(
      3.9984*msq2(2,2)*msu2(2,2) + 7.2828*Sqr(msq2(2,2)) + 4.802*Sqr(msu2(2,2))))*
      TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)) + TDelta(0.99*Sqr(
      mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*(1.9992*msq2(2,2)*msu2(2,2)*(1.02*
      msq2(2,2) - 0.98*msu2(2,2) + 0.99*Sqr(mAInput))*Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2)) + (3.183624*Cube(msq2(2,2)) + 0.9411919999999999*Cube(msu2(2,2))
      + 7.137144*msu2(2,2)*Sqr(msq2(2,2)) - 0.99*Sqr(mAInput)*Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2)) + 12.734903999999998*msq2(2,2)*Sqr(msu2(2,2)))*TDelta(0.99*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))))/(msq2(2,2)*msu2(2,2)*Quad(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),
      0.98*msu2(2,2))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) +
      Log(0.99*Sqr(mAInput))*((-3*Log(1.02*msq2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2
      (2,2) + 0.99*Sqr(mAInput)))/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (3*Log(
      0.98*msu2(2,2))*(-1.02*msq2(2,2) + 0.98*msu2(2,2) + 0.99*Sqr(mAInput)))/Quad
      (1.02*msq2(2,2) - 0.98*msu2(2,2)) + (0.5104082449306253*((1.02*msq2(2,2) -
      0.98*msu2(2,2))*(1.061208*Cube(msq2(2,2)) + 2.8235759999999996*Cube(msu2(2,2
      )) - 0.9702989999999999*Power6(mAInput) + 2.9402999999999997*(1.02*msq2(2,2)
      + 0.98*msu2(2,2))*Quad(mAInput) - 3.058776*msu2(2,2)*Sqr(msq2(2,2)) -
      0.9796079999999999*msq2(2,2)*Sqr(msu2(2,2)) - 0.99*Sqr(mAInput)*(3.1212*Sqr(
      msq2(2,2)) + 4.802*Sqr(msu2(2,2))))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),
      1.02*msq2(2,2)) + TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*(-
      0.9996*msq2(2,2)*msu2(2,2)*(-0.9801*Quad(mAInput) + Sqr(1.02*msq2(2,2) -
      0.98*msu2(2,2))) + (2.9988*msq2(2,2)*msu2(2,2) + 0.99*(1.02*msq2(2,2) - 0.98
      *msu2(2,2))*Sqr(mAInput) - 1.0404*Sqr(msq2(2,2)) + 2.8811999999999998*Sqr(
      msu2(2,2)))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))))/(msq2
      (2,2)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2
      ,2),1.02*msq2(2,2)))) + Log(1.02*msq2(2,2))*((6.12*Log(Sqr(SCALE))*msq2(2,2)
      )/Quad(1.02*msq2(2,2) - 0.98*msu2(2,2)) + (Log(0.98*msu2(2,2))*(13.26*msq2(2
      ,2) + 8.82*msu2(2,2) - 10.89*Sqr(mAInput)))/Quad(1.02*msq2(2,2) - 0.98*msu2(
      2,2)) + (0.5206164098292378*(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-1.061208
      *Cube(msq2(2,2)) - 8.470728*Cube(msu2(2,2)) + 0.9702989999999999*Power6(
      mAInput) - 0.9801*(3.06*msq2(2,2) + 4.9*msu2(2,2))*Quad(mAInput) + 1.019592*
      msu2(2,2)*Sqr(msq2(2,2)) + 8.816472*msq2(2,2)*Sqr(msu2(2,2)) + 0.99*Sqr(
      mAInput)*(3.9984*msq2(2,2)*msu2(2,2) + 3.1212*Sqr(msq2(2,2)) + 8.6436*Sqr(
      msu2(2,2))))*TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)) +
      TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*(0.98*msu2(2,2)*Sqr(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*(-0.9801*Quad(mAInput) + 1.9404*msu2(2,2)*
      Sqr(mAInput) + 1.0404*Sqr(msq2(2,2)) - 0.9603999999999999*Sqr(msu2(2,2))) +
      (1.061208*Cube(msq2(2,2)) + 3.7647679999999997*Cube(msu2(2,2)) + 2.039184*
      msu2(2,2)*Sqr(msq2(2,2)) - 0.99*Sqr(mAInput)*Sqr(1.02*msq2(2,2) - 0.98*msu2(
      2,2)) + 16.653336*msq2(2,2)*Sqr(msu2(2,2)))*TDelta(0.99*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2)))))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msu2
      (2,2))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))*TDelta(0.99*
      Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))) + (0.9611687812379854*(0.99*(
      1.02*msq2(2,2) - 0.98*msu2(2,2))*(-7.140000000000001*msq2(2,2) + 0.99*Sqr(
      mAInput))*Sqr(mAInput) + (-8.16*msq2(2,2) + 0.98*msu2(2,2))*TDelta(0.99*Sqr(
      mAInput),1.02*msq2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(mAInput),1.02*msq2(2,
      2),1.02*msq2(2,2)))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(msq2(2,2))) +
      (0.5312412345196305*(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2))*(-
      3.8811959999999996*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Power6(mAInput) +
      0.96059601*Power8(mAInput) - 3.96*(1.02*msq2(2,2) + 0.98*msu2(2,2))*Sqr(
      mAInput)*Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)) + Sqr(1.02*msq2(2,2) - 0.98*
      msu2(2,2))*(-1.9992*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2)) -
      2.8811999999999998*Sqr(msu2(2,2))) + 0.9801*Quad(mAInput)*(3.9984*msq2(2,2)*
      msu2(2,2) + 6.2424*Sqr(msq2(2,2)) + 5.7623999999999995*Sqr(msu2(2,2)))) + 2*
      (-4.998*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(msq2(2,2)) + 11.524799999999999*Sqr
      (msu2(2,2)))*Sqr(TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2))) -
      (1.02*msq2(2,2) - 0.98*msu2(2,2))*(3.183624*Cube(msq2(2,2)) -
      16.000263999999998*Cube(msu2(2,2)) + 2.9402999999999997*(1.02*msq2(2,2) -
      2.94*msu2(2,2))*Quad(mAInput) - 15.29388*msu2(2,2)*Sqr(msq2(2,2)) -
      5.9399999999999995*Sqr(mAInput)*(-1.9992*msq2(2,2)*msu2(2,2) + 1.0404*Sqr(
      msq2(2,2)) - 2.8811999999999998*Sqr(msu2(2,2))) + 28.408631999999997*msq2(2,
      2)*Sqr(msu2(2,2)))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2)))*
      TPhi(0.99*Sqr(mAInput),1.02*msq2(2,2),0.98*msu2(2,2)))/(Cube(msu2(2,2))*Quad
      (1.02*msq2(2,2) - 0.98*msu2(2,2))*TDelta(0.99*Sqr(mAInput),1.02*msq2(2,2),
      0.98*msu2(2,2))) - (0.49019607843137253*(Sqr(1.02*msq2(2,2) - 0.98*msu2(2,2)
      )*Sqr(0.99*Sqr(mAInput) + 1.02*msq2(2,2) - 0.98*msu2(2,2)) - 6*Sqr(TDelta(
      0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (1.02*msq2(2,2) - 0.98*
      msu2(2,2))*(3.06*msq2(2,2) - 2.94*msu2(2,2) + 3.96*Sqr(mAInput))*TDelta(0.99
      *Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2)))*TPhi(0.99*Sqr(mAInput),0.98*
      msu2(2,2),1.02*msq2(2,2)))/(msq2(2,2)*Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*
      TDelta(0.99*Sqr(mAInput),0.98*msu2(2,2),1.02*msq2(2,2))) + (
      1.0412328196584757*(0.99*(-1.02*msq2(2,2) + 0.98*msu2(2,2))*(-5.88*msu2(2,2)
      + 0.99*Sqr(mAInput))*Sqr(mAInput) + (1.02*msq2(2,2) - 4.9*msu2(2,2))*TDelta(
      0.99*Sqr(mAInput),0.98*msu2(2,2),0.98*msu2(2,2)))*TPhi(0.99*Sqr(mAInput),
      0.98*msu2(2,2),0.98*msu2(2,2)))/(Quad(1.02*msq2(2,2) - 0.98*msu2(2,2))*Sqr(
      msu2(2,2))))))/(1 + Sqr(TanBeta))))/Sqr(TanBeta)))/Quad(3.141592653589793)),
      0) + IF(TwoLoopAtauAtau >= 1, (0.01171875*(((AtauInput - MuInput*TanBeta)*
      Quad(Yd(2,2))*((AbInput - MuInput*TanBeta)*(Log(1.03*msq2(2,2))*(4/(0.97*
      msd2(2,2) - 1.03*msq2(2,2)) + (4.04*Log(1.01*mse2(2,2))*mse2(2,2))/((-1.01*
      mse2(2,2) + 0.99*msl2(2,2))*(0.97*msd2(2,2) - 1.03*msq2(2,2)))) + (4*Log(
      1.03*msq2(2,2))*Log(Sqr(SCALE)))/(0.97*msd2(2,2) - 1.03*msq2(2,2)) + (3.96*
      Log(0.99*msl2(2,2))*Log(1.03*msq2(2,2))*msl2(2,2))/((1.01*mse2(2,2) - 0.99*
      msl2(2,2))*(0.97*msd2(2,2) - 1.03*msq2(2,2))) + Log(0.97*msd2(2,2))*((4.04*
      Log(1.01*mse2(2,2))*mse2(2,2))/((1.01*mse2(2,2) - 0.99*msl2(2,2))*(0.97*msd2
      (2,2) - 1.03*msq2(2,2))) + (3.96*Log(0.99*msl2(2,2))*msl2(2,2))/((-1.01*mse2
      (2,2) + 0.99*msl2(2,2))*(0.97*msd2(2,2) - 1.03*msq2(2,2))) + 4/(-0.97*msd2(2
      ,2) + 1.03*msq2(2,2)) + (4*Log(Sqr(SCALE)))/(-0.97*msd2(2,2) + 1.03*msq2(2,2
      )))) + Cube(AbInput - MuInput*TanBeta)*(Log(1.03*msq2(2,2))*((-4*(0.97*msd2(
      2,2) + 1.03*msq2(2,2)))/Cube(0.97*msd2(2,2) - 1.03*msq2(2,2)) - (4.04*Log(
      1.01*mse2(2,2))*mse2(2,2)*(0.97*msd2(2,2) + 1.03*msq2(2,2)))/(Cube(0.97*msd2
      (2,2) - 1.03*msq2(2,2))*(-1.01*mse2(2,2) + 0.99*msl2(2,2)))) + Log(0.97*msd2
      (2,2))*((4*(0.97*msd2(2,2) + 1.03*msq2(2,2)))/Cube(0.97*msd2(2,2) - 1.03*
      msq2(2,2)) + (4*Log(Sqr(SCALE))*(0.97*msd2(2,2) + 1.03*msq2(2,2)))/Cube(0.97
      *msd2(2,2) - 1.03*msq2(2,2)) + (4.04*Log(1.01*mse2(2,2))*mse2(2,2)*(0.97*
      msd2(2,2) + 1.03*msq2(2,2)))/(Cube(0.97*msd2(2,2) - 1.03*msq2(2,2))*(-1.01*
      mse2(2,2) + 0.99*msl2(2,2))) + (3.96*Log(0.99*msl2(2,2))*msl2(2,2)*(0.97*
      msd2(2,2) + 1.03*msq2(2,2)))/(Cube(0.97*msd2(2,2) - 1.03*msq2(2,2))*(1.01*
      mse2(2,2) - 0.99*msl2(2,2)))) + Log(Sqr(SCALE))*((-4*Log(1.03*msq2(2,2))*(
      0.97*msd2(2,2) + 1.03*msq2(2,2)))/Cube(0.97*msd2(2,2) - 1.03*msq2(2,2)) - 8/
      Sqr(0.97*msd2(2,2) - 1.03*msq2(2,2))) + Log(0.99*msl2(2,2))*((3.96*Log(1.03*
      msq2(2,2))*msl2(2,2)*(0.97*msd2(2,2) + 1.03*msq2(2,2)))/(Cube(0.97*msd2(2,2)
      - 1.03*msq2(2,2))*(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + (7.92*msl2(2,2))/((-
      1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.97*msd2(2,2) - 1.03*msq2(2,2)))) - 8/
      Sqr(0.97*msd2(2,2) - 1.03*msq2(2,2)) + (8.08*Log(1.01*mse2(2,2))*mse2(2,2))/
      ((1.01*mse2(2,2) - 0.99*msl2(2,2))*Sqr(0.97*msd2(2,2) - 1.03*msq2(2,2)))))*
      Sqr(Ye(2,2)))/(Quad(1 + (0.0625*(1 + Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/
      Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta
      )) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*
      Sqr(Yu(2,2)))/(Sqr(3.141592653589793)*Sqr(TanBeta)) + 0.5*((-0.03125*Sqr(
      AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,2))
      )))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(3.141592653589793)) - (0.03125*Sqr(
      AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2))
      )))/(Sqrt(Abs(msq2(2,2)*msu2(2,2)))*Sqr(3.141592653589793))) + (0.0625*(1 +
      Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.375*(-1 +
      2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6
      )(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2))/MuInput)))/Sqr(
      3.141592653589793) + (0.08333333333333333*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(
      SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input)
      - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2)
      )/M3Input))/M3Input))/Sqr(3.141592653589793) + (0.0625*(1 + Sqr(TanBeta))*
      Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/
      TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/
      MuInput))/(Sqr(3.141592653589793)*Sqr(TanBeta)))*Sqrt(1/(1 + Sqr(TanBeta)))*
      Sqrt(Sqr(TanBeta)/(1 + Sqr(TanBeta)))) + (Quad(Ye(2,2))*((AbInput - MuInput*
      TanBeta)*(AtauInput - MuInput*TanBeta)*((4*Log(1.01*mse2(2,2)))/(-1.01*mse2(
      2,2) + 0.99*msl2(2,2)) + (4*Log(1.01*mse2(2,2))*Log(Sqr(SCALE)))/(-1.01*mse2
      (2,2) + 0.99*msl2(2,2)) + Log(0.97*msd2(2,2))*((3.88*Log(1.01*mse2(2,2))*
      msd2(2,2))/((1.01*mse2(2,2) - 0.99*msl2(2,2))*(0.97*msd2(2,2) - 1.03*msq2(2,
      2))) + (3.88*Log(0.99*msl2(2,2))*msd2(2,2))/((-1.01*mse2(2,2) + 0.99*msl2(2,
      2))*(0.97*msd2(2,2) - 1.03*msq2(2,2)))) + (4.12*Log(1.01*mse2(2,2))*Log(1.03
      *msq2(2,2))*msq2(2,2))/((-1.01*mse2(2,2) + 0.99*msl2(2,2))*(0.97*msd2(2,2) -
      1.03*msq2(2,2))) + Log(0.99*msl2(2,2))*(4/(1.01*mse2(2,2) - 0.99*msl2(2,2))
      + (4*Log(Sqr(SCALE)))/(1.01*mse2(2,2) - 0.99*msl2(2,2)) + (4.12*Log(1.03*
      msq2(2,2))*msq2(2,2))/((-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-0.97*msd2(2,2) +
      1.03*msq2(2,2))))) + (AbInput - MuInput*TanBeta)*Cube(AtauInput - MuInput*
      TanBeta)*((-4*Log(1.01*mse2(2,2))*(1.01*mse2(2,2) + 0.99*msl2(2,2)))/Cube(-
      1.01*mse2(2,2) + 0.99*msl2(2,2)) + Log(0.99*msl2(2,2))*((4*(1.01*mse2(2,2) +
      0.99*msl2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (4*Log(Sqr(SCALE))
      *(1.01*mse2(2,2) + 0.99*msl2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) +
      (4.12*Log(1.03*msq2(2,2))*(1.01*mse2(2,2) + 0.99*msl2(2,2))*msq2(2,2))/(Cube
      (-1.01*mse2(2,2) + 0.99*msl2(2,2))*(0.97*msd2(2,2) - 1.03*msq2(2,2)))) + Log
      (Sqr(SCALE))*((-4*Log(1.01*mse2(2,2))*(1.01*mse2(2,2) + 0.99*msl2(2,2)))/
      Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - 8/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2
      ,2))) + Log(0.97*msd2(2,2))*((3.88*Log(1.01*mse2(2,2))*msd2(2,2)*(1.01*mse2(
      2,2) + 0.99*msl2(2,2)))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(0.97*msd2(2
      ,2) - 1.03*msq2(2,2))) - (3.88*Log(0.99*msl2(2,2))*msd2(2,2)*(1.01*mse2(2,2)
      + 0.99*msl2(2,2)))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(0.97*msd2(2,2) -
      1.03*msq2(2,2))) + (7.76*msd2(2,2))/((0.97*msd2(2,2) - 1.03*msq2(2,2))*Sqr(-
      1.01*mse2(2,2) + 0.99*msl2(2,2)))) + Log(1.03*msq2(2,2))*((-4.12*Log(1.01*
      mse2(2,2))*(1.01*mse2(2,2) + 0.99*msl2(2,2))*msq2(2,2))/(Cube(-1.01*mse2(2,2
      ) + 0.99*msl2(2,2))*(0.97*msd2(2,2) - 1.03*msq2(2,2))) + (8.24*msq2(2,2))/((
      -0.97*msd2(2,2) + 1.03*msq2(2,2))*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)))) -
      8/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))))*Sqr(Yd(2,2)))/(Sqrt(1/(1 + Sqr(
      TanBeta)))*Sqrt(Sqr(TanBeta)/(1 + Sqr(TanBeta)))*Sqr(1 + (
      0.006332573977646111*(1 + Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) +
      (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 +
      Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2))
      )/Sqr(TanBeta) + 0.5*((-0.0031662869888230555*Sqr(AbInput - MuInput*TanBeta)
      *Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*
      msq2(2,2))) - (0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2
      ,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2))))
      + 0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput
      )/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(
      1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2
      ,2))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(
      SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input)
      - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2)
      )/M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      *(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*
      TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/
      Sqr(TanBeta)))))/Quad(3.141592653589793) + (0.00390625*Power6(Ye(2,2))*(5 -
      4*Log(1.01*mse2(2,2)) + Log(0.99*msl2(2,2))*(-6 - 6*Log(Sqr(SCALE))) + (10 -
      4*Log(1.01*mse2(2,2)))*Log(Sqr(SCALE)) + 5*Sqr(Log(Sqr(SCALE))) + 2*Sqr(Log(
      1.01*mse2(2,2))) + 3*Sqr(Log(0.99*msl2(2,2))) + Power6(AtauInput - MuInput*
      TanBeta)*((4*PolyLog(2,(1.0101010101010102*(-1.01*mse2(2,2) + 0.99*msl2(2,2)
      ))/msl2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (2*PolyLog(2,1 - (
      0.9801980198019802*msl2(2,2))/mse2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2
      ,2)) + ((-5.05*mse2(2,2) - 2.9699999999999998*msl2(2,2))*Sqr(Log(1.01*mse2(2
      ,2))))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + ((-5.05*mse2(2,2) - 8.91*
      msl2(2,2))*Sqr(Log(0.99*msl2(2,2))))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))
      + Log(Sqr(SCALE))*((-5.9399999999999995*Log(1.01*mse2(2,2))*msl2(2,2))/Quad(
      -1.01*mse2(2,2) + 0.99*msl2(2,2)) + (1.000100010001*(-4.9995*mse2(2,2)*msl2(
      2,2) + 1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2))))/(Cube(-1.01*mse2(2,2)
      + 0.99*msl2(2,2))*mse2(2,2)*msl2(2,2))) + (1.0101010101010102*Log(1.01*mse2(
      2,2))*(-2.9997*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) - 9.801*Sqr(msl2(
      2,2))))/(msl2(2,2)*Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + (1.000100010001
      *(-10.9989*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2
      ))))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*mse2(2,2)*msl2(2,2)) + Log(0.99
      *msl2(2,2))*((5.9399999999999995*Log(Sqr(SCALE))*msl2(2,2))/Quad(-1.01*mse2(
      2,2) + 0.99*msl2(2,2)) + (2*Log(1.01*mse2(2,2))*(5.05*mse2(2,2) +
      5.9399999999999995*msl2(2,2)))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (
      0.9900990099009901*(12.9987*mse2(2,2)*msl2(2,2) - 3.0603*Sqr(mse2(2,2)) +
      1.9602*Sqr(msl2(2,2))))/(mse2(2,2)*Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)))))
      + Quad(AtauInput - MuInput*TanBeta)*(((5.05*mse2(2,2) + 2.9699999999999998*
      msl2(2,2))*Sqr(Log(1.01*mse2(2,2))))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))
      + ((-9.09*mse2(2,2) - 6.93*msl2(2,2))*Sqr(Log(0.99*msl2(2,2))))/Cube(-1.01*
      mse2(2,2) + 0.99*msl2(2,2)) - (6*PolyLog(2,1 - (0.9801980198019802*msl2(2,2)
      )/mse2(2,2)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + Log(0.99*msl2(2,2))*((
      4*Log(1.01*mse2(2,2))*(1.01*mse2(2,2) + 0.99*msl2(2,2)))/Cube(-1.01*mse2(2,2
      ) + 0.99*msl2(2,2)) + (2*Log(Sqr(SCALE))*(7.07*mse2(2,2) + 4.95*msl2(2,2)))/
      Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (3.9603960396039604*(-9.999*mse2(2,
      2)*msl2(2,2) - 3.0603*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))))/(Cube(-1.01*
      mse2(2,2) + 0.99*msl2(2,2))*mse2(2,2))) + (1.0101010101010102*Log(1.01*mse2(
      2,2))*(-43.9956*mse2(2,2)*msl2(2,2) + 2.0402*Sqr(mse2(2,2)) -
      5.880599999999999*Sqr(msl2(2,2))))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*
      msl2(2,2)) + (2.000200020002*(-26.9973*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2
      (2,2)) + 1.9602*Sqr(msl2(2,2))))/(mse2(2,2)*msl2(2,2)*Sqr(-1.01*mse2(2,2) +
      0.99*msl2(2,2))) + Log(Sqr(SCALE))*((-2*Log(1.01*mse2(2,2))*(7.07*mse2(2,2)
      + 4.95*msl2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (2.000200020002*
      (-14.9985*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) + 1.9602*Sqr(msl2(2,2)
      )))/(mse2(2,2)*msl2(2,2)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))))) + Sqr(
      AtauInput - MuInput*TanBeta)*((-6*PolyLog(2,1 - (0.9801980198019802*msl2(2,2
      ))/mse2(2,2)))/(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (9*Sqr(Log(1.01*mse2(2,2
      ))))/(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + ((-9.09*mse2(2,2) + 6.93*msl2(2,2)
      )*Sqr(Log(0.99*msl2(2,2))))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + Log(Sqr(
      SCALE))*((2*Log(1.01*mse2(2,2))*(-9.09*mse2(2,2) + 7.92*msl2(2,2)))/Sqr(-
      1.01*mse2(2,2) + 0.99*msl2(2,2)) + (1.000100010001*(-0.9999*mse2(2,2)*msl2(2
      ,2) + 1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2))))/(mse2(2,2)*(-1.01*mse2
      (2,2) + 0.99*msl2(2,2))*msl2(2,2))) + (1.000100010001*(-2.9997*mse2(2,2)*
      msl2(2,2) + 1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2))))/(mse2(2,2)*(-
      1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)) + (1.0101010101010102*Log(1.01*
      mse2(2,2))*(-20.997899999999998*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2))
      + 15.6816*Sqr(msl2(2,2))))/(msl2(2,2)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)))
      + Log(0.99*msl2(2,2))*((1.98*Log(1.01*mse2(2,2))*msl2(2,2))/Sqr(-1.01*mse2(2
      ,2) + 0.99*msl2(2,2)) - (2*Log(Sqr(SCALE))*(-9.09*mse2(2,2) + 7.92*msl2(2,2)
      ))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (0.9900990099009901*(-
      16.998299999999997*mse2(2,2)*msl2(2,2) + 19.3819*Sqr(mse2(2,2)) + 1.9602*Sqr
      (msl2(2,2))))/(mse2(2,2)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))))) + (1 + Sqr
      (TanBeta))*(Log(Sqr(SCALE))*(Log(1.01*mse2(2,2)) - (1.000100010001*(1.01*
      mse2(2,2) + 1.98*msl2(2,2))*(1.01*mse2(2,2) + 0.99*msl2(2,2) - 1.9*Sqr(
      MuInput)))/(mse2(2,2)*msl2(2,2))) - 2*Sqr(Log(Sqr(SCALE))) + (
      0.5000500050005*(3.4294999999999995*(1.01*mse2(2,2) + 1.98*msl2(2,2))*Power6
      (MuInput) - 0.9999*mse2(2,2)*msl2(2,2)*(2.9997*mse2(2,2)*msl2(2,2) + 2.0402*
      Sqr(mse2(2,2)) + 3.9204*Sqr(msl2(2,2))) - 2.7075*Quad(MuInput)*(6.9993*mse2(
      2,2)*msl2(2,2) + 2.0402*Sqr(mse2(2,2)) + 3.9204*Sqr(msl2(2,2))) + 0.95*Sqr(
      MuInput)*(2.060602*Cube(mse2(2,2)) + 3.8811959999999996*Cube(msl2(2,2)) +
      13.128687000000001*msl2(2,2)*Sqr(mse2(2,2)) + 16.828317000000002*mse2(2,2)*
      Sqr(msl2(2,2)))))/(mse2(2,2)*msl2(2,2)*(-1.01*mse2(2,2) + 0.95*Sqr(MuInput))
      *(-0.99*msl2(2,2) + 0.95*Sqr(MuInput))) + (3.61*PolyLog(2,(
      1.0526315789473684*(-0.99*msl2(2,2) + 0.95*Sqr(MuInput)))/Sqr(MuInput))*Quad
      (MuInput))/Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput)) - (1.98*msl2(2,2)*(0.99*
      msl2(2,2) - 1.9*Sqr(MuInput))*Sqr(Log(0.99*msl2(2,2))))/Sqr(0.99*msl2(2,2) -
      0.95*Sqr(MuInput)) + (2*PolyLog(2,1 - (1.0421052631578946*msl2(2,2))/Sqr(
      MuInput))*(0.9025*Quad(MuInput) + 1.881*msl2(2,2)*Sqr(MuInput) - 0.9801*Sqr(
      msl2(2,2))))/Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput)) + Sqr(AtauInput -
      MuInput*TanBeta)*((4*PolyLog(2,1 - (1.063157894736842*mse2(2,2))/Sqr(MuInput
      ))*(0.99*msl2(2,2) - 0.95*Sqr(MuInput)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2
      )) + (4*PolyLog(2,(1.0526315789473684*(-0.99*msl2(2,2) + 0.95*Sqr(MuInput)))
      /Sqr(MuInput))*(-0.99*msl2(2,2) + 0.95*Sqr(MuInput)))/Sqr(-1.01*mse2(2,2) +
      0.99*msl2(2,2)) + (5.9994*mse2(2,2)*msl2(2,2) + 3.8379999999999996*mse2(2,2)
      *Sqr(MuInput) - 7.524*msl2(2,2)*Sqr(MuInput) - 2.0402*Sqr(mse2(2,2)) +
      3.9204*Sqr(msl2(2,2)))/(-1.009899*msl2(2,2)*Sqr(mse2(2,2)) + 0.989901*mse2(2
      ,2)*Sqr(msl2(2,2))) + (1.0101010101010102*Log(1.01*mse2(2,2))*(
      0.9702989999999999*Cube(msl2(2,2))*(-11.11*mse2(2,2) + 4.75*Sqr(MuInput)) +
      1.9381899999999999*(-1.01*mse2(2,2) + 0.95*Sqr(MuInput))*Sqr(MuInput)*Sqr(
      mse2(2,2)) + 0.99*msl2(2,2)*(2.060602*Cube(mse2(2,2)) + 3.4294999999999995*
      Power6(MuInput) - 11.849825*mse2(2,2)*Quad(MuInput) + 0.9690949999999999*Sqr
      (MuInput)*Sqr(mse2(2,2))) + 0.9801*(-4.5125*Quad(MuInput) +
      15.351999999999999*mse2(2,2)*Sqr(MuInput) + 1.0201*Sqr(mse2(2,2)))*Sqr(msl2(
      2,2))))/(msl2(2,2)*(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*(-1.01*mse2(2,2) +
      0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + Log(Sqr(SCALE))*
      ((Log(1.01*mse2(2,2))*(7.92*msl2(2,2) - 2*(1.01*mse2(2,2) + 1.9*Sqr(MuInput)
      )))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (3.9996*mse2(2,2)*msl2(2,2) +
      3.8379999999999996*mse2(2,2)*Sqr(MuInput) - 7.524*msl2(2,2)*Sqr(MuInput) -
      2.0402*Sqr(mse2(2,2)) + 3.9204*Sqr(msl2(2,2)))/(-1.009899*msl2(2,2)*Sqr(mse2
      (2,2)) + 0.989901*mse2(2,2)*Sqr(msl2(2,2)))) + (1.98*msl2(2,2)*Sqr(Log(0.99*
      msl2(2,2)))*(1.9998*mse2(2,2)*msl2(2,2) + 4.5125*Quad(MuInput) -
      3.8379999999999996*mse2(2,2)*Sqr(MuInput) - 5.642999999999999*msl2(2,2)*Sqr(
      MuInput) + 2.9402999999999997*Sqr(msl2(2,2))))/(Sqr(-1.01*mse2(2,2) + 0.99*
      msl2(2,2))*Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput))) + Log(0.95*Sqr(MuInput))
      *((3.8003800380037998*Sqr(MuInput)*(0.9405*msl2(2,2)*(1.01*mse2(2,2) - 1.9*
      Sqr(MuInput))*Sqr(MuInput) + 0.9594999999999999*mse2(2,2)*(-1.01*mse2(2,2) +
      0.95*Sqr(MuInput))*Sqr(MuInput) + 0.9801*(-1.01*mse2(2,2) + 1.9*Sqr(MuInput)
      )*Sqr(msl2(2,2))))/(mse2(2,2)*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*(
      0.99*msl2(2,2) - 0.95*Sqr(MuInput))*(-1.01*mse2(2,2) + 0.95*Sqr(MuInput))) +
      (2*Log(1.01*mse2(2,2))*(1.547561875*Power10(MuInput) - 0.81450625*(7.07*mse2
      (2,2) + 2.9699999999999998*msl2(2,2))*Power8(MuInput) - 1.9796040197999998*
      Cube(msl2(2,2))*Sqr(mse2(2,2)) + 1.8808118999999999*mse2(2,2)*(
      3.0300000000000002*mse2(2,2) + 1.98*msl2(2,2))*Sqr(MuInput)*Sqr(msl2(2,2)) +
      1.7147499999999998*Power6(MuInput)*(4.9995*mse2(2,2)*msl2(2,2) + 3.0603*Sqr(
      mse2(2,2)) + 1.9602*Sqr(msl2(2,2))) - 0.9025*Quad(MuInput)*(2.060602*Cube(
      mse2(2,2)) + 0.9702989999999999*Cube(msl2(2,2)) + 4.039596*msl2(2,2)*Sqr(
      mse2(2,2)) + 12.868713*mse2(2,2)*Sqr(msl2(2,2)))))/(Sqr(-1.01*mse2(2,2) +
      0.99*msl2(2,2))*Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2)
      + 0.95*Sqr(MuInput)))) + Log(0.99*msl2(2,2))*((Log(Sqr(SCALE))*(2.02*mse2(2,
      2) - 7.92*msl2(2,2) + 3.8*Sqr(MuInput)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2
      )) + (0.9900990099009901*(-3.8811959999999996*Cube(msl2(2,2))*(-1.01*mse2(2,
      2) + 0.95*Sqr(MuInput)) + 0.9999*mse2(2,2)*msl2(2,2)*(4.5125*Quad(MuInput) -
      1.9189999999999998*mse2(2,2)*Sqr(MuInput) - 7.1407*Sqr(mse2(2,2))) -
      0.9594999999999999*mse2(2,2)*Sqr(MuInput)*(3.61*Quad(MuInput) - 6.7165*mse2(
      2,2)*Sqr(MuInput) + 1.0201*Sqr(mse2(2,2))) + 0.9801*(3.61*Quad(MuInput) -
      12.4735*mse2(2,2)*Sqr(MuInput) + 11.2211*Sqr(mse2(2,2)))*Sqr(msl2(2,2))))/(
      mse2(2,2)*(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*(-1.01*mse2(2,2) + 0.95*Sqr(
      MuInput))*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + Log(1.01*mse2(2,2))*((-
      4.04*mse2(2,2))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (4*(1.01*mse2(2,2) +
      0.99*msl2(2,2)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (5.9399999999999995
      *msl2(2,2))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (3.96*msl2(2,2)*(0.99*
      msl2(2,2) - 1.9*Sqr(MuInput)))/((-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.99*
      msl2(2,2) - 0.95*Sqr(MuInput))) + (2.02*mse2(2,2)*(1.01*mse2(2,2) - 1.9*Sqr(
      MuInput)))/((1.01*mse2(2,2) - 0.99*msl2(2,2))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr
      (MuInput)))) + (2*Log(0.95*Sqr(MuInput))*(-1.547561875*Power10(MuInput) +
      0.81450625*(7.07*mse2(2,2) + 2.9699999999999998*msl2(2,2))*Power8(MuInput) +
      1.9796040197999998*Cube(msl2(2,2))*Sqr(mse2(2,2)) - 1.8808118999999999*mse2(
      2,2)*(3.0300000000000002*mse2(2,2) + 1.98*msl2(2,2))*Sqr(MuInput)*Sqr(msl2(2
      ,2)) - 1.7147499999999998*Power6(MuInput)*(4.9995*mse2(2,2)*msl2(2,2) +
      3.0603*Sqr(mse2(2,2)) + 1.9602*Sqr(msl2(2,2))) + 0.9025*Quad(MuInput)*(
      2.060602*Cube(mse2(2,2)) + 0.9702989999999999*Cube(msl2(2,2)) + 4.039596*
      msl2(2,2)*Sqr(mse2(2,2)) + 12.868713*mse2(2,2)*Sqr(msl2(2,2)))))/(Sqr(-1.01*
      mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*Sqr(-
      1.01*mse2(2,2) + 0.95*Sqr(MuInput)))) + (2.02*mse2(2,2)*Sqr(Log(1.01*mse2(2,
      2)))*(3.61*Quad(MuInput) + 0.99*msl2(2,2)*(1.01*mse2(2,2) - 1.9*Sqr(MuInput)
      ) - 5.757*mse2(2,2)*Sqr(MuInput) + 3.0603*Sqr(mse2(2,2))))/(Sqr(-1.01*mse2(2
      ,2) + 0.99*msl2(2,2))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)))) + Quad(
      AtauInput - MuInput*TanBeta)*((2*PolyLog(2,1 - (1.063157894736842*mse2(2,2))
      /Sqr(MuInput))*(-1.9998*mse2(2,2)*msl2(2,2) - 2.7075*Quad(MuInput) + 5.643*
      msl2(2,2)*Sqr(MuInput) + 1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2))))/
      Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (PolyLog(2,(1.0526315789473684*(-
      0.99*msl2(2,2) + 0.95*Sqr(MuInput)))/Sqr(MuInput))*(5.415*Quad(MuInput) +
      3.96*msl2(2,2)*(1.01*mse2(2,2) - 2.8499999999999996*Sqr(MuInput)) - 2.0402*
      Sqr(mse2(2,2)) + 3.9204*Sqr(msl2(2,2))))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,
      2)) + (1.000100010001*(-1.92119202*Quad(msl2(2,2))*(-1.01*mse2(2,2) + 0.95*
      Sqr(MuInput)) + 0.9690949999999999*Sqr(MuInput)*Sqr(mse2(2,2))*(1.805*Quad(
      MuInput) - 2.8785*mse2(2,2)*Sqr(MuInput) + 1.0201*Sqr(mse2(2,2))) +
      0.9702989999999999*Cube(msl2(2,2))*(5.415*Quad(MuInput) - 34.541999999999994
      *mse2(2,2)*Sqr(MuInput) + 32.6432*Sqr(mse2(2,2))) - 0.9999*mse2(2,2)*msl2(2,
      2)*(1.030301*Cube(mse2(2,2)) + 13.717999999999998*Power6(MuInput) - 4.557625
      *mse2(2,2)*Quad(MuInput) - 13.567329999999998*Sqr(MuInput)*Sqr(mse2(2,2))) -
      0.9801*(15.454514999999999*Cube(mse2(2,2)) + 3.4294999999999995*Power6(
      MuInput) - 41.93015*mse2(2,2)*Quad(MuInput) + 30.041945*Sqr(MuInput)*Sqr(
      mse2(2,2)))*Sqr(msl2(2,2))))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*mse2(2,
      2)*msl2(2,2)*(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*(-1.01*mse2(2,2) + 0.95*
      Sqr(MuInput))) + Log(Sqr(SCALE))*((Log(1.01*mse2(2,2))*(-5.9994*mse2(2,2)*
      msl2(2,2) + 11.286*msl2(2,2)*Sqr(MuInput) + 7.1407*Sqr(mse2(2,2)) -
      12.741299999999999*Sqr(msl2(2,2))))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) +
      (1.000100010001*(-1.9405979999999998*Cube(msl2(2,2)) + 9.999*mse2(2,2)*msl2(
      2,2)*(1.01*mse2(2,2) + 0.95*Sqr(MuInput)) + 1.0201*(1.01*mse2(2,2) - 1.9*Sqr
      (MuInput))*Sqr(mse2(2,2)) + 0.9801*(-21.21*mse2(2,2) + 3.8*Sqr(MuInput))*Sqr
      (msl2(2,2))))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*mse2(2,2)*msl2(2,2)))
      - (1.98*msl2(2,2)*Sqr(Log(0.99*msl2(2,2)))*(4.851495*Cube(msl2(2,2)) +
      0.9594999999999999*mse2(2,2)*Sqr(MuInput)*(-2.02*mse2(2,2) + 4.75*Sqr(
      MuInput)) + 0.99*msl2(2,2)*(5.415*Quad(MuInput) - 9.595*mse2(2,2)*Sqr(
      MuInput) + 1.0201*Sqr(mse2(2,2))) + 4.9005*(1.01*mse2(2,2) - 1.9*Sqr(MuInput
      ))*Sqr(msl2(2,2))))/(Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.99*msl2(2,
      2) - 0.95*Sqr(MuInput))) + Log(0.99*msl2(2,2))*((Log(Sqr(SCALE))*(
      5.9399999999999995*msl2(2,2)*(1.01*mse2(2,2) - 1.9*Sqr(MuInput)) - 7.1407*
      Sqr(mse2(2,2)) + 12.741299999999999*Sqr(msl2(2,2))))/Quad(-1.01*mse2(2,2) +
      0.99*msl2(2,2)) + (0.49504950495049505*(0.9298466524999999*Cube(mse2(2,2))*
      Quad(MuInput)*(15.15*mse2(2,2) - 12.35*Sqr(MuInput)) + 3.8039601995999996*
      Power5(msl2(2,2))*(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)) + 0.96059601*Quad(
      msl2(2,2))*(-7.22*Quad(MuInput) + 73.88149999999999*mse2(2,2)*Sqr(MuInput) -
      72.4271*Sqr(mse2(2,2))) + 1.89981*mse2(2,2)*msl2(2,2)*Sqr(MuInput)*(-
      25.757524999999998*Cube(mse2(2,2)) - 15.432749999999999*Power6(MuInput) +
      3.6461*mse2(2,2)*Quad(MuInput) + 35.856514999999995*Sqr(MuInput)*Sqr(mse2(2,
      2))) + 1.9405979999999998*Cube(msl2(2,2))*(6.181806*Cube(mse2(2,2)) +
      1.7147499999999998*Power6(MuInput) - 79.302675*mse2(2,2)*Quad(MuInput) +
      78.49669499999999*Sqr(MuInput)*Sqr(mse2(2,2))) + 0.989901*mse2(2,2)*(
      27.818126999999997*Cube(mse2(2,2)) + 124.31937499999998*Power6(MuInput) -
      112.11757499999999*mse2(2,2)*Quad(MuInput) - 47.485654999999994*Sqr(MuInput)
      *Sqr(mse2(2,2)))*Sqr(msl2(2,2))))/(mse2(2,2)*Quad(-1.01*mse2(2,2) + 0.99*
      msl2(2,2))*(-1.01*mse2(2,2) + 0.95*Sqr(MuInput))*Sqr(0.99*msl2(2,2) - 0.95*
      Sqr(MuInput))) + Log(1.01*mse2(2,2))*((6.0600000000000005*mse2(2,2))/Cube(
      1.01*mse2(2,2) - 0.99*msl2(2,2)) + (9.9*msl2(2,2))/Cube(-1.01*mse2(2,2) +
      0.99*msl2(2,2)) + (21.9978*mse2(2,2)*msl2(2,2))/Quad(-1.01*mse2(2,2) + 0.99*
      msl2(2,2)) - 2/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (4*Sqr(1.01*mse2(2,2)
      + 0.99*msl2(2,2)))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (1.98*(1.01*mse2
      (2,2) + 0.99*msl2(2,2))*msl2(2,2)*(0.99*msl2(2,2) - 1.9*Sqr(MuInput)))/(Cube
      (-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput)))
      + (1.01*mse2(2,2)*(1.01*mse2(2,2) + 0.99*msl2(2,2))*(-1.01*mse2(2,2) + 1.9*
      Sqr(MuInput)))/(Cube(1.01*mse2(2,2) - 0.99*msl2(2,2))*Sqr(-1.01*mse2(2,2) +
      0.95*Sqr(MuInput)))) + (Log(0.95*Sqr(MuInput))*(9.192517537499999*msl2(2,2)*
      Power10(MuInput) + 1.99960002*Sqr(mse2(2,2))*(-1.9998*mse2(2,2)*msl2(2,2) +
      1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2)))*Sqr(msl2(2,2)) +
      5.092807499999999*msl2(2,2)*Power6(MuInput)*(9.999*mse2(2,2)*msl2(2,2) +
      3.0603*Sqr(mse2(2,2)) + 2.9402999999999997*Sqr(msl2(2,2))) - 0.81450625*
      Power8(MuInput)*(27.9972*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) +
      24.502499999999998*Sqr(msl2(2,2))) + 3.79962*mse2(2,2)*msl2(2,2)*Sqr(MuInput
      )*(-1.030301*Cube(mse2(2,2)) + 1.9405979999999998*Cube(msl2(2,2)) + 1.009899
      *msl2(2,2)*Sqr(mse2(2,2)) + 6.9293070000000005*mse2(2,2)*Sqr(msl2(2,2))) -
      0.8934749999999999*msl2(2,2)*Quad(MuInput)*(-4.121204*Cube(mse2(2,2)) +
      2.910897*Cube(msl2(2,2)) + 41.40585900000001*msl2(2,2)*Sqr(mse2(2,2)) +
      43.555644*mse2(2,2)*Sqr(msl2(2,2)))))/(Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)
      )*Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(
      MuInput)))) + Log(0.95*Sqr(MuInput))*((-1.9001900190018999*Sqr(MuInput)*(-
      0.9206402499999999*Quad(MuInput)*Sqr(0.95*Sqr(MuInput) - 1.01*mse2(2,2))*Sqr
      (mse2(2,2)) + 0.9024097499999999*mse2(2,2)*msl2(2,2)*Quad(MuInput)*(1.805*
      Quad(MuInput) - 4.7975*mse2(2,2)*Sqr(MuInput) + 2.0402*Sqr(mse2(2,2))) +
      0.96059601*Quad(msl2(2,2))*(1.805*Quad(MuInput) - 2.8785*mse2(2,2)*Sqr(
      MuInput) + 2.0402*Sqr(mse2(2,2))) + 0.9702989999999999*Cube(msl2(2,2))*(
      2.060602*Cube(mse2(2,2)) - 3.4294999999999995*Power6(MuInput) + 7.2922*mse2(
      2,2)*Quad(MuInput) - 8.721855*Sqr(MuInput)*Sqr(mse2(2,2))) + 0.9801*(-
      4.32974375*mse2(2,2)*Power6(MuInput) + 1.6290125*Power8(MuInput) -
      1.04060401*Quad(mse2(2,2)) + 6.44448175*Quad(MuInput)*Sqr(mse2(2,2)))*Sqr(
      msl2(2,2))))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*mse2(2,2)*msl2(2,2)*Sqr
      (0.99*msl2(2,2) - 0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)
      )) + (Log(1.01*mse2(2,2))*(-9.1925175375*msl2(2,2)*Power10(MuInput) +
      1.99960002*Sqr(mse2(2,2))*Sqr(msl2(2,2))*(1.9998*mse2(2,2)*msl2(2,2) -
      1.0201*Sqr(mse2(2,2)) + 1.9602*Sqr(msl2(2,2))) - 5.092807499999999*msl2(2,2)
      *Power6(MuInput)*(9.999*mse2(2,2)*msl2(2,2) + 3.0603*Sqr(mse2(2,2)) +
      2.9402999999999997*Sqr(msl2(2,2))) + 0.81450625*Power8(MuInput)*(27.9972*
      mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) + 24.502499999999998*Sqr(msl2(2,
      2))) - 3.79962*mse2(2,2)*msl2(2,2)*Sqr(MuInput)*(-1.030301*Cube(mse2(2,2)) +
      1.9405979999999998*Cube(msl2(2,2)) + 1.009899*msl2(2,2)*Sqr(mse2(2,2)) +
      6.9293070000000005*mse2(2,2)*Sqr(msl2(2,2))) + 0.8934749999999999*msl2(2,2)*
      Quad(MuInput)*(-4.121204*Cube(mse2(2,2)) + 2.910897*Cube(msl2(2,2)) +
      41.40585900000001*msl2(2,2)*Sqr(mse2(2,2)) + 43.555644*mse2(2,2)*Sqr(msl2(2,
      2)))))/(Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.99*msl2(2,2) - 0.95*Sqr
      (MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)))) - (1.01*mse2(2,2)*(
      1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(Log(1.01*mse2(2,2)))*(7.22*Quad(MuInput
      ) + 0.99*msl2(2,2)*(1.01*mse2(2,2) - 1.9*Sqr(MuInput)) - 13.433*mse2(2,2)*
      Sqr(MuInput) + 7.1407*Sqr(mse2(2,2))))/(Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2
      ))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput))) + (0.5050505050505051*Log(1.01*
      mse2(2,2))*(-1.9575718999999998*Cube(mse2(2,2))*Sqr(MuInput)*Sqr(0.95*Sqr(
      MuInput) - 1.01*mse2(2,2)) + 0.96059601*Quad(msl2(2,2))*(-24.3675*Quad(
      MuInput) + 61.407999999999994*mse2(2,2)*Sqr(MuInput) - 33.6633*Sqr(mse2(2,2)
      )) + 0.9999*mse2(2,2)*msl2(2,2)*(-16.453026249999997*mse2(2,2)*Power6(
      MuInput) - 9.774075*Power8(MuInput) + 2.08120802*Quad(mse2(2,2)) -
      44.04536775*Cube(mse2(2,2))*Sqr(MuInput) + 64.4448175*Quad(MuInput)*Sqr(mse2
      (2,2))) + 0.9702989999999999*Cube(msl2(2,2))*(-51.515049999999995*Cube(mse2(
      2,2)) + 40.29662499999999*Power6(MuInput) - 147.66705*mse2(2,2)*Quad(MuInput
      ) + 148.271535*Sqr(MuInput)*Sqr(mse2(2,2))) + 0.9801*(102.1819525*mse2(2,2)*
      Power6(MuInput) - 19.54815*Power8(MuInput) + 46.82718045*Quad(mse2(2,2)) -
      25.448434699999996*Cube(mse2(2,2))*Sqr(MuInput) - 92.98466525*Quad(MuInput)*
      Sqr(mse2(2,2)))*Sqr(msl2(2,2))))/(msl2(2,2)*Quad(-1.01*mse2(2,2) + 0.99*msl2
      (2,2))*(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(
      MuInput)))) + Log(0.99*msl2(2,2))*(3*Log(Sqr(SCALE)) + (0.49504950495049505*
      (3.8811959999999996*Cube(msl2(2,2))*(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)) +
      0.911525*mse2(2,2)*Quad(MuInput)*(-13.13*mse2(2,2) + 14.25*Sqr(MuInput)) +
      1.881*msl2(2,2)*Sqr(MuInput)*(1.805*Quad(MuInput) - 6.7165*mse2(2,2)*Sqr(
      MuInput) + 3.0603*Sqr(mse2(2,2))) + 0.9801*(-7.22*Quad(MuInput) +
      14.392499999999998*mse2(2,2)*Sqr(MuInput) - 5.1005*Sqr(mse2(2,2)))*Sqr(msl2(
      2,2))))/(mse2(2,2)*(-1.01*mse2(2,2) + 0.95*Sqr(MuInput))*Sqr(0.99*msl2(2,2)
      - 0.95*Sqr(MuInput))) + (Log(1.01*mse2(2,2))*(-1.7318974999999999*mse2(2,2)*
      Power6(MuInput) + 1.6290125*Power8(MuInput) + 0.911525*mse2(2,2)*(1.01*mse2(
      2,2) - 3.96*msl2(2,2))*Quad(MuInput) + 1.89981*mse2(2,2)*(1.01*mse2(2,2) +
      0.99*msl2(2,2))*msl2(2,2)*Sqr(MuInput) - 0.99980001*Sqr(mse2(2,2))*Sqr(msl2(
      2,2))))/(Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*
      Sqr(MuInput))) + (Log(0.95*Sqr(MuInput))*(-1.7147499999999998*(-4.04*mse2(2,
      2) + 0.99*msl2(2,2))*Power6(MuInput) - 4.07253125*Power8(MuInput) - 3.79962*
      mse2(2,2)*(1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*Sqr(MuInput) + 0.9025*
      Quad(MuInput)*(7.9992*mse2(2,2)*msl2(2,2) - 4.0804*Sqr(mse2(2,2)) + 0.9801*
      Sqr(msl2(2,2))) + 1.99960002*Sqr(mse2(2,2))*Sqr(msl2(2,2))))/(Sqr(0.99*msl2(
      2,2) - 0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)))) + Log(
      0.95*Sqr(MuInput))*((Log(1.01*mse2(2,2))*(1.6976024999999997*msl2(2,2)*
      Power6(MuInput) - 2.44351875*Power8(MuInput) - 0.8934749999999999*(-8.08*
      mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*Quad(MuInput) - 3.79962*mse2(2,2)*(
      1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*Sqr(MuInput) + 1.99960002*Sqr(
      mse2(2,2))*Sqr(msl2(2,2))))/(Sqr(0.99*msl2(2,2) - 0.95*Sqr(MuInput))*Sqr(-
      1.01*mse2(2,2) + 0.95*Sqr(MuInput))) + (0.9500950095009499*Sqr(MuInput)*(-
      1.6290125*(1.01*mse2(2,2) + 1.98*msl2(2,2))*Power8(MuInput) - 1.99960002*(
      1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(mse2(2,2))*Sqr(msl2(2,2)) + 0.949905*
      mse2(2,2)*msl2(2,2)*Sqr(MuInput)*(7.9992*mse2(2,2)*msl2(2,2) - 2.0402*Sqr(
      mse2(2,2)) + 2.9402999999999997*Sqr(msl2(2,2))) + 0.8573749999999999*Power6(
      MuInput)*(0.9999*mse2(2,2)*msl2(2,2) + 4.0804*Sqr(mse2(2,2)) + 7.8408*Sqr(
      msl2(2,2))) - 1.805*Quad(MuInput)*(1.030301*Cube(mse2(2,2)) +
      1.9405979999999998*Cube(msl2(2,2)) - 1.009899*msl2(2,2)*Sqr(mse2(2,2)) +
      3.959604*mse2(2,2)*Sqr(msl2(2,2)))))/(mse2(2,2)*msl2(2,2)*Sqr(0.99*msl2(2,2)
      - 0.95*Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)))) + (1.01*mse2
      (2,2)*(-1.01*mse2(2,2) + 1.9*Sqr(MuInput))*Sqr(Log(1.01*mse2(2,2))))/Sqr(-
      1.01*mse2(2,2) + 0.95*Sqr(MuInput)) + (2*PolyLog(2,1 - (1.063157894736842*
      mse2(2,2))/Sqr(MuInput))*(0.9025*Quad(MuInput) + 1.9189999999999998*mse2(2,2
      )*Sqr(MuInput) - 1.0201*Sqr(mse2(2,2))))/Sqr(-1.01*mse2(2,2) + 0.95*Sqr(
      MuInput)) + (0.5050505050505051*Log(1.01*mse2(2,2))*(-1.9189999999999998*
      mse2(2,2)*Sqr(MuInput)*Sqr(0.95*Sqr(MuInput) - 1.01*mse2(2,2)) + 0.99*msl2(2
      ,2)*(2.060602*Cube(mse2(2,2)) - 11.145874999999998*Power6(MuInput) +
      12.76135*mse2(2,2)*Quad(MuInput) - 8.721855*Sqr(MuInput)*Sqr(mse2(2,2))) +
      0.9801*(8.1225*Quad(MuInput) - 3.8379999999999996*mse2(2,2)*Sqr(MuInput) +
      1.0201*Sqr(mse2(2,2)))*Sqr(msl2(2,2))))/(msl2(2,2)*(0.99*msl2(2,2) - 0.95*
      Sqr(MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput))) + (Sqr(Log(0.95*Sqr(
      MuInput)))*(-3.4637949999999997*mse2(2,2)*Power6(MuInput) + 3.258025*Power8(
      MuInput) + 1.82305*mse2(2,2)*(1.01*mse2(2,2) - 3.96*msl2(2,2))*Quad(MuInput)
      + 3.79962*mse2(2,2)*(1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*Sqr(MuInput)
      - 1.99960002*Sqr(mse2(2,2))*Sqr(msl2(2,2))))/(Sqr(0.99*msl2(2,2) - 0.95*Sqr(
      MuInput))*Sqr(-1.01*mse2(2,2) + 0.95*Sqr(MuInput)))) + Sqr(TanBeta)*(1.5 -
      2.5*Log(1.01*mse2(2,2)) + Log(0.99*msl2(2,2))*(-4.5 - 3*Log(Sqr(SCALE))) +
      Sqr(3.141592653589793) - (1.9504950495049505*Sqr(mAInput))/mse2(2,2) - (
      0.994949494949495*Sqr(mAInput))/msl2(2,2) + Log(0.985*Sqr(mAInput))*(7 - 3*
      Log(1.01*mse2(2,2)) - 3*Log(0.99*msl2(2,2)) + 0.985*(1.9801980198019802/mse2
      (2,2) + 1.0101010101010102/msl2(2,2))*Sqr(mAInput)) + Log(Sqr(SCALE))*(-Log(
      1.01*mse2(2,2)) - (0.985098509850985*(1.01*mse2(2,2) + 1.98*msl2(2,2))*Sqr(
      mAInput))/(mse2(2,2)*msl2(2,2))) + 3*Sqr(Log(0.985*Sqr(mAInput))) + 2*Sqr(
      Log(Sqr(SCALE))) + 2*Sqr(Log(1.01*mse2(2,2))) + 3*Sqr(Log(0.99*msl2(2,2))) +
      (0.9802960494069208*(0.970225*Quad(mAInput) - 5.9691*mse2(2,2)*Sqr(mAInput)
      - TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*mse2(2,2)))*TPhi(0.985*Sqr(
      mAInput),1.01*mse2(2,2),1.01*mse2(2,2)))/Sqr(mse2(2,2)) + Sqr(AtauInput +
      MuInput/TanBeta)*(-1.9801980198019802/mse2(2,2) + Log(Sqr(SCALE))*(-
      1.9801980198019802/mse2(2,2) - 1.0101010101010102/msl2(2,2)) -
      1.0101010101010102/msl2(2,2) + (0.4901480247034604*Log(0.99*msl2(2,2))*((-
      Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + 0.955671625*Power6(mAInput) -
      2.910675*(1.01*mse2(2,2) + 0.99*msl2(2,2))*Quad(mAInput) + 2.955*Sqr(mAInput
      )*(-1.0201*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))))*TDelta(0.985*Sqr(mAInput
      ),1.01*mse2(2,2),0.99*msl2(2,2)) + (-1.01*mse2(2,2)*(0.970225*Quad(mAInput)
      - 1.9897*mse2(2,2)*Sqr(mAInput) + 1.0201*Sqr(mse2(2,2)) - 0.9801*Sqr(msl2(2,
      2))) + (2.02*mse2(2,2) + 0.99*msl2(2,2) - 0.985*Sqr(mAInput))*TDelta(0.985*
      Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))*TDelta(0.985*Sqr(mAInput),0.99*
      msl2(2,2),1.01*mse2(2,2))))/(Sqr(mse2(2,2))*TDelta(0.985*Sqr(mAInput),1.01*
      mse2(2,2),0.99*msl2(2,2))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2
      (2,2))) + (0.5000500050005*Log(1.01*mse2(2,2))*((Cube(-1.01*mse2(2,2) + 0.99
      *msl2(2,2)) + 0.955671625*Power6(mAInput) - 0.970225*(5.05*mse2(2,2) + 0.99*
      msl2(2,2))*Quad(mAInput) - 0.985*Sqr(mAInput)*(3.9996*mse2(2,2)*msl2(2,2) -
      5.1005*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))))*TDelta(0.985*Sqr(mAInput),
      1.01*mse2(2,2),0.99*msl2(2,2)) - (1.9998*mse2(2,2)*msl2(2,2)*(-1.01*mse2(2,2
      ) + 0.99*msl2(2,2) + 0.985*Sqr(mAInput)) + (-3.0300000000000002*mse2(2,2) +
      0.99*msl2(2,2) + 0.985*Sqr(mAInput))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2
      ),0.99*msl2(2,2)))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2)))
      )/(mse2(2,2)*msl2(2,2)*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,
      2))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2))) + (
      0.4950990148519802*Log(0.985*Sqr(mAInput))*((-0.955671625*(1.01*mse2(2,2) +
      0.99*msl2(2,2))*Power6(mAInput) + Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) -
      0.985*(5.05*mse2(2,2) + 2.9699999999999998*msl2(2,2))*Sqr(mAInput)*Sqr(-1.01
      *mse2(2,2) + 0.99*msl2(2,2)) + 0.970225*Quad(mAInput)*(3.9996*mse2(2,2)*msl2
      (2,2) + 5.1005*Sqr(mse2(2,2)) + 2.9402999999999997*Sqr(msl2(2,2))))*TDelta(
      0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)) + (-0.9999*mse2(2,2)*msl2(
      2,2)*(-0.970225*Quad(mAInput) + Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + (
      2.9997*mse2(2,2)*msl2(2,2) + 0.985*(1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(
      mAInput) - 1.0201*Sqr(mse2(2,2)) - 0.9801*Sqr(msl2(2,2)))*TDelta(0.985*Sqr(
      mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))*TDelta(0.985*Sqr(mAInput),0.99*msl2
      (2,2),1.01*mse2(2,2))))/(msl2(2,2)*Sqr(mse2(2,2))*TDelta(0.985*Sqr(mAInput),
      1.01*mse2(2,2),0.99*msl2(2,2))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01
      *mse2(2,2))) + (0.5050505050505051*(-Sqr(0.985*Sqr(mAInput) - 1.01*mse2(2,2)
      + 0.99*msl2(2,2)) + TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))
      )*TPhi(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))/(msl2(2,2)*TDelta(
      0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))) + (0.48529507396382227*(-
      3.8226865*(1.01*mse2(2,2) + 0.99*msl2(2,2))*Power6(mAInput) + 0.941336550625
      *Power8(mAInput) + Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - 3.94*(1.01*mse2(
      2,2) + 0.99*msl2(2,2))*Sqr(mAInput)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) +
      1.94045*Quad(mAInput)*(1.9998*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) +
      2.9402999999999997*Sqr(msl2(2,2))) + 2*Sqr(TDelta(0.985*Sqr(mAInput),0.99*
      msl2(2,2),1.01*mse2(2,2))) - 3*(0.970225*Quad(mAInput) - 1.97*(1.01*mse2(2,2
      ) + 0.99*msl2(2,2))*Sqr(mAInput) + Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)))*
      TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2)))*TPhi(0.985*Sqr(
      mAInput),0.99*msl2(2,2),1.01*mse2(2,2)))/(Cube(mse2(2,2))*TDelta(0.985*Sqr(
      mAInput),0.99*msl2(2,2),1.01*mse2(2,2)))) + (1.0203040506070808*(0.970225*
      Quad(mAInput) - 6.8260499999999995*msl2(2,2)*Sqr(mAInput) - TDelta(0.985*Sqr
      (mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))*TPhi(0.985*Sqr(mAInput),0.99*msl2(
      2,2),0.99*msl2(2,2)))/Sqr(msl2(2,2)) + (AtauInput + MuInput/TanBeta)*(
      AtauInput - MuInput*TanBeta)*(Log(0.985*Sqr(mAInput))*((-2*Log(1.01*mse2(2,2
      )))/(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (2*Log(0.99*msl2(2,2)))/(-1.01*mse2
      (2,2) + 0.99*msl2(2,2))) + Log(0.99*msl2(2,2))*(-12/(-1.01*mse2(2,2) + 0.99*
      msl2(2,2)) - (2*Log(1.01*mse2(2,2)))/(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (
      12*Log(Sqr(SCALE)))/(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + (12*Log(1.01*mse2(
      2,2)))/(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (12*Log(1.01*mse2(2,2))*Log(Sqr(
      SCALE)))/(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (4*Sqr(Log(1.01*mse2(2,2))))/(
      -1.01*mse2(2,2) + 0.99*msl2(2,2)) + (6*Sqr(Log(0.99*msl2(2,2))))/(-1.01*mse2
      (2,2) + 0.99*msl2(2,2)) + (1.9605920988138417*(-0.985*(-6.0600000000000005*
      mse2(2,2) + 0.985*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.985*Sqr(mAInput),
      1.01*mse2(2,2),1.01*mse2(2,2)))*TPhi(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*
      mse2(2,2)))/((-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(mse2(2,2))) + (
      2.0202020202020203*(-1.01*mse2(2,2) + 0.99*msl2(2,2) + 0.985*Sqr(mAInput))*
      TPhi(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))/((-1.01*mse2(2,2) +
      0.99*msl2(2,2))*msl2(2,2)) - (2.0406081012141617*(-0.985*(-6.93*msl2(2,2) +
      0.985*Sqr(mAInput))*Sqr(mAInput) + TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),
      0.99*msl2(2,2)))*TPhi(0.985*Sqr(mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))/((-
      1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(msl2(2,2)))) + (AtauInput + MuInput/
      TanBeta)*Cube(AtauInput - MuInput*TanBeta)*((-12*Log(1.01*mse2(2,2))*(
      3.0300000000000002*mse2(2,2) + 0.99*msl2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*
      msl2(2,2)) + Log(0.99*msl2(2,2))*((12*Log(Sqr(SCALE))*(1.01*mse2(2,2) + 0.99
      *msl2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (12*(1.01*mse2(2,2) +
      2.9699999999999998*msl2(2,2)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (2*
      Log(1.01*mse2(2,2))*(1.01*mse2(2,2) + 2.9699999999999998*msl2(2,2) - 1.97*
      Sqr(mAInput)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + Log(0.985*Sqr(
      mAInput))*((2*Log(0.99*msl2(2,2))*(-1.01*mse2(2,2) + 0.99*msl2(2,2) - 5.91*
      Sqr(mAInput)))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (2*Log(1.01*mse2(2,2
      ))*(1.01*mse2(2,2) - 0.99*msl2(2,2) + 5.91*Sqr(mAInput)))/Cube(-1.01*mse2(2,
      2) + 0.99*msl2(2,2))) + (4*(1.01*mse2(2,2) + 0.99*msl2(2,2) - 0.985*Sqr(
      mAInput))*Sqr(Log(1.01*mse2(2,2))))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) +
      (2*(-3.0300000000000002*mse2(2,2) - 4.95*msl2(2,2) + 3.94*Sqr(mAInput))*Sqr(
      Log(0.99*msl2(2,2))))/Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + Log(Sqr(SCALE
      ))*((-12*Log(1.01*mse2(2,2))*(1.01*mse2(2,2) + 0.99*msl2(2,2)))/Cube(-1.01*
      mse2(2,2) + 0.99*msl2(2,2)) - 24/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))) - 48
      /Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (1.9605920988138417*(0.985*(-1.01*
      mse2(2,2) + 0.99*msl2(2,2))*(-6.0600000000000005*mse2(2,2) + 0.985*Sqr(
      mAInput))*Sqr(mAInput) - (-3.0300000000000002*mse2(2,2) + 0.99*msl2(2,2))*
      TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*mse2(2,2)))*TPhi(0.985*Sqr(
      mAInput),1.01*mse2(2,2),1.01*mse2(2,2)))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2
      ,2))*Sqr(mse2(2,2))) - (2.0202020202020203*((-1.01*mse2(2,2) + 0.99*msl2(2,2
      ))*(-1.01*mse2(2,2) + 0.99*msl2(2,2) + 0.985*Sqr(mAInput)) - 2*TDelta(0.985*
      Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))*TPhi(0.985*Sqr(mAInput),1.01*
      mse2(2,2),0.99*msl2(2,2)))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)
      ) + (2.0406081012141617*(0.985*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-6.93*
      msl2(2,2) + 0.985*Sqr(mAInput))*Sqr(mAInput) + (1.01*mse2(2,2) - 4.95*msl2(2
      ,2))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))*TPhi(0.985*
      Sqr(mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))/(Cube(-1.01*mse2(2,2) + 0.99*
      msl2(2,2))*Sqr(msl2(2,2)))) + Sqr(AtauInput - MuInput*TanBeta)*((-2*Sqr(Log(
      1.01*mse2(2,2))))/(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + Log(0.985*Sqr(mAInput
      ))*((1.97019701970197*(1.01*mse2(2,2) - 1.98*msl2(2,2))*Sqr(mAInput))/(mse2(
      2,2)*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)) - (Log(1.01*mse2(2,2))*(
      5.05*mse2(2,2) - 4.95*msl2(2,2) + 0.985*Sqr(mAInput)))/Sqr(-1.01*mse2(2,2) +
      0.99*msl2(2,2)) + (Log(0.99*msl2(2,2))*(5.05*mse2(2,2) - 4.95*msl2(2,2) +
      0.985*Sqr(mAInput)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + Log(0.99*msl2(
      2,2))*((3.0300000000000002*mse2(2,2) + 0.99*msl2(2,2) - 1.97*Sqr(mAInput))/
      Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (Log(1.01*mse2(2,2))*(1.01*mse2(2,2)
      + 0.99*msl2(2,2) - 0.985*Sqr(mAInput)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)
      ) - (2*Log(Sqr(SCALE))*(-0.99*msl2(2,2) + 0.985*Sqr(mAInput)))/Sqr(-1.01*
      mse2(2,2) + 0.99*msl2(2,2))) + (Log(1.01*mse2(2,2))*(-5.05*mse2(2,2) + 0.99*
      msl2(2,2) + 1.97*Sqr(mAInput)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + ((-
      3.0300000000000002*mse2(2,2) + 0.99*msl2(2,2) + 0.985*Sqr(mAInput))*Sqr(Log(
      0.99*msl2(2,2))))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (-3.9996*mse2(2,2)
      *msl2(2,2) - 1.9897*mse2(2,2)*Sqr(mAInput) + 3.9006*msl2(2,2)*Sqr(mAInput))/
      (-1.009899*msl2(2,2)*Sqr(mse2(2,2)) + 0.989901*mse2(2,2)*Sqr(msl2(2,2))) +
      Log(Sqr(SCALE))*((2*Log(1.01*mse2(2,2))*(-0.99*msl2(2,2) + 0.985*Sqr(mAInput
      )))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (-1.9998*mse2(2,2)*msl2(2,2) -
      1.9897*mse2(2,2)*Sqr(mAInput) + 3.9006*msl2(2,2)*Sqr(mAInput))/(-1.009899*
      msl2(2,2)*Sqr(mse2(2,2)) + 0.989901*mse2(2,2)*Sqr(msl2(2,2)))) + (
      0.9802960494069208*(-0.985*(-6.0600000000000005*mse2(2,2) + 0.985*Sqr(
      mAInput))*Sqr(mAInput) + TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*mse2(
      2,2)))*TPhi(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*mse2(2,2)))/((-1.01*mse2(
      2,2) + 0.99*msl2(2,2))*Sqr(mse2(2,2))) + (1.0101010101010102*TDelta(0.985*
      Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))*TPhi(0.985*Sqr(mAInput),1.01*
      mse2(2,2),0.99*msl2(2,2)))/(msl2(2,2)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)))
      + (1.0203040506070808*(0.985*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-6.93*msl2(
      2,2) + 0.985*Sqr(mAInput))*Sqr(mAInput) + (1.01*mse2(2,2) - 1.98*msl2(2,2))*
      TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))*TPhi(0.985*Sqr(
      mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))/(Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,
      2))*Sqr(msl2(2,2))) + Sqr(AtauInput + MuInput/TanBeta)*((2*Sqr(Log(1.01*mse2
      (2,2))))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (3*Sqr(Log(0.99*msl2(2,2)))
      )/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (-2.02*mse2(2,2) + 3.96*msl2(2,2))
      /(-1.009899*msl2(2,2)*Sqr(mse2(2,2)) + 0.989901*mse2(2,2)*Sqr(msl2(2,2))) +
      Log(Sqr(SCALE))*((2*Log(1.01*mse2(2,2)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2
      )) + (-2.02*mse2(2,2) + 3.96*msl2(2,2))/(-1.009899*msl2(2,2)*Sqr(mse2(2,2))
      + 0.989901*mse2(2,2)*Sqr(msl2(2,2)))) + (1.000100010001*Log(1.01*mse2(2,2))*
      (-((-1.01*mse2(2,2) + 0.99*msl2(2,2))*(0.970225*(1.01*mse2(2,2) + 1.98*msl2(
      2,2))*Quad(mAInput) - 3.9006*(2.02*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*Sqr
      (mAInput) + (-1.01*mse2(2,2) + 1.98*msl2(2,2))*Sqr(-1.01*mse2(2,2) + 0.99*
      msl2(2,2)))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))) + (
      1.9998*mse2(2,2)*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*(-1.01*mse2(2,
      2) + 0.99*msl2(2,2) + 0.985*Sqr(mAInput)) + (0.9999*mse2(2,2)*msl2(2,2) -
      1.0201*Sqr(mse2(2,2)) + 1.9602*Sqr(msl2(2,2)))*TDelta(0.985*Sqr(mAInput),
      1.01*mse2(2,2),0.99*msl2(2,2)))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),
      1.01*mse2(2,2))))/(mse2(2,2)*msl2(2,2)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))
      *TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))*TDelta(0.985*Sqr(
      mAInput),0.99*msl2(2,2),1.01*mse2(2,2))) + Log(0.99*msl2(2,2))*((-5*Log(1.01
      *mse2(2,2)))/Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (2*Log(Sqr(SCALE)))/Sqr
      (-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (0.9802960494069208*((-1.01*mse2(2,2) +
      0.99*msl2(2,2))*(-0.955671625*Power6(mAInput) + 0.970225*(4.04*mse2(2,2) +
      2.9699999999999998*msl2(2,2))*Quad(mAInput) + 0.99*msl2(2,2)*Sqr(-1.01*mse2(
      2,2) + 0.99*msl2(2,2)) - 0.985*Sqr(mAInput)*(1.9998*mse2(2,2)*msl2(2,2) +
      7.1407*Sqr(mse2(2,2)) + 2.9402999999999997*Sqr(msl2(2,2))))*TDelta(0.985*Sqr
      (mAInput),1.01*mse2(2,2),0.99*msl2(2,2)) + (-1.01*mse2(2,2)*(-1.01*mse2(2,2)
      + 0.99*msl2(2,2))*(-0.970225*Quad(mAInput) + 1.9897*mse2(2,2)*Sqr(mAInput) -
      1.0201*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))) + (-1.9998*mse2(2,2)*msl2(2,2
      ) + 0.985*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(mAInput) + 1.0201*Sqr(mse2(
      2,2)) - 0.9801*Sqr(msl2(2,2)))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99
      *msl2(2,2)))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2))))/(Sqr
      (mse2(2,2))*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))*TDelta(0.985*Sqr(mAInput),
      1.01*mse2(2,2),0.99*msl2(2,2))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01
      *mse2(2,2)))) + Log(0.985*Sqr(mAInput))*(-(Log(1.01*mse2(2,2))/Sqr(-1.01*
      mse2(2,2) + 0.99*msl2(2,2))) + Log(0.99*msl2(2,2))/Sqr(-1.01*mse2(2,2) +
      0.99*msl2(2,2)) + (0.9901980297039604*((0.94611490875*msl2(2,2)*Power6(
      mAInput) - Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + 0.970225*Quad(mAInput)*(
      -1.9998*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) - 2.9402999999999997*Sqr
      (msl2(2,2))) + 0.97515*msl2(2,2)*Sqr(mAInput)*(-1.9998*mse2(2,2)*msl2(2,2) -
      1.0201*Sqr(mse2(2,2)) + 2.9402999999999997*Sqr(msl2(2,2))))*TDelta(0.985*Sqr
      (mAInput),1.01*mse2(2,2),0.99*msl2(2,2)) + (0.9999*mse2(2,2)*msl2(2,2)*(-
      0.970225*Quad(mAInput) + Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + (-2.9997*
      mse2(2,2)*msl2(2,2) - 0.97515*msl2(2,2)*Sqr(mAInput) + 1.0201*Sqr(mse2(2,2))
      + 0.9801*Sqr(msl2(2,2)))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(
      2,2)))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2))))/((-1.01*
      mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2)*Sqr(mse2(2,2))*TDelta(0.985*Sqr(
      mAInput),1.01*mse2(2,2),0.99*msl2(2,2))*TDelta(0.985*Sqr(mAInput),0.99*msl2(
      2,2),1.01*mse2(2,2)))) + (0.9802960494069208*(0.985*(-6.0600000000000005*
      mse2(2,2) + 0.985*Sqr(mAInput))*Sqr(mAInput) - TDelta(0.985*Sqr(mAInput),
      1.01*mse2(2,2),1.01*mse2(2,2)))*TPhi(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*
      mse2(2,2)))/(Sqr(mse2(2,2))*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + (
      1.0101010101010102*((-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.985*Sqr(mAInput
      ) - 1.01*mse2(2,2) + 0.99*msl2(2,2)) + 0.985*Sqr(mAInput)*TDelta(0.985*Sqr(
      mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))*TPhi(0.985*Sqr(mAInput),1.01*mse2(2
      ,2),0.99*msl2(2,2)))/(msl2(2,2)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))*TDelta
      (0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))) + (0.9705901479276445*(-
      ((-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-3.8226865*(1.01*mse2(2,2) + 0.99*msl2(
      2,2))*Power6(mAInput) + 0.941336550625*Power8(mAInput) + Quad(-1.01*mse2(2,2
      ) + 0.99*msl2(2,2)) - 3.8615939999999997*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*
      Sqr(mAInput)*Sqr(msl2(2,2)) + 0.970225*Quad(mAInput)*(3.9996*mse2(2,2)*msl2(
      2,2) + 6.1206*Sqr(mse2(2,2)) + 5.880599999999999*Sqr(msl2(2,2))))) - 2*(-
      2.02*mse2(2,2) + 0.99*msl2(2,2))*Sqr(TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2
      ),1.01*mse2(2,2))) + (0.970225*(-5.05*mse2(2,2) + 2.9699999999999998*msl2(2,
      2))*Quad(mAInput) + (-5.05*mse2(2,2) + 2.9699999999999998*msl2(2,2))*Sqr(-
      1.01*mse2(2,2) + 0.99*msl2(2,2)) + 0.985*Sqr(mAInput)*(3.9996*mse2(2,2)*msl2
      (2,2) + 14.2814*Sqr(mse2(2,2)) - 5.880599999999999*Sqr(msl2(2,2))))*TDelta(
      0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2)))*TPhi(0.985*Sqr(mAInput),
      0.99*msl2(2,2),1.01*mse2(2,2)))/(Cube(mse2(2,2))*Sqr(-1.01*mse2(2,2) + 0.99*
      msl2(2,2))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2))) + (
      1.0203040506070808*(0.985*(-6.93*msl2(2,2) + 0.985*Sqr(mAInput))*Sqr(mAInput
      ) - TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))*TPhi(0.985*Sqr
      (mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))/(Sqr(-1.01*mse2(2,2) + 0.99*msl2(2
      ,2))*Sqr(msl2(2,2))))) + Quad(AtauInput - MuInput*TanBeta)*((1.000100010001*
      (2.9997*mse2(2,2)*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*msl2(2,2) + 0.985*Sqr(
      mAInput)*(-4.9995*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(
      msl2(2,2)))))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2))*mse2(2,2)*msl2(2,2)) +
      Log(0.99*msl2(2,2))*((3*Log(Sqr(SCALE))*(1.9503*msl2(2,2)*Sqr(mAInput) +
      1.0201*Sqr(mse2(2,2)) - 0.9801*Sqr(msl2(2,2))))/Quad(-1.01*mse2(2,2) + 0.99*
      msl2(2,2)) + (1.5*(3.9006*msl2(2,2)*Sqr(mAInput) + 1.0201*Sqr(mse2(2,2)) -
      0.9801*Sqr(msl2(2,2))))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))) + Log(Sqr(
      SCALE))*((1.000100010001*(5.9994*mse2(2,2)*(-1.01*mse2(2,2) + 0.99*msl2(2,2)
      )*msl2(2,2) + 0.985*Sqr(mAInput)*(-4.9995*mse2(2,2)*msl2(2,2) + 1.0201*Sqr(
      mse2(2,2)) - 1.9602*Sqr(msl2(2,2)))))/(Cube(-1.01*mse2(2,2) + 0.99*msl2(2,2)
      )*mse2(2,2)*msl2(2,2)) + (3*Log(1.01*mse2(2,2))*(-1.9503*msl2(2,2)*Sqr(
      mAInput) - 1.0201*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))))/Quad(-1.01*mse2(2
      ,2) + 0.99*msl2(2,2))) + (1.5*Log(1.01*mse2(2,2))*(-3.9006*msl2(2,2)*Sqr(
      mAInput) - 1.0201*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))))/Quad(-1.01*mse2(2
      ,2) + 0.99*msl2(2,2)) + Log(0.985*Sqr(mAInput))*((3*Log(1.01*mse2(2,2))*(
      1.9503*msl2(2,2)*Sqr(mAInput) + 1.0201*Sqr(mse2(2,2)) - 0.9801*Sqr(msl2(2,2)
      )))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (3*Log(0.99*msl2(2,2))*(-1.9503
      *msl2(2,2)*Sqr(mAInput) - 1.0201*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))))/
      Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (1.000100010001*(5.9994*mse2(2,2)*(
      1.01*mse2(2,2) - 0.99*msl2(2,2))*msl2(2,2) + 0.985*Sqr(mAInput)*(4.9995*mse2
      (2,2)*msl2(2,2) - 1.0201*Sqr(mse2(2,2)) + 1.9602*Sqr(msl2(2,2)))))/(Cube(-
      1.01*mse2(2,2) + 0.99*msl2(2,2))*mse2(2,2)*msl2(2,2))) + Sqr(AtauInput +
      MuInput/TanBeta)*((-2*(3.0300000000000002*mse2(2,2) + 0.99*msl2(2,2) - 1.97*
      Sqr(mAInput))*Sqr(Log(1.01*mse2(2,2))))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2
      )) + ((-3.0300000000000002*mse2(2,2) - 10.89*msl2(2,2) + 6.895*Sqr(mAInput))
      *Sqr(Log(0.99*msl2(2,2))))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + Log(Sqr(
      SCALE))*((-5.9399999999999995*Log(1.01*mse2(2,2))*msl2(2,2))/Quad(-1.01*mse2
      (2,2) + 0.99*msl2(2,2)) + (1.000100010001*(-4.9995*mse2(2,2)*msl2(2,2) +
      1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2))))/(Cube(-1.01*mse2(2,2) + 0.99
      *msl2(2,2))*mse2(2,2)*msl2(2,2))) + (1.000100010001*(-10.9989*mse2(2,2)*msl2
      (2,2) + 1.0201*Sqr(mse2(2,2)) - 1.9602*Sqr(msl2(2,2))))/(Cube(-1.01*mse2(2,2
      ) + 0.99*msl2(2,2))*mse2(2,2)*msl2(2,2)) - (0.5000500050005*Log(1.01*mse2(2,
      2))*(Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-3.090903*Cube(mse2(2,2)) -
      2.910897*Cube(msl2(2,2)) + 0.955671625*Power6(mAInput) - 0.970225*(
      3.0300000000000002*mse2(2,2) + 4.95*msl2(2,2))*Quad(mAInput) - 5.049495*msl2
      (2,2)*Sqr(mse2(2,2)) + 10.888911*mse2(2,2)*Sqr(msl2(2,2)) + 0.985*Sqr(
      mAInput)*(3.9996*mse2(2,2)*msl2(2,2) + 5.1005*Sqr(mse2(2,2)) + 6.8607*Sqr(
      msl2(2,2))))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)) + (
      1.9998*mse2(2,2)*msl2(2,2)*(-1.01*mse2(2,2) + 0.99*msl2(2,2) + 0.985*Sqr(
      mAInput))*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (1.030301*Cube(mse2(2,2))
      + 2.910897*Cube(msl2(2,2)) + 13.128687000000001*msl2(2,2)*Sqr(mse2(2,2)) -
      0.985*Sqr(mAInput)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) +
      6.9293070000000005*mse2(2,2)*Sqr(msl2(2,2)))*TDelta(0.985*Sqr(mAInput),1.01*
      mse2(2,2),0.99*msl2(2,2)))*TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*
      mse2(2,2))))/(mse2(2,2)*msl2(2,2)*Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))*
      TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))*TDelta(0.985*Sqr(
      mAInput),0.99*msl2(2,2),1.01*mse2(2,2))) + Log(0.985*Sqr(mAInput))*((3*Log(
      1.01*mse2(2,2))*(1.01*mse2(2,2) - 0.99*msl2(2,2) + 0.985*Sqr(mAInput)))/Quad
      (-1.01*mse2(2,2) + 0.99*msl2(2,2)) - (3*Log(0.99*msl2(2,2))*(1.01*mse2(2,2)
      - 0.99*msl2(2,2) + 0.985*Sqr(mAInput)))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2
      )) + (0.4950990148519802*((-1.01*mse2(2,2) + 0.99*msl2(2,2))*(3.090903*Cube(
      mse2(2,2)) + 0.9702989999999999*Cube(msl2(2,2)) - 0.955671625*Power6(mAInput
      ) + 2.910675*(1.01*mse2(2,2) + 0.99*msl2(2,2))*Quad(mAInput) - 1.009899*msl2
      (2,2)*Sqr(mse2(2,2)) - 2.9697029999999995*mse2(2,2)*Sqr(msl2(2,2)) - 0.985*
      Sqr(mAInput)*(5.1005*Sqr(mse2(2,2)) + 2.9402999999999997*Sqr(msl2(2,2))))*
      TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)) + (-0.9999*mse2(2,2
      )*msl2(2,2)*(-0.970225*Quad(mAInput) + Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))
      ) + (2.9997*mse2(2,2)*msl2(2,2) + 0.985*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*
      Sqr(mAInput) + 3.0603*Sqr(mse2(2,2)) - 0.9801*Sqr(msl2(2,2)))*TDelta(0.985*
      Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))*TDelta(0.985*Sqr(mAInput),0.99*
      msl2(2,2),1.01*mse2(2,2))))/(msl2(2,2)*Sqr(mse2(2,2))*Sqr(-1.01*mse2(2,2) +
      0.99*msl2(2,2))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))*
      TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2)))) + Log(0.99*msl2(2
      ,2))*((5.9399999999999995*Log(Sqr(SCALE))*msl2(2,2))/Quad(-1.01*mse2(2,2) +
      0.99*msl2(2,2)) + (Log(1.01*mse2(2,2))*(9.09*mse2(2,2) + 12.87*msl2(2,2) -
      10.834999999999999*Sqr(mAInput)))/Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + (
      0.4901480247034604*(Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-
      9.272708999999999*Cube(mse2(2,2)) - 0.9702989999999999*Cube(msl2(2,2)) +
      0.955671625*Power6(mAInput) - 0.970225*(5.05*mse2(2,2) + 2.9699999999999998*
      msl2(2,2))*Quad(mAInput) + 9.089091000000002*msl2(2,2)*Sqr(mse2(2,2)) +
      0.989901*mse2(2,2)*Sqr(msl2(2,2)) + 0.985*Sqr(mAInput)*(3.9996*mse2(2,2)*
      msl2(2,2) + 9.1809*Sqr(mse2(2,2)) + 2.9402999999999997*Sqr(msl2(2,2))))*
      TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)) + (1.01*mse2(2,2)*
      Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-0.970225*Quad(mAInput) + 1.9897*mse2
      (2,2)*Sqr(mAInput) - 1.0201*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2))) + (
      4.121204*Cube(mse2(2,2)) + 0.9702989999999999*Cube(msl2(2,2)) +
      17.168283000000002*msl2(2,2)*Sqr(mse2(2,2)) - 0.985*Sqr(mAInput)*Sqr(-1.01*
      mse2(2,2) + 0.99*msl2(2,2)) + 1.979802*mse2(2,2)*Sqr(msl2(2,2)))*TDelta(
      0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2)))*TDelta(0.985*Sqr(mAInput)
      ,0.99*msl2(2,2),1.01*mse2(2,2))))/(Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))*
      Sqr(mse2(2,2))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))*
      TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2)))) + (
      0.9802960494069208*(0.985*(1.01*mse2(2,2) - 0.99*msl2(2,2))*(-
      6.0600000000000005*mse2(2,2) + 0.985*Sqr(mAInput))*Sqr(mAInput) + (-5.05*
      mse2(2,2) + 0.99*msl2(2,2))*TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*
      mse2(2,2)))*TPhi(0.985*Sqr(mAInput),1.01*mse2(2,2),1.01*mse2(2,2)))/(Quad(-
      1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(mse2(2,2))) - (0.5050505050505051*(Sqr(
      -1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(0.985*Sqr(mAInput) - 1.01*mse2(2,2) +
      0.99*msl2(2,2)) - 6*Sqr(TDelta(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2
      ,2))) + (-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-3.0300000000000002*mse2(2,2) +
      2.9699999999999998*msl2(2,2) + 3.94*Sqr(mAInput))*TDelta(0.985*Sqr(mAInput),
      1.01*mse2(2,2),0.99*msl2(2,2)))*TPhi(0.985*Sqr(mAInput),1.01*mse2(2,2),0.99*
      msl2(2,2)))/(msl2(2,2)*Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))*TDelta(0.985*
      Sqr(mAInput),1.01*mse2(2,2),0.99*msl2(2,2))) + (0.48529507396382227*(Sqr(-
      1.01*mse2(2,2) + 0.99*msl2(2,2))*(-3.8226865*(1.01*mse2(2,2) + 0.99*msl2(2,2
      ))*Power6(mAInput) + 0.941336550625*Power8(mAInput) - 3.94*(1.01*mse2(2,2) +
      0.99*msl2(2,2))*Sqr(mAInput)*Sqr(-1.01*mse2(2,2) + 0.99*msl2(2,2)) + Sqr(-
      1.01*mse2(2,2) + 0.99*msl2(2,2))*(-1.9998*mse2(2,2)*msl2(2,2) - 3.0603*Sqr(
      mse2(2,2)) + 0.9801*Sqr(msl2(2,2))) + 0.970225*Quad(mAInput)*(3.9996*mse2(2,
      2)*msl2(2,2) + 6.1206*Sqr(mse2(2,2)) + 5.880599999999999*Sqr(msl2(2,2)))) +
      2*(-4.9995*mse2(2,2)*msl2(2,2) + 12.2412*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,
      2)))*Sqr(TDelta(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2))) - (-1.01*
      mse2(2,2) + 0.99*msl2(2,2))*(-17.515117*Cube(mse2(2,2)) + 2.910897*Cube(msl2
      (2,2)) + 2.910675*(-3.0300000000000002*mse2(2,2) + 0.99*msl2(2,2))*Quad(
      mAInput) + 29.287071000000005*msl2(2,2)*Sqr(mse2(2,2)) - 5.91*Sqr(mAInput)*(
      -1.9998*mse2(2,2)*msl2(2,2) - 3.0603*Sqr(mse2(2,2)) + 0.9801*Sqr(msl2(2,2)))
      - 14.848514999999999*mse2(2,2)*Sqr(msl2(2,2)))*TDelta(0.985*Sqr(mAInput),
      0.99*msl2(2,2),1.01*mse2(2,2)))*TPhi(0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*
      mse2(2,2)))/(Cube(mse2(2,2))*Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))*TDelta(
      0.985*Sqr(mAInput),0.99*msl2(2,2),1.01*mse2(2,2))) + (1.0203040506070808*(
      0.985*(-1.01*mse2(2,2) + 0.99*msl2(2,2))*(-6.93*msl2(2,2) + 0.985*Sqr(
      mAInput))*Sqr(mAInput) + (1.01*mse2(2,2) - 7.92*msl2(2,2))*TDelta(0.985*Sqr(
      mAInput),0.99*msl2(2,2),0.99*msl2(2,2)))*TPhi(0.985*Sqr(mAInput),0.99*msl2(2
      ,2),0.99*msl2(2,2)))/(Quad(-1.01*mse2(2,2) + 0.99*msl2(2,2))*Sqr(msl2(2,2)))
      )))))/Quad(3.141592653589793), 0))*UnitStep(-2 + LambdaLoopOrder) + UnitStep
      (-1 + LambdaLoopOrder)*(0.006332573977646111*(-0.09*Quad(g1) - 0.3*Sqr(g1)*
      Sqr(g2) - Quad(g2)*(0.75 - 0.16666666666666666*Sqr(Cos(2*ArcTan(TanBeta)))))
      - 0.0010554289962743518*(2*Log(Sqr(M2Input)/Sqr(SCALE))*Quad(g2)*(1 + (
      DeltaEFT*Sqr(v))/Sqr(M2Input)) + Log(Sqr(MuInput)/Sqr(SCALE))*(0.36*Quad(g1)
      + Quad(g2))*(1 + (DeltaEFT*Sqr(v))/Sqr(MuInput)))*Sqr(Cos(2*ArcTan(TanBeta))
      ) + 0.006332573977646111*(0.5*Log(Sqr(MuInput)/Sqr(SCALE))*(1 + (DeltaEFT*
      Sqr(v))/Sqr(MuInput))*(-2*(Sqr(g2)/(1 + Sqr(TanBeta)) + (0.6*Sqr(g1)*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)))*((0.6*Sqr(g1))/(1 + Sqr(TanBeta)) + (Sqr(g2)*
      Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.5*(0.6*Sqr(g1) + Sqr(g2))*((0.6*Sqr(g1
      ))/(1 + Sqr(TanBeta)) + (3*Sqr(g2))/(1 + Sqr(TanBeta)) + (0.6*Sqr(g1)*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)) + (3*Sqr(g2)*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*
      Sqr(Cos(2*ArcTan(TanBeta))) - (0.36*Quad(g1))/Sqr(1 + Sqr(TanBeta)) - (5*
      Quad(g2))/Sqr(1 + Sqr(TanBeta)) - (0.36*Quad(g1)*Quad(TanBeta))/Sqr(1 + Sqr(
      TanBeta)) - (5*Quad(g2)*Quad(TanBeta))/Sqr(1 + Sqr(TanBeta)) - (2.4*Sqr(g1)*
      Sqr(g2)*Sqr(TanBeta))/Sqr(1 + Sqr(TanBeta))) + (0.4*TanBeta*Sqr(g1)*(1 + (
      DeltaEFT*Sqr(v))/Min(Sqr(M1Input),Sqr(MuInput)))*(-2*((0.6*Sqr(g1))/(1 + Sqr
      (TanBeta)) + (0.6*Sqr(g1)*Sqr(TanBeta))/(1 + Sqr(TanBeta))) + 0.25*(0.6*Sqr(
      g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta))))*TCf0(M1Input/MuInput))/(1 + Sqr(
      TanBeta)) + (2*TanBeta*Sqr(g2)*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(M2Input),Sqr(
      MuInput)))*(-2*(Sqr(g2)/(1 + Sqr(TanBeta)) + (Sqr(g2)*Sqr(TanBeta))/(1 + Sqr
      (TanBeta))) + 0.25*(0.6*Sqr(g1) + Sqr(g2))*Sqr(Cos(2*ArcTan(TanBeta))))*TCf0
      (M2Input/MuInput))/(1 + Sqr(TanBeta)) + 0.08333333333333333*(0.6*Sqr(g1) +
      Sqr(g2))*((0.6*Sqr(g1))/(1 + Sqr(TanBeta)) + (0.6*Sqr(g1)*Sqr(TanBeta))/(1 +
      Sqr(TanBeta)))*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(M1Input),Sqr(MuInput)))*Sqr(
      Cos(2*ArcTan(TanBeta)))*TCg0(M1Input/MuInput) + 0.25*(0.6*Sqr(g1) + Sqr(g2))
      *(Sqr(g2)/(1 + Sqr(TanBeta)) + (Sqr(g2)*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*(1
       + (DeltaEFT*Sqr(v))/Min(Sqr(M2Input),Sqr(MuInput)))*Sqr(Cos(2*ArcTan(
      TanBeta)))*TCg0(M2Input/MuInput) - 0.5833333333333334*(1 + (DeltaEFT*Sqr(v))
      /Min(Sqr(M1Input),Sqr(MuInput)))*((0.36*Quad(g1))/Sqr(1 + Sqr(TanBeta)) + (
      0.36*Quad(g1)*Quad(TanBeta))/Sqr(1 + Sqr(TanBeta)))*TCf(1)(M1Input/MuInput)
      - 2.25*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(M2Input),Sqr(MuInput)))*(Quad(g2)/Sqr(
      1 + Sqr(TanBeta)) + (Quad(g2)*Quad(TanBeta))/Sqr(1 + Sqr(TanBeta)))*TCf(2)(
      M2Input/MuInput) - (0.54*Quad(g1)*Sqr(TanBeta)*(1 + (DeltaEFT*Sqr(v))/Min(
      Sqr(M1Input),Sqr(MuInput)))*TCf(3)(M1Input/MuInput))/Sqr(1 + Sqr(TanBeta)) -
      (3.5*Quad(g2)*Sqr(TanBeta)*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(M2Input),Sqr(
      MuInput)))*TCf(4)(M2Input/MuInput))/Sqr(1 + Sqr(TanBeta)) - (1.6*Sqr(g1)*Sqr
      (g2)*Sqr(TanBeta)*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(M1Input),Sqr(M2Input),Sqr(
      MuInput)))*TCf(5)(M1Input/MuInput,M2Input/MuInput))/Sqr(1 + Sqr(TanBeta)) -
      1.1666666666666667*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(M1Input),Sqr(M2Input),Sqr(
      MuInput)))*((0.6*Sqr(g1)*Sqr(g2))/Sqr(1 + Sqr(TanBeta)) + (0.6*Quad(TanBeta)
      *Sqr(g1)*Sqr(g2))/Sqr(1 + Sqr(TanBeta)))*TCf(6)(M1Input/MuInput,M2Input/
      MuInput) - (0.2*Sqr(g1)*Sqr(g2)*Sqr(TanBeta)*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(
      M1Input),Sqr(M2Input),Sqr(MuInput)))*TCf(7)(M1Input/MuInput,M2Input/MuInput)
      )/Sqr(1 + Sqr(TanBeta)) - (2.065591117977289*g1*g2*TanBeta*((
      0.7745966692414834*g1*g2)/(1 + Sqr(TanBeta)) + (0.7745966692414834*g1*g2*Sqr
      (TanBeta))/(1 + Sqr(TanBeta)))*(1 + (DeltaEFT*Sqr(v))/Min(Sqr(M1Input),Sqr(
      M2Input),Sqr(MuInput)))*TCf(8)(M1Input/MuInput,M2Input/MuInput))/(1 + Sqr(
      TanBeta))) + 0.006332573977646111*(1 + (DeltaEFT*Sqr(v))/Min(Abs(mse2(2,2)),
      Abs(msl2(2,2)),Sqr(MuInput)))*(Log(mse2(2,2)/Sqr(SCALE))*Sqr(Ye(2,2))*(-0.6*
      Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(Ye(2,2))) + Log(msl2(2,2)/Sqr(SCALE))*
      Sqr(Ye(2,2))*(-0.5*Cos(2*ArcTan(TanBeta))*(-0.6*Sqr(g1) + Sqr(g2)) + Sqr(Ye(
      2,2))) + (2*Quad(Ye(2,2))*Sqr(AtauInput - MuInput*TanBeta)*(TCF(1)(Sqrt(Abs(
      msl2(2,2)/mse2(2,2)))) - (0.08333333333333333*Sqr(AtauInput - MuInput*
      TanBeta)*TCF(2)(Sqrt(Abs(msl2(2,2)/mse2(2,2)))))/Sqrt(Abs(mse2(2,2)*msl2(2,2
      )))))/Sqrt(Abs(mse2(2,2)*msl2(2,2))) + (0.25*Cos(2*ArcTan(TanBeta))*Sqr(
      AtauInput - MuInput*TanBeta)*Sqr(Ye(2,2))*(-0.9*Sqr(g1)*TCF(3)(Sqrt(Abs(msl2
      (2,2)/mse2(2,2)))) + (0.3*Sqr(g1) - Sqr(g2))*TCF(4)(Sqrt(Abs(msl2(2,2)/mse2(
      2,2))))))/Sqrt(Abs(mse2(2,2)*msl2(2,2))) - (0.08333333333333333*(0.6*Sqr(g1)
      + Sqr(g2))*Sqr(AtauInput - MuInput*TanBeta)*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(
      Ye(2,2))*TCF(5)(Sqrt(Abs(msl2(2,2)/mse2(2,2)))))/Sqrt(Abs(mse2(2,2)*msl2(2,2
      )))) + 0.006332573977646111*(1 + (DeltaEFT*Sqr(v))/Min(Abs(msd2(2,2)),Abs(
      msq2(2,2)),Abs(msu2(2,2)),Sqr(M3Input),Sqr(MuInput)))*((3*Log(msd2(2,2)/Sqr(
      SCALE))*Sqr(Yd(2,2))*(-0.2*Cos(2*ArcTan(TanBeta))*Sqr(g1) + Sqr(Yd(2,2))/Sqr
      (1 + (0.006332573977646111*(1 + Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(
      SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) +
      ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(
      Yu(2,2)))/Sqr(TanBeta) + 0.5*((-0.0031662869888230555*Sqr(AbInput - MuInput*
      TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(
      2,2)*msq2(2,2))) - (0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr
      (Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,
      2)))) + 0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(
      MuInput)/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(
      TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(
      Sqrt(msq2(2,2))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(
      M3Input)/Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,
      2))/M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,
      Sqrt(msd2(2,2))/M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta)
      )*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput -
      MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/
      MuInput))/MuInput))/Sqr(TanBeta))))/Sqr(1 + (0.006332573977646111*(1 + Sqr(
      TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(
      SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/Sqr(TanBeta) + 0.5*
      ((-0.0031662869888230555*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*msq2(2,2))) - (
      0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))) +
      0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/
      Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1
      + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2
      ))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)
      ) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((
      AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2))/
      M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(
      0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*
      TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/
      Sqr(TanBeta)) + (3*Log(msq2(2,2)/Sqr(SCALE))*Sqr(Yd(2,2))*(-0.5*Cos(2*ArcTan
      (TanBeta))*(0.2*Sqr(g1) + Sqr(g2)) + Sqr(Yd(2,2))/Sqr(1 + (
      0.006332573977646111*(1 + Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) +
      (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 +
      Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2))
      )/Sqr(TanBeta) + 0.5*((-0.0031662869888230555*Sqr(AbInput - MuInput*TanBeta)
      *Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*
      msq2(2,2))) - (0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2
      ,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2))))
      + 0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput
      )/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(
      1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2
      ,2))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(
      SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input)
      - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2)
      )/M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))
      *(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*
      TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/
      Sqr(TanBeta))))/Sqr(1 + (0.006332573977646111*(1 + Sqr(TanBeta))*(0.25*Log(
      Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1
      + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 +
      Sqr(TanBeta)))*Sqr(Yu(2,2)))/Sqr(TanBeta) + 0.5*((-0.0031662869888230555*Sqr
      (AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,2)
      ))))/Sqrt(Abs(msd2(2,2)*msq2(2,2))) - (0.0031662869888230555*Sqr(AtInput -
      MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(
      Abs(msq2(2,2)*msu2(2,2)))) + 0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(
      2,2))*(0.75*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/
      Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/
      MuInput) + TCF(6)(Sqrt(msq2(2,2))/MuInput)) + 0.008443431970194815*Sqr(g3)*(
      1 + Log(Sqr(M3Input)/Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(
      Sqrt(msq2(2,2))/M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2
      ))/M3Input,Sqrt(msd2(2,2))/M3Input))/M3Input) + (0.006332573977646111*(1 +
      Sqr(TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(
      AtInput - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(
      2,2))/MuInput))/MuInput))/Sqr(TanBeta)) + (6*Quad(Yd(2,2))*Sqr(AbInput -
      MuInput*TanBeta)*(TCF(1)(Sqrt(Abs(msq2(2,2)/msd2(2,2)))) - (
      0.08333333333333333*Sqr(AbInput - MuInput*TanBeta)*TCF(2)(Sqrt(Abs(msq2(2,2)
      /msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*msq2(2,2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2
      )))*Quad(1 + (0.006332573977646111*(1 + Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)
      /Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE))))/(1 + Sqr(
      TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(
      TanBeta)))*Sqr(Yu(2,2)))/Sqr(TanBeta) + 0.5*((-0.0031662869888230555*Sqr(
      AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,2))
      )))/Sqrt(Abs(msd2(2,2)*msq2(2,2))) - (0.0031662869888230555*Sqr(AtInput -
      MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(
      Abs(msq2(2,2)*msu2(2,2)))) + 0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(
      2,2))*(0.75*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/
      Sqr(SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/
      MuInput) + TCF(6)(Sqrt(msq2(2,2))/MuInput)) + 0.008443431970194815*Sqr(g3)*(
      1 + Log(Sqr(M3Input)/Sqr(SCALE)) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(
      Sqrt(msq2(2,2))/M3Input) - ((AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2
      ))/M3Input,Sqrt(msd2(2,2))/M3Input))/M3Input) + (0.006332573977646111*(1 +
      Sqr(TanBeta))*Sqr(Yu(2,2))*(0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(
      AtInput - MuInput/TanBeta)*TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(
      2,2))/MuInput))/MuInput))/Sqr(TanBeta))) + (0.75*Cos(2*ArcTan(TanBeta))*Sqr(
      AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*(-0.3*Sqr(g1)*TCF(3)(Sqrt(Abs(msq2(2
      ,2)/msd2(2,2)))) + (-0.3*Sqr(g1) - Sqr(g2))*TCF(4)(Sqrt(Abs(msq2(2,2)/msd2(2
      ,2))))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(1 + (0.006332573977646111*(1 +
      Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(
      SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/Sqr(TanBeta) + 0.5*
      ((-0.0031662869888230555*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*msq2(2,2))) - (
      0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))) +
      0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/
      Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1
      + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2
      ))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)
      ) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((
      AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2))/
      M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(
      0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*
      TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/
      Sqr(TanBeta))) - (0.25*(0.6*Sqr(g1) + Sqr(g2))*Sqr(AbInput - MuInput*TanBeta
      )*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Yd(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/msd2(2,
      2)))))/(Sqrt(Abs(msd2(2,2)*msq2(2,2)))*Sqr(1 + (0.006332573977646111*(1 +
      Sqr(TanBeta))*(0.25*Log(Sqr(MuInput)/Sqr(SCALE)) + (0.125*(-1 + 2*Log(Sqr(
      mAInput)/Sqr(SCALE))))/(1 + Sqr(TanBeta)) + ((-1 + Log(Sqr(mAInput)/Sqr(
      SCALE)))*Sqr(TanBeta))/(1 + Sqr(TanBeta)))*Sqr(Yu(2,2)))/Sqr(TanBeta) + 0.5*
      ((-0.0031662869888230555*Sqr(AbInput - MuInput*TanBeta)*Sqr(Yd(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msd2(2,2)))))/Sqrt(Abs(msd2(2,2)*msq2(2,2))) - (
      0.0031662869888230555*Sqr(AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*TCF(5)(
      Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2)))) +
      0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yd(2,2))*(0.75*Log(Sqr(MuInput)/
      Sqr(SCALE)) + (0.375*(-1 + 2*Log(Sqr(mAInput)/Sqr(SCALE)))*Sqr(TanBeta))/(1
      + Sqr(TanBeta)) + 0.5*TCF(6)(Sqrt(msd2(2,2))/MuInput) + TCF(6)(Sqrt(msq2(2,2
      ))/MuInput)) + 0.008443431970194815*Sqr(g3)*(1 + Log(Sqr(M3Input)/Sqr(SCALE)
      ) + TCF(6)(Sqrt(msd2(2,2))/M3Input) + TCF(6)(Sqrt(msq2(2,2))/M3Input) - ((
      AbInput - MuInput*TanBeta)*TCF(9)(Sqrt(msq2(2,2))/M3Input,Sqrt(msd2(2,2))/
      M3Input))/M3Input) + (0.006332573977646111*(1 + Sqr(TanBeta))*Sqr(Yu(2,2))*(
      0.5*TCF(6)(Sqrt(msu2(2,2))/MuInput) + (0.5*(AtInput - MuInput/TanBeta)*
      TanBeta*TCF(9)(Sqrt(msq2(2,2))/MuInput,Sqrt(msu2(2,2))/MuInput))/MuInput))/
      Sqr(TanBeta)))) + 0.006332573977646111*(0.0033333333333333335*(6*Quad(g1)*(
      Log(msd2(0,0)/Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msd2(0,0))) + Log(msd2(
      1,1)/Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msd2(1,1))) + Log(msd2(2,2)/Sqr(
      SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msd2(2,2)))) + 18*Quad(g1)*(Log(mse2(0,0)
      /Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(mse2(0,0))) + Log(mse2(1,1)/Sqr(
      SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(mse2(1,1))) + Log(mse2(2,2)/Sqr(SCALE))*(
      1 + (DeltaEFT*Sqr(v))/Abs(mse2(2,2)))) + (9*Quad(g1) + 25*Quad(g2))*(Log(
      msl2(0,0)/Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msl2(0,0))) + Log(msl2(1,1)
      /Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msl2(1,1))) + Log(msl2(2,2)/Sqr(
      SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msl2(2,2)))) + 3*(Quad(g1) + 25*Quad(g2))
      *(Log(msq2(0,0)/Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msq2(0,0))) + Log(
      msq2(1,1)/Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msq2(1,1))) + Log(msq2(2,2)
      /Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msq2(2,2)))) + 24*Quad(g1)*(Log(msu2
      (0,0)/Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msu2(0,0))) + Log(msu2(1,1)/Sqr
      (SCALE))*(1 + (DeltaEFT*Sqr(v))/Abs(msu2(1,1))) + Log(msu2(2,2)/Sqr(SCALE))*
      (1 + (DeltaEFT*Sqr(v))/Abs(msu2(2,2)))))*Sqr(Cos(2*ArcTan(TanBeta))) + (1 +
      (DeltaEFT*Sqr(v))/Sqr(mAInput))*(0.00020833333333333335*Log(Sqr(mAInput)/Sqr
      (SCALE))*(261*Quad(g1) + 1325*Quad(g2) + 630*Sqr(g1)*Sqr(g2) - 4*Cos(4*
      ArcTan(TanBeta))*(9*Quad(g1) + 175*Quad(g2) + 90*Sqr(g1)*Sqr(g2)) - 9*Cos(8*
      ArcTan(TanBeta))*Sqr(3*Sqr(g1) + 5*Sqr(g2))) - 0.1875*Sqr(0.6*Sqr(g1) + Sqr(
      g2))*Sqr(Sin(4*ArcTan(TanBeta)))) + 3*Log(msu2(2,2)/Sqr(SCALE))*(1 + (
      DeltaEFT*Sqr(v))/Abs(msu2(2,2)))*Sqr(Yu(2,2))*(0.4*Cos(2*ArcTan(TanBeta))*
      Sqr(g1) + Sqr(Yu(2,2))) + 3*Log(msq2(2,2)/Sqr(SCALE))*(1 + (DeltaEFT*Sqr(v))
      /Abs(msq2(2,2)))*Sqr(Yu(2,2))*(0.5*Cos(2*ArcTan(TanBeta))*(-0.2*Sqr(g1) +
      Sqr(g2)) + Sqr(Yu(2,2))) + (1 + (DeltaEFT*Sqr(v))/Min(Abs(msq2(2,2)),Abs(
      msu2(2,2))))*((6*Quad(Yu(2,2))*Sqr(AtInput - MuInput/TanBeta)*(TCF(1)(Sqrt(
      Abs(msq2(2,2)/msu2(2,2)))) - (0.08333333333333333*Sqr(AtInput - MuInput/
      TanBeta)*TCF(2)(Sqrt(Abs(msq2(2,2)/msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2
      )))))/Sqrt(Abs(msq2(2,2)*msu2(2,2))) + (0.75*Cos(2*ArcTan(TanBeta))*Sqr(
      AtInput - MuInput/TanBeta)*Sqr(Yu(2,2))*(0.6*Sqr(g1)*TCF(3)(Sqrt(Abs(msq2(2,
      2)/msu2(2,2)))) + Sqr(g2)*TCF(4)(Sqrt(Abs(msq2(2,2)/msu2(2,2))))))/Sqrt(Abs(
      msq2(2,2)*msu2(2,2))) - (0.25*(0.6*Sqr(g1) + Sqr(g2))*Sqr(AtInput - MuInput/
      TanBeta)*Sqr(Cos(2*ArcTan(TanBeta)))*Sqr(Yu(2,2))*TCF(5)(Sqrt(Abs(msq2(2,2)/
      msu2(2,2)))))/Sqrt(Abs(msq2(2,2)*msu2(2,2))))))));


   check_non_perturbative();
}

bool HSSUSY_high_scale_constraint<Two_scale>::check_non_perturbative()
{
   bool problem = false;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);
   const auto Lambdax = MODELPARAMETER(Lambdax);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Ye = MODELPARAMETER(Ye);

   if (MaxAbsValue(g1) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::g1, MaxAbsValue(g1), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::g1);
   }
   if (MaxAbsValue(g2) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::g2, MaxAbsValue(g2), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::g2);
   }
   if (MaxAbsValue(g3) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::g3, MaxAbsValue(g3), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::g3);
   }
   if (MaxAbsValue(Lambdax) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Lambdax, MaxAbsValue(Lambdax), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Lambdax);
   }
   if (MaxAbsValue(Yu(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu0_0, MaxAbsValue(Yu(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu0_0);
   }

   if (MaxAbsValue(Yu(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu0_1, MaxAbsValue(Yu(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu0_1);
   }

   if (MaxAbsValue(Yu(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu0_2, MaxAbsValue(Yu(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu0_2);
   }

   if (MaxAbsValue(Yu(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu1_0, MaxAbsValue(Yu(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu1_0);
   }

   if (MaxAbsValue(Yu(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu1_1, MaxAbsValue(Yu(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu1_1);
   }

   if (MaxAbsValue(Yu(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu1_2, MaxAbsValue(Yu(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu1_2);
   }

   if (MaxAbsValue(Yu(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu2_0, MaxAbsValue(Yu(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu2_0);
   }

   if (MaxAbsValue(Yu(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu2_1, MaxAbsValue(Yu(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu2_1);
   }

   if (MaxAbsValue(Yu(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yu2_2, MaxAbsValue(Yu(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yu2_2);
   }
   if (MaxAbsValue(Yd(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd0_0, MaxAbsValue(Yd(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd0_0);
   }

   if (MaxAbsValue(Yd(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd0_1, MaxAbsValue(Yd(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd0_1);
   }

   if (MaxAbsValue(Yd(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd0_2, MaxAbsValue(Yd(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd0_2);
   }

   if (MaxAbsValue(Yd(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd1_0, MaxAbsValue(Yd(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd1_0);
   }

   if (MaxAbsValue(Yd(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd1_1, MaxAbsValue(Yd(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd1_1);
   }

   if (MaxAbsValue(Yd(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd1_2, MaxAbsValue(Yd(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd1_2);
   }

   if (MaxAbsValue(Yd(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd2_0, MaxAbsValue(Yd(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd2_0);
   }

   if (MaxAbsValue(Yd(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd2_1, MaxAbsValue(Yd(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd2_1);
   }

   if (MaxAbsValue(Yd(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Yd2_2, MaxAbsValue(Yd(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Yd2_2);
   }
   if (MaxAbsValue(Ye(0,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye0_0, MaxAbsValue(Ye(0,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye0_0);
   }

   if (MaxAbsValue(Ye(0,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye0_1, MaxAbsValue(Ye(0,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye0_1);
   }

   if (MaxAbsValue(Ye(0,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye0_2, MaxAbsValue(Ye(0,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye0_2);
   }

   if (MaxAbsValue(Ye(1,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye1_0, MaxAbsValue(Ye(1,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye1_0);
   }

   if (MaxAbsValue(Ye(1,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye1_1, MaxAbsValue(Ye(1,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye1_1);
   }

   if (MaxAbsValue(Ye(1,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye1_2, MaxAbsValue(Ye(1,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye1_2);
   }

   if (MaxAbsValue(Ye(2,0)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye2_0, MaxAbsValue(Ye(2,0)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye2_0);
   }

   if (MaxAbsValue(Ye(2,1)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye2_1, MaxAbsValue(Ye(2,1)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye2_1);
   }

   if (MaxAbsValue(Ye(2,2)) > 3.5449077018110318) {
      problem = true;
      model->get_problems().flag_non_perturbative_parameter(HSSUSY_info::Ye2_2, MaxAbsValue(Ye(2,2)), model->get_scale(), 3.5449077018110318);
   } else {
      model->get_problems().unflag_non_perturbative_parameter(HSSUSY_info::Ye2_2);
   }


   return problem;
}

double HSSUSY_high_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

double HSSUSY_high_scale_constraint<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

const HSSUSY_input_parameters& HSSUSY_high_scale_constraint<Two_scale>::get_input_parameters() const
{
   return model->get_input();
}

HSSUSY<Two_scale>* HSSUSY_high_scale_constraint<Two_scale>::get_model() const
{
   return model;
}

void HSSUSY_high_scale_constraint<Two_scale>::set_model(Model* model_)
{
   model = cast_model<HSSUSY<Two_scale>*>(model_);
}

void HSSUSY_high_scale_constraint<Two_scale>::set_scale(double s)
{
   scale = s;
}

void HSSUSY_high_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   initial_scale_guess = 0.;
   model = nullptr;
}

void HSSUSY_high_scale_constraint<Two_scale>::initialize()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto msq2 = INPUTPARAMETER(msq2);
   const auto msu2 = INPUTPARAMETER(msu2);

   initial_scale_guess = IF(MSUSY != 0, MSUSY, Sqrt(Sqrt((51200*AtInput*MuInput*
      TanBeta - 25600*Sqr(MuInput) + ((25600 + msq2(2,2))*(25600 + msu2(2,2)) -
      25600*Sqr(AtInput))*Sqr(TanBeta))/Sqr(TanBeta))));

   scale = initial_scale_guess;
}

void HSSUSY_high_scale_constraint<Two_scale>::update_scale()
{
   check_model_ptr();

   const auto MSUSY = INPUTPARAMETER(MSUSY);
   const auto TanBeta = INPUTPARAMETER(TanBeta);
   const auto MuInput = INPUTPARAMETER(MuInput);
   const auto AtInput = INPUTPARAMETER(AtInput);
   const auto msq2 = INPUTPARAMETER(msq2);
   const auto msu2 = INPUTPARAMETER(msu2);
   const auto v = MODELPARAMETER(v);
   const auto Yu = MODELPARAMETER(Yu);

   scale = IF(MSUSY != 0, MSUSY, 0.7071067811865476*Sqrt(Sqrt((2*msq2(2,2)*(1 +
      Sqr(TanBeta))*(2*msu2(2,2)*(1 + Sqr(TanBeta)) + Sqr(TanBeta)*Sqr(v)*Sqr(Yu(2
      ,2))) + Sqr(v)*Sqr(Yu(2,2))*(4*AtInput*MuInput*TanBeta*(1 + Sqr(TanBeta)) -
      2*Sqr(MuInput)*(1 + Sqr(TanBeta)) + Sqr(TanBeta)*(2*msu2(2,2)*(1 + Sqr(
      TanBeta)) - 2*Sqr(AtInput)*(1 + Sqr(TanBeta)) + Sqr(TanBeta)*Sqr(v)*Sqr(Yu(2
      ,2)))))/Sqr(1 + Sqr(TanBeta)))));


}

void HSSUSY_high_scale_constraint<Two_scale>::check_model_ptr() const
{
   if (!model)
      throw SetupError("HSSUSY_high_scale_constraint<Two_scale>: "
                       "model pointer is zero!");
}

} // namespace flexiblesusy
