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

// File generated at Sun 28 Aug 2016 15:02:32

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda4.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda4;

   beta_Lambda4 = Re(oneOver16PiSqr*(-4*g1d*g1dp*g2u*g2up + 4*Lambda1*
      Lambda4 + 4*Lambda2*Lambda4 + 8*Lambda3*Lambda4 + 4*Lambda6*Lambda7 + 6*
      Lambda4*traceYdAdjYd + 12*traceYdAdjYuYuAdjYd + 2*Lambda4*traceYeAdjYe +
      6*Lambda4*traceYuAdjYu - 1.8*Lambda4*Sqr(g1) + 3*Lambda4*Sqr(g1d) +
      Lambda4*Sqr(g1dp) - 9*Lambda4*Sqr(g2) + 1.8*Sqr(g1)*Sqr(g2) + 3*Lambda4*
      Sqr(g2u) + 4*Sqr(g1d)*Sqr(g2u) + Lambda4*Sqr(g2up) + 4*Sqr(Lambda4) + 8*
      Sqr(Lambda5) + 10*Sqr(Lambda6) + 10*Sqr(Lambda7)));


   return beta_Lambda4;
}

/**
 * Calculates the two-loop beta function of Lambda4.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd;


   double beta_Lambda4;

   const double beta_Lambda4_1 = Re(-0.005*twoLoop*(3*Power(g1,4)*(-511*
      Lambda4 + 940*Sqr(g2)) - 15*Sqr(g1)*(-176*Power(g2,4) + 64*Lambda1*
      Lambda4 + 64*Lambda2*Lambda4 + 32*Lambda3*Lambda4 + 64*Lambda6*Lambda7 +
      Sqr(g1dp)*(5*Lambda4 - 4*Sqr(g2)) + 80*Lambda1*Sqr(g2) + 80*Lambda2*Sqr(
      g2) + 16*Lambda3*Sqr(g2) + 102*Lambda4*Sqr(g2) + 3*Sqr(g1d)*(5*Lambda4 +
      28*Sqr(g2)) + 15*Lambda4*Sqr(g2u) + 84*Sqr(g2)*Sqr(g2u) + 5*Lambda4*Sqr(
      g2up) - 4*Sqr(g2)*Sqr(g2up) + 64*Sqr(Lambda4) + 128*Sqr(Lambda5) + 112*
      Sqr(Lambda6) + 112*Sqr(Lambda7)) + 25*(-24*Power(g1d,3)*g1dp*g2u*g2up + 9
      *Power(g1dp,4)*Lambda4 + 111*Power(g2,4)*Lambda4 + 45*Power(g2u,4)*
      Lambda4 + 9*Power(g2up,4)*Lambda4 + 640*Lambda1*Lambda3*Lambda4 + 640*
      Lambda2*Lambda3*Lambda4 + 320*Lambda1*Lambda6*Lambda7 + 320*Lambda2*
      Lambda6*Lambda7 + 640*Lambda3*Lambda6*Lambda7 + 1280*Lambda4*Lambda6*
      Lambda7 + 768*Lambda5*Lambda6*Lambda7 - 288*Lambda3*Lambda4*Sqr(g2) + 96*
      Lambda2*Lambda4*Sqr(g2u) + 96*Lambda3*Lambda4*Sqr(g2u) + 48*Lambda6*
      Lambda7*Sqr(g2u) - 165*Lambda4*Sqr(g2)*Sqr(g2u) + 9*Power(g1d,4)*(5*
      Lambda4 + 8*Sqr(g2u)) + 32*Lambda2*Lambda4*Sqr(g2up) + 32*Lambda3*Lambda4
      *Sqr(g2up) + 16*Lambda6*Lambda7*Sqr(g2up) - 15*Lambda4*Sqr(g2)*Sqr(g2up)
      + 18*Lambda4*Sqr(g2u)*Sqr(g2up) - 8*g1d*g1dp*g2u*g2up*(8*Lambda3 + 12*
      Lambda4 + 5*Sqr(g1dp) - 8*Sqr(g2) + 3*Sqr(g2u) + 5*Sqr(g2up)) + 224*
      Lambda4*Sqr(Lambda1) + 224*Lambda4*Sqr(Lambda2) + 224*Lambda4*Sqr(Lambda3
      ) + 320*Lambda1*Sqr(Lambda4) + 320*Lambda2*Sqr(Lambda4) + 224*Lambda3*Sqr
      (Lambda4) - 144*Sqr(g2)*Sqr(Lambda4) + 48*Sqr(g2u)*Sqr(Lambda4) + 16*Sqr(
      g2up)*Sqr(Lambda4) + 384*Lambda1*Sqr(Lambda5) + 384*Lambda2*Sqr(Lambda5)
      + 384*Lambda3*Sqr(Lambda5) + 208*Lambda4*Sqr(Lambda5) - 432*Sqr(g2)*Sqr(
      Lambda5) + 96*Sqr(g2u)*Sqr(Lambda5) + 32*Sqr(g2up)*Sqr(Lambda5) + 1184*
      Lambda1*Sqr(Lambda6) + 160*Lambda2*Sqr(Lambda6) + 576*Lambda3*Sqr(Lambda6
      ) + 544*Lambda4*Sqr(Lambda6) + 640*Lambda5*Sqr(Lambda6) - 432*Sqr(g2)*Sqr
      (Lambda6) + Sqr(g1d)*(72*Power(g2u,4) + 96*Lambda1*Lambda4 + 96*Lambda3*
      Lambda4 + 48*Lambda6*Lambda7 + Sqr(g1dp)*(18*Lambda4 - 8*Sqr(g2u)) + 64*
      Lambda3*Sqr(g2u) + 84*Lambda4*Sqr(g2u) - Sqr(g2)*(165*Lambda4 + 128*Sqr(
      g2u)) - 8*Sqr(g2u)*Sqr(g2up) + 48*Sqr(Lambda4) + 96*Sqr(Lambda5) + 240*
      Sqr(Lambda6)) + Sqr(g1dp)*(-15*Lambda4*Sqr(g2) + 4*(-(Lambda4*Sqr(g2up))
      + 4*(2*Lambda1*Lambda4 + 2*Lambda3*Lambda4 + Lambda6*Lambda7 + Sqr(
      Lambda4) + 2*Sqr(Lambda5) + 5*Sqr(Lambda6)))) - 432*Sqr(g2)*Sqr(Lambda7)
      + 240*Sqr(g2u)*Sqr(Lambda7) + 80*Sqr(g2up)*Sqr(Lambda7))));
   const double beta_Lambda4_2 = Re(twoLoop*(-24*Lambda3*Lambda4*
      traceYdAdjYd - 12*Lambda6*Lambda7*traceYdAdjYd - 13.5*Lambda4*
      traceYdAdjYdYdAdjYd - 12*traceYdAdjYdYdAdjYuYuAdjYd - 24*Lambda3*
      traceYdAdjYuYuAdjYd - 33*Lambda4*traceYdAdjYuYuAdjYd - 12*
      traceYdAdjYuYuAdjYdYdAdjYd - 24*traceYdAdjYuYuAdjYuYuAdjYd - 8*Lambda3*
      Lambda4*traceYeAdjYe - 4*Lambda6*Lambda7*traceYeAdjYe - 4.5*Lambda4*
      traceYeAdjYeYeAdjYe - 24*Lambda3*Lambda4*traceYuAdjYu - 12*Lambda6*
      Lambda7*traceYuAdjYu - 13.5*Lambda4*traceYuAdjYuYuAdjYu + 1.25*Lambda4*
      traceYdAdjYd*Sqr(g1) + 0.8*traceYdAdjYuYuAdjYd*Sqr(g1) + 3.75*Lambda4*
      traceYeAdjYe*Sqr(g1) + 4.25*Lambda4*traceYuAdjYu*Sqr(g1) + 11.25*Lambda4*
      traceYdAdjYd*Sqr(g2) + 3.75*Lambda4*traceYeAdjYe*Sqr(g2) + 11.25*Lambda4*
      traceYuAdjYu*Sqr(g2) + 5.4*traceYdAdjYd*Sqr(g1)*Sqr(g2) + 6.6*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) + 12.6*traceYuAdjYu*Sqr(g1)*Sqr(g2) + 40*
      Lambda4*traceYdAdjYd*Sqr(g3) + 64*traceYdAdjYuYuAdjYd*Sqr(g3) + 40*
      Lambda4*traceYuAdjYu*Sqr(g3) - 12*traceYdAdjYd*Sqr(Lambda4) - 4*
      traceYeAdjYe*Sqr(Lambda4) - 12*traceYuAdjYu*Sqr(Lambda4) - 24*
      traceYdAdjYd*Sqr(Lambda5) - 8*traceYeAdjYe*Sqr(Lambda5) - 24*traceYuAdjYu
      *Sqr(Lambda5) - 60*traceYdAdjYd*Sqr(Lambda6) - 20*traceYeAdjYe*Sqr(
      Lambda6) - 72*Lambda3*Sqr(Lambda7) - 68*Lambda4*Sqr(Lambda7) - 80*Lambda5
      *Sqr(Lambda7) - 60*traceYuAdjYu*Sqr(Lambda7) - 4*Lambda1*(2*Lambda4*(3*
      traceYdAdjYd + traceYeAdjYe) + 5*Sqr(Lambda7)) - 4*Lambda2*(6*Lambda4*
      traceYuAdjYu + 37*Sqr(Lambda7))));

   beta_Lambda4 = beta_Lambda4_1 + beta_Lambda4_2;


   return beta_Lambda4;
}

/**
 * Calculates the three-loop beta function of Lambda4.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return beta_Lambda4;
}

} // namespace flexiblesusy
