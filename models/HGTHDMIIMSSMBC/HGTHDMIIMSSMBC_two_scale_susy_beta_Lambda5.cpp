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

// File generated at Tue 12 Jul 2016 10:37:30

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda5.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda5;

   beta_Lambda5 = Re(oneOver16PiSqr*(4*Lambda1*Lambda5 + 4*Lambda2*
      Lambda5 + 8*Lambda3*Lambda5 + 12*Lambda4*Lambda5 + 4*Lambda6*Lambda7 + 6*
      Lambda5*traceYdAdjYd + 2*Lambda5*traceYeAdjYe + 6*Lambda5*traceYuAdjYu -
      1.8*Lambda5*Sqr(g1) + 3*Lambda5*Sqr(g1d) + Lambda5*Sqr(g1dp) - 9*Lambda5*
      Sqr(g2) + 3*Lambda5*Sqr(g2u) + Lambda5*Sqr(g2up) + 10*Sqr(Lambda6) + 10*
      Sqr(Lambda7)));


   return beta_Lambda5;
}

/**
 * Calculates the two-loop beta function of Lambda5.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda5;

   const double beta_Lambda5_1 = Re(0.005*twoLoop*(1533*Power(g1,4)*
      Lambda5 + 5*Sqr(g1)*(-96*Lambda1*Lambda5 - 96*Lambda2*Lambda5 + 384*
      Lambda3*Lambda5 + 576*Lambda4*Lambda5 - 96*Lambda6*Lambda7 + 50*Lambda5*
      traceYdAdjYd + 150*Lambda5*traceYeAdjYe + 170*Lambda5*traceYuAdjYu + 45*
      Lambda5*Sqr(g1d) + 15*Lambda5*Sqr(g1dp) + 114*Lambda5*Sqr(g2) + 45*
      Lambda5*Sqr(g2u) + 15*Lambda5*Sqr(g2up) + 480*Sqr(Lambda6) + 480*Sqr(
      Lambda7)) - 25*(5*Power(g1d,4)*Lambda5 + Power(g1dp,4)*Lambda5 + 111*
      Power(g2,4)*Lambda5 + 5*Power(g2u,4)*Lambda5 - 128*g1d*g1dp*g2u*g2up*
      Lambda5 + Power(g2up,4)*Lambda5 + 640*Lambda1*Lambda3*Lambda5 + 640*
      Lambda2*Lambda3*Lambda5 + 704*Lambda1*Lambda4*Lambda5 + 704*Lambda2*
      Lambda4*Lambda5 + 608*Lambda3*Lambda4*Lambda5 - 48*Power(Lambda5,3) + 320
      *Lambda1*Lambda6*Lambda7 + 320*Lambda2*Lambda6*Lambda7 + 640*Lambda3*
      Lambda6*Lambda7 + 704*Lambda4*Lambda6*Lambda7 + 1344*Lambda5*Lambda6*
      Lambda7 + 192*Lambda1*Lambda5*traceYdAdjYd + 192*Lambda3*Lambda5*
      traceYdAdjYd + 288*Lambda4*Lambda5*traceYdAdjYd + 96*Lambda6*Lambda7*
      traceYdAdjYd + 12*Lambda5*traceYdAdjYdYdAdjYd + 264*Lambda5*
      traceYdAdjYuYuAdjYd + 64*Lambda1*Lambda5*traceYeAdjYe + 64*Lambda3*
      Lambda5*traceYeAdjYe + 96*Lambda4*Lambda5*traceYeAdjYe + 32*Lambda6*
      Lambda7*traceYeAdjYe + 4*Lambda5*traceYeAdjYeYeAdjYe + 192*Lambda2*
      Lambda5*traceYuAdjYu + 192*Lambda3*Lambda5*traceYuAdjYu + 288*Lambda4*
      Lambda5*traceYuAdjYu - 288*Lambda3*Lambda5*Sqr(g2) - 576*Lambda4*Lambda5*
      Sqr(g2) - 90*Lambda5*traceYdAdjYd*Sqr(g2) - 30*Lambda5*traceYeAdjYe*Sqr(
      g2) - 90*Lambda5*traceYuAdjYu*Sqr(g2) + 96*Lambda2*Lambda5*Sqr(g2u) + 96*
      Lambda3*Lambda5*Sqr(g2u) + 144*Lambda4*Lambda5*Sqr(g2u) + 48*Lambda6*
      Lambda7*Sqr(g2u) - 165*Lambda5*Sqr(g2)*Sqr(g2u) + 32*Lambda2*Lambda5*Sqr(
      g2up) + 32*Lambda3*Lambda5*Sqr(g2up) + 48*Lambda4*Lambda5*Sqr(g2up) + 16*
      Lambda6*Lambda7*Sqr(g2up) - 15*Lambda5*Sqr(g2)*Sqr(g2up) + 2*Lambda5*Sqr(
      g2u)*Sqr(g2up) - 320*Lambda5*traceYdAdjYd*Sqr(g3) - 320*Lambda5*
      traceYuAdjYu*Sqr(g3) + 224*Lambda5*Sqr(Lambda1) + 224*Lambda5*Sqr(Lambda2
      ) + 224*Lambda5*Sqr(Lambda3) + 256*Lambda5*Sqr(Lambda4) + 1184*Lambda1*
      Sqr(Lambda6) + 160*Lambda2*Sqr(Lambda6) + 576*Lambda3*Sqr(Lambda6) + 608*
      Lambda4*Sqr(Lambda6) + 576*Lambda5*Sqr(Lambda6) + 480*traceYdAdjYd*Sqr(
      Lambda6) + 160*traceYeAdjYe*Sqr(Lambda6) - 432*Sqr(g2)*Sqr(Lambda6) + Sqr
      (g1d)*(96*Lambda1*Lambda5 + 96*Lambda3*Lambda5 + 144*Lambda4*Lambda5 + 48
      *Lambda6*Lambda7 + 2*Lambda5*Sqr(g1dp) - 165*Lambda5*Sqr(g2) + 68*Lambda5
      *Sqr(g2u) + 240*Sqr(Lambda6)) + Sqr(g1dp)*(-15*Lambda5*Sqr(g2) + 4*(8*
      Lambda1*Lambda5 + 8*Lambda3*Lambda5 + 12*Lambda4*Lambda5 + 4*Lambda6*
      Lambda7 - 5*Lambda5*Sqr(g2up) + 20*Sqr(Lambda6))) + 160*Lambda1*Sqr(
      Lambda7) + 1184*Lambda2*Sqr(Lambda7) + 576*Lambda3*Sqr(Lambda7) + 608*
      Lambda4*Sqr(Lambda7) + 576*Lambda5*Sqr(Lambda7) - 432*Sqr(g2)*Sqr(Lambda7
      ) + 240*Sqr(g2u)*Sqr(Lambda7) + 80*Sqr(g2up)*Sqr(Lambda7))));
   const double beta_Lambda5_2 = Re(-1.5*twoLoop*(8*Lambda6*Lambda7*
      traceYuAdjYu + Lambda5*traceYuAdjYuYuAdjYu + 40*traceYuAdjYu*Sqr(Lambda7)
      ));

   beta_Lambda5 = beta_Lambda5_1 + beta_Lambda5_2;


   return beta_Lambda5;
}

/**
 * Calculates the three-loop beta function of Lambda5.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

} // namespace flexiblesusy
