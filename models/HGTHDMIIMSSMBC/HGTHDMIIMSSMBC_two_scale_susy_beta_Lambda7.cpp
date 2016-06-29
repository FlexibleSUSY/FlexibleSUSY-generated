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

// File generated at Wed 29 Jun 2016 11:27:56

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda7.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(oneOver16PiSqr*(6*Lambda3*Lambda6 + 4*Lambda4*
      Lambda6 + 2*Lambda5*Lambda6 + 24*Lambda2*Lambda7 + 6*Lambda3*Lambda7 + 8*
      Lambda4*Lambda7 + 10*Lambda5*Lambda7 + 3*Lambda7*traceYdAdjYd + Lambda7*
      traceYeAdjYe + 9*Lambda7*traceYuAdjYu - 1.8*Lambda7*Sqr(g1) + 1.5*Lambda7
      *Sqr(g1d) + 0.5*Lambda7*Sqr(g1dp) - 9*Lambda7*Sqr(g2) + 4.5*Lambda7*Sqr(
      g2u) + 1.5*Lambda7*Sqr(g2up)));


   return beta_Lambda7;
}

/**
 * Calculates the two-loop beta function of Lambda7.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda7;

   const double beta_Lambda7_1 = Re(0.0025*twoLoop*(6*Power(g1,4)*(90*
      Lambda6 + 601*Lambda7) + 5*Sqr(g1)*(384*Lambda4*Lambda6 - 96*Lambda5*
      Lambda6 + 1728*Lambda2*Lambda7 + 480*Lambda4*Lambda7 + 960*Lambda5*
      Lambda7 + 288*Lambda3*(2*Lambda6 + Lambda7) + 50*Lambda7*traceYdAdjYd +
      45*Lambda7*Sqr(g1d) + 15*Lambda7*Sqr(g1dp) + 12*(10*Lambda6 + 29*Lambda7)
      *Sqr(g2) + 135*Lambda7*Sqr(g2u) + 45*Lambda7*Sqr(g2up)) - 25*(576*Lambda1
      *Lambda3*Lambda6 + 576*Lambda2*Lambda3*Lambda6 + 448*Lambda1*Lambda4*
      Lambda6 + 448*Lambda2*Lambda4*Lambda6 + 896*Lambda3*Lambda4*Lambda6 + 320
      *Lambda1*Lambda5*Lambda6 + 320*Lambda2*Lambda5*Lambda6 + 640*Lambda3*
      Lambda5*Lambda6 + 704*Lambda4*Lambda5*Lambda6 + 672*Power(Lambda6,3) - 6*
      Power(g2,4)*(30*Lambda6 - 7*Lambda7) + 45*Power(g1d,4)*Lambda7 + 9*Power(
      g1dp,4)*Lambda7 + 55*Power(g2u,4)*Lambda7 - 128*g1d*g1dp*g2u*g2up*Lambda7
      + 11*Power(g2up,4)*Lambda7 + 576*Lambda1*Lambda3*Lambda7 + 2112*Lambda2*
      Lambda3*Lambda7 + 448*Lambda1*Lambda4*Lambda7 + 2240*Lambda2*Lambda4*
      Lambda7 + 1088*Lambda3*Lambda4*Lambda7 + 320*Lambda1*Lambda5*Lambda7 +
      2368*Lambda2*Lambda5*Lambda7 + 1152*Lambda3*Lambda5*Lambda7 + 1216*
      Lambda4*Lambda5*Lambda7 + 1776*Power(Lambda7,3) + 576*Lambda3*Lambda6*
      traceYdAdjYd + 384*Lambda4*Lambda6*traceYdAdjYd + 192*Lambda5*Lambda6*
      traceYdAdjYd + 288*Lambda3*Lambda7*traceYdAdjYd + 384*Lambda4*Lambda7*
      traceYdAdjYd + 480*Lambda5*Lambda7*traceYdAdjYd + 108*Lambda7*
      traceYdAdjYdYdAdjYd + 336*Lambda7*traceYdAdjYuYuAdjYd + 192*Lambda3*
      Lambda6*traceYeAdjYe + 128*Lambda4*Lambda6*traceYeAdjYe + 64*Lambda5*
      Lambda6*traceYeAdjYe + 96*Lambda3*Lambda6*Sqr(g1dp) + 64*Lambda4*Lambda6*
      Sqr(g1dp) + 32*Lambda5*Lambda6*Sqr(g1dp) + 48*Lambda3*Lambda7*Sqr(g1dp) +
      64*Lambda4*Lambda7*Sqr(g1dp) + 80*Lambda5*Lambda7*Sqr(g1dp) + 1152*
      Lambda2*Lambda7*Sqr(g2u) + 144*Lambda3*Lambda7*Sqr(g2u) + 192*Lambda4*
      Lambda7*Sqr(g2u) + 240*Lambda5*Lambda7*Sqr(g2u) + 2*Sqr(g1d)*(48*Lambda5*
      Lambda6 + 120*Lambda5*Lambda7 + 96*Lambda4*(Lambda6 + Lambda7) + 72*
      Lambda3*(2*Lambda6 + Lambda7) + 9*Lambda7*Sqr(g1dp) + 52*Lambda7*Sqr(g2u)
      ) + 384*Lambda2*Lambda7*Sqr(g2up) + 48*Lambda3*Lambda7*Sqr(g2up) + 64*
      Lambda4*Lambda7*Sqr(g2up) + 80*Lambda5*Lambda7*Sqr(g2up) - 8*Lambda7*Sqr(
      g1dp)*Sqr(g2up) + 22*Lambda7*Sqr(g2u)*Sqr(g2up) - 3*Sqr(g2)*(96*Lambda3*(
      2*Lambda6 + Lambda7) + 96*Lambda4*(Lambda6 + 2*Lambda7) + Lambda7*(576*
      Lambda2 + 288*Lambda5 + 30*traceYdAdjYd + 55*Sqr(g1d) + 5*Sqr(g1dp) + 165
      *Sqr(g2u) + 15*Sqr(g2up))) - 320*Lambda7*traceYdAdjYd*Sqr(g3) - 96*
      Lambda7*Sqr(Lambda1) + 5088*Lambda7*Sqr(Lambda2) + 576*Lambda6*Sqr(
      Lambda3) + 512*Lambda7*Sqr(Lambda3) + 544*Lambda6*Sqr(Lambda4) + 544*
      Lambda7*Sqr(Lambda4) + 672*Lambda6*Sqr(Lambda5) + 576*Lambda7*Sqr(Lambda5
      ) + 528*Lambda7*Sqr(Lambda6) + 2016*Lambda6*Sqr(Lambda7))));
   const double beta_Lambda7_2 = Re(0.125*Lambda7*twoLoop*(3*(5*
      traceYeAdjYe + 17*traceYuAdjYu)*Sqr(g1) + 15*(traceYeAdjYe + 9*
      traceYuAdjYu)*Sqr(g2) - 2*(40*Lambda5*traceYeAdjYe + 9*
      traceYeAdjYeYeAdjYe + 576*Lambda2*traceYuAdjYu + 120*Lambda5*traceYuAdjYu
      + 24*Lambda3*(traceYeAdjYe + 3*traceYuAdjYu) + 32*Lambda4*(traceYeAdjYe
      + 3*traceYuAdjYu) + 33*traceYuAdjYuYuAdjYu - 240*traceYuAdjYu*Sqr(g3))));

   beta_Lambda7 = beta_Lambda7_1 + beta_Lambda7_2;


   return beta_Lambda7;
}

/**
 * Calculates the three-loop beta function of Lambda7.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

} // namespace flexiblesusy
