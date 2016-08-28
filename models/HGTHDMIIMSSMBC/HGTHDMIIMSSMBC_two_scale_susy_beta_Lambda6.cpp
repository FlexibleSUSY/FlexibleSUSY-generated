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

// File generated at Sun 28 Aug 2016 15:02:24

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda6.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(oneOver16PiSqr*(24*Lambda1*Lambda6 + 6*Lambda3*
      Lambda6 + 8*Lambda4*Lambda6 + 10*Lambda5*Lambda6 + 6*Lambda3*Lambda7 + 4*
      Lambda4*Lambda7 + 2*Lambda5*Lambda7 + 9*Lambda6*traceYdAdjYd + 3*Lambda6*
      traceYeAdjYe + 3*Lambda6*traceYuAdjYu - 1.8*Lambda6*Sqr(g1) + 4.5*Lambda6
      *Sqr(g1d) + 1.5*Lambda6*Sqr(g1dp) - 9*Lambda6*Sqr(g2) + 1.5*Lambda6*Sqr(
      g2u) + 0.5*Lambda6*Sqr(g2up)));


   return beta_Lambda6;
}

/**
 * Calculates the two-loop beta function of Lambda6.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda6;

   const double beta_Lambda6_1 = Re(0.0025*twoLoop*(6*Power(g1,4)*(601*
      Lambda6 + 90*Lambda7) + 15*Sqr(g1)*(576*Lambda1*Lambda6 + 96*Lambda3*
      Lambda6 + 160*Lambda4*Lambda6 + 320*Lambda5*Lambda6 + 192*Lambda3*Lambda7
      + 128*Lambda4*Lambda7 - 32*Lambda5*Lambda7 + 50*Lambda6*traceYdAdjYd +
      150*Lambda6*traceYeAdjYe + 45*Lambda6*Sqr(g1d) + 15*Lambda6*Sqr(g1dp) +
      116*Lambda6*Sqr(g2) + 40*Lambda7*Sqr(g2) + 15*Lambda6*Sqr(g2u) + 5*
      Lambda6*Sqr(g2up)) - 25*(55*Power(g1d,4)*Lambda6 + 11*Power(g1dp,4)*
      Lambda6 + 42*Power(g2,4)*Lambda6 + 45*Power(g2u,4)*Lambda6 - 128*g1d*g1dp
      *g2u*g2up*Lambda6 + 9*Power(g2up,4)*Lambda6 + 2112*Lambda1*Lambda3*
      Lambda6 + 576*Lambda2*Lambda3*Lambda6 + 2240*Lambda1*Lambda4*Lambda6 +
      448*Lambda2*Lambda4*Lambda6 + 1088*Lambda3*Lambda4*Lambda6 + 2368*Lambda1
      *Lambda5*Lambda6 + 320*Lambda2*Lambda5*Lambda6 + 1152*Lambda3*Lambda5*
      Lambda6 + 1216*Lambda4*Lambda5*Lambda6 + 1776*Power(Lambda6,3) - 180*
      Power(g2,4)*Lambda7 + 576*Lambda1*Lambda3*Lambda7 + 576*Lambda2*Lambda3*
      Lambda7 + 448*Lambda1*Lambda4*Lambda7 + 448*Lambda2*Lambda4*Lambda7 + 896
      *Lambda3*Lambda4*Lambda7 + 320*Lambda1*Lambda5*Lambda7 + 320*Lambda2*
      Lambda5*Lambda7 + 640*Lambda3*Lambda5*Lambda7 + 704*Lambda4*Lambda5*
      Lambda7 + 672*Power(Lambda7,3) + 2304*Lambda1*Lambda6*traceYdAdjYd + 288*
      Lambda3*Lambda6*traceYdAdjYd + 384*Lambda4*Lambda6*traceYdAdjYd + 480*
      Lambda5*Lambda6*traceYdAdjYd + 132*Lambda6*traceYdAdjYdYdAdjYd + 336*
      Lambda6*traceYdAdjYuYuAdjYd + 768*Lambda1*Lambda6*traceYeAdjYe + 96*
      Lambda3*Lambda6*traceYeAdjYe + 128*Lambda4*Lambda6*traceYeAdjYe - 1728*
      Lambda1*Lambda6*Sqr(g2) - 288*Lambda3*Lambda6*Sqr(g2) - 576*Lambda4*
      Lambda6*Sqr(g2) - 864*Lambda5*Lambda6*Sqr(g2) - 576*Lambda3*Lambda7*Sqr(
      g2) - 288*Lambda4*Lambda7*Sqr(g2) - 270*Lambda6*traceYdAdjYd*Sqr(g2) - 90
      *Lambda6*traceYeAdjYe*Sqr(g2) + 144*Lambda3*Lambda6*Sqr(g2u) + 192*
      Lambda4*Lambda6*Sqr(g2u) + 240*Lambda5*Lambda6*Sqr(g2u) + 288*Lambda3*
      Lambda7*Sqr(g2u) + 192*Lambda4*Lambda7*Sqr(g2u) + 96*Lambda5*Lambda7*Sqr(
      g2u) - 165*Lambda6*Sqr(g2)*Sqr(g2u) + Lambda6*Sqr(g1d)*(22*Sqr(g1dp) -
      495*Sqr(g2) + 8*(6*(24*Lambda1 + 3*Lambda3 + 4*Lambda4 + 5*Lambda5) + 13*
      Sqr(g2u))) + Lambda6*Sqr(g1dp)*(-45*Sqr(g2) + 8*(48*Lambda1 + 6*Lambda3 +
      8*Lambda4 + 10*Lambda5 - Sqr(g2up))) + 48*Lambda3*Lambda6*Sqr(g2up) + 64
      *Lambda4*Lambda6*Sqr(g2up) + 80*Lambda5*Lambda6*Sqr(g2up) + 96*Lambda3*
      Lambda7*Sqr(g2up) + 64*Lambda4*Lambda7*Sqr(g2up) + 32*Lambda5*Lambda7*Sqr
      (g2up) - 15*Lambda6*Sqr(g2)*Sqr(g2up) + 18*Lambda6*Sqr(g2u)*Sqr(g2up) -
      960*Lambda6*traceYdAdjYd*Sqr(g3) + 5088*Lambda6*Sqr(Lambda1) - 96*Lambda6
      *Sqr(Lambda2) + 512*Lambda6*Sqr(Lambda3) + 576*Lambda7*Sqr(Lambda3) + 544
      *Lambda6*Sqr(Lambda4) + 544*Lambda7*Sqr(Lambda4) + 576*Lambda6*Sqr(
      Lambda5) + 672*Lambda7*Sqr(Lambda5) + 2016*Lambda7*Sqr(Lambda6) + 528*
      Lambda6*Sqr(Lambda7))));
   const double beta_Lambda6_2 = Re(-0.125*twoLoop*(96*(3*Lambda3 + 2*
      Lambda4)*Lambda7*traceYuAdjYu + 16*Lambda5*(6*Lambda7*traceYuAdjYu + 5*
      Lambda6*(traceYeAdjYe + 3*traceYuAdjYu)) + Lambda6*(22*
      traceYeAdjYeYeAdjYe + 144*Lambda3*traceYuAdjYu + 192*Lambda4*traceYuAdjYu
      + 54*traceYuAdjYuYuAdjYu - 17*traceYuAdjYu*Sqr(g1) - 45*traceYuAdjYu*Sqr
      (g2) - 160*traceYuAdjYu*Sqr(g3))));

   beta_Lambda6 = beta_Lambda6_1 + beta_Lambda6_2;


   return beta_Lambda6;
}

/**
 * Calculates the three-loop beta function of Lambda6.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

} // namespace flexiblesusy
