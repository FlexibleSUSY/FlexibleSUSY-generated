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

// File generated at Fri 20 Oct 2017 08:36:22

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda3.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda3;

   beta_Lambda3 = Re(oneOver16PiSqr*(0.27*Quad(g1) + 2.25*Quad(g2) - 9*
      Lambda3*Sqr(g2) - 0.9*Sqr(g1)*(2*Lambda3 + Sqr(g2)) + 2*(6*Lambda1*
      Lambda3 + 6*Lambda2*Lambda3 + 2*Lambda1*Lambda4 + 2*Lambda2*Lambda4 + 8*
      Lambda6*Lambda7 + 3*Lambda3*traceYdAdjYd - 6*traceYdAdjYuYuAdjYd +
      Lambda3*traceYeAdjYe + 3*Lambda3*traceYuAdjYu + 2*Sqr(Lambda3) + Sqr(
      Lambda4) + Sqr(Lambda5) + 2*Sqr(Lambda6) + 2*Sqr(Lambda7))));


   return beta_Lambda3;
}

/**
 * Calculates the 2-loop beta function of Lambda3.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_2_loop(const Susy_traces& susy_traces) const
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


   double beta_Lambda3;

   const double beta_Lambda3_1 = Re(-0.005*twoLoop*(765*Power6(g1) - 3*
      Quad(g1)*(180*Lambda1 + 180*Lambda2 + 631*Lambda3 + 60*Lambda4 + 30*
      traceYdAdjYd - 150*traceYeAdjYe + 335*Sqr(g2)) - 5*Sqr(g1)*(65*Quad(g2) -
      6*(20*Lambda1 + 20*Lambda2 - 11*Lambda3 + 12*Lambda4 + 18*traceYdAdjYd +
      22*traceYeAdjYe)*Sqr(g2) + 2*(96*Lambda1*(3*Lambda3 + Lambda4) + 96*
      Lambda2*(3*Lambda3 + Lambda4) + 384*Lambda6*Lambda7 + 25*Lambda3*
      traceYdAdjYd - 16*traceYdAdjYuYuAdjYd + 75*Lambda3*traceYeAdjYe + 24*Sqr(
      Lambda3) - 24*Sqr(Lambda4) + 48*Sqr(Lambda5) + 24*Sqr(Lambda6) + 24*Sqr(
      Lambda7))) - 25*(259*Power6(g2) + (180*Lambda1 + 180*Lambda2 - 71*Lambda3
      + 60*Lambda4 - 18*traceYdAdjYd - 6*traceYeAdjYe)*Quad(g2) + 6*Sqr(g2)*(
      -16*Lambda3*Lambda4 + 48*Lambda1*(2*Lambda3 + Lambda4) + 48*Lambda2*(2*
      Lambda3 + Lambda4) + 144*Lambda6*Lambda7 + 15*Lambda3*traceYdAdjYd + 5*
      Lambda3*traceYeAdjYe + 8*Sqr(Lambda3) + 8*Sqr(Lambda4)) - 4*(352*Lambda3*
      Lambda6*Lambda7 + 176*Lambda4*Lambda6*Lambda7 + 144*Lambda5*Lambda6*
      Lambda7 + 96*Lambda6*Lambda7*traceYdAdjYd + 27*Lambda3*
      traceYdAdjYdYdAdjYd - 24*traceYdAdjYdYdAdjYuYuAdjYd - 30*Lambda3*
      traceYdAdjYuYuAdjYd - 48*traceYdAdjYuYuAdjYdYdAdjYd - 72*
      traceYdAdjYuYuAdjYuYuAdjYd + 32*Lambda6*Lambda7*traceYeAdjYe + 9*Lambda3*
      traceYeAdjYeYeAdjYe + 24*Cube(Lambda3) + 24*Cube(Lambda4) - 80*Lambda3*
      traceYdAdjYd*Sqr(g3) + 128*traceYdAdjYuYuAdjYd*Sqr(g3) + 8*(15*Lambda3 +
      4*Lambda4)*Sqr(Lambda1) + 8*(15*Lambda3 + 4*Lambda4)*Sqr(Lambda2) + 8*
      Lambda4*Sqr(Lambda3) + 24*traceYdAdjYd*Sqr(Lambda3) + 8*traceYeAdjYe*Sqr(
      Lambda3) + 32*Lambda3*Sqr(Lambda4) + 12*traceYdAdjYd*Sqr(Lambda4) + 4*
      traceYeAdjYe*Sqr(Lambda4) + 36*Lambda3*Sqr(Lambda5) + 88*Lambda4*Sqr(
      Lambda5) + 12*traceYdAdjYd*Sqr(Lambda5) + 4*traceYeAdjYe*Sqr(Lambda5) +
      120*Lambda3*Sqr(Lambda6) + 136*Lambda4*Sqr(Lambda6) + 136*Lambda5*Sqr(
      Lambda6) + 48*traceYdAdjYd*Sqr(Lambda6) + 16*traceYeAdjYe*Sqr(Lambda6) +
      120*Lambda3*Sqr(Lambda7) + 136*Lambda4*Sqr(Lambda7) + 136*Lambda5*Sqr(
      Lambda7) + 8*Lambda1*(22*Lambda6*Lambda7 + 6*Lambda4*traceYdAdjYd + 2*
      Lambda4*traceYeAdjYe + 2*Lambda3*(4*Lambda4 + 9*traceYdAdjYd + 3*
      traceYeAdjYe) + 18*Sqr(Lambda3) + 7*Sqr(Lambda4) + 9*Sqr(Lambda5) + 31*
      Sqr(Lambda6) + 11*Sqr(Lambda7)) + 8*Lambda2*(8*Lambda3*Lambda4 + 22*
      Lambda6*Lambda7 + 18*Sqr(Lambda3) + 7*Sqr(Lambda4) + 9*Sqr(Lambda5) + 11*
      Sqr(Lambda6) + 31*Sqr(Lambda7))))));
   const double beta_Lambda3_2 = Re(-0.01*twoLoop*(171*traceYuAdjYu*Quad(
      g1) + 5*traceYuAdjYu*Sqr(g1)*(-85*Lambda3 + 126*Sqr(g2)) + 25*(9*
      traceYuAdjYu*Quad(g2) - 45*Lambda3*traceYuAdjYu*Sqr(g2) - 160*Lambda3*
      traceYuAdjYu*Sqr(g3) + 6*(16*Lambda2*(3*Lambda3 + Lambda4)*traceYuAdjYu +
      9*Lambda3*traceYuAdjYuYuAdjYu + 8*traceYuAdjYu*Sqr(Lambda3) + 4*
      traceYuAdjYu*(4*Lambda7*(2*Lambda6 + Lambda7) + Sqr(Lambda4) + Sqr(
      Lambda5))))));

   beta_Lambda3 = beta_Lambda3_1 + beta_Lambda3_2;


   return beta_Lambda3;
}

/**
 * Calculates the 3-loop beta function of Lambda3.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

} // namespace flexiblesusy
