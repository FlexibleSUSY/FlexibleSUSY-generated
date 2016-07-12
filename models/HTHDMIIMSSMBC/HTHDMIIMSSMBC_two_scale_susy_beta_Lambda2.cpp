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

// File generated at Tue 12 Jul 2016 10:32:40

#include "HTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda2.
 *
 * @return one-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(oneOver16PiSqr*(0.135*Power(g1,4) + 1.125*Power(g2,4
      ) + 2*Lambda3*Lambda4 + 12*Lambda2*traceYuAdjYu - 6*traceYuAdjYuYuAdjYu -
      1.8*Lambda2*Sqr(g1) - 9*Lambda2*Sqr(g2) + 0.45*Sqr(g1)*Sqr(g2) + 24*Sqr(
      Lambda2) + 2*Sqr(Lambda3) + Sqr(Lambda4) + Sqr(Lambda5) + 12*Sqr(Lambda7)
      ));


   return beta_Lambda2;
}

/**
 * Calculates the two-loop beta function of Lambda2.
 *
 * @return two-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYuYuAdjYu =
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(twoLoop*(-1.9125*Power(g1,6) + 16.1875*Power(g2,6) +
      10.365*Power(g1,4)*Lambda2 - 1.375*Power(g2,4)*Lambda2 - 312*Power(
      Lambda2,3) + 0.9*Power(g1,4)*Lambda3 + 7.5*Power(g2,4)*Lambda3 - 8*Power(
      Lambda3,3) + 0.45*Power(g1,4)*Lambda4 + 3.75*Power(g2,4)*Lambda4 - 20*
      Lambda2*Lambda3*Lambda4 - 6*Power(Lambda4,3) - 36*Lambda3*Lambda6*Lambda7
      - 28*Lambda4*Lambda6*Lambda7 - 20*Lambda5*Lambda6*Lambda7 - 9*Lambda2*
      traceYdAdjYuYuAdjYd + 6*traceYdAdjYuYuAdjYuYuAdjYd - 1.71*Power(g1,4)*
      traceYuAdjYu - 2.25*Power(g2,4)*traceYuAdjYu - 3*Lambda2*
      traceYuAdjYuYuAdjYu + 30*traceYuAdjYuYuAdjYuYuAdjYu - 4.1875*Power(g2,4)*
      Sqr(g1) + 2.4*Lambda3*Lambda4*Sqr(g1) + 8.5*Lambda2*traceYuAdjYu*Sqr(g1)
      - 1.6*traceYuAdjYuYuAdjYu*Sqr(g1) - 4.5375*Power(g1,4)*Sqr(g2) + 12*
      Lambda3*Lambda4*Sqr(g2) + 22.5*Lambda2*traceYuAdjYu*Sqr(g2) + 5.85*
      Lambda2*Sqr(g1)*Sqr(g2) + 1.5*Lambda4*Sqr(g1)*Sqr(g2) + 6.3*traceYuAdjYu*
      Sqr(g1)*Sqr(g2) + 80*Lambda2*traceYuAdjYu*Sqr(g3) - 32*
      traceYuAdjYuYuAdjYu*Sqr(g3) - 144*traceYuAdjYu*Sqr(Lambda2) + 21.6*Sqr(g1
      )*Sqr(Lambda2) + 108*Sqr(g2)*Sqr(Lambda2) - 20*Lambda2*Sqr(Lambda3) - 12*
      Lambda4*Sqr(Lambda3) + 2.4*Sqr(g1)*Sqr(Lambda3) + 12*Sqr(g2)*Sqr(Lambda3)
      - 12*Lambda2*Sqr(Lambda4) - 16*Lambda3*Sqr(Lambda4) + 1.2*Sqr(g1)*Sqr(
      Lambda4) + 3*Sqr(g2)*Sqr(Lambda4) - 14*Lambda2*Sqr(Lambda5) - 20*Lambda3*
      Sqr(Lambda5) - 22*Lambda4*Sqr(Lambda5) - 0.6*Sqr(g1)*Sqr(Lambda5) + 6*
      Lambda2*Sqr(Lambda6) - 18*Lambda3*Sqr(Lambda6) - 14*Lambda4*Sqr(Lambda6)
      - 10*Lambda5*Sqr(Lambda6) - 318*Lambda2*Sqr(Lambda7) - 66*Lambda3*Sqr(
      Lambda7) - 70*Lambda4*Sqr(Lambda7) - 74*Lambda5*Sqr(Lambda7) - 36*
      traceYuAdjYu*Sqr(Lambda7) + 10.8*Sqr(g1)*Sqr(Lambda7) + 54*Sqr(g2)*Sqr(
      Lambda7) - 6*traceYdAdjYd*(2*Lambda3*Lambda4 + 2*Sqr(Lambda3) + Sqr(
      Lambda4) + Sqr(Lambda5) + 6*Sqr(Lambda7)) - 2*traceYeAdjYe*(2*Lambda3*
      Lambda4 + 2*Sqr(Lambda3) + Sqr(Lambda4) + Sqr(Lambda5) + 6*Sqr(Lambda7)))
      );


   return beta_Lambda2;
}

/**
 * Calculates the three-loop beta function of Lambda2.
 *
 * @return three-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

} // namespace flexiblesusy
