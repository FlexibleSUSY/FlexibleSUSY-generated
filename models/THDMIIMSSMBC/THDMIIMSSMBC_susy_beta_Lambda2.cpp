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

// File generated at Tue 10 Oct 2017 21:12:30

#include "THDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda2.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(oneOver16PiSqr*(2*Lambda3*Lambda4 + 12*Lambda2*
      traceYuAdjYu - 6*traceYuAdjYuYuAdjYu + 0.135*Quad(g1) + 1.125*Quad(g2) -
      9*Lambda2*Sqr(g2) + 0.45*Sqr(g1)*(-4*Lambda2 + Sqr(g2)) + 24*Sqr(Lambda2)
      + 2*Sqr(Lambda3) + Sqr(Lambda4) + Sqr(Lambda5) + 12*Sqr(Lambda7)));


   return beta_Lambda2;
}

/**
 * Calculates the 2-loop beta function of Lambda2.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_2_loop(const Susy_traces& susy_traces) const
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

   beta_Lambda2 = Re(twoLoop*(-20*Lambda2*Lambda3*Lambda4 - 36*Lambda3*
      Lambda6*Lambda7 - 28*Lambda4*Lambda6*Lambda7 - 20*Lambda5*Lambda6*Lambda7
      - 12*Lambda3*Lambda4*traceYdAdjYd - 9*Lambda2*traceYdAdjYuYuAdjYd + 6*
      traceYdAdjYuYuAdjYuYuAdjYd - 4*Lambda3*Lambda4*traceYeAdjYe - 3*Lambda2*
      traceYuAdjYuYuAdjYu + 30*traceYuAdjYuYuAdjYuYuAdjYu - 312*Cube(Lambda2) -
      8*Cube(Lambda3) - 6*Cube(Lambda4) - 1.7685*Power6(g1) + 18.1875*Power6(
      g2) - 0.375*(17*Lambda2 - 20*Lambda3 - 10*Lambda4 + 6*traceYuAdjYu)*Quad(
      g2) - 0.0225*Quad(g1)*(-434*Lambda2 - 40*Lambda3 - 20*Lambda4 + 76*
      traceYuAdjYu + 191*Sqr(g2)) + 80*Lambda2*traceYuAdjYu*Sqr(g3) - 32*
      traceYuAdjYuYuAdjYu*Sqr(g3) - 144*traceYuAdjYu*Sqr(Lambda2) - 20*Lambda2*
      Sqr(Lambda3) - 12*Lambda4*Sqr(Lambda3) - 12*traceYdAdjYd*Sqr(Lambda3) - 4
      *traceYeAdjYe*Sqr(Lambda3) - 12*Lambda2*Sqr(Lambda4) - 16*Lambda3*Sqr(
      Lambda4) - 6*traceYdAdjYd*Sqr(Lambda4) - 2*traceYeAdjYe*Sqr(Lambda4) - 14
      *Lambda2*Sqr(Lambda5) - 20*Lambda3*Sqr(Lambda5) - 22*Lambda4*Sqr(Lambda5)
      - 6*traceYdAdjYd*Sqr(Lambda5) - 2*traceYeAdjYe*Sqr(Lambda5) + 6*Lambda2*
      Sqr(Lambda6) - 18*Lambda3*Sqr(Lambda6) - 14*Lambda4*Sqr(Lambda6) - 10*
      Lambda5*Sqr(Lambda6) - 318*Lambda2*Sqr(Lambda7) - 66*Lambda3*Sqr(Lambda7)
      - 70*Lambda4*Sqr(Lambda7) - 74*Lambda5*Sqr(Lambda7) - 36*traceYdAdjYd*
      Sqr(Lambda7) - 12*traceYeAdjYe*Sqr(Lambda7) - 36*traceYuAdjYu*Sqr(Lambda7
      ) + 1.5*Sqr(g2)*(8*Lambda3*Lambda4 + 15*Lambda2*traceYuAdjYu + 72*Sqr(
      Lambda2) + 8*Sqr(Lambda3) + 2*Sqr(Lambda4) + 36*Sqr(Lambda7)) - 0.0125*
      Sqr(g1)*(303*Quad(g2) - 12*(39*Lambda2 + 10*Lambda4 + 42*traceYuAdjYu)*
      Sqr(g2) - 8*(24*Lambda3*Lambda4 + 85*Lambda2*traceYuAdjYu - 16*
      traceYuAdjYuYuAdjYu + 216*Sqr(Lambda2) + 24*Sqr(Lambda3) + 12*Sqr(Lambda4
      ) - 6*Sqr(Lambda5) + 108*Sqr(Lambda7)))));


   return beta_Lambda2;
}

/**
 * Calculates the 3-loop beta function of Lambda2.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

} // namespace flexiblesusy