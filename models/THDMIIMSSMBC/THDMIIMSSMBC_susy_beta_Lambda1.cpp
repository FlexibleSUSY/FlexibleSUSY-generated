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

// File generated at Tue 10 Oct 2017 21:12:23

#include "THDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda1.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(oneOver16PiSqr*(2*Lambda3*Lambda4 + 12*Lambda1*
      traceYdAdjYd - 6*traceYdAdjYdYdAdjYd + 4*Lambda1*traceYeAdjYe - 2*
      traceYeAdjYeYeAdjYe + 0.135*Quad(g1) + 1.125*Quad(g2) - 9*Lambda1*Sqr(g2)
      + 0.45*Sqr(g1)*(-4*Lambda1 + Sqr(g2)) + 24*Sqr(Lambda1) + 2*Sqr(Lambda3)
      + Sqr(Lambda4) + Sqr(Lambda5) + 12*Sqr(Lambda6)));


   return beta_Lambda1;
}

/**
 * Calculates the 2-loop beta function of Lambda1.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYdAdjYdYdAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe =
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(twoLoop*(-20*Lambda1*Lambda3*Lambda4 - 36*Lambda3*
      Lambda6*Lambda7 - 28*Lambda4*Lambda6*Lambda7 - 20*Lambda5*Lambda6*Lambda7
      - 3*Lambda1*traceYdAdjYdYdAdjYd + 30*traceYdAdjYdYdAdjYdYdAdjYd - 9*
      Lambda1*traceYdAdjYuYuAdjYd + 6*traceYdAdjYuYuAdjYdYdAdjYd - Lambda1*
      traceYeAdjYeYeAdjYe + 10*traceYeAdjYeYeAdjYeYeAdjYe - 12*Lambda3*Lambda4*
      traceYuAdjYu - 312*Cube(Lambda1) - 8*Cube(Lambda3) - 6*Cube(Lambda4) -
      1.7685*Power6(g1) + 18.1875*Power6(g2) - 0.375*(17*Lambda1 + 2*(-10*
      Lambda3 - 5*Lambda4 + 3*traceYdAdjYd + traceYeAdjYe))*Quad(g2) - 0.0225*
      Quad(g1)*(-2*(217*Lambda1 + 10*(2*Lambda3 + Lambda4 + traceYdAdjYd - 5*
      traceYeAdjYe)) + 191*Sqr(g2)) + 80*Lambda1*traceYdAdjYd*Sqr(g3) - 32*
      traceYdAdjYdYdAdjYd*Sqr(g3) - 144*traceYdAdjYd*Sqr(Lambda1) - 48*
      traceYeAdjYe*Sqr(Lambda1) - 20*Lambda1*Sqr(Lambda3) - 12*Lambda4*Sqr(
      Lambda3) - 12*traceYuAdjYu*Sqr(Lambda3) - 12*Lambda1*Sqr(Lambda4) - 16*
      Lambda3*Sqr(Lambda4) - 6*traceYuAdjYu*Sqr(Lambda4) - 14*Lambda1*Sqr(
      Lambda5) - 20*Lambda3*Sqr(Lambda5) - 22*Lambda4*Sqr(Lambda5) - 6*
      traceYuAdjYu*Sqr(Lambda5) - 318*Lambda1*Sqr(Lambda6) - 66*Lambda3*Sqr(
      Lambda6) - 70*Lambda4*Sqr(Lambda6) - 74*Lambda5*Sqr(Lambda6) - 36*
      traceYdAdjYd*Sqr(Lambda6) - 12*traceYeAdjYe*Sqr(Lambda6) - 36*
      traceYuAdjYu*Sqr(Lambda6) + 1.5*Sqr(g2)*(5*Lambda1*(3*traceYdAdjYd +
      traceYeAdjYe) + 72*Sqr(Lambda1) + 2*(4*Lambda3*Lambda4 + 4*Sqr(Lambda3) +
      Sqr(Lambda4) + 18*Sqr(Lambda6))) - 0.0125*Sqr(g1)*(303*Quad(g2) - 12*(39
      *Lambda1 + 10*Lambda4 + 18*traceYdAdjYd + 22*traceYeAdjYe)*Sqr(g2) - 8*(
      25*Lambda1*(traceYdAdjYd + 3*traceYeAdjYe) + 216*Sqr(Lambda1) + 2*(12*
      Lambda3*Lambda4 + 4*traceYdAdjYdYdAdjYd - 12*traceYeAdjYeYeAdjYe + 12*Sqr
      (Lambda3) + 6*Sqr(Lambda4) - 3*Sqr(Lambda5) + 54*Sqr(Lambda6)))) + 6*
      Lambda1*Sqr(Lambda7) - 18*Lambda3*Sqr(Lambda7) - 14*Lambda4*Sqr(Lambda7)
      - 10*Lambda5*Sqr(Lambda7)));


   return beta_Lambda1;
}

/**
 * Calculates the 3-loop beta function of Lambda1.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

} // namespace flexiblesusy
