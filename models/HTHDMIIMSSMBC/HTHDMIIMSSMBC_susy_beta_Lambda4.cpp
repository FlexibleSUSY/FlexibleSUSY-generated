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

// File generated at Fri 20 Oct 2017 08:36:19

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda4.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda4;

   beta_Lambda4 = Re(oneOver16PiSqr*(-9*Lambda4*Sqr(g2) + 1.8*Sqr(g1)*(
      -Lambda4 + Sqr(g2)) + 2*(2*Lambda1*Lambda4 + 2*Lambda2*Lambda4 + 4*
      Lambda3*Lambda4 + 2*Lambda6*Lambda7 + 3*Lambda4*traceYdAdjYd + 6*
      traceYdAdjYuYuAdjYd + Lambda4*traceYeAdjYe + 3*Lambda4*traceYuAdjYu + 2*
      Sqr(Lambda4) + 4*Sqr(Lambda5) + 5*Sqr(Lambda6) + 5*Sqr(Lambda7))));


   return beta_Lambda4;
}

/**
 * Calculates the 2-loop beta function of Lambda4.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_2_loop(const Susy_traces& susy_traces) const
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

   beta_Lambda4 = Re(twoLoop*(-80*Lambda1*Lambda3*Lambda4 - 80*Lambda2*
      Lambda3*Lambda4 - 40*Lambda1*Lambda6*Lambda7 - 40*Lambda2*Lambda6*Lambda7
      - 80*Lambda3*Lambda6*Lambda7 - 160*Lambda4*Lambda6*Lambda7 - 96*Lambda5*
      Lambda6*Lambda7 - 24*Lambda1*Lambda4*traceYdAdjYd - 24*Lambda3*Lambda4*
      traceYdAdjYd - 12*Lambda6*Lambda7*traceYdAdjYd - 13.5*Lambda4*
      traceYdAdjYdYdAdjYd - 12*traceYdAdjYdYdAdjYuYuAdjYd - 24*Lambda3*
      traceYdAdjYuYuAdjYd - 33*Lambda4*traceYdAdjYuYuAdjYd - 12*
      traceYdAdjYuYuAdjYdYdAdjYd - 24*traceYdAdjYuYuAdjYuYuAdjYd - 8*Lambda1*
      Lambda4*traceYeAdjYe - 8*Lambda3*Lambda4*traceYeAdjYe - 4*Lambda6*Lambda7
      *traceYeAdjYe - 4.5*Lambda4*traceYeAdjYeYeAdjYe - 24*Lambda2*Lambda4*
      traceYuAdjYu - 24*Lambda3*Lambda4*traceYuAdjYu - 12*Lambda6*Lambda7*
      traceYuAdjYu - 13.5*Lambda4*traceYuAdjYuYuAdjYu - 23.875*Lambda4*Quad(g2)
      - 0.015*Quad(g1)*(-511*Lambda4 + 940*Sqr(g2)) + 40*Lambda4*traceYdAdjYd*
      Sqr(g3) + 64*traceYdAdjYuYuAdjYd*Sqr(g3) + 40*Lambda4*traceYuAdjYu*Sqr(g3
      ) - 28*Lambda4*Sqr(Lambda1) - 28*Lambda4*Sqr(Lambda2) - 28*Lambda4*Sqr(
      Lambda3) - 40*Lambda1*Sqr(Lambda4) - 40*Lambda2*Sqr(Lambda4) - 28*Lambda3
      *Sqr(Lambda4) - 12*traceYdAdjYd*Sqr(Lambda4) - 4*traceYeAdjYe*Sqr(Lambda4
      ) - 12*traceYuAdjYu*Sqr(Lambda4) - 48*Lambda1*Sqr(Lambda5) - 48*Lambda2*
      Sqr(Lambda5) - 48*Lambda3*Sqr(Lambda5) - 26*Lambda4*Sqr(Lambda5) - 24*
      traceYdAdjYd*Sqr(Lambda5) - 8*traceYeAdjYe*Sqr(Lambda5) - 24*traceYuAdjYu
      *Sqr(Lambda5) - 148*Lambda1*Sqr(Lambda6) - 20*Lambda2*Sqr(Lambda6) - 72*
      Lambda3*Sqr(Lambda6) - 68*Lambda4*Sqr(Lambda6) - 80*Lambda5*Sqr(Lambda6)
      - 60*traceYdAdjYd*Sqr(Lambda6) - 20*traceYeAdjYe*Sqr(Lambda6) - 0.05*Sqr(
      g1)*(-96*Lambda1*Lambda4 - 96*Lambda2*Lambda4 - 48*Lambda3*Lambda4 - 96*
      Lambda6*Lambda7 - 25*Lambda4*traceYdAdjYd - 16*traceYdAdjYuYuAdjYd - 75*
      Lambda4*traceYeAdjYe - 85*Lambda4*traceYuAdjYu + 200*Quad(g2) - 3*(40*
      Lambda1 + 40*Lambda2 + 8*Lambda3 + 51*Lambda4 + 36*traceYdAdjYd + 44*
      traceYeAdjYe + 84*traceYuAdjYu)*Sqr(g2) - 96*Sqr(Lambda4) - 192*Sqr(
      Lambda5) - 168*Sqr(Lambda6) - 168*Sqr(Lambda7)) - 20*Lambda1*Sqr(Lambda7)
      - 148*Lambda2*Sqr(Lambda7) - 72*Lambda3*Sqr(Lambda7) - 68*Lambda4*Sqr(
      Lambda7) - 80*Lambda5*Sqr(Lambda7) - 60*traceYuAdjYu*Sqr(Lambda7) + 0.75*
      Sqr(g2)*(48*Lambda3*Lambda4 + 5*Lambda4*(3*traceYdAdjYd + traceYeAdjYe +
      3*traceYuAdjYu) + 24*Sqr(Lambda4) + 72*(Sqr(Lambda5) + Sqr(Lambda6) + Sqr
      (Lambda7)))));


   return beta_Lambda4;
}

/**
 * Calculates the 3-loop beta function of Lambda4.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return beta_Lambda4;
}

} // namespace flexiblesusy
