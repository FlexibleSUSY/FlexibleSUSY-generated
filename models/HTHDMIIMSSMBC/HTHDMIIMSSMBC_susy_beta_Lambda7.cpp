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

// File generated at Fri 20 Oct 2017 08:36:16

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda7.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(oneOver16PiSqr*(2*Lambda5*Lambda6 + 24*Lambda2*
      Lambda7 + 10*Lambda5*Lambda7 + 6*Lambda3*(Lambda6 + Lambda7) + 4*Lambda4*
      (Lambda6 + 2*Lambda7) + 3*Lambda7*traceYdAdjYd + Lambda7*traceYeAdjYe + 9
      *Lambda7*traceYuAdjYu - 1.8*Lambda7*Sqr(g1) - 9*Lambda7*Sqr(g2)));


   return beta_Lambda7;
}

/**
 * Calculates the 2-loop beta function of Lambda7.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(twoLoop*(-36*Lambda1*Lambda3*Lambda6 - 36*Lambda2*
      Lambda3*Lambda6 - 28*Lambda1*Lambda4*Lambda6 - 28*Lambda2*Lambda4*Lambda6
      - 56*Lambda3*Lambda4*Lambda6 - 20*Lambda1*Lambda5*Lambda6 - 20*Lambda2*
      Lambda5*Lambda6 - 40*Lambda3*Lambda5*Lambda6 - 44*Lambda4*Lambda5*Lambda6
      - 36*Lambda1*Lambda3*Lambda7 - 132*Lambda2*Lambda3*Lambda7 - 28*Lambda1*
      Lambda4*Lambda7 - 140*Lambda2*Lambda4*Lambda7 - 68*Lambda3*Lambda4*
      Lambda7 - 20*Lambda1*Lambda5*Lambda7 - 148*Lambda2*Lambda5*Lambda7 - 72*
      Lambda3*Lambda5*Lambda7 - 76*Lambda4*Lambda5*Lambda7 - 36*Lambda3*Lambda6
      *traceYdAdjYd - 24*Lambda4*Lambda6*traceYdAdjYd - 12*Lambda5*Lambda6*
      traceYdAdjYd - 18*Lambda3*Lambda7*traceYdAdjYd - 24*Lambda4*Lambda7*
      traceYdAdjYd - 30*Lambda5*Lambda7*traceYdAdjYd - 6.75*Lambda7*
      traceYdAdjYdYdAdjYd - 21*Lambda7*traceYdAdjYuYuAdjYd - 12*Lambda3*Lambda6
      *traceYeAdjYe - 8*Lambda4*Lambda6*traceYeAdjYe - 4*Lambda5*Lambda6*
      traceYeAdjYe - 6*Lambda3*Lambda7*traceYeAdjYe - 8*Lambda4*Lambda7*
      traceYeAdjYe - 10*Lambda5*Lambda7*traceYeAdjYe - 2.25*Lambda7*
      traceYeAdjYeYeAdjYe - 144*Lambda2*Lambda7*traceYuAdjYu - 18*Lambda3*
      Lambda7*traceYuAdjYu - 24*Lambda4*Lambda7*traceYuAdjYu - 30*Lambda5*
      Lambda7*traceYuAdjYu - 8.25*Lambda7*traceYuAdjYuYuAdjYu - 42*Cube(Lambda6
      ) - 111*Cube(Lambda7) + 0.015*(90*Lambda6 + 601*Lambda7)*Quad(g1) + 0.125
      *(90*Lambda6 - 101*Lambda7)*Quad(g2) + 0.375*(48*Lambda3*(2*Lambda6 +
      Lambda7) + 48*Lambda4*(Lambda6 + 2*Lambda7) + Lambda7*(288*Lambda2 + 144*
      Lambda5 + 5*(3*traceYdAdjYd + traceYeAdjYe + 9*traceYuAdjYu)))*Sqr(g2) +
      0.025*Sqr(g1)*(192*Lambda4*Lambda6 - 48*Lambda5*Lambda6 + 864*Lambda2*
      Lambda7 + 240*Lambda4*Lambda7 + 480*Lambda5*Lambda7 + 144*Lambda3*(2*
      Lambda6 + Lambda7) + 25*Lambda7*traceYdAdjYd + 75*Lambda7*traceYeAdjYe +
      255*Lambda7*traceYuAdjYu + 6*(10*Lambda6 + 29*Lambda7)*Sqr(g2)) + 20*
      Lambda7*traceYdAdjYd*Sqr(g3) + 60*Lambda7*traceYuAdjYu*Sqr(g3) + 6*
      Lambda7*Sqr(Lambda1) - 318*Lambda7*Sqr(Lambda2) - 36*Lambda6*Sqr(Lambda3)
      - 32*Lambda7*Sqr(Lambda3) - 34*Lambda6*Sqr(Lambda4) - 34*Lambda7*Sqr(
      Lambda4) - 42*Lambda6*Sqr(Lambda5) - 36*Lambda7*Sqr(Lambda5) - 33*Lambda7
      *Sqr(Lambda6) - 126*Lambda6*Sqr(Lambda7)));


   return beta_Lambda7;
}

/**
 * Calculates the 3-loop beta function of Lambda7.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

} // namespace flexiblesusy
