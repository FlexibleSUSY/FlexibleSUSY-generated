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

// File generated at Tue 5 Sep 2017 10:18:00

#include "HTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda6.
 *
 * @return one-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(oneOver16PiSqr*(24*Lambda1*Lambda6 + 6*Lambda3*
      Lambda6 + 8*Lambda4*Lambda6 + 10*Lambda5*Lambda6 + 6*Lambda3*Lambda7 + 4*
      Lambda4*Lambda7 + 2*Lambda5*Lambda7 + 9*Lambda6*traceYdAdjYd + 3*Lambda6*
      traceYeAdjYe + 3*Lambda6*traceYuAdjYu - 1.8*Lambda6*Sqr(g1) - 9*Lambda6*
      Sqr(g2)));


   return beta_Lambda6;
}

/**
 * Calculates the two-loop beta function of Lambda6.
 *
 * @return two-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(twoLoop*(-132*Lambda1*Lambda3*Lambda6 - 36*Lambda2*
      Lambda3*Lambda6 - 140*Lambda1*Lambda4*Lambda6 - 28*Lambda2*Lambda4*
      Lambda6 - 68*Lambda3*Lambda4*Lambda6 - 148*Lambda1*Lambda5*Lambda6 - 20*
      Lambda2*Lambda5*Lambda6 - 72*Lambda3*Lambda5*Lambda6 - 76*Lambda4*Lambda5
      *Lambda6 - 111*Power(Lambda6,3) - 0.125*Power(g2,4)*(101*Lambda6 - 90*
      Lambda7) - 36*Lambda1*Lambda3*Lambda7 - 36*Lambda2*Lambda3*Lambda7 - 28*
      Lambda1*Lambda4*Lambda7 - 28*Lambda2*Lambda4*Lambda7 - 56*Lambda3*Lambda4
      *Lambda7 - 20*Lambda1*Lambda5*Lambda7 - 20*Lambda2*Lambda5*Lambda7 - 40*
      Lambda3*Lambda5*Lambda7 - 44*Lambda4*Lambda5*Lambda7 - 42*Power(Lambda7,3
      ) + 0.015*Power(g1,4)*(601*Lambda6 + 90*Lambda7) - 144*Lambda1*Lambda6*
      traceYdAdjYd - 18*Lambda3*Lambda6*traceYdAdjYd - 24*Lambda4*Lambda6*
      traceYdAdjYd - 30*Lambda5*Lambda6*traceYdAdjYd - 8.25*Lambda6*
      traceYdAdjYdYdAdjYd - 21*Lambda6*traceYdAdjYuYuAdjYd - 48*Lambda1*Lambda6
      *traceYeAdjYe - 6*Lambda3*Lambda6*traceYeAdjYe - 8*Lambda4*Lambda6*
      traceYeAdjYe - 10*Lambda5*Lambda6*traceYeAdjYe - 2.75*Lambda6*
      traceYeAdjYeYeAdjYe - 18*Lambda3*Lambda6*traceYuAdjYu - 24*Lambda4*
      Lambda6*traceYuAdjYu - 30*Lambda5*Lambda6*traceYuAdjYu - 36*Lambda3*
      Lambda7*traceYuAdjYu - 24*Lambda4*Lambda7*traceYuAdjYu - 12*Lambda5*
      Lambda7*traceYuAdjYu - 6.75*Lambda6*traceYuAdjYuYuAdjYu + 1.125*(96*
      Lambda1*Lambda6 + 32*Lambda4*Lambda6 + 48*Lambda5*Lambda6 + 16*Lambda4*
      Lambda7 + 16*Lambda3*(Lambda6 + 2*Lambda7) + 15*Lambda6*traceYdAdjYd + 5*
      Lambda6*traceYeAdjYe + 5*Lambda6*traceYuAdjYu)*Sqr(g2) + 0.025*Sqr(g1)*(
      864*Lambda1*Lambda6 + 144*Lambda3*Lambda6 + 240*Lambda4*Lambda6 + 480*
      Lambda5*Lambda6 + 288*Lambda3*Lambda7 + 192*Lambda4*Lambda7 - 48*Lambda5*
      Lambda7 + 75*Lambda6*traceYdAdjYd + 225*Lambda6*traceYeAdjYe + 85*Lambda6
      *traceYuAdjYu + 6*(29*Lambda6 + 10*Lambda7)*Sqr(g2)) + 60*Lambda6*
      traceYdAdjYd*Sqr(g3) + 20*Lambda6*traceYuAdjYu*Sqr(g3) - 318*Lambda6*Sqr(
      Lambda1) + 6*Lambda6*Sqr(Lambda2) - 32*Lambda6*Sqr(Lambda3) - 36*Lambda7*
      Sqr(Lambda3) - 34*Lambda6*Sqr(Lambda4) - 34*Lambda7*Sqr(Lambda4) - 36*
      Lambda6*Sqr(Lambda5) - 42*Lambda7*Sqr(Lambda5) - 126*Lambda7*Sqr(Lambda6)
      - 33*Lambda6*Sqr(Lambda7)));


   return beta_Lambda6;
}

/**
 * Calculates the three-loop beta function of Lambda6.
 *
 * @return three-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

} // namespace flexiblesusy
