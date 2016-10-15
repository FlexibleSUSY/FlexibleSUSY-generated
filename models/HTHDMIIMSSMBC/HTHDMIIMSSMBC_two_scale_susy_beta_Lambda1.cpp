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

// File generated at Sat 15 Oct 2016 15:19:33

#include "HTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda1.
 *
 * @return one-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(oneOver16PiSqr*(0.135*Power(g1,4) + 1.125*Power(g2,4
      ) + 2*Lambda3*Lambda4 + 12*Lambda1*traceYdAdjYd - 6*traceYdAdjYdYdAdjYd +
      4*Lambda1*traceYeAdjYe - 2*traceYeAdjYeYeAdjYe - 1.8*Lambda1*Sqr(g1) - 9
      *Lambda1*Sqr(g2) + 0.45*Sqr(g1)*Sqr(g2) + 24*Sqr(Lambda1) + 2*Sqr(Lambda3
      ) + Sqr(Lambda4) + Sqr(Lambda5) + 12*Sqr(Lambda6)));


   return beta_Lambda1;
}

/**
 * Calculates the two-loop beta function of Lambda1.
 *
 * @return two-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_two_loop(const Susy_traces& susy_traces) const
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

   beta_Lambda1 = Re(twoLoop*(-1.9125*Power(g1,6) + 16.1875*Power(g2,6) +
      10.365*Power(g1,4)*Lambda1 - 1.375*Power(g2,4)*Lambda1 - 312*Power(
      Lambda1,3) + 0.9*Power(g1,4)*Lambda3 + 7.5*Power(g2,4)*Lambda3 - 8*Power(
      Lambda3,3) + 0.45*Power(g1,4)*Lambda4 + 3.75*Power(g2,4)*Lambda4 - 20*
      Lambda1*Lambda3*Lambda4 - 6*Power(Lambda4,3) - 36*Lambda3*Lambda6*Lambda7
      - 28*Lambda4*Lambda6*Lambda7 - 20*Lambda5*Lambda6*Lambda7 - 3*Lambda1*
      traceYdAdjYdYdAdjYd + 30*traceYdAdjYdYdAdjYdYdAdjYd - 9*Lambda1*
      traceYdAdjYuYuAdjYd + 6*traceYdAdjYuYuAdjYdYdAdjYd - Lambda1*
      traceYeAdjYeYeAdjYe + 10*traceYeAdjYeYeAdjYeYeAdjYe - 12*Lambda3*Lambda4*
      traceYuAdjYu - 4.1875*Power(g2,4)*Sqr(g1) + 2.4*Lambda3*Lambda4*Sqr(g1) +
      0.8*traceYdAdjYdYdAdjYd*Sqr(g1) - 2.4*traceYeAdjYeYeAdjYe*Sqr(g1) -
      4.5375*Power(g1,4)*Sqr(g2) + 12*Lambda3*Lambda4*Sqr(g2) + 5.85*Lambda1*
      Sqr(g1)*Sqr(g2) + 1.5*Lambda4*Sqr(g1)*Sqr(g2) - 32*traceYdAdjYdYdAdjYd*
      Sqr(g3) + 21.6*Sqr(g1)*Sqr(Lambda1) + 108*Sqr(g2)*Sqr(Lambda1) - 20*
      Lambda1*Sqr(Lambda3) - 12*Lambda4*Sqr(Lambda3) - 12*traceYuAdjYu*Sqr(
      Lambda3) + 2.4*Sqr(g1)*Sqr(Lambda3) + 12*Sqr(g2)*Sqr(Lambda3) - 12*
      Lambda1*Sqr(Lambda4) - 16*Lambda3*Sqr(Lambda4) - 6*traceYuAdjYu*Sqr(
      Lambda4) + 1.2*Sqr(g1)*Sqr(Lambda4) + 3*Sqr(g2)*Sqr(Lambda4) - 14*Lambda1
      *Sqr(Lambda5) - 20*Lambda3*Sqr(Lambda5) - 22*Lambda4*Sqr(Lambda5) - 6*
      traceYuAdjYu*Sqr(Lambda5) - 0.6*Sqr(g1)*Sqr(Lambda5) - 318*Lambda1*Sqr(
      Lambda6) - 66*Lambda3*Sqr(Lambda6) - 70*Lambda4*Sqr(Lambda6) - 74*Lambda5
      *Sqr(Lambda6) - 36*traceYuAdjYu*Sqr(Lambda6) + 10.8*Sqr(g1)*Sqr(Lambda6)
      + 54*Sqr(g2)*Sqr(Lambda6) - 0.15*traceYeAdjYe*(15*Power(g1,4) - 2*Sqr(g1)
      *(25*Lambda1 + 11*Sqr(g2)) + 5*(Power(g2,4) - 10*Lambda1*Sqr(g2) + 64*Sqr
      (Lambda1) + 16*Sqr(Lambda6))) + 0.05*traceYdAdjYd*(9*Power(g1,4) + Sqr(g1
      )*(50*Lambda1 + 54*Sqr(g2)) - 5*(9*Power(g2,4) - 90*Lambda1*Sqr(g2) - 320
      *Lambda1*Sqr(g3) + 576*Sqr(Lambda1) + 144*Sqr(Lambda6))) + 6*Lambda1*Sqr(
      Lambda7) - 18*Lambda3*Sqr(Lambda7) - 14*Lambda4*Sqr(Lambda7) - 10*Lambda5
      *Sqr(Lambda7)));


   return beta_Lambda1;
}

/**
 * Calculates the three-loop beta function of Lambda1.
 *
 * @return three-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

} // namespace flexiblesusy
