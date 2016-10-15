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

// File generated at Sat 15 Oct 2016 15:25:11

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g1dp.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_g1dp;

   beta_g1dp = Re(0.05*g1dp*oneOver16PiSqr*(60*traceYdAdjYd + 20*
      traceYeAdjYe - 9*Sqr(g1) + 45*Sqr(g1d) + 25*Sqr(g1dp) - 45*Sqr(g2) + 10*
      Sqr(g2up)));


   return beta_g1dp;
}

/**
 * Calculates the two-loop beta function of g1dp.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_g1dp;

   beta_g1dp = Re(twoLoop*(0.645*Power(g1,4)*g1dp - 6.1875*Power(g1d,4)*
      g1dp - 0.75*Power(g1dp,5) - 3.75*g1dp*Power(g2,4) - 0.4375*g1dp*Power(
      g2up,4) - 6*Power(g1dp,3)*Lambda1 - 3*g1d*g2u*g2up*Lambda4 + g1dp*Lambda3
      *Lambda4 - 6.75*g1dp*traceYdAdjYdYdAdjYd - 2.25*g1dp*traceYdAdjYuYuAdjYd
      - 2.25*g1dp*traceYeAdjYeYeAdjYe + 1.93125*Power(g1dp,3)*Sqr(g1) - 0.5625*
      Power(g1dp,3)*Sqr(g1d) - 6*g1dp*Lambda1*Sqr(g1d) + 1.18125*g1dp*Sqr(g1)*
      Sqr(g1d) + 5.15625*Power(g1dp,3)*Sqr(g2) - 4.5*g1d*g2u*g2up*Sqr(g2) -
      1.35*g1dp*Sqr(g1)*Sqr(g2) + 17.15625*g1dp*Sqr(g1d)*Sqr(g2) + 0.375*g1dp*
      traceYeAdjYe*(5*Sqr(g1) - 3*Sqr(g1d) - 3*Sqr(g1dp) + 5*Sqr(g2)) - 1.3125*
      g1dp*Sqr(g1d)*Sqr(g2u) - 0.4375*Power(g1dp,3)*Sqr(g2up) - 2*g1dp*Lambda3*
      Sqr(g2up) - g1dp*Lambda4*Sqr(g2up) - 2.25*g1dp*traceYuAdjYu*Sqr(g2up) -
      0.2625*g1dp*Sqr(g1)*Sqr(g2up) + 3.1875*g1dp*Sqr(g2)*Sqr(g2up) - 1.3125*
      g1dp*Sqr(g2u)*Sqr(g2up) + 0.125*g1dp*traceYdAdjYd*(5*Sqr(g1) - 27*Sqr(g1d
      ) - 27*Sqr(g1dp) + 45*Sqr(g2) + 160*Sqr(g3)) + 6*g1dp*Sqr(Lambda1) + g1dp
      *Sqr(Lambda3) + g1dp*Sqr(Lambda4) + 1.5*g1dp*Sqr(Lambda5) + 4.5*g1dp*Sqr(
      Lambda6) + 1.5*g1dp*Sqr(Lambda7)));


   return beta_g1dp;
}

/**
 * Calculates the three-loop beta function of g1dp.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1dp;

   beta_g1dp = 0;


   return beta_g1dp;
}

} // namespace flexiblesusy
