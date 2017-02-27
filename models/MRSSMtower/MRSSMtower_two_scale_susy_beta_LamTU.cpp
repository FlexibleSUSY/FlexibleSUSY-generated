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

// File generated at Mon 27 Feb 2017 13:21:56

#include "MRSSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of LamTU.
 *
 * @return one-loop beta function
 */
double MRSSMtower_susy_parameters::calc_beta_LamTU_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_LamTU;

   beta_LamTU = Re(oneOver16PiSqr*(3*LamTU*traceYuAdjYu + 2*LamTU*AbsSqr(
      LamSU) + LamTU*AbsSqr(LamTD) - 0.6*LamTU*Sqr(g1) - 7*LamTU*Sqr(g2) + 4*
      Conj(LamTU)*Sqr(LamTU)));


   return beta_LamTU;
}

/**
 * Calculates the two-loop beta function of LamTU.
 *
 * @return two-loop beta function
 */
double MRSSMtower_susy_parameters::calc_beta_LamTU_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_LamTU;

   beta_LamTU = Re(-0.1*LamTU*twoLoop*(-45*Power(g1,4) - 485*Power(g2,4)
      + 30*traceYdAdjYuYuAdjYd + 90*traceYuAdjYuYuAdjYu + 75*traceYuAdjYu*
      AbsSqr(LamTU) + 10*AbsSqr(LamSU)*(3*traceYuAdjYu + 4*AbsSqr(LamSD) + 8*
      AbsSqr(LamTU)) - 8*traceYuAdjYu*Sqr(g1) - 6*AbsSqr(LamTU)*Sqr(g1) - 110*
      AbsSqr(LamTU)*Sqr(g2) - 18*Sqr(g1)*Sqr(g2) + 2*AbsSqr(LamTD)*(15*
      traceYdAdjYd + 5*traceYeAdjYe + 10*AbsSqr(LamSD) + 15*AbsSqr(LamTU) - 3*
      Sqr(g1) + 5*Sqr(g2)) - 160*traceYuAdjYu*Sqr(g3) + 60*Sqr(LamSU)*Sqr(Conj(
      LamSU)) + 30*Sqr(LamTD)*Sqr(Conj(LamTD)) + 105*Sqr(LamTU)*Sqr(Conj(LamTU)
      )));


   return beta_LamTU;
}

/**
 * Calculates the three-loop beta function of LamTU.
 *
 * @return three-loop beta function
 */
double MRSSMtower_susy_parameters::calc_beta_LamTU_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamTU;

   beta_LamTU = 0;


   return beta_LamTU;
}

} // namespace flexiblesusy
