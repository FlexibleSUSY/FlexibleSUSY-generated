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

// File generated at Fri 8 Jan 2016 12:05:05

#include "MRSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of MuD.
 *
 * @return one-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_MuD;

   beta_MuD = Re(oneOver16PiSqr*(3*MuD*traceYdAdjYd + MuD*traceYeAdjYe +
      2*MuD*AbsSqr(LamSD) + 3*MuD*AbsSqr(LamTD) - 0.6*MuD*Sqr(g1) - 3*MuD*Sqr(
      g2)));


   return beta_MuD;
}

/**
 * Calculates the two-loop beta function of MuD.
 *
 * @return two-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_MuD;

   beta_MuD = Re(0.1*MuD*twoLoop*(45*Power(g1,4) + 165*Power(g2,4) - 90*
      traceYdAdjYdYdAdjYd - 30*traceYdAdjYuYuAdjYd - 30*traceYeAdjYeYeAdjYe -
      10*AbsSqr(LamSD)*(3*traceYdAdjYd + traceYeAdjYe + 4*AbsSqr(LamSU) + 6*
      AbsSqr(LamTD)) - 4*traceYdAdjYd*Sqr(g1) + 12*traceYeAdjYe*Sqr(g1) + 18*
      Sqr(g1)*Sqr(g2) + 15*AbsSqr(LamTD)*(-3*traceYdAdjYd - traceYeAdjYe - 2*
      AbsSqr(LamTU) + 8*Sqr(g2)) + 160*traceYdAdjYd*Sqr(g3) - 60*Sqr(LamSD)*Sqr
      (Conj(LamSD)) - 75*Sqr(LamTD)*Sqr(Conj(LamTD))));


   return beta_MuD;
}

/**
 * Calculates the three-loop beta function of MuD.
 *
 * @return three-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuD;

   beta_MuD = 0;


   return beta_MuD;
}

} // namespace flexiblesusy
