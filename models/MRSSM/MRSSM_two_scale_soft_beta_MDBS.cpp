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

// File generated at Fri 26 Jun 2015 18:59:15

#include "MRSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MDBS.
 *
 * @return one-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_MDBS_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MDBS;

   beta_MDBS = Re(0.4*MDBS*oneOver16PiSqr*(5*AbsSqr(LamSD) + 5*AbsSqr(
      LamSU) + 18*Sqr(g1)));


   return beta_MDBS;
}

/**
 * Calculates the two-loop beta function of MDBS.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_MDBS_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MDBS;

   beta_MDBS = Re(-0.04*MDBS*twoLoop*(-150*AbsSqr(LamSU)*(-traceYuAdjYu -
      AbsSqr(LamTU) + Sqr(g2)) - 50*AbsSqr(LamSD)*(-3*traceYdAdjYd -
      traceYeAdjYe - 3*AbsSqr(LamTD) + 3*Sqr(g2)) + Sqr(g1)*(70*traceYdAdjYd +
      90*traceYeAdjYe + 130*traceYuAdjYu + 45*AbsSqr(LamTD) + 45*AbsSqr(LamTU)
      - 208*Sqr(g1) - 180*Sqr(g2) - 440*Sqr(g3)) + 100*Sqr(LamSD)*Sqr(Conj(
      LamSD)) + 100*Sqr(LamSU)*Sqr(Conj(LamSU))));


   return beta_MDBS;
}

/**
 * Calculates the three-loop beta function of MDBS.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_MDBS_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDBS;

   beta_MDBS = 0;


   return beta_MDBS;
}

} // namespace flexiblesusy
