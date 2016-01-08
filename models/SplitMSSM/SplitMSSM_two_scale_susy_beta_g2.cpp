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

// File generated at Fri 8 Jan 2016 15:06:54

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g2.
 *
 * @return one-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(-1.1666666666666667*Power(g2,3)*oneOver16PiSqr);


   return beta_g2;
}

/**
 * Calculates the two-loop beta function of g2.
 *
 * @return two-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2;

   beta_g2 = Re(0.016666666666666666*Power(g2,3)*twoLoop*(-90*
      traceYdAdjYd - 30*traceYeAdjYe - 90*traceYuAdjYu + 72*Sqr(g1) + 2120*Sqr(
      g2) - 165*Sqr(g2d) - 165*Sqr(g2u) + 720*Sqr(g3) - 15*Sqr(gYd) - 15*Sqr(
      gYu)));


   return beta_g2;
}

/**
 * Calculates the three-loop beta function of g2.
 *
 * @return three-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

} // namespace flexiblesusy