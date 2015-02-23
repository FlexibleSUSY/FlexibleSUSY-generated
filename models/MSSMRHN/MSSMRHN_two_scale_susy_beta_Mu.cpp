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

// File generated at Mon 23 Feb 2015 13:39:24

#include "MSSMRHN_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Mu.
 *
 * @return one-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_Mu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_Mu;

   beta_Mu = oneOver16PiSqr*(3*traceYdAdjYd*Mu + traceYeAdjYe*Mu + 3*
      traceYuAdjYu*Mu + traceYvAdjYv*Mu - 0.6*Mu*Sqr(g1) - 3*Mu*Sqr(g2));


   return beta_Mu;
}

/**
 * Calculates the two-loop beta function of Mu.
 *
 * @return two-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_Mu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;


   double beta_Mu;

   beta_Mu = 0.02*twoLoop*Mu*(207*Power(g1,4) + 375*Power(g2,4) - 450*
      traceYdAdjYdYdAdjYd - 300*traceYdAdjYuYuAdjYd - 150*traceYeAdjYeYeAdjYe -
      100*traceYeAdjYvYvAdjYe - 450*traceYuAdjYuYuAdjYu - 150*
      traceYvAdjYvYvAdjYv + 60*traceYeAdjYe*Sqr(g1) + 40*traceYuAdjYu*Sqr(g1) +
      90*Sqr(g1)*Sqr(g2) - 20*traceYdAdjYd*(Sqr(g1) - 40*Sqr(g3)) + 800*
      traceYuAdjYu*Sqr(g3));


   return beta_Mu;
}

} // namespace flexiblesusy
