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

// File generated at Tue 5 Sep 2017 12:33:49

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

   beta_Mu = Re(oneOver16PiSqr*(-0.6*Mu*Sqr(g1) + Mu*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + traceYvAdjYv - 3*Sqr(g2))));


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

   beta_Mu = Re(0.02*twoLoop*Mu*(207*Power(g1,4) + 10*Sqr(g1)*(-2*
      traceYdAdjYd + 6*traceYeAdjYe + 4*traceYuAdjYu + 9*Sqr(g2)) + 25*(15*
      Power(g2,4) + 2*(-9*traceYdAdjYdYdAdjYd - 6*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - 2*traceYeAdjYvYvAdjYe - 9*traceYuAdjYuYuAdjYu - 3*
      traceYvAdjYvYvAdjYv + 16*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3)))));


   return beta_Mu;
}

/**
 * Calculates the three-loop beta function of Mu.
 *
 * @return three-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_Mu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
