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

// File generated at Mon 5 Mar 2018 18:45:56

#include "MSSMRHN_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Mu.
 *
 * @return 1-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_Mu_1_loop(const Susy_traces& susy_traces) const
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
 * Calculates the 2-loop beta function of Mu.
 *
 * @return 2-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_Mu_2_loop(const Susy_traces& susy_traces) const
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

   beta_Mu = Re(0.02*twoLoop*Mu*(207*Quad(g1) + 10*Sqr(g1)*(-2*
      traceYdAdjYd + 6*traceYeAdjYe + 4*traceYuAdjYu + 9*Sqr(g2)) + 25*(15*Quad
      (g2) + 2*(-9*traceYdAdjYdYdAdjYd - 6*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - 2*traceYeAdjYvYvAdjYe - 9*traceYuAdjYuYuAdjYu - 3*
      traceYvAdjYvYvAdjYv + 16*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3)))));


   return beta_Mu;
}

/**
 * Calculates the 3-loop beta function of Mu.
 *
 * @return 3-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_Mu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

/**
 * Calculates the 4-loop beta function of Mu.
 *
 * @return 4-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_Mu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
