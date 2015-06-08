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

// File generated at Mon 8 Jun 2015 17:42:23

#include "TMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vT.
 *
 * @return one-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_vT_one_loop(const Susy_traces& susy_traces) const
{


   double beta_vT;

   beta_vT = Re(oneOver16PiSqr*vT*(-AbsSqr(Lambdax) + 4*Sqr(g2)));


   return beta_vT;
}

/**
 * Calculates the two-loop beta function of vT.
 *
 * @return two-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_vT_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vT;

   beta_vT = Re(twoLoop*(-12.666666666666666*Power(g2,4)*vT + 0.2*vT*
      AbsSqr(Lambdax)*(15*traceYdAdjYd + 5*traceYeAdjYe - 3*(-5*traceYuAdjYu +
      Sqr(g1) + 5*Sqr(g2))) + 3*vT*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vT;
}

} // namespace flexiblesusy
