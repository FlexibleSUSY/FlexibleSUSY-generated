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

// File generated at Fri 8 Jan 2016 15:08:10

#include "TMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of MT.
 *
 * @return one-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_one_loop(const Susy_traces& susy_traces) const
{


   double beta_MT;

   beta_MT = Re(oneOver16PiSqr*(2*MT*AbsSqr(Lambdax) - 8*MT*Sqr(g2)));


   return beta_MT;
}

/**
 * Calculates the two-loop beta function of MT.
 *
 * @return two-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MT;

   beta_MT = Re(0.4*MT*twoLoop*(140*Power(g2,4) + AbsSqr(Lambdax)*(-15*
      traceYdAdjYd - 5*traceYeAdjYe - 15*traceYuAdjYu + 3*Sqr(g1) - 5*Sqr(g2))
      - 15*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_MT;
}

/**
 * Calculates the three-loop beta function of MT.
 *
 * @return three-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MT;

   beta_MT = 0;


   return beta_MT;
}

} // namespace flexiblesusy
