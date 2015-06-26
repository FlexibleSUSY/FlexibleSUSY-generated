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

// File generated at Fri 26 Jun 2015 19:09:49

#include "SMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of MS.
 *
 * @return one-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_one_loop(const Susy_traces& susy_traces) const
{


   double beta_MS;

   beta_MS = Re(4*MS*oneOver16PiSqr*(AbsSqr(Kappa) + AbsSqr(Lambdax)));


   return beta_MS;
}

/**
 * Calculates the two-loop beta function of MS.
 *
 * @return two-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MS;

   beta_MS = Re(-0.8*MS*twoLoop*(20*AbsSqr(Kappa)*AbsSqr(Lambdax) +
      AbsSqr(Lambdax)*(15*traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu + 10*
      AbsSqr(Lambdax) - 3*Sqr(g1) - 15*Sqr(g2)) + 20*Sqr(Conj(Kappa))*Sqr(Kappa
      )));


   return beta_MS;
}

/**
 * Calculates the three-loop beta function of MS.
 *
 * @return three-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MS;

   beta_MS = 0;


   return beta_MS;
}

} // namespace flexiblesusy
