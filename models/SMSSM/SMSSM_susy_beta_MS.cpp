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

// File generated at Tue 22 Jan 2019 17:17:28

#include "SMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of MS.
 *
 * @return 1-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_1_loop(const Susy_traces& susy_traces) const
{


   double beta_MS;

   beta_MS = Re(4*MS*oneOver16PiSqr*(AbsSqr(Kappa) + AbsSqr(Lambdax)));


   return beta_MS;
}

/**
 * Calculates the 2-loop beta function of MS.
 *
 * @return 2-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MS;

   beta_MS = Re(-0.8*MS*twoLoop*(15*traceYdAdjYd*AbsSqr(Lambdax) + 5*
      traceYeAdjYe*AbsSqr(Lambdax) + 15*traceYuAdjYu*AbsSqr(Lambdax) + 20*
      AbsSqr(Kappa)*AbsSqr(Lambdax) - 3*AbsSqr(Lambdax)*Sqr(g1) - 15*AbsSqr(
      Lambdax)*Sqr(g2) + 20*Sqr(Conj(Kappa))*Sqr(Kappa) + 10*Sqr(Conj(Lambdax))
      *Sqr(Lambdax)));


   return beta_MS;
}

/**
 * Calculates the 3-loop beta function of MS.
 *
 * @return 3-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MS;

   beta_MS = 0;


   return beta_MS;
}

/**
 * Calculates the 4-loop beta function of MS.
 *
 * @return 4-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MS;

   beta_MS = 0;


   return beta_MS;
}

/**
 * Calculates the 5-loop beta function of MS.
 *
 * @return 5-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_MS_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MS;

   beta_MS = 0;


   return beta_MS;
}

} // namespace flexiblesusy
