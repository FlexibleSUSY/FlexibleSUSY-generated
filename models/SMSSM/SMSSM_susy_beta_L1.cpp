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

// File generated at Sun 26 Aug 2018 14:31:01

#include "SMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of L1.
 *
 * @return 1-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_L1_1_loop(const Susy_traces& susy_traces) const
{


   double beta_L1;

   beta_L1 = Re(2*L1*oneOver16PiSqr*(AbsSqr(Kappa) + AbsSqr(Lambdax)));


   return beta_L1;
}

/**
 * Calculates the 2-loop beta function of L1.
 *
 * @return 2-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_L1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_L1;

   beta_L1 = Re(-0.4*L1*twoLoop*(20*AbsSqr(Kappa)*AbsSqr(Lambdax) + AbsSqr(
      Lambdax)*(10*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu - 3*Sqr(g2))) + 20*Sqr(Conj(Kappa))*Sqr(
      Kappa)));


   return beta_L1;
}

/**
 * Calculates the 3-loop beta function of L1.
 *
 * @return 3-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_L1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_L1;

   beta_L1 = 0;


   return beta_L1;
}

/**
 * Calculates the 4-loop beta function of L1.
 *
 * @return 4-loop beta function
 */
double SMSSM_susy_parameters::calc_beta_L1_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_L1;

   beta_L1 = 0;


   return beta_L1;
}

} // namespace flexiblesusy
