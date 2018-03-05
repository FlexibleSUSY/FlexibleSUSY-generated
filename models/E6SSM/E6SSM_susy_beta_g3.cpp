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

// File generated at Mon 5 Mar 2018 17:57:41

#include "E6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g3.
 *
 * @return 1-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_g3_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g3;

   beta_g3 = Re(0);


   return beta_g3;
}

/**
 * Calculates the 2-loop beta function of g3.
 *
 * @return 2-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_g3_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;


   double beta_g3;

   beta_g3 = Re(twoLoop*Cube(g3)*(-2*traceKappaAdjKappa - 4*traceYdAdjYd
      - 4*traceYuAdjYu + 3*Sqr(g1) + 9*Sqr(g2) + 48*Sqr(g3) + 3*Sqr(gN)));


   return beta_g3;
}

/**
 * Calculates the 3-loop beta function of g3.
 *
 * @return 3-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_g3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = 0;


   return beta_g3;
}

/**
 * Calculates the 4-loop beta function of g3.
 *
 * @return 4-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_g3_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = 0;


   return beta_g3;
}

} // namespace flexiblesusy
