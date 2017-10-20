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

// File generated at Fri 20 Oct 2017 08:37:53

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g1.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g1_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = Re(4.5*oneOver16PiSqr*Cube(g1));


   return beta_g1;
}

/**
 * Calculates the 2-loop beta function of g1.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g1;

   beta_g1 = Re(0.01*twoLoop*Cube(g1)*(416*Sqr(g1) + 5*(-10*traceYdAdjYd
      - 30*traceYeAdjYe - 34*traceYuAdjYu + 72*Sqr(g2) - 9*Sqr(g2d) - 9*Sqr(g2u
      ) + 176*Sqr(g3) - 3*Sqr(gYd) - 3*Sqr(gYu))));


   return beta_g1;
}

/**
 * Calculates the 3-loop beta function of g1.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = 0;


   return beta_g1;
}

} // namespace flexiblesusy
