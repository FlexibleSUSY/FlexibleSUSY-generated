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

// File generated at Tue 10 Oct 2017 21:35:03

#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g1.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_g1_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = Re(7.2*oneOver16PiSqr*Cube(g1));


   return beta_g1;
}

/**
 * Calculates the 2-loop beta function of g1.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_g1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g1;

   beta_g1 = Re(0.04*twoLoop*Cube(g1)*(-70*traceYdAdjYd - 90*traceYeAdjYe
      - 130*traceYuAdjYu - 30*AbsSqr(LamSD) - 30*AbsSqr(LamSU) - 45*AbsSqr(
      LamTD) - 45*AbsSqr(LamTU) + 208*Sqr(g1) + 180*Sqr(g2) + 440*Sqr(g3)));


   return beta_g1;
}

/**
 * Calculates the 3-loop beta function of g1.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_g1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = 0;


   return beta_g1;
}

} // namespace flexiblesusy