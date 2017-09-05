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

// File generated at Tue 5 Sep 2017 12:33:51

#include "MSSMRHN_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vu.
 *
 * @return one-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_vu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_vu;

   beta_vu = Re(0.1*oneOver16PiSqr*vu*(3*Sqr(g1) + 5*(-2*(3*traceYuAdjYu
      + traceYvAdjYv) + 3*Sqr(g2))));


   return beta_vu;
}

/**
 * Calculates the two-loop beta function of vu.
 *
 * @return two-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_vu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;


   double beta_vu;

   beta_vu = Re(-0.005*twoLoop*vu*(207*Power(g1,4) + 10*Sqr(g1)*(34*
      traceYuAdjYu + 6*traceYvAdjYv + 9*Sqr(g2)) + 25*(11*Power(g2,4) + 12*(3*
      traceYuAdjYu + traceYvAdjYv)*Sqr(g2) - 8*(3*traceYdAdjYuYuAdjYd +
      traceYeAdjYvYvAdjYe + 9*traceYuAdjYuYuAdjYu + 3*traceYvAdjYvYvAdjYv - 16*
      traceYuAdjYu*Sqr(g3)))));


   return beta_vu;
}

/**
 * Calculates the three-loop beta function of vu.
 *
 * @return three-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_vu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return beta_vu;
}

} // namespace flexiblesusy
