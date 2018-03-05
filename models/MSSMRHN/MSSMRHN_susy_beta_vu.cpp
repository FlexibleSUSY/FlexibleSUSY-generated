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

// File generated at Mon 5 Mar 2018 18:45:58

#include "MSSMRHN_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vu.
 *
 * @return 1-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_vu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_vu;

   beta_vu = Re(0.1*oneOver16PiSqr*vu*(3*Sqr(g1) + 5*(-2*(3*traceYuAdjYu
      + traceYvAdjYv) + 3*Sqr(g2))));


   return beta_vu;
}

/**
 * Calculates the 2-loop beta function of vu.
 *
 * @return 2-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_vu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;


   double beta_vu;

   beta_vu = Re(-0.005*twoLoop*vu*(207*Quad(g1) + 10*Sqr(g1)*(34*
      traceYuAdjYu + 6*traceYvAdjYv + 9*Sqr(g2)) + 25*(11*Quad(g2) + 12*(3*
      traceYuAdjYu + traceYvAdjYv)*Sqr(g2) - 8*(3*traceYdAdjYuYuAdjYd +
      traceYeAdjYvYvAdjYe + 9*traceYuAdjYuYuAdjYu + 3*traceYvAdjYvYvAdjYv - 16*
      traceYuAdjYu*Sqr(g3)))));


   return beta_vu;
}

/**
 * Calculates the 3-loop beta function of vu.
 *
 * @return 3-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_vu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return beta_vu;
}

/**
 * Calculates the 4-loop beta function of vu.
 *
 * @return 4-loop beta function
 */
double MSSMRHN_susy_parameters::calc_beta_vu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return beta_vu;
}

} // namespace flexiblesusy
