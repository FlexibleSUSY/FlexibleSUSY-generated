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

// File generated at Fri 20 Oct 2017 08:35:27

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2up.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2up;

   beta_g2up = Re(0.05*g2up*oneOver16PiSqr*(-9*Sqr(g1) + 5*(12*
      traceYuAdjYu + 2*Sqr(g1dp) - 9*Sqr(g2) + 9*Sqr(g2u) + 5*Sqr(g2up))));


   return beta_g2up;
}

/**
 * Calculates the 2-loop beta function of g2up.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2up;

   beta_g2up = Re(twoLoop*(-1.5*g1d*g1dp*g2u*(2*Lambda4 + 3*Sqr(g2)) -
      1.3125*g2up*Sqr(g1d)*(Sqr(g1dp) + Sqr(g2u)) - 0.00125*g2up*(-516*Quad(g1)
      + 5*Sqr(g1)*(-340*traceYuAdjYu + 42*Sqr(g1dp) + 216*Sqr(g2) - 189*Sqr(
      g2u) - 309*Sqr(g2up)) + 25*(14*Quad(g1dp) + 120*Quad(g2) + Sqr(g1dp)*(8*(
      8*Lambda3 + 4*Lambda4 + 9*traceYdAdjYd + 3*traceYeAdjYe) - 102*Sqr(g2) +
      14*Sqr(g2up)) - 3*Sqr(g2)*(60*traceYuAdjYu + 183*Sqr(g2u) + 55*Sqr(g2up))
      + 2*(99*Quad(g2u) + Sqr(g2u)*(96*Lambda2 + 54*traceYuAdjYu + 9*Sqr(g2up)
      ) + 2*(6*Quad(g2up) + 3*(16*Lambda2 + 9*traceYuAdjYu)*Sqr(g2up) - 2*(4*
      Lambda3*Lambda4 - 9*traceYdAdjYuYuAdjYd - 27*traceYuAdjYuYuAdjYu + 80*
      traceYuAdjYu*Sqr(g3) + 24*Sqr(Lambda2) + 4*Sqr(Lambda3) + 4*Sqr(Lambda4)
      + 6*Sqr(Lambda5) + 6*Sqr(Lambda6) + 18*Sqr(Lambda7))))))));


   return beta_g2up;
}

/**
 * Calculates the 3-loop beta function of g2up.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2up;

   beta_g2up = 0;


   return beta_g2up;
}

} // namespace flexiblesusy
