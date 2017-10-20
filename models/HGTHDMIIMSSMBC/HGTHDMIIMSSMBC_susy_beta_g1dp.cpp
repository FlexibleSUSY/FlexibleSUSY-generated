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

// File generated at Fri 20 Oct 2017 08:35:24

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g1dp.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_g1dp;

   beta_g1dp = Re(0.05*g1dp*oneOver16PiSqr*(-9*Sqr(g1) + 5*(12*
      traceYdAdjYd + 4*traceYeAdjYe + 9*Sqr(g1d) + 5*Sqr(g1dp) - 9*Sqr(g2) + 2*
      Sqr(g2up))));


   return beta_g1dp;
}

/**
 * Calculates the 2-loop beta function of g1dp.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_g1dp;

   beta_g1dp = Re(0.00125*twoLoop*(516*g1dp*Quad(g1) + 5*g1dp*Sqr(g1)*(
      100*traceYdAdjYd + 300*traceYeAdjYe + 189*Sqr(g1d) + 309*Sqr(g1dp) - 216*
      Sqr(g2) - 42*Sqr(g2up)) - 25*(198*g1dp*Quad(g1d) + 48*g1d*g2u*g2up*(2*
      Lambda4 + 3*Sqr(g2)) + 3*g1dp*Sqr(g1d)*(64*Lambda1 + 36*traceYdAdjYd + 12
      *traceYeAdjYe + 6*Sqr(g1dp) - 183*Sqr(g2) + 14*Sqr(g2u)) + g1dp*(24*Quad(
      g1dp) + Sqr(g1dp)*(-165*Sqr(g2) + 2*(96*Lambda1 + 54*traceYdAdjYd + 18*
      traceYeAdjYe + 7*Sqr(g2up))) + 2*(-16*Lambda3*Lambda4 + 108*
      traceYdAdjYdYdAdjYd + 36*traceYdAdjYuYuAdjYd + 36*traceYeAdjYeYeAdjYe +
      60*Quad(g2) + 7*Quad(g2up) + 32*Lambda3*Sqr(g2up) + 16*Lambda4*Sqr(g2up)
      + 36*traceYuAdjYu*Sqr(g2up) + 21*Sqr(g2u)*Sqr(g2up) - 3*Sqr(g2)*(30*
      traceYdAdjYd + 10*traceYeAdjYe + 17*Sqr(g2up)) - 320*traceYdAdjYd*Sqr(g3)
      - 96*Sqr(Lambda1) - 16*Sqr(Lambda3) - 16*Sqr(Lambda4) - 24*Sqr(Lambda5)
      - 72*Sqr(Lambda6) - 24*Sqr(Lambda7))))));


   return beta_g1dp;
}

/**
 * Calculates the 3-loop beta function of g1dp.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1dp;

   beta_g1dp = 0;


   return beta_g1dp;
}

} // namespace flexiblesusy
