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

// File generated at Fri 20 Oct 2017 08:35:28

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2u.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(0.05*g2u*oneOver16PiSqr*(-9*Sqr(g1) + 5*(12*traceYuAdjYu
      + 2*Sqr(g1d) - 33*Sqr(g2) + 11*Sqr(g2u) + 3*Sqr(g2up))));


   return beta_g2u;
}

/**
 * Calculates the 2-loop beta function of g2u.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(0.00125*twoLoop*(516*g2u*Quad(g1) + 5*g2u*Sqr(g1)*(340*
      traceYuAdjYu - 42*Sqr(g1d) + 72*Sqr(g2) + 435*Sqr(g2u) + 63*Sqr(g2up)) -
      25*(42*g2u*Quad(g1d) + 16*g1d*g1dp*g2up*(2*Lambda4 + 3*Sqr(g2)) + 2*g2u*
      Sqr(g1d)*(32*Lambda3 - 16*Lambda4 + 36*traceYdAdjYd + 12*traceYeAdjYe + 7
      *Sqr(g1dp) - 43*Sqr(g2) + 21*Sqr(g2u)) + g2u*(1016*Quad(g2) - Sqr(g2)*(
      180*traceYuAdjYu + 875*Sqr(g2u) + 111*Sqr(g2up)) + 2*(-16*Lambda3*Lambda4
      + 36*traceYdAdjYuYuAdjYd + 108*traceYuAdjYuYuAdjYu + 56*Quad(g2u) + 5*
      Quad(g2up) + 32*Lambda2*Sqr(g2up) + 18*traceYuAdjYu*Sqr(g2up) + 7*Sqr(
      g1dp)*Sqr(g2up) + Sqr(g2u)*(160*Lambda2 + 90*traceYuAdjYu + 59*Sqr(g2up))
      - 320*traceYuAdjYu*Sqr(g3) - 96*Sqr(Lambda2) - 16*Sqr(Lambda3) - 16*Sqr(
      Lambda4) - 24*Sqr(Lambda5) - 24*Sqr(Lambda6) - 72*Sqr(Lambda7))))));


   return beta_g2u;
}

/**
 * Calculates the 3-loop beta function of g2u.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

} // namespace flexiblesusy
