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

// File generated at Tue 8 Mar 2016 19:13:20

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g2u.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(0.05*g2u*oneOver16PiSqr*(60*traceYuAdjYu - 9*Sqr(g1) + 5
      *(2*Sqr(g1d) - 33*Sqr(g2) + 11*Sqr(g2u) + 3*Sqr(g2up))));


   return beta_g2u;
}

/**
 * Calculates the two-loop beta function of g2u.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(twoLoop*(0.645*Power(g1,4)*g2u - 1.3125*Power(g1d,4)*g2u
      - 31.75*Power(g2,4)*g2u - 3.5*Power(g2u,5) - 0.3125*g2u*Power(g2up,4) -
      10*Power(g2u,3)*Lambda2 - g1d*g1dp*g2up*Lambda4 + g2u*Lambda3*Lambda4 -
      2.25*g2u*traceYdAdjYuYuAdjYd - 5.625*Power(g2u,3)*traceYuAdjYu - 6.75*g2u
      *traceYuAdjYuYuAdjYu + 2.71875*Power(g2u,3)*Sqr(g1) + 2.125*g2u*
      traceYuAdjYu*Sqr(g1) - 1.3125*Power(g2u,3)*Sqr(g1d) - 2*g2u*Lambda3*Sqr(
      g1d) + g2u*Lambda4*Sqr(g1d) - 2.25*g2u*traceYdAdjYd*Sqr(g1d) - 0.75*g2u*
      traceYeAdjYe*Sqr(g1d) - 0.2625*g2u*Sqr(g1)*Sqr(g1d) - 0.4375*g2u*Sqr(g1d)
      *Sqr(g1dp) + 27.34375*Power(g2u,3)*Sqr(g2) - 1.5*g1d*g1dp*g2up*Sqr(g2) +
      5.625*g2u*traceYuAdjYu*Sqr(g2) + 0.45*g2u*Sqr(g1)*Sqr(g2) + 2.6875*g2u*
      Sqr(g1d)*Sqr(g2) - 3.6875*Power(g2u,3)*Sqr(g2up) - 2*g2u*Lambda2*Sqr(g2up
      ) - 1.125*g2u*traceYuAdjYu*Sqr(g2up) + 0.39375*g2u*Sqr(g1)*Sqr(g2up) -
      0.4375*g2u*Sqr(g1dp)*Sqr(g2up) + 3.46875*g2u*Sqr(g2)*Sqr(g2up) + 20*g2u*
      traceYuAdjYu*Sqr(g3) + 6*g2u*Sqr(Lambda2) + g2u*Sqr(Lambda3) + g2u*Sqr(
      Lambda4) + 1.5*g2u*Sqr(Lambda5) + 1.5*g2u*Sqr(Lambda6) + 4.5*g2u*Sqr(
      Lambda7)));


   return beta_g2u;
}

/**
 * Calculates the three-loop beta function of g2u.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

} // namespace flexiblesusy
