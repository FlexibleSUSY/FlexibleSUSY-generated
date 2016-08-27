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

// File generated at Sat 27 Aug 2016 11:51:29

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g1d.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_g1d;

   beta_g1d = Re(0.05*g1d*oneOver16PiSqr*(60*traceYdAdjYd + 20*
      traceYeAdjYe - 9*Sqr(g1) + 55*Sqr(g1d) + 15*Sqr(g1dp) - 165*Sqr(g2) + 10*
      Sqr(g2u)));


   return beta_g1d;
}

/**
 * Calculates the two-loop beta function of g1d.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_g1d;

   beta_g1d = Re(twoLoop*(0.645*Power(g1,4)*g1d - 3.5*Power(g1d,5) -
      0.3125*g1d*Power(g1dp,4) - 31.75*g1d*Power(g2,4) - 1.3125*g1d*Power(g2u,4
      ) - 10*Power(g1d,3)*Lambda1 - g1dp*g2u*g2up*Lambda4 + g1d*Lambda3*Lambda4
      - 6.75*g1d*traceYdAdjYdYdAdjYd - 2.25*g1d*traceYdAdjYuYuAdjYd - 2.25*g1d
      *traceYeAdjYeYeAdjYe + 2.71875*Power(g1d,3)*Sqr(g1) - 3.6875*Power(g1d,3)
      *Sqr(g1dp) - 2*g1d*Lambda1*Sqr(g1dp) + 0.39375*g1d*Sqr(g1)*Sqr(g1dp) -
      0.375*g1d*traceYeAdjYe*(-5*Sqr(g1) + 5*Sqr(g1d) + Sqr(g1dp) - 5*Sqr(g2))
      + 27.34375*Power(g1d,3)*Sqr(g2) - 1.5*g1dp*g2u*g2up*Sqr(g2) + 0.45*g1d*
      Sqr(g1)*Sqr(g2) + 3.46875*g1d*Sqr(g1dp)*Sqr(g2) - 1.3125*Power(g1d,3)*Sqr
      (g2u) - 2*g1d*Lambda3*Sqr(g2u) + g1d*Lambda4*Sqr(g2u) - 2.25*g1d*
      traceYuAdjYu*Sqr(g2u) - 0.2625*g1d*Sqr(g1)*Sqr(g2u) + 2.6875*g1d*Sqr(g2)*
      Sqr(g2u) - 0.4375*g1d*Sqr(g1dp)*Sqr(g2up) - 0.4375*g1d*Sqr(g2u)*Sqr(g2up)
      + 0.125*g1d*traceYdAdjYd*(5*Sqr(g1) - 45*Sqr(g1d) - 9*Sqr(g1dp) + 45*Sqr
      (g2) + 160*Sqr(g3)) + 6*g1d*Sqr(Lambda1) + g1d*Sqr(Lambda3) + g1d*Sqr(
      Lambda4) + 1.5*g1d*Sqr(Lambda5) + 4.5*g1d*Sqr(Lambda6) + 1.5*g1d*Sqr(
      Lambda7)));


   return beta_g1d;
}

/**
 * Calculates the three-loop beta function of g1d.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1d;

   beta_g1d = 0;


   return beta_g1d;
}

} // namespace flexiblesusy
