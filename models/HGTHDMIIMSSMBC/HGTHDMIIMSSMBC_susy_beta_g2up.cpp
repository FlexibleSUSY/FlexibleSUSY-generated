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

// File generated at Sun 26 Aug 2018 14:06:36

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

   beta_g2up = Re(0.05*g2up*oneOver16PiSqr*(-9*Sqr(g1) + 5*(12*traceYuAdjYu + 2
      *Sqr(g1dp) - 9*Sqr(g2) + 9*Sqr(g2u) + 5*Sqr(g2up))));


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

   beta_g2up = Re(twoLoop*(-3*g1d*g1dp*g2u*Lambda4 + g2up*Lambda3*Lambda4 -
      2.25*g2up*traceYdAdjYuYuAdjYd - 6.75*g2up*traceYuAdjYuYuAdjYu + 1.5*g2up*
      AbsSqr(Lambda5) + 1.5*g2up*AbsSqr(Lambda6) + 4.5*g2up*AbsSqr(Lambda7) - 6
      *Lambda2*Cube(g2up) - 3.375*traceYuAdjYu*Cube(g2up) - 0.75*Power5(g2up) +
      0.645*g2up*Quad(g1) - 0.4375*g2up*Quad(g1dp) - 3.75*g2up*Quad(g2) -
      6.1875*g2up*Quad(g2u) + 2.125*g2up*traceYuAdjYu*Sqr(g1) + 1.93125*Cube(
      g2up)*Sqr(g1) - 2*g2up*Lambda3*Sqr(g1dp) - g2up*Lambda4*Sqr(g1dp) - 2.25*
      g2up*traceYdAdjYd*Sqr(g1dp) - 0.75*g2up*traceYeAdjYe*Sqr(g1dp) - 0.4375*
      Cube(g2up)*Sqr(g1dp) - 0.2625*g2up*Sqr(g1)*Sqr(g1dp) - 1.3125*g2up*Sqr(
      g1d)*Sqr(g1dp) - 4.5*g1d*g1dp*g2u*Sqr(g2) + 5.625*g2up*traceYuAdjYu*Sqr(
      g2) + 5.15625*Cube(g2up)*Sqr(g2) - 1.35*g2up*Sqr(g1)*Sqr(g2) + 3.1875*
      g2up*Sqr(g1dp)*Sqr(g2) - 6*g2up*Lambda2*Sqr(g2u) - 3.375*g2up*
      traceYuAdjYu*Sqr(g2u) - 0.5625*Cube(g2up)*Sqr(g2u) + 1.18125*g2up*Sqr(g1)
      *Sqr(g2u) - 1.3125*g2up*Sqr(g1d)*Sqr(g2u) + 17.15625*g2up*Sqr(g2)*Sqr(g2u
      ) + 20*g2up*traceYuAdjYu*Sqr(g3) + 6*g2up*Sqr(Lambda2) + g2up*Sqr(Lambda3
      ) + g2up*Sqr(Lambda4)));


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

/**
 * Calculates the 4-loop beta function of g2up.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2up;

   beta_g2up = 0;


   return beta_g2up;
}

} // namespace flexiblesusy
