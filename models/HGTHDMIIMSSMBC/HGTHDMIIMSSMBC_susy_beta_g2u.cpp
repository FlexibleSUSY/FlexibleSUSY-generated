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

// File generated at Sun 26 Aug 2018 14:06:37

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

   beta_g2u = Re(0.05*g2u*oneOver16PiSqr*(-9*Sqr(g1) + 5*(12*traceYuAdjYu + 2*
      Sqr(g1d) - 33*Sqr(g2) + 11*Sqr(g2u) + 3*Sqr(g2up))));


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

   beta_g2u = Re(twoLoop*(-(g1d*g1dp*g2up*Lambda4) + g2u*Lambda3*Lambda4 - 2.25
      *g2u*traceYdAdjYuYuAdjYd - 6.75*g2u*traceYuAdjYuYuAdjYu + 1.5*g2u*AbsSqr(
      Lambda5) + 1.5*g2u*AbsSqr(Lambda6) + 4.5*g2u*AbsSqr(Lambda7) - 10*Lambda2
      *Cube(g2u) - 5.625*traceYuAdjYu*Cube(g2u) - 3.5*Power5(g2u) + 0.645*g2u*
      Quad(g1) - 1.3125*g2u*Quad(g1d) - 31.75*g2u*Quad(g2) - 0.3125*g2u*Quad(
      g2up) + 2.125*g2u*traceYuAdjYu*Sqr(g1) + 2.71875*Cube(g2u)*Sqr(g1) - 2*
      g2u*Lambda3*Sqr(g1d) + g2u*Lambda4*Sqr(g1d) - 2.25*g2u*traceYdAdjYd*Sqr(
      g1d) - 0.75*g2u*traceYeAdjYe*Sqr(g1d) - 1.3125*Cube(g2u)*Sqr(g1d) -
      0.2625*g2u*Sqr(g1)*Sqr(g1d) - 0.4375*g2u*Sqr(g1d)*Sqr(g1dp) - 1.5*g1d*
      g1dp*g2up*Sqr(g2) + 5.625*g2u*traceYuAdjYu*Sqr(g2) + 27.34375*Cube(g2u)*
      Sqr(g2) + 0.45*g2u*Sqr(g1)*Sqr(g2) + 2.6875*g2u*Sqr(g1d)*Sqr(g2) - 2*g2u*
      Lambda2*Sqr(g2up) - 1.125*g2u*traceYuAdjYu*Sqr(g2up) - 3.6875*Cube(g2u)*
      Sqr(g2up) + 0.39375*g2u*Sqr(g1)*Sqr(g2up) - 0.4375*g2u*Sqr(g1dp)*Sqr(g2up
      ) + 3.46875*g2u*Sqr(g2)*Sqr(g2up) + 20*g2u*traceYuAdjYu*Sqr(g3) + 6*g2u*
      Sqr(Lambda2) + g2u*Sqr(Lambda3) + g2u*Sqr(Lambda4)));


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

/**
 * Calculates the 4-loop beta function of g2u.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

} // namespace flexiblesusy
