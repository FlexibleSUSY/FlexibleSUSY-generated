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

// File generated at Fri 10 Apr 2020 19:25:42

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

   beta_g2u = Re(0.05*g2u*oneOver16PiSqr*(60*traceYuAdjYu - 9*Sqr(g1) + 10*Sqr(
      g1d) - 165*Sqr(g2) + 55*Sqr(g2u) + 15*Sqr(g2up)));


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

   beta_g2u = Re(0.00125*twoLoop*(-800*g1d*g1dp*g2up*Lambda4 + 800*g2u*Lambda3*
      Lambda4 - 1800*g2u*traceYdAdjYuYuAdjYd - 5400*g2u*traceYuAdjYuYuAdjYu +
      1200*g2u*AbsSqr(Lambda5) + 1200*g2u*AbsSqr(Lambda6) + 3600*g2u*AbsSqr(
      Lambda7) - 8000*Lambda2*Cube(g2u) - 4500*traceYuAdjYu*Cube(g2u) - 2800*
      Power5(g2u) + 516*g2u*Quad(g1) - 1050*g2u*Quad(g1d) - 25400*g2u*Quad(g2)
      - 250*g2u*Quad(g2up) + 1700*g2u*traceYuAdjYu*Sqr(g1) + 2175*Cube(g2u)*Sqr
      (g1) - 1600*g2u*Lambda3*Sqr(g1d) + 800*g2u*Lambda4*Sqr(g1d) - 1800*g2u*
      traceYdAdjYd*Sqr(g1d) - 600*g2u*traceYeAdjYe*Sqr(g1d) - 1050*Cube(g2u)*
      Sqr(g1d) - 210*g2u*Sqr(g1)*Sqr(g1d) - 350*g2u*Sqr(g1d)*Sqr(g1dp) - 1200*
      g1d*g1dp*g2up*Sqr(g2) + 4500*g2u*traceYuAdjYu*Sqr(g2) + 21875*Cube(g2u)*
      Sqr(g2) + 360*g2u*Sqr(g1)*Sqr(g2) + 2150*g2u*Sqr(g1d)*Sqr(g2) - 1600*g2u*
      Lambda2*Sqr(g2up) - 900*g2u*traceYuAdjYu*Sqr(g2up) - 2950*Cube(g2u)*Sqr(
      g2up) + 315*g2u*Sqr(g1)*Sqr(g2up) - 350*g2u*Sqr(g1dp)*Sqr(g2up) + 2775*
      g2u*Sqr(g2)*Sqr(g2up) + 16000*g2u*traceYuAdjYu*Sqr(g3) + 4800*g2u*Sqr(
      Lambda2) + 800*g2u*Sqr(Lambda3) + 800*g2u*Sqr(Lambda4)));


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

/**
 * Calculates the 5-loop beta function of g2u.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2u_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

} // namespace flexiblesusy
