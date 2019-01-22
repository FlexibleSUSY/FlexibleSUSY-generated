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

// File generated at Tue 22 Jan 2019 16:22:21

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g1d.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_g1d;

   beta_g1d = Re(0.05*g1d*oneOver16PiSqr*(60*traceYdAdjYd + 20*traceYeAdjYe - 9
      *Sqr(g1) + 55*Sqr(g1d) + 15*Sqr(g1dp) - 165*Sqr(g2) + 10*Sqr(g2u)));


   return beta_g1d;
}

/**
 * Calculates the 2-loop beta function of g1d.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_g1d;

   beta_g1d = Re(0.00125*twoLoop*(-800*g1dp*g2u*g2up*Lambda4 + 800*g1d*Lambda3*
      Lambda4 - 5400*g1d*traceYdAdjYdYdAdjYd - 1800*g1d*traceYdAdjYuYuAdjYd -
      1800*g1d*traceYeAdjYeYeAdjYe + 1200*g1d*AbsSqr(Lambda5) + 3600*g1d*AbsSqr
      (Lambda6) + 1200*g1d*AbsSqr(Lambda7) - 8000*Lambda1*Cube(g1d) - 4500*
      traceYdAdjYd*Cube(g1d) - 1500*traceYeAdjYe*Cube(g1d) - 2800*Power5(g1d) +
      516*g1d*Quad(g1) - 250*g1d*Quad(g1dp) - 25400*g1d*Quad(g2) - 1050*g1d*
      Quad(g2u) + 500*g1d*traceYdAdjYd*Sqr(g1) + 1500*g1d*traceYeAdjYe*Sqr(g1)
      + 2175*Cube(g1d)*Sqr(g1) - 1600*g1d*Lambda1*Sqr(g1dp) - 900*g1d*
      traceYdAdjYd*Sqr(g1dp) - 300*g1d*traceYeAdjYe*Sqr(g1dp) - 2950*Cube(g1d)*
      Sqr(g1dp) + 315*g1d*Sqr(g1)*Sqr(g1dp) - 1200*g1dp*g2u*g2up*Sqr(g2) + 4500
      *g1d*traceYdAdjYd*Sqr(g2) + 1500*g1d*traceYeAdjYe*Sqr(g2) + 21875*Cube(
      g1d)*Sqr(g2) + 360*g1d*Sqr(g1)*Sqr(g2) + 2775*g1d*Sqr(g1dp)*Sqr(g2) -
      1600*g1d*Lambda3*Sqr(g2u) + 800*g1d*Lambda4*Sqr(g2u) - 1800*g1d*
      traceYuAdjYu*Sqr(g2u) - 1050*Cube(g1d)*Sqr(g2u) - 210*g1d*Sqr(g1)*Sqr(g2u
      ) + 2150*g1d*Sqr(g2)*Sqr(g2u) - 350*g1d*Sqr(g1dp)*Sqr(g2up) - 350*g1d*Sqr
      (g2u)*Sqr(g2up) + 16000*g1d*traceYdAdjYd*Sqr(g3) + 4800*g1d*Sqr(Lambda1)
      + 800*g1d*Sqr(Lambda3) + 800*g1d*Sqr(Lambda4)));


   return beta_g1d;
}

/**
 * Calculates the 3-loop beta function of g1d.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1d;

   beta_g1d = 0;


   return beta_g1d;
}

/**
 * Calculates the 4-loop beta function of g1d.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1d;

   beta_g1d = 0;


   return beta_g1d;
}

/**
 * Calculates the 5-loop beta function of g1d.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1d;

   beta_g1d = 0;


   return beta_g1d;
}

} // namespace flexiblesusy
