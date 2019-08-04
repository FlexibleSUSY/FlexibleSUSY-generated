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

// File generated at Sun 4 Aug 2019 19:38:27

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

   beta_g1dp = Re(0.05*g1dp*oneOver16PiSqr*(60*traceYdAdjYd + 20*traceYeAdjYe -
      9*Sqr(g1) + 45*Sqr(g1d) + 25*Sqr(g1dp) - 45*Sqr(g2) + 10*Sqr(g2up)));


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

   beta_g1dp = Re(0.00125*twoLoop*(-2400*g1d*g2u*g2up*Lambda4 + 800*g1dp*
      Lambda3*Lambda4 - 5400*g1dp*traceYdAdjYdYdAdjYd - 1800*g1dp*
      traceYdAdjYuYuAdjYd - 1800*g1dp*traceYeAdjYeYeAdjYe + 1200*g1dp*AbsSqr(
      Lambda5) + 3600*g1dp*AbsSqr(Lambda6) + 1200*g1dp*AbsSqr(Lambda7) - 4800*
      Lambda1*Cube(g1dp) - 2700*traceYdAdjYd*Cube(g1dp) - 900*traceYeAdjYe*Cube
      (g1dp) - 600*Power5(g1dp) + 516*g1dp*Quad(g1) - 4950*g1dp*Quad(g1d) -
      3000*g1dp*Quad(g2) - 350*g1dp*Quad(g2up) + 500*g1dp*traceYdAdjYd*Sqr(g1)
      + 1500*g1dp*traceYeAdjYe*Sqr(g1) + 1545*Cube(g1dp)*Sqr(g1) - 4800*g1dp*
      Lambda1*Sqr(g1d) - 2700*g1dp*traceYdAdjYd*Sqr(g1d) - 900*g1dp*
      traceYeAdjYe*Sqr(g1d) - 450*Cube(g1dp)*Sqr(g1d) + 945*g1dp*Sqr(g1)*Sqr(
      g1d) - 3600*g1d*g2u*g2up*Sqr(g2) + 4500*g1dp*traceYdAdjYd*Sqr(g2) + 1500*
      g1dp*traceYeAdjYe*Sqr(g2) + 4125*Cube(g1dp)*Sqr(g2) - 1080*g1dp*Sqr(g1)*
      Sqr(g2) + 13725*g1dp*Sqr(g1d)*Sqr(g2) - 1050*g1dp*Sqr(g1d)*Sqr(g2u) -
      1600*g1dp*Lambda3*Sqr(g2up) - 800*g1dp*Lambda4*Sqr(g2up) - 1800*g1dp*
      traceYuAdjYu*Sqr(g2up) - 350*Cube(g1dp)*Sqr(g2up) - 210*g1dp*Sqr(g1)*Sqr(
      g2up) + 2550*g1dp*Sqr(g2)*Sqr(g2up) - 1050*g1dp*Sqr(g2u)*Sqr(g2up) +
      16000*g1dp*traceYdAdjYd*Sqr(g3) + 4800*g1dp*Sqr(Lambda1) + 800*g1dp*Sqr(
      Lambda3) + 800*g1dp*Sqr(Lambda4)));


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

/**
 * Calculates the 4-loop beta function of g1dp.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1dp;

   beta_g1dp = 0;


   return beta_g1dp;
}

/**
 * Calculates the 5-loop beta function of g1dp.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1dp_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1dp;

   beta_g1dp = 0;


   return beta_g1dp;
}

} // namespace flexiblesusy
