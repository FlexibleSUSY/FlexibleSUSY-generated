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

// File generated at Sun 26 Aug 2018 14:06:35

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

   beta_g1d = Re(0.05*g1d*oneOver16PiSqr*(-9*Sqr(g1) + 5*(12*traceYdAdjYd + 4*
      traceYeAdjYe + 11*Sqr(g1d) + 3*Sqr(g1dp) - 33*Sqr(g2) + 2*Sqr(g2u))));


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

   beta_g1d = Re(twoLoop*(-(g1dp*g2u*g2up*Lambda4) + g1d*Lambda3*Lambda4 - 6.75
      *g1d*traceYdAdjYdYdAdjYd - 2.25*g1d*traceYdAdjYuYuAdjYd - 2.25*g1d*
      traceYeAdjYeYeAdjYe + 1.5*g1d*AbsSqr(Lambda5) + 4.5*g1d*AbsSqr(Lambda6) +
      1.5*g1d*AbsSqr(Lambda7) - 10*Lambda1*Cube(g1d) - 5.625*traceYdAdjYd*Cube(
      g1d) - 1.875*traceYeAdjYe*Cube(g1d) - 3.5*Power5(g1d) + 0.645*g1d*Quad(g1
      ) - 0.3125*g1d*Quad(g1dp) - 31.75*g1d*Quad(g2) - 1.3125*g1d*Quad(g2u) +
      0.625*g1d*traceYdAdjYd*Sqr(g1) + 1.875*g1d*traceYeAdjYe*Sqr(g1) + 2.71875
      *Cube(g1d)*Sqr(g1) - 2*g1d*Lambda1*Sqr(g1dp) - 1.125*g1d*traceYdAdjYd*Sqr
      (g1dp) - 0.375*g1d*traceYeAdjYe*Sqr(g1dp) - 3.6875*Cube(g1d)*Sqr(g1dp) +
      0.39375*g1d*Sqr(g1)*Sqr(g1dp) - 1.5*g1dp*g2u*g2up*Sqr(g2) + 5.625*g1d*
      traceYdAdjYd*Sqr(g2) + 1.875*g1d*traceYeAdjYe*Sqr(g2) + 27.34375*Cube(g1d
      )*Sqr(g2) + 0.45*g1d*Sqr(g1)*Sqr(g2) + 3.46875*g1d*Sqr(g1dp)*Sqr(g2) - 2*
      g1d*Lambda3*Sqr(g2u) + g1d*Lambda4*Sqr(g2u) - 2.25*g1d*traceYuAdjYu*Sqr(
      g2u) - 1.3125*Cube(g1d)*Sqr(g2u) - 0.2625*g1d*Sqr(g1)*Sqr(g2u) + 2.6875*
      g1d*Sqr(g2)*Sqr(g2u) - 0.4375*g1d*Sqr(g1dp)*Sqr(g2up) - 0.4375*g1d*Sqr(
      g2u)*Sqr(g2up) + 20*g1d*traceYdAdjYd*Sqr(g3) + 6*g1d*Sqr(Lambda1) + g1d*
      Sqr(Lambda3) + g1d*Sqr(Lambda4)));


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

} // namespace flexiblesusy
