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

// File generated at Tue 22 Jan 2019 16:46:37

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2d.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2d;

   beta_g2d = Re(0.05*oneOver16PiSqr*(20*g2u*gYd*gYu + 60*g2d*traceYdAdjYd + 20
      *g2d*traceYeAdjYe + 60*g2d*traceYuAdjYu + 55*Cube(g2d) - 9*g2d*Sqr(g1) -
      165*g2d*Sqr(g2) + 20*g2d*Sqr(g2u) + 15*g2d*Sqr(gYd) + 10*g2d*Sqr(gYu)));


   return beta_g2d;
}

/**
 * Calculates the 2-loop beta function of g2d.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2d;

   beta_g2d = Re(0.0004166666666666667*twoLoop*(-7200*g2u*gYd*gYu*traceYdAdjYd
      - 16200*g2d*traceYdAdjYdYdAdjYd + 3600*g2d*traceYdAdjYuYuAdjYd - 2400*g2u
      *gYd*gYu*traceYeAdjYe - 5400*g2d*traceYeAdjYeYeAdjYe - 7200*g2u*gYd*gYu*
      traceYuAdjYu - 16200*g2d*traceYuAdjYuYuAdjYu - 13500*traceYdAdjYd*Cube(
      g2d) - 4500*traceYeAdjYe*Cube(g2d) - 13500*traceYuAdjYu*Cube(g2d) - 5400*
      gYd*gYu*Cube(g2u) - 3600*g2u*gYu*Cube(gYd) - 3000*g2u*gYd*Cube(gYu) -
      2400*g2u*gYd*gYu*Lambdax - 12000*Cube(g2d)*Lambdax - 8400*Power5(g2d) +
      1404*g2d*Quad(g1) - 81800*g2d*Quad(g2) - 3300*g2d*Quad(g2u) - 750*g2d*
      Quad(gYd) - 1350*g2d*Quad(gYu) + 360*g2u*gYd*gYu*Sqr(g1) + 1500*g2d*
      traceYdAdjYd*Sqr(g1) + 4500*g2d*traceYeAdjYe*Sqr(g1) + 5100*g2d*
      traceYuAdjYu*Sqr(g1) + 6525*Cube(g2d)*Sqr(g1) + 5400*g2u*gYd*gYu*Sqr(g2)
      + 13500*g2d*traceYdAdjYd*Sqr(g2) + 4500*g2d*traceYeAdjYe*Sqr(g2) + 13500*
      g2d*traceYuAdjYu*Sqr(g2) + 65625*Cube(g2d)*Sqr(g2) + 1080*g2d*Sqr(g1)*Sqr
      (g2) - 9600*g2u*gYd*gYu*Sqr(g2d) + 1800*g2d*traceYdAdjYd*Sqr(g2u) + 600*
      g2d*traceYeAdjYe*Sqr(g2u) + 1800*g2d*traceYuAdjYu*Sqr(g2u) - 8100*Cube(
      g2d)*Sqr(g2u) - 2400*g2d*Lambdax*Sqr(g2u) + 360*g2d*Sqr(g1)*Sqr(g2u) +
      10200*g2d*Sqr(g2)*Sqr(g2u) + 48000*g2d*traceYdAdjYd*Sqr(g3) + 48000*g2d*
      traceYuAdjYu*Sqr(g3) - 2700*g2d*traceYdAdjYd*Sqr(gYd) - 900*g2d*
      traceYeAdjYe*Sqr(gYd) - 2700*g2d*traceYuAdjYu*Sqr(gYd) - 8850*Cube(g2d)*
      Sqr(gYd) - 2400*g2d*Lambdax*Sqr(gYd) + 945*g2d*Sqr(g1)*Sqr(gYd) + 8325*
      g2d*Sqr(g2)*Sqr(gYd) - 4650*g2d*Sqr(g2u)*Sqr(gYd) - 2250*Cube(g2d)*Sqr(
      gYu) + 450*g2d*Sqr(g1)*Sqr(gYu) + 2250*g2d*Sqr(g2)*Sqr(gYu) - 1950*g2d*
      Sqr(g2u)*Sqr(gYu) - 3600*g2d*Sqr(gYd)*Sqr(gYu) + 3600*g2d*Sqr(Lambdax)));


   return beta_g2d;
}

/**
 * Calculates the 3-loop beta function of g2d.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2d;

   beta_g2d = 0;


   return beta_g2d;
}

/**
 * Calculates the 4-loop beta function of g2d.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2d;

   beta_g2d = 0;


   return beta_g2d;
}

/**
 * Calculates the 5-loop beta function of g2d.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2d;

   beta_g2d = 0;


   return beta_g2d;
}

} // namespace flexiblesusy
