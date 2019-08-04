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

// File generated at Sun 4 Aug 2019 19:24:27

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2u.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(0.05*oneOver16PiSqr*(20*g2d*gYd*gYu + 60*g2u*traceYdAdjYd + 20
      *g2u*traceYeAdjYe + 60*g2u*traceYuAdjYu + 55*Cube(g2u) - 9*g2u*Sqr(g1) -
      165*g2u*Sqr(g2) + 20*g2u*Sqr(g2d) + 10*g2u*Sqr(gYd) + 15*g2u*Sqr(gYu)));


   return beta_g2u;
}

/**
 * Calculates the 2-loop beta function of g2u.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(0.0004166666666666667*twoLoop*(-7200*g2d*gYd*gYu*traceYdAdjYd
      - 16200*g2u*traceYdAdjYdYdAdjYd + 3600*g2u*traceYdAdjYuYuAdjYd - 2400*g2d
      *gYd*gYu*traceYeAdjYe - 5400*g2u*traceYeAdjYeYeAdjYe - 7200*g2d*gYd*gYu*
      traceYuAdjYu - 16200*g2u*traceYuAdjYuYuAdjYu - 5400*gYd*gYu*Cube(g2d) -
      13500*traceYdAdjYd*Cube(g2u) - 4500*traceYeAdjYe*Cube(g2u) - 13500*
      traceYuAdjYu*Cube(g2u) - 3000*g2d*gYu*Cube(gYd) - 3600*g2d*gYd*Cube(gYu)
      - 2400*g2d*gYd*gYu*Lambdax - 12000*Cube(g2u)*Lambdax - 8400*Power5(g2u) +
      1404*g2u*Quad(g1) - 81800*g2u*Quad(g2) - 3300*g2u*Quad(g2d) - 1350*g2u*
      Quad(gYd) - 750*g2u*Quad(gYu) + 360*g2d*gYd*gYu*Sqr(g1) + 1500*g2u*
      traceYdAdjYd*Sqr(g1) + 4500*g2u*traceYeAdjYe*Sqr(g1) + 5100*g2u*
      traceYuAdjYu*Sqr(g1) + 6525*Cube(g2u)*Sqr(g1) + 5400*g2d*gYd*gYu*Sqr(g2)
      + 13500*g2u*traceYdAdjYd*Sqr(g2) + 4500*g2u*traceYeAdjYe*Sqr(g2) + 13500*
      g2u*traceYuAdjYu*Sqr(g2) + 65625*Cube(g2u)*Sqr(g2) + 1080*g2u*Sqr(g1)*Sqr
      (g2) + 1800*g2u*traceYdAdjYd*Sqr(g2d) + 600*g2u*traceYeAdjYe*Sqr(g2d) +
      1800*g2u*traceYuAdjYu*Sqr(g2d) - 8100*Cube(g2u)*Sqr(g2d) - 2400*g2u*
      Lambdax*Sqr(g2d) + 360*g2u*Sqr(g1)*Sqr(g2d) + 10200*g2u*Sqr(g2)*Sqr(g2d)
      - 9600*g2d*gYd*gYu*Sqr(g2u) + 48000*g2u*traceYdAdjYd*Sqr(g3) + 48000*g2u*
      traceYuAdjYu*Sqr(g3) - 2250*Cube(g2u)*Sqr(gYd) + 450*g2u*Sqr(g1)*Sqr(gYd)
      + 2250*g2u*Sqr(g2)*Sqr(gYd) - 1950*g2u*Sqr(g2d)*Sqr(gYd) - 2700*g2u*
      traceYdAdjYd*Sqr(gYu) - 900*g2u*traceYeAdjYe*Sqr(gYu) - 2700*g2u*
      traceYuAdjYu*Sqr(gYu) - 8850*Cube(g2u)*Sqr(gYu) - 2400*g2u*Lambdax*Sqr(
      gYu) + 945*g2u*Sqr(g1)*Sqr(gYu) + 8325*g2u*Sqr(g2)*Sqr(gYu) - 4650*g2u*
      Sqr(g2d)*Sqr(gYu) - 3600*g2u*Sqr(gYd)*Sqr(gYu) + 3600*g2u*Sqr(Lambdax)));


   return beta_g2u;
}

/**
 * Calculates the 3-loop beta function of g2u.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_3_loop(const Susy_traces& susy_traces) const
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
double SplitMSSM_susy_parameters::calc_beta_g2u_4_loop(const Susy_traces& susy_traces) const
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
double SplitMSSM_susy_parameters::calc_beta_g2u_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

} // namespace flexiblesusy
