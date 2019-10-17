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

// File generated at Wed 16 Oct 2019 21:52:12

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of gYd.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_gYd;

   beta_gYd = Re(0.05*oneOver16PiSqr*(60*g2d*g2u*gYu + 60*gYd*traceYdAdjYd + 20
      *gYd*traceYeAdjYe + 60*gYd*traceYuAdjYu + 25*Cube(gYd) - 9*gYd*Sqr(g1) -
      45*gYd*Sqr(g2) + 45*gYd*Sqr(g2d) + 30*gYd*Sqr(g2u) + 40*gYd*Sqr(gYu)));


   return beta_gYd;
}

/**
 * Calculates the 2-loop beta function of gYd.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_gYd;

   beta_gYd = Re(0.00125*twoLoop*(-7200*g2d*g2u*gYu*traceYdAdjYd - 5400*gYd*
      traceYdAdjYdYdAdjYd + 1200*gYd*traceYdAdjYuYuAdjYd - 2400*g2d*g2u*gYu*
      traceYeAdjYe - 1800*gYd*traceYeAdjYeYeAdjYe - 7200*g2d*g2u*gYu*
      traceYuAdjYu - 5400*gYd*traceYuAdjYuYuAdjYu - 3600*g2u*gYu*Cube(g2d) -
      6600*g2d*gYu*Cube(g2u) - 2700*traceYdAdjYd*Cube(gYd) - 900*traceYeAdjYe*
      Cube(gYd) - 2700*traceYuAdjYu*Cube(gYd) - 1800*g2d*g2u*Cube(gYu) - 2400*
      g2d*g2u*gYu*Lambdax - 2400*Cube(gYd)*Lambdax - 600*Power5(gYd) + 468*gYd*
      Quad(g1) - 3400*gYd*Quad(g2) - 4950*gYd*Quad(g2d) - 2250*gYd*Quad(g2u) -
      1800*gYd*Quad(gYu) + 360*g2d*g2u*gYu*Sqr(g1) + 500*gYd*traceYdAdjYd*Sqr(
      g1) + 1500*gYd*traceYeAdjYe*Sqr(g1) + 1700*gYd*traceYuAdjYu*Sqr(g1) +
      1545*Cube(gYd)*Sqr(g1) + 10200*g2d*g2u*gYu*Sqr(g2) + 4500*gYd*
      traceYdAdjYd*Sqr(g2) + 1500*gYd*traceYeAdjYe*Sqr(g2) + 4500*gYd*
      traceYuAdjYu*Sqr(g2) + 4125*Cube(gYd)*Sqr(g2) - 1080*gYd*Sqr(g1)*Sqr(g2)
      - 2700*gYd*traceYdAdjYd*Sqr(g2d) - 900*gYd*traceYeAdjYe*Sqr(g2d) - 2700*
      gYd*traceYuAdjYu*Sqr(g2d) - 450*Cube(gYd)*Sqr(g2d) - 2400*gYd*Lambdax*Sqr
      (g2d) + 945*gYd*Sqr(g1)*Sqr(g2d) + 13725*gYd*Sqr(g2)*Sqr(g2d) - 1350*Cube
      (gYd)*Sqr(g2u) + 450*gYd*Sqr(g1)*Sqr(g2u) + 8250*gYd*Sqr(g2)*Sqr(g2u) -
      2100*gYd*Sqr(g2d)*Sqr(g2u) + 16000*gYd*traceYdAdjYd*Sqr(g3) + 16000*gYd*
      traceYuAdjYu*Sqr(g3) - 4800*g2d*g2u*gYu*Sqr(gYd) - 4200*gYd*traceYdAdjYd*
      Sqr(gYu) - 1400*gYd*traceYeAdjYe*Sqr(gYu) - 4200*gYd*traceYuAdjYu*Sqr(gYu
      ) - 3000*Cube(gYd)*Sqr(gYu) - 2400*gYd*Lambdax*Sqr(gYu) + 60*gYd*Sqr(g1)*
      Sqr(gYu) + 3900*gYd*Sqr(g2)*Sqr(gYu) - 3750*gYd*Sqr(g2d)*Sqr(gYu) - 3750*
      gYd*Sqr(g2u)*Sqr(gYu) + 1200*gYd*Sqr(Lambdax)));


   return beta_gYd;
}

/**
 * Calculates the 3-loop beta function of gYd.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYd;

   beta_gYd = 0;


   return beta_gYd;
}

/**
 * Calculates the 4-loop beta function of gYd.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYd;

   beta_gYd = 0;


   return beta_gYd;
}

/**
 * Calculates the 5-loop beta function of gYd.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYd;

   beta_gYd = 0;


   return beta_gYd;
}

} // namespace flexiblesusy
