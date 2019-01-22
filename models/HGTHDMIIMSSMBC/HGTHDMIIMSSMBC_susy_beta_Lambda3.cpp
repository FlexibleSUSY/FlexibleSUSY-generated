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

// File generated at Tue 22 Jan 2019 16:22:17

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda3.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda3;

   beta_Lambda3 = Re(0.01*oneOver16PiSqr*(200*g1d*g1dp*g2u*g2up + 1200*Lambda1*
      Lambda3 + 1200*Lambda2*Lambda3 + 400*Lambda1*Lambda4 + 400*Lambda2*
      Lambda4 + 600*Lambda3*traceYdAdjYd - 1200*traceYdAdjYuYuAdjYd + 200*
      Lambda3*traceYeAdjYe + 600*Lambda3*traceYuAdjYu + 200*AbsSqr(Lambda5) +
      400*AbsSqr(Lambda6) + 400*AbsSqr(Lambda7) + 800*Lambda7*Conj(Lambda6) +
      800*Lambda6*Conj(Lambda7) + 27*Quad(g1) + 225*Quad(g2) - 180*Lambda3*Sqr(
      g1) + 300*Lambda3*Sqr(g1d) + 100*Lambda3*Sqr(g1dp) - 900*Lambda3*Sqr(g2)
      - 90*Sqr(g1)*Sqr(g2) + 300*Lambda3*Sqr(g2u) - 500*Sqr(g1d)*Sqr(g2u) + 100
      *Lambda3*Sqr(g2up) - 100*Sqr(g1dp)*Sqr(g2up) + 400*Sqr(Lambda3) + 200*Sqr
      (Lambda4)));


   return beta_Lambda3;
}

/**
 * Calculates the 2-loop beta function of Lambda3.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYuYuAdjYd;


   double beta_Lambda3;

   const double beta_Lambda3_1 = Re(-0.005*twoLoop*(800*g1d*g1dp*g2u*g2up*
      Lambda3 + 6400*Lambda1*Lambda3*Lambda4 + 6400*Lambda2*Lambda3*Lambda4 +
      300*g1dp*g2u*g2up*Cube(g1d) + 500*g1d*g2u*g2up*Cube(g1dp) + 300*g1d*g1dp*
      g2up*Cube(g2u) + 500*g1d*g1dp*g2u*Cube(g2up) + 2400*Cube(Lambda3) + 765*
      Power6(g1) - 4875*Power6(g2) - 540*Lambda1*Quad(g1) - 540*Lambda2*Quad(g1
      ) - 1893*Lambda3*Quad(g1) - 180*Lambda4*Quad(g1) + 1125*Lambda3*Quad(g1d)
      + 225*Lambda3*Quad(g1dp) - 4500*Lambda1*Quad(g2) - 4500*Lambda2*Quad(g2)
      - 225*Lambda3*Quad(g2) - 1500*Lambda4*Quad(g2) + 1125*Lambda3*Quad(g2u) +
      225*Lambda3*Quad(g2up) - 2880*Lambda1*Lambda3*Sqr(g1) - 2880*Lambda2*
      Lambda3*Sqr(g1) - 960*Lambda1*Lambda4*Sqr(g1) - 960*Lambda2*Lambda4*Sqr(
      g1) - 645*Quad(g2)*Sqr(g1) + 7200*Lambda1*Lambda3*Sqr(g1d) + 2400*Lambda1
      *Lambda4*Sqr(g1d) + 27*Quad(g1)*Sqr(g1d) + 3825*Quad(g2)*Sqr(g1d) - 2850*
      Quad(g2u)*Sqr(g1d) - 225*Lambda3*Sqr(g1)*Sqr(g1d) + 2400*Lambda1*Lambda3*
      Sqr(g1dp) + 800*Lambda1*Lambda4*Sqr(g1dp) + 9*Quad(g1)*Sqr(g1dp) + 75*
      Quad(g2)*Sqr(g1dp) - 350*Quad(g2up)*Sqr(g1dp) - 75*Lambda3*Sqr(g1)*Sqr(
      g1dp) + 450*Lambda3*Sqr(g1d)*Sqr(g1dp) - 800*g1d*g1dp*g2u*g2up*Sqr(g2) -
      14400*Lambda1*Lambda3*Sqr(g2) - 14400*Lambda2*Lambda3*Sqr(g2) - 7200*
      Lambda1*Lambda4*Sqr(g2) - 7200*Lambda2*Lambda4*Sqr(g2) + 2400*Lambda3*
      Lambda4*Sqr(g2) - 1005*Quad(g1)*Sqr(g2) + 600*Lambda1*Sqr(g1)*Sqr(g2) +
      600*Lambda2*Sqr(g1)*Sqr(g2) - 330*Lambda3*Sqr(g1)*Sqr(g2) + 360*Lambda4*
      Sqr(g1)*Sqr(g2) - 4125*Lambda3*Sqr(g1d)*Sqr(g2) + 630*Sqr(g1)*Sqr(g1d)*
      Sqr(g2) - 375*Lambda3*Sqr(g1dp)*Sqr(g2) - 30*Sqr(g1)*Sqr(g1dp)*Sqr(g2) +
      7200*Lambda2*Lambda3*Sqr(g2u) + 2400*Lambda2*Lambda4*Sqr(g2u) + 27*Quad(
      g1)*Sqr(g2u) - 2850*Quad(g1d)*Sqr(g2u) + 3825*Quad(g2)*Sqr(g2u) - 225*
      Lambda3*Sqr(g1)*Sqr(g2u) - 1100*Lambda3*Sqr(g1d)*Sqr(g2u) - 350*Sqr(g1d)*
      Sqr(g1dp)*Sqr(g2u) - 4125*Lambda3*Sqr(g2)*Sqr(g2u) + 630*Sqr(g1)*Sqr(g2)*
      Sqr(g2u) + 4000*Sqr(g1d)*Sqr(g2)*Sqr(g2u) + 2400*Lambda2*Lambda3*Sqr(g2up
      ) + 800*Lambda2*Lambda4*Sqr(g2up) + 9*Quad(g1)*Sqr(g2up) - 350*Quad(g1dp)
      *Sqr(g2up) + 75*Quad(g2)*Sqr(g2up) - 75*Lambda3*Sqr(g1)*Sqr(g2up) - 100*
      Lambda3*Sqr(g1dp)*Sqr(g2up) - 450*Sqr(g1d)*Sqr(g1dp)*Sqr(g2up) - 375*
      Lambda3*Sqr(g2)*Sqr(g2up) - 30*Sqr(g1)*Sqr(g2)*Sqr(g2up) + 450*Lambda3*
      Sqr(g2u)*Sqr(g2up) - 350*Sqr(g1d)*Sqr(g2u)*Sqr(g2up) - 450*Sqr(g1dp)*Sqr(
      g2u)*Sqr(g2up) + 12000*Lambda3*Sqr(Lambda1) + 3200*Lambda4*Sqr(Lambda1) +
      12000*Lambda3*Sqr(Lambda2) + 3200*Lambda4*Sqr(Lambda2) + 14400*Lambda1*
      Sqr(Lambda3) + 14400*Lambda2*Sqr(Lambda3) + 800*Lambda4*Sqr(Lambda3) -
      240*Sqr(g1)*Sqr(Lambda3) + 1200*Sqr(g1d)*Sqr(Lambda3) + 400*Sqr(g1dp)*Sqr
      (Lambda3) - 1200*Sqr(g2)*Sqr(Lambda3) + 1200*Sqr(g2u)*Sqr(Lambda3) + 400*
      Sqr(g2up)*Sqr(Lambda3) + 5600*Lambda1*Sqr(Lambda4) + 5600*Lambda2*Sqr(
      Lambda4) + 240*Sqr(g1)*Sqr(Lambda4) + 600*Sqr(g1d)*Sqr(Lambda4) + 200*Sqr
      (g1dp)*Sqr(Lambda4) - 1200*Sqr(g2)*Sqr(Lambda4) + 600*Sqr(g2u)*Sqr(
      Lambda4) + 200*Sqr(g2up)*Sqr(Lambda4)));
   const double beta_Lambda3_2 = Re(-0.01*twoLoop*(3300*Lambda5*Lambda6*Lambda7
       + 7200*Lambda1*Lambda3*traceYdAdjYd + 2400*Lambda1*Lambda4*traceYdAdjYd
      + 1350*Lambda3*traceYdAdjYdYdAdjYd - 1200*traceYdAdjYdYdAdjYuYuAdjYd -
      1500*Lambda3*traceYdAdjYuYuAdjYd - 2400*traceYdAdjYuYuAdjYdYdAdjYd - 3600
      *traceYdAdjYuYuAdjYuYuAdjYd + 2400*Lambda1*Lambda3*traceYeAdjYe + 800*
      Lambda1*Lambda4*traceYeAdjYe + 450*Lambda3*traceYeAdjYeYeAdjYe + 7200*
      Lambda2*Lambda3*traceYuAdjYu + 2400*Lambda2*Lambda4*traceYuAdjYu + 1350*
      Lambda3*traceYuAdjYuYuAdjYu + 3600*Lambda1*AbsSqr(Lambda5) + 3600*Lambda2
      *AbsSqr(Lambda5) + 1800*Lambda3*AbsSqr(Lambda5) + 4400*Lambda4*AbsSqr(
      Lambda5) + 600*traceYdAdjYd*AbsSqr(Lambda5) + 200*traceYeAdjYe*AbsSqr(
      Lambda5) + 600*traceYuAdjYu*AbsSqr(Lambda5) + 11800*Lambda1*AbsSqr(
      Lambda6) + 4400*Lambda2*AbsSqr(Lambda6) + 5700*Lambda3*AbsSqr(Lambda6) +
      6500*Lambda4*AbsSqr(Lambda6) + 2400*traceYdAdjYd*AbsSqr(Lambda6) + 800*
      traceYeAdjYe*AbsSqr(Lambda6) + 4400*Lambda1*AbsSqr(Lambda7) + 11800*
      Lambda2*AbsSqr(Lambda7) + 5700*Lambda3*AbsSqr(Lambda7) + 6500*Lambda4*
      AbsSqr(Lambda7) + 2400*traceYuAdjYu*AbsSqr(Lambda7) + 4100*Lambda1*
      Lambda7*Conj(Lambda6) + 4100*Lambda2*Lambda7*Conj(Lambda6) + 8500*Lambda3
      *Lambda7*Conj(Lambda6) + 4100*Lambda4*Lambda7*Conj(Lambda6) + 2400*
      Lambda7*traceYdAdjYd*Conj(Lambda6) + 800*Lambda7*traceYeAdjYe*Conj(
      Lambda6) + 2400*Lambda7*traceYuAdjYu*Conj(Lambda6) + 4100*Lambda1*Lambda6
      *Conj(Lambda7) + 4100*Lambda2*Lambda6*Conj(Lambda7) + 8500*Lambda3*
      Lambda6*Conj(Lambda7) + 4100*Lambda4*Lambda6*Conj(Lambda7) + 2400*Lambda6
      *traceYdAdjYd*Conj(Lambda7) + 800*Lambda6*traceYeAdjYe*Conj(Lambda7) +
      2400*Lambda6*traceYuAdjYu*Conj(Lambda7) + 3300*Conj(Lambda5)*Conj(Lambda6
      )*Conj(Lambda7) + 1200*Cube(Lambda4) - 45*traceYdAdjYd*Quad(g1) + 225*
      traceYeAdjYe*Quad(g1) + 171*traceYuAdjYu*Quad(g1) + 225*traceYdAdjYd*Quad
      (g2) + 75*traceYeAdjYe*Quad(g2) + 225*traceYuAdjYu*Quad(g2) - 125*Lambda3
      *traceYdAdjYd*Sqr(g1) + 80*traceYdAdjYuYuAdjYd*Sqr(g1) - 375*Lambda3*
      traceYeAdjYe*Sqr(g1) - 425*Lambda3*traceYuAdjYu*Sqr(g1) - 240*AbsSqr(
      Lambda5)*Sqr(g1) - 120*AbsSqr(Lambda6)*Sqr(g1) - 120*AbsSqr(Lambda7)*Sqr(
      g1) - 960*Lambda7*Conj(Lambda6)*Sqr(g1) - 960*Lambda6*Conj(Lambda7)*Sqr(
      g1) + 300*AbsSqr(Lambda5)*Sqr(g1d) + 1200*AbsSqr(Lambda6)*Sqr(g1d) + 1200
      *Lambda7*Conj(Lambda6)*Sqr(g1d) + 1200*Lambda6*Conj(Lambda7)*Sqr(g1d) +
      100*AbsSqr(Lambda5)*Sqr(g1dp) + 400*AbsSqr(Lambda6)*Sqr(g1dp) + 400*
      Lambda7*Conj(Lambda6)*Sqr(g1dp) + 400*Lambda6*Conj(Lambda7)*Sqr(g1dp) -
      1125*Lambda3*traceYdAdjYd*Sqr(g2) - 375*Lambda3*traceYeAdjYe*Sqr(g2) -
      1125*Lambda3*traceYuAdjYu*Sqr(g2) - 5400*Lambda7*Conj(Lambda6)*Sqr(g2) -
      5400*Lambda6*Conj(Lambda7)*Sqr(g2) + 270*traceYdAdjYd*Sqr(g1)*Sqr(g2) +
      330*traceYeAdjYe*Sqr(g1)*Sqr(g2) + 630*traceYuAdjYu*Sqr(g1)*Sqr(g2) + 300
      *AbsSqr(Lambda5)*Sqr(g2u) + 1200*AbsSqr(Lambda7)*Sqr(g2u) + 1200*Lambda7*
      Conj(Lambda6)*Sqr(g2u) + 1200*Lambda6*Conj(Lambda7)*Sqr(g2u) + 100*AbsSqr
      (Lambda5)*Sqr(g2up) + 400*AbsSqr(Lambda7)*Sqr(g2up) + 400*Lambda7*Conj(
      Lambda6)*Sqr(g2up) + 400*Lambda6*Conj(Lambda7)*Sqr(g2up) - 4000*Lambda3*
      traceYdAdjYd*Sqr(g3) + 6400*traceYdAdjYuYuAdjYd*Sqr(g3) - 4000*Lambda3*
      traceYuAdjYu*Sqr(g3) + 1200*traceYdAdjYd*Sqr(Lambda3) + 400*traceYeAdjYe*
      Sqr(Lambda3) + 1200*traceYuAdjYu*Sqr(Lambda3) + 1600*Lambda3*Sqr(Lambda4)
      + 600*traceYdAdjYd*Sqr(Lambda4) + 200*traceYeAdjYe*Sqr(Lambda4) + 600*
      traceYuAdjYu*Sqr(Lambda4) + 3250*Lambda5*Sqr(Lambda6) + 3250*Lambda5*Sqr(
      Lambda7) + 3250*Conj(Lambda5)*Sqr(Conj(Lambda6))));
   const double beta_Lambda3_3 = Re(-32.5*twoLoop*Conj(Lambda5)*Sqr(Conj(
      Lambda7)));

   beta_Lambda3 = beta_Lambda3_1 + beta_Lambda3_2 + beta_Lambda3_3;


   return beta_Lambda3;
}

/**
 * Calculates the 3-loop beta function of Lambda3.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

/**
 * Calculates the 4-loop beta function of Lambda3.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

/**
 * Calculates the 5-loop beta function of Lambda3.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

} // namespace flexiblesusy
