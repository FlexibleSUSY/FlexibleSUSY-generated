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

// File generated at Sun 4 Aug 2019 19:38:24

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda2.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(0.005*oneOver16PiSqr*(400*Lambda3*Lambda4 + 2400*Lambda2*
      traceYuAdjYu - 1200*traceYuAdjYuYuAdjYu + 200*AbsSqr(Lambda5) + 2400*
      AbsSqr(Lambda7) + 27*Quad(g1) + 225*Quad(g2) - 500*Quad(g2u) - 100*Quad(
      g2up) - 360*Lambda2*Sqr(g1) - 1800*Lambda2*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) +
      1200*Lambda2*Sqr(g2u) + 400*Lambda2*Sqr(g2up) - 200*Sqr(g2u)*Sqr(g2up) +
      4800*Sqr(Lambda2) + 400*Sqr(Lambda3) + 200*Sqr(Lambda4)));


   return beta_Lambda2;
}

/**
 * Calculates the 2-loop beta function of Lambda2.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYuYuAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYuYuAdjYu = TRACE_STRUCT.
      traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambda2;

   const double beta_Lambda2_1 = Re(-0.0025*twoLoop*(8000*Lambda2*Lambda3*
      Lambda4 + 3400*Lambda5*Lambda6*Lambda7 + 4800*Lambda3*Lambda4*
      traceYdAdjYd + 3600*Lambda2*traceYdAdjYuYuAdjYd - 2400*
      traceYdAdjYuYuAdjYuYuAdjYd + 1600*Lambda3*Lambda4*traceYeAdjYe + 1200*
      Lambda2*traceYuAdjYuYuAdjYu - 12000*traceYuAdjYuYuAdjYuYuAdjYu + 5600*
      Lambda2*AbsSqr(Lambda5) + 8000*Lambda3*AbsSqr(Lambda5) + 8800*Lambda4*
      AbsSqr(Lambda5) + 2400*traceYdAdjYd*AbsSqr(Lambda5) + 800*traceYeAdjYe*
      AbsSqr(Lambda5) - 2400*Lambda2*AbsSqr(Lambda6) + 7200*Lambda3*AbsSqr(
      Lambda6) + 5600*Lambda4*AbsSqr(Lambda6) - 1200*Lambda1*Lambda7*Conj(
      Lambda6) + 6600*Lambda3*Lambda7*Conj(Lambda6) + 5000*Lambda4*Lambda7*Conj
      (Lambda6) - 1200*Lambda1*Lambda6*Conj(Lambda7) + 6600*Lambda3*Lambda6*
      Conj(Lambda7) + 124800*Cube(Lambda2) + 3200*Cube(Lambda3) + 2400*Cube(
      Lambda4) + 765*Power6(g1) - 4875*Power6(g2) - 4700*Power6(g2u) - 500*
      Power6(g2up) - 4146*Lambda2*Quad(g1) - 360*Lambda3*Quad(g1) - 180*Lambda4
      *Quad(g1) + 684*traceYuAdjYu*Quad(g1) - 3450*Lambda2*Quad(g2) - 3000*
      Lambda3*Quad(g2) - 1500*Lambda4*Quad(g2) + 900*traceYuAdjYu*Quad(g2) +
      500*Lambda2*Quad(g2u) + 100*Lambda2*Quad(g2up) - 960*Lambda3*Lambda4*Sqr(
      g1) - 3400*Lambda2*traceYuAdjYu*Sqr(g1) + 640*traceYuAdjYuYuAdjYu*Sqr(g1)
      + 240*AbsSqr(Lambda5)*Sqr(g1) + 1995*Quad(g2)*Sqr(g1) + 2400*Lambda3*
      Lambda4*Sqr(g1d) + 1200*AbsSqr(Lambda5)*Sqr(g1d) - 1000*Quad(g2u)*Sqr(g1d
      ) + 800*Lambda3*Lambda4*Sqr(g1dp) + 400*AbsSqr(Lambda5)*Sqr(g1dp) - 200*
      Quad(g2up)*Sqr(g1dp) - 4800*Lambda3*Lambda4*Sqr(g2) - 9000*Lambda2*
      traceYuAdjYu*Sqr(g2) + 1815*Quad(g1)*Sqr(g2) + 4000*Quad(g2u)*Sqr(g2) -
      2340*Lambda2*Sqr(g1)*Sqr(g2) - 600*Lambda4*Sqr(g1)*Sqr(g2) - 2520*
      traceYuAdjYu*Sqr(g1)*Sqr(g2) + 54*Quad(g1)*Sqr(g2u) + 7650*Quad(g2)*Sqr(
      g2u) - 1700*Quad(g2up)*Sqr(g2u) - 900*Lambda2*Sqr(g1)*Sqr(g2u) + 1800*
      Lambda2*Sqr(g1d)*Sqr(g2u) - 16500*Lambda2*Sqr(g2)*Sqr(g2u) - 1260*Sqr(g1)
      *Sqr(g2)*Sqr(g2u) + 18*Quad(g1)*Sqr(g2up) + 150*Quad(g2)*Sqr(g2up) - 1100
      *Quad(g2u)*Sqr(g2up) - 300*Lambda2*Sqr(g1)*Sqr(g2up) + 600*Lambda2*Sqr(
      g1dp)*Sqr(g2up) - 1500*Lambda2*Sqr(g2)*Sqr(g2up) + 60*Sqr(g1)*Sqr(g2)*Sqr
      (g2up) + 200*Lambda2*Sqr(g2u)*Sqr(g2up) - 200*Sqr(g1d)*Sqr(g2u)*Sqr(g2up)
      - 200*Sqr(g1dp)*Sqr(g2u)*Sqr(g2up) + 800*Sqr(g2)*Sqr(g2u)*Sqr(g2up) -
      32000*Lambda2*traceYuAdjYu*Sqr(g3) + 12800*traceYuAdjYuYuAdjYu*Sqr(g3) +
      57600*traceYuAdjYu*Sqr(Lambda2) - 8640*Sqr(g1)*Sqr(Lambda2) - 43200*Sqr(
      g2)*Sqr(Lambda2) + 28800*Sqr(g2u)*Sqr(Lambda2) + 9600*Sqr(g2up)*Sqr(
      Lambda2) + 8000*Lambda2*Sqr(Lambda3) + 4800*Lambda4*Sqr(Lambda3) + 4800*
      traceYdAdjYd*Sqr(Lambda3) + 1600*traceYeAdjYe*Sqr(Lambda3) - 960*Sqr(g1)*
      Sqr(Lambda3) + 2400*Sqr(g1d)*Sqr(Lambda3) + 800*Sqr(g1dp)*Sqr(Lambda3) -
      4800*Sqr(g2)*Sqr(Lambda3) + 4800*Lambda2*Sqr(Lambda4) + 6400*Lambda3*Sqr(
      Lambda4) + 2400*traceYdAdjYd*Sqr(Lambda4) + 800*traceYeAdjYe*Sqr(Lambda4)
      - 480*Sqr(g1)*Sqr(Lambda4) + 1200*Sqr(g1d)*Sqr(Lambda4) + 400*Sqr(g1dp)*
      Sqr(Lambda4) - 1200*Sqr(g2)*Sqr(Lambda4) + 2000*Lambda5*Sqr(Lambda6) +
      14200*Lambda5*Sqr(Lambda7) + 2000*Conj(Lambda5)*Sqr(Conj(Lambda6))));
   const double beta_Lambda2_2 = Re(-0.1*twoLoop*Conj(Lambda7)*(125*Lambda4*
      Lambda6 + 3120*Lambda2*Lambda7 + 630*Lambda3*Lambda7 + 670*Lambda4*
      Lambda7 + 360*Lambda7*traceYdAdjYd + 120*Lambda7*traceYeAdjYe + 360*
      Lambda7*traceYuAdjYu + 85*Conj(Lambda5)*Conj(Lambda6) + 355*Conj(Lambda5)
      *Conj(Lambda7) - 108*Lambda7*Sqr(g1) + 180*Lambda7*Sqr(g1d) + 60*Lambda7*
      Sqr(g1dp) - 540*Lambda7*Sqr(g2) + 180*Lambda7*Sqr(g2u) + 60*Lambda7*Sqr(
      g2up)));

   beta_Lambda2 = beta_Lambda2_1 + beta_Lambda2_2;


   return beta_Lambda2;
}

/**
 * Calculates the 3-loop beta function of Lambda2.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

/**
 * Calculates the 4-loop beta function of Lambda2.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

/**
 * Calculates the 5-loop beta function of Lambda2.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

} // namespace flexiblesusy
