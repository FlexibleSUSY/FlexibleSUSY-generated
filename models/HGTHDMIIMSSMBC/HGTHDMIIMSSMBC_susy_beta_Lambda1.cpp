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

// File generated at Sun 4 Aug 2019 19:38:20

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda1.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(0.005*oneOver16PiSqr*(400*Lambda3*Lambda4 + 2400*Lambda1*
      traceYdAdjYd - 1200*traceYdAdjYdYdAdjYd + 800*Lambda1*traceYeAdjYe - 400*
      traceYeAdjYeYeAdjYe + 200*AbsSqr(Lambda5) + 2400*AbsSqr(Lambda6) + 27*
      Quad(g1) - 500*Quad(g1d) - 100*Quad(g1dp) + 225*Quad(g2) - 360*Lambda1*
      Sqr(g1) + 1200*Lambda1*Sqr(g1d) + 400*Lambda1*Sqr(g1dp) - 200*Sqr(g1d)*
      Sqr(g1dp) - 1800*Lambda1*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 4800*Sqr(Lambda1)
      + 400*Sqr(Lambda3) + 200*Sqr(Lambda4)));


   return beta_Lambda1;
}

/**
 * Calculates the 2-loop beta function of Lambda1.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYdAdjYdYdAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe = TRACE_STRUCT.
      traceYeAdjYeYeAdjYeYeAdjYe;


   double beta_Lambda1;

   const double beta_Lambda1_1 = Re(-0.0025*twoLoop*(8000*Lambda1*Lambda3*
      Lambda4 + 3400*Lambda5*Lambda6*Lambda7 + 1200*Lambda1*traceYdAdjYdYdAdjYd
       - 12000*traceYdAdjYdYdAdjYdYdAdjYd + 3600*Lambda1*traceYdAdjYuYuAdjYd -
      2400*traceYdAdjYuYuAdjYdYdAdjYd + 400*Lambda1*traceYeAdjYeYeAdjYe - 4000*
      traceYeAdjYeYeAdjYeYeAdjYe + 4800*Lambda3*Lambda4*traceYuAdjYu + 5600*
      Lambda1*AbsSqr(Lambda5) + 8000*Lambda3*AbsSqr(Lambda5) + 8800*Lambda4*
      AbsSqr(Lambda5) + 2400*traceYuAdjYu*AbsSqr(Lambda5) + 124800*Cube(Lambda1
      ) + 3200*Cube(Lambda3) + 2400*Cube(Lambda4) + 765*Power6(g1) - 4700*
      Power6(g1d) - 500*Power6(g1dp) - 4875*Power6(g2) - 4146*Lambda1*Quad(g1)
      - 360*Lambda3*Quad(g1) - 180*Lambda4*Quad(g1) - 180*traceYdAdjYd*Quad(g1)
      + 900*traceYeAdjYe*Quad(g1) + 500*Lambda1*Quad(g1d) + 100*Lambda1*Quad(
      g1dp) - 3450*Lambda1*Quad(g2) - 3000*Lambda3*Quad(g2) - 1500*Lambda4*Quad
      (g2) + 900*traceYdAdjYd*Quad(g2) + 300*traceYeAdjYe*Quad(g2) - 960*
      Lambda3*Lambda4*Sqr(g1) - 1000*Lambda1*traceYdAdjYd*Sqr(g1) - 320*
      traceYdAdjYdYdAdjYd*Sqr(g1) - 3000*Lambda1*traceYeAdjYe*Sqr(g1) + 960*
      traceYeAdjYeYeAdjYe*Sqr(g1) + 240*AbsSqr(Lambda5)*Sqr(g1) - 4320*AbsSqr(
      Lambda6)*Sqr(g1) + 1995*Quad(g2)*Sqr(g1) + 7200*AbsSqr(Lambda6)*Sqr(g1d)
      + 54*Quad(g1)*Sqr(g1d) - 1700*Quad(g1dp)*Sqr(g1d) + 7650*Quad(g2)*Sqr(g1d
      ) - 900*Lambda1*Sqr(g1)*Sqr(g1d) + 2400*AbsSqr(Lambda6)*Sqr(g1dp) + 18*
      Quad(g1)*Sqr(g1dp) - 1100*Quad(g1d)*Sqr(g1dp) + 150*Quad(g2)*Sqr(g1dp) -
      300*Lambda1*Sqr(g1)*Sqr(g1dp) + 200*Lambda1*Sqr(g1d)*Sqr(g1dp) - 4800*
      Lambda3*Lambda4*Sqr(g2) - 9000*Lambda1*traceYdAdjYd*Sqr(g2) - 3000*
      Lambda1*traceYeAdjYe*Sqr(g2) - 21600*AbsSqr(Lambda6)*Sqr(g2) + 1815*Quad(
      g1)*Sqr(g2) + 4000*Quad(g1d)*Sqr(g2) - 2340*Lambda1*Sqr(g1)*Sqr(g2) - 600
      *Lambda4*Sqr(g1)*Sqr(g2) - 1080*traceYdAdjYd*Sqr(g1)*Sqr(g2) - 1320*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) - 16500*Lambda1*Sqr(g1d)*Sqr(g2) - 1260*Sqr(
      g1)*Sqr(g1d)*Sqr(g2) - 1500*Lambda1*Sqr(g1dp)*Sqr(g2) + 60*Sqr(g1)*Sqr(
      g1dp)*Sqr(g2) + 800*Sqr(g1d)*Sqr(g1dp)*Sqr(g2) + 2400*Lambda3*Lambda4*Sqr
      (g2u) + 1200*AbsSqr(Lambda5)*Sqr(g2u) - 1000*Quad(g1d)*Sqr(g2u) + 1800*
      Lambda1*Sqr(g1d)*Sqr(g2u) - 200*Sqr(g1d)*Sqr(g1dp)*Sqr(g2u) + 800*Lambda3
      *Lambda4*Sqr(g2up) + 400*AbsSqr(Lambda5)*Sqr(g2up) - 200*Quad(g1dp)*Sqr(
      g2up) + 600*Lambda1*Sqr(g1dp)*Sqr(g2up) - 200*Sqr(g1d)*Sqr(g1dp)*Sqr(g2up
      ) - 32000*Lambda1*traceYdAdjYd*Sqr(g3) + 12800*traceYdAdjYdYdAdjYd*Sqr(g3
      ) + 57600*traceYdAdjYd*Sqr(Lambda1) + 19200*traceYeAdjYe*Sqr(Lambda1) -
      8640*Sqr(g1)*Sqr(Lambda1) + 28800*Sqr(g1d)*Sqr(Lambda1) + 9600*Sqr(g1dp)*
      Sqr(Lambda1) - 43200*Sqr(g2)*Sqr(Lambda1) + 8000*Lambda1*Sqr(Lambda3) +
      4800*Lambda4*Sqr(Lambda3) + 4800*traceYuAdjYu*Sqr(Lambda3) - 960*Sqr(g1)*
      Sqr(Lambda3) - 4800*Sqr(g2)*Sqr(Lambda3) + 2400*Sqr(g2u)*Sqr(Lambda3) +
      800*Sqr(g2up)*Sqr(Lambda3) + 4800*Lambda1*Sqr(Lambda4) + 6400*Lambda3*Sqr
      (Lambda4) + 2400*traceYuAdjYu*Sqr(Lambda4) - 480*Sqr(g1)*Sqr(Lambda4) -
      1200*Sqr(g2)*Sqr(Lambda4) + 1200*Sqr(g2u)*Sqr(Lambda4) + 400*Sqr(g2up)*
      Sqr(Lambda4) + 14200*Lambda5*Sqr(Lambda6) + 2000*Lambda5*Sqr(Lambda7)));
   const double beta_Lambda1_2 = Re(-0.5*twoLoop*(624*Lambda1*AbsSqr(Lambda6) +
      126*Lambda3*AbsSqr(Lambda6) + 134*Lambda4*AbsSqr(Lambda6) + 72*
      traceYdAdjYd*AbsSqr(Lambda6) + 24*traceYeAdjYe*AbsSqr(Lambda6) + 72*
      traceYuAdjYu*AbsSqr(Lambda6) - 12*Lambda1*AbsSqr(Lambda7) + 36*Lambda3*
      AbsSqr(Lambda7) + 28*Lambda4*AbsSqr(Lambda7) - 6*Lambda2*Lambda7*Conj(
      Lambda6) + 33*Lambda3*Lambda7*Conj(Lambda6) + 25*Lambda4*Lambda7*Conj(
      Lambda6) - 6*Lambda2*Lambda6*Conj(Lambda7) + 33*Lambda3*Lambda6*Conj(
      Lambda7) + 25*Lambda4*Lambda6*Conj(Lambda7) + 17*Conj(Lambda5)*Conj(
      Lambda6)*Conj(Lambda7) + 36*AbsSqr(Lambda6)*Sqr(g2u) + 12*AbsSqr(Lambda6)
      *Sqr(g2up) + 71*Conj(Lambda5)*Sqr(Conj(Lambda6)) + 10*Conj(Lambda5)*Sqr(
      Conj(Lambda7))));

   beta_Lambda1 = beta_Lambda1_1 + beta_Lambda1_2;


   return beta_Lambda1;
}

/**
 * Calculates the 3-loop beta function of Lambda1.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

/**
 * Calculates the 4-loop beta function of Lambda1.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

/**
 * Calculates the 5-loop beta function of Lambda1.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

} // namespace flexiblesusy
