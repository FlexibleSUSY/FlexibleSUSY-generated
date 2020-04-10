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

// File generated at Fri 10 Apr 2020 19:22:06

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda3.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda3;

   beta_Lambda3 = Re(0.01*oneOver16PiSqr*(1200*Lambda1*Lambda3 + 1200*Lambda2*
      Lambda3 + 400*Lambda1*Lambda4 + 400*Lambda2*Lambda4 + 600*Lambda3*
      traceYdAdjYd - 1200*traceYdAdjYuYuAdjYd + 200*Lambda3*traceYeAdjYe + 600*
      Lambda3*traceYuAdjYu + 200*AbsSqr(Lambda5) + 400*AbsSqr(Lambda6) + 400*
      AbsSqr(Lambda7) + 800*Lambda7*Conj(Lambda6) + 800*Lambda6*Conj(Lambda7) +
      27*Quad(g1) + 225*Quad(g2) - 180*Lambda3*Sqr(g1) - 900*Lambda3*Sqr(g2) -
      90*Sqr(g1)*Sqr(g2) + 400*Sqr(Lambda3) + 200*Sqr(Lambda4)));


   return beta_Lambda3;
}

/**
 * Calculates the 2-loop beta function of Lambda3.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_2_loop(const Susy_traces& susy_traces) const
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

   const double beta_Lambda3_1 = Re(-0.005*twoLoop*(6400*Lambda1*Lambda3*
      Lambda4 + 6400*Lambda2*Lambda3*Lambda4 + 6600*Lambda5*Lambda6*Lambda7 +
      14400*Lambda1*Lambda3*traceYdAdjYd + 4800*Lambda1*Lambda4*traceYdAdjYd +
      2700*Lambda3*traceYdAdjYdYdAdjYd - 2400*traceYdAdjYdYdAdjYuYuAdjYd - 3000
      *Lambda3*traceYdAdjYuYuAdjYd - 4800*traceYdAdjYuYuAdjYdYdAdjYd - 7200*
      traceYdAdjYuYuAdjYuYuAdjYd + 4800*Lambda1*Lambda3*traceYeAdjYe + 1600*
      Lambda1*Lambda4*traceYeAdjYe + 900*Lambda3*traceYeAdjYeYeAdjYe + 14400*
      Lambda2*Lambda3*traceYuAdjYu + 4800*Lambda2*Lambda4*traceYuAdjYu + 2700*
      Lambda3*traceYuAdjYuYuAdjYu + 7200*Lambda1*AbsSqr(Lambda5) + 7200*Lambda2
      *AbsSqr(Lambda5) + 3600*Lambda3*AbsSqr(Lambda5) + 8800*Lambda4*AbsSqr(
      Lambda5) + 1200*traceYdAdjYd*AbsSqr(Lambda5) + 400*traceYeAdjYe*AbsSqr(
      Lambda5) + 1200*traceYuAdjYu*AbsSqr(Lambda5) + 23600*Lambda1*AbsSqr(
      Lambda6) + 8800*Lambda2*AbsSqr(Lambda6) + 11400*Lambda3*AbsSqr(Lambda6) +
      13000*Lambda4*AbsSqr(Lambda6) + 8200*Lambda1*Lambda7*Conj(Lambda6) + 2400
      *Cube(Lambda3) + 2400*Cube(Lambda4) + 765*Power6(g1) - 6475*Power6(g2) -
      540*Lambda1*Quad(g1) - 540*Lambda2*Quad(g1) - 1893*Lambda3*Quad(g1) - 180
      *Lambda4*Quad(g1) - 90*traceYdAdjYd*Quad(g1) + 450*traceYeAdjYe*Quad(g1)
      + 342*traceYuAdjYu*Quad(g1) - 4500*Lambda1*Quad(g2) - 4500*Lambda2*Quad(
      g2) + 1775*Lambda3*Quad(g2) - 1500*Lambda4*Quad(g2) + 450*traceYdAdjYd*
      Quad(g2) + 150*traceYeAdjYe*Quad(g2) + 450*traceYuAdjYu*Quad(g2) - 2880*
      Lambda1*Lambda3*Sqr(g1) - 2880*Lambda2*Lambda3*Sqr(g1) - 960*Lambda1*
      Lambda4*Sqr(g1) - 960*Lambda2*Lambda4*Sqr(g1) - 250*Lambda3*traceYdAdjYd*
      Sqr(g1) + 160*traceYdAdjYuYuAdjYd*Sqr(g1) - 750*Lambda3*traceYeAdjYe*Sqr(
      g1) - 850*Lambda3*traceYuAdjYu*Sqr(g1) - 480*AbsSqr(Lambda5)*Sqr(g1) -
      240*AbsSqr(Lambda6)*Sqr(g1) - 1920*Lambda7*Conj(Lambda6)*Sqr(g1) - 325*
      Quad(g2)*Sqr(g1) - 14400*Lambda1*Lambda3*Sqr(g2) - 14400*Lambda2*Lambda3*
      Sqr(g2) - 7200*Lambda1*Lambda4*Sqr(g2) - 7200*Lambda2*Lambda4*Sqr(g2) +
      2400*Lambda3*Lambda4*Sqr(g2) - 2250*Lambda3*traceYdAdjYd*Sqr(g2) - 750*
      Lambda3*traceYeAdjYe*Sqr(g2) - 2250*Lambda3*traceYuAdjYu*Sqr(g2) - 10800*
      Lambda7*Conj(Lambda6)*Sqr(g2) - 1005*Quad(g1)*Sqr(g2) + 600*Lambda1*Sqr(
      g1)*Sqr(g2) + 600*Lambda2*Sqr(g1)*Sqr(g2) - 330*Lambda3*Sqr(g1)*Sqr(g2) +
      360*Lambda4*Sqr(g1)*Sqr(g2) + 540*traceYdAdjYd*Sqr(g1)*Sqr(g2) + 660*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) + 1260*traceYuAdjYu*Sqr(g1)*Sqr(g2) - 8000*
      Lambda3*traceYdAdjYd*Sqr(g3) + 12800*traceYdAdjYuYuAdjYd*Sqr(g3) - 8000*
      Lambda3*traceYuAdjYu*Sqr(g3) + 12000*Lambda3*Sqr(Lambda1) + 3200*Lambda4*
      Sqr(Lambda1) + 12000*Lambda3*Sqr(Lambda2) + 3200*Lambda4*Sqr(Lambda2) +
      14400*Lambda1*Sqr(Lambda3) + 14400*Lambda2*Sqr(Lambda3) + 800*Lambda4*Sqr
      (Lambda3) + 2400*traceYdAdjYd*Sqr(Lambda3) + 800*traceYeAdjYe*Sqr(Lambda3
      ) + 2400*traceYuAdjYu*Sqr(Lambda3) - 240*Sqr(g1)*Sqr(Lambda3) - 1200*Sqr(
      g2)*Sqr(Lambda3) + 5600*Lambda1*Sqr(Lambda4) + 5600*Lambda2*Sqr(Lambda4)
      + 3200*Lambda3*Sqr(Lambda4) + 1200*traceYdAdjYd*Sqr(Lambda4) + 400*
      traceYeAdjYe*Sqr(Lambda4) + 1200*traceYuAdjYu*Sqr(Lambda4) + 240*Sqr(g1)*
      Sqr(Lambda4) - 1200*Sqr(g2)*Sqr(Lambda4) + 6500*Lambda5*Sqr(Lambda6) +
      6500*Lambda5*Sqr(Lambda7)));
   const double beta_Lambda3_2 = Re(-0.1*twoLoop*(240*traceYdAdjYd*AbsSqr(
      Lambda6) + 80*traceYeAdjYe*AbsSqr(Lambda6) + 440*Lambda1*AbsSqr(Lambda7)
      + 1180*Lambda2*AbsSqr(Lambda7) + 570*Lambda3*AbsSqr(Lambda7) + 650*
      Lambda4*AbsSqr(Lambda7) + 240*traceYuAdjYu*AbsSqr(Lambda7) + 410*Lambda2*
      Lambda7*Conj(Lambda6) + 850*Lambda3*Lambda7*Conj(Lambda6) + 410*Lambda4*
      Lambda7*Conj(Lambda6) + 240*Lambda7*traceYdAdjYd*Conj(Lambda6) + 80*
      Lambda7*traceYeAdjYe*Conj(Lambda6) + 240*Lambda7*traceYuAdjYu*Conj(
      Lambda6) + 410*Lambda1*Lambda6*Conj(Lambda7) + 410*Lambda2*Lambda6*Conj(
      Lambda7) + 850*Lambda3*Lambda6*Conj(Lambda7) + 410*Lambda4*Lambda6*Conj(
      Lambda7) + 240*Lambda6*traceYdAdjYd*Conj(Lambda7) + 80*Lambda6*
      traceYeAdjYe*Conj(Lambda7) + 240*Lambda6*traceYuAdjYu*Conj(Lambda7) + 330
      *Conj(Lambda5)*Conj(Lambda6)*Conj(Lambda7) - 12*AbsSqr(Lambda7)*Sqr(g1) -
      96*Lambda6*Conj(Lambda7)*Sqr(g1) - 540*Lambda6*Conj(Lambda7)*Sqr(g2) +
      325*Conj(Lambda5)*Sqr(Conj(Lambda6)) + 325*Conj(Lambda5)*Sqr(Conj(Lambda7
      ))));

   beta_Lambda3 = beta_Lambda3_1 + beta_Lambda3_2;


   return beta_Lambda3;
}

/**
 * Calculates the 3-loop beta function of Lambda3.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_3_loop(const Susy_traces& susy_traces) const
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
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_4_loop(const Susy_traces& susy_traces) const
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
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

} // namespace flexiblesusy
