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

// File generated at Sun 4 Aug 2019 18:57:52

#include "THDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda1.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(0.005*oneOver16PiSqr*(400*Lambda3*Lambda4 + 2400*Lambda1*
      traceYdAdjYd - 1200*traceYdAdjYdYdAdjYd + 800*Lambda1*traceYeAdjYe - 400*
      traceYeAdjYeYeAdjYe + 200*AbsSqr(Lambda5) + 2400*AbsSqr(Lambda6) + 27*
      Quad(g1) + 225*Quad(g2) - 360*Lambda1*Sqr(g1) - 1800*Lambda1*Sqr(g2) + 90
      *Sqr(g1)*Sqr(g2) + 4800*Sqr(Lambda1) + 400*Sqr(Lambda3) + 200*Sqr(Lambda4
      )));


   return beta_Lambda1;
}

/**
 * Calculates the 2-loop beta function of Lambda1.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_2_loop(const Susy_traces& susy_traces) const
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

   beta_Lambda1 = Re(0.0005*twoLoop*(-40000*Lambda1*Lambda3*Lambda4 - 17000*
      Lambda5*Lambda6*Lambda7 - 6000*Lambda1*traceYdAdjYdYdAdjYd + 60000*
      traceYdAdjYdYdAdjYdYdAdjYd - 18000*Lambda1*traceYdAdjYuYuAdjYd + 12000*
      traceYdAdjYuYuAdjYdYdAdjYd - 2000*Lambda1*traceYeAdjYeYeAdjYe + 20000*
      traceYeAdjYeYeAdjYeYeAdjYe - 24000*Lambda3*Lambda4*traceYuAdjYu - 28000*
      Lambda1*AbsSqr(Lambda5) - 40000*Lambda3*AbsSqr(Lambda5) - 44000*Lambda4*
      AbsSqr(Lambda5) - 12000*traceYuAdjYu*AbsSqr(Lambda5) - 624000*Lambda1*
      AbsSqr(Lambda6) - 126000*Lambda3*AbsSqr(Lambda6) - 134000*Lambda4*AbsSqr(
      Lambda6) - 72000*traceYdAdjYd*AbsSqr(Lambda6) - 24000*traceYeAdjYe*AbsSqr
      (Lambda6) - 72000*traceYuAdjYu*AbsSqr(Lambda6) + 12000*Lambda1*AbsSqr(
      Lambda7) - 36000*Lambda3*AbsSqr(Lambda7) - 28000*Lambda4*AbsSqr(Lambda7)
      + 6000*Lambda2*Lambda7*Conj(Lambda6) - 33000*Lambda3*Lambda7*Conj(Lambda6
      ) - 25000*Lambda4*Lambda7*Conj(Lambda6) + 6000*Lambda2*Lambda6*Conj(
      Lambda7) - 33000*Lambda3*Lambda6*Conj(Lambda7) - 25000*Lambda4*Lambda6*
      Conj(Lambda7) - 17000*Conj(Lambda5)*Conj(Lambda6)*Conj(Lambda7) - 624000*
      Cube(Lambda1) - 16000*Cube(Lambda3) - 12000*Cube(Lambda4) - 3537*Power6(
      g1) + 36375*Power6(g2) + 19530*Lambda1*Quad(g1) + 1800*Lambda3*Quad(g1) +
      900*Lambda4*Quad(g1) + 900*traceYdAdjYd*Quad(g1) - 4500*traceYeAdjYe*Quad
      (g1) - 12750*Lambda1*Quad(g2) + 15000*Lambda3*Quad(g2) + 7500*Lambda4*
      Quad(g2) - 4500*traceYdAdjYd*Quad(g2) - 1500*traceYeAdjYe*Quad(g2) + 4800
      *Lambda3*Lambda4*Sqr(g1) + 5000*Lambda1*traceYdAdjYd*Sqr(g1) + 1600*
      traceYdAdjYdYdAdjYd*Sqr(g1) + 15000*Lambda1*traceYeAdjYe*Sqr(g1) - 4800*
      traceYeAdjYeYeAdjYe*Sqr(g1) - 1200*AbsSqr(Lambda5)*Sqr(g1) + 21600*AbsSqr
      (Lambda6)*Sqr(g1) - 7575*Quad(g2)*Sqr(g1) + 24000*Lambda3*Lambda4*Sqr(g2)
      + 45000*Lambda1*traceYdAdjYd*Sqr(g2) + 15000*Lambda1*traceYeAdjYe*Sqr(g2)
      + 108000*AbsSqr(Lambda6)*Sqr(g2) - 8595*Quad(g1)*Sqr(g2) + 11700*Lambda1*
      Sqr(g1)*Sqr(g2) + 3000*Lambda4*Sqr(g1)*Sqr(g2) + 5400*traceYdAdjYd*Sqr(g1
      )*Sqr(g2) + 6600*traceYeAdjYe*Sqr(g1)*Sqr(g2) + 160000*Lambda1*
      traceYdAdjYd*Sqr(g3) - 64000*traceYdAdjYdYdAdjYd*Sqr(g3) - 288000*
      traceYdAdjYd*Sqr(Lambda1) - 96000*traceYeAdjYe*Sqr(Lambda1) + 43200*Sqr(
      g1)*Sqr(Lambda1) + 216000*Sqr(g2)*Sqr(Lambda1) - 40000*Lambda1*Sqr(
      Lambda3) - 24000*Lambda4*Sqr(Lambda3) - 24000*traceYuAdjYu*Sqr(Lambda3) +
      4800*Sqr(g1)*Sqr(Lambda3) + 24000*Sqr(g2)*Sqr(Lambda3) - 24000*Lambda1*
      Sqr(Lambda4) - 32000*Lambda3*Sqr(Lambda4) - 12000*traceYuAdjYu*Sqr(
      Lambda4) + 2400*Sqr(g1)*Sqr(Lambda4) + 6000*Sqr(g2)*Sqr(Lambda4) - 71000*
      Lambda5*Sqr(Lambda6) - 10000*Lambda5*Sqr(Lambda7) - 71000*Conj(Lambda5)*
      Sqr(Conj(Lambda6)) - 10000*Conj(Lambda5)*Sqr(Conj(Lambda7))));


   return beta_Lambda1;
}

/**
 * Calculates the 3-loop beta function of Lambda1.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_3_loop(const Susy_traces& susy_traces) const
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
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_4_loop(const Susy_traces& susy_traces) const
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
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

} // namespace flexiblesusy
