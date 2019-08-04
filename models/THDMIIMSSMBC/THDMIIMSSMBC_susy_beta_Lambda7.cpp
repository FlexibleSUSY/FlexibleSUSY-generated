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

// File generated at Sun 4 Aug 2019 18:57:51

#include "THDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda7.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(0.2*oneOver16PiSqr*(30*Lambda3*Lambda6 + 20*Lambda4*
      Lambda6 + 120*Lambda2*Lambda7 + 30*Lambda3*Lambda7 + 40*Lambda4*Lambda7 +
      15*Lambda7*traceYdAdjYd + 5*Lambda7*traceYeAdjYe + 45*Lambda7*
      traceYuAdjYu + 10*Conj(Lambda5)*Conj(Lambda6) + 50*Conj(Lambda5)*Conj(
      Lambda7) - 9*Lambda7*Sqr(g1) - 45*Lambda7*Sqr(g2)));


   return beta_Lambda7;
}

/**
 * Calculates the 2-loop beta function of Lambda7.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(0.005*twoLoop*(1200*Lambda1*Lambda2*Lambda6 - 6600*Lambda1
      *Lambda3*Lambda6 - 6600*Lambda2*Lambda3*Lambda6 - 5000*Lambda1*Lambda4*
      Lambda6 - 5000*Lambda2*Lambda4*Lambda6 - 10600*Lambda3*Lambda4*Lambda6 -
      7200*Lambda1*Lambda3*Lambda7 - 25200*Lambda2*Lambda3*Lambda7 - 5600*
      Lambda1*Lambda4*Lambda7 - 26800*Lambda2*Lambda4*Lambda7 - 13000*Lambda3*
      Lambda4*Lambda7 - 7200*Lambda3*Lambda6*traceYdAdjYd - 4800*Lambda4*
      Lambda6*traceYdAdjYd - 3600*Lambda3*Lambda7*traceYdAdjYd - 4800*Lambda4*
      Lambda7*traceYdAdjYd - 1350*Lambda7*traceYdAdjYdYdAdjYd - 4200*Lambda7*
      traceYdAdjYuYuAdjYd - 2400*Lambda3*Lambda6*traceYeAdjYe - 1600*Lambda4*
      Lambda6*traceYeAdjYe - 1200*Lambda3*Lambda7*traceYeAdjYe - 1600*Lambda4*
      Lambda7*traceYeAdjYe - 450*Lambda7*traceYeAdjYeYeAdjYe - 28800*Lambda2*
      Lambda7*traceYuAdjYu - 3600*Lambda3*Lambda7*traceYuAdjYu - 4800*Lambda4*
      Lambda7*traceYuAdjYu - 450*Lambda7*traceYuAdjYuYuAdjYu - 8100*Lambda6*
      AbsSqr(Lambda5) - 6900*Lambda7*AbsSqr(Lambda5) - 2200*Lambda7*AbsSqr(
      Lambda6) - 16800*Lambda6*AbsSqr(Lambda7) - 3400*Lambda1*Conj(Lambda5)*
      Conj(Lambda6) - 3400*Lambda2*Conj(Lambda5)*Conj(Lambda6) - 7400*Lambda3*
      Conj(Lambda5)*Conj(Lambda6) - 8200*Lambda4*Conj(Lambda5)*Conj(Lambda6) -
      2400*traceYdAdjYd*Conj(Lambda5)*Conj(Lambda6) - 800*traceYeAdjYe*Conj(
      Lambda5)*Conj(Lambda6) - 4000*Lambda1*Conj(Lambda5)*Conj(Lambda7) - 28400
      *Lambda2*Conj(Lambda5)*Conj(Lambda7) - 13800*Lambda3*Conj(Lambda5)*Conj(
      Lambda7) - 14600*Lambda4*Conj(Lambda5)*Conj(Lambda7) - 6000*traceYdAdjYd*
      Conj(Lambda5)*Conj(Lambda7) - 2000*traceYeAdjYe*Conj(Lambda5)*Conj(
      Lambda7) - 6000*traceYuAdjYu*Conj(Lambda5)*Conj(Lambda7) + 270*Lambda6*
      Quad(g1) + 1683*Lambda7*Quad(g1) + 2250*Lambda6*Quad(g2) - 3525*Lambda7*
      Quad(g2) + 1440*Lambda3*Lambda6*Sqr(g1) + 960*Lambda4*Lambda6*Sqr(g1) +
      4320*Lambda2*Lambda7*Sqr(g1) + 720*Lambda3*Lambda7*Sqr(g1) + 1200*Lambda4
      *Lambda7*Sqr(g1) + 125*Lambda7*traceYdAdjYd*Sqr(g1) + 375*Lambda7*
      traceYeAdjYe*Sqr(g1) + 1275*Lambda7*traceYuAdjYu*Sqr(g1) - 240*Conj(
      Lambda5)*Conj(Lambda6)*Sqr(g1) + 2400*Conj(Lambda5)*Conj(Lambda7)*Sqr(g1)
      + 7200*Lambda3*Lambda6*Sqr(g2) + 3600*Lambda4*Lambda6*Sqr(g2) + 21600*
      Lambda2*Lambda7*Sqr(g2) + 3600*Lambda3*Lambda7*Sqr(g2) + 7200*Lambda4*
      Lambda7*Sqr(g2) + 1125*Lambda7*traceYdAdjYd*Sqr(g2) + 375*Lambda7*
      traceYeAdjYe*Sqr(g2) + 3375*Lambda7*traceYuAdjYu*Sqr(g2) + 10800*Conj(
      Lambda5)*Conj(Lambda7)*Sqr(g2) + 300*Lambda6*Sqr(g1)*Sqr(g2) + 870*
      Lambda7*Sqr(g1)*Sqr(g2) + 4000*Lambda7*traceYdAdjYd*Sqr(g3) + 12000*
      Lambda7*traceYuAdjYu*Sqr(g3) + 1200*Lambda7*Sqr(Lambda1) - 62400*Lambda7*
      Sqr(Lambda2) - 6900*Lambda6*Sqr(Lambda3) - 6100*Lambda7*Sqr(Lambda3) -
      6500*Lambda6*Sqr(Lambda4) - 6500*Lambda7*Sqr(Lambda4) - 8400*Conj(Lambda6
      )*Sqr(Lambda6) - 4400*Conj(Lambda7)*Sqr(Lambda6) - 8400*Conj(Lambda6)*Sqr
      (Lambda7) - 22200*Conj(Lambda7)*Sqr(Lambda7)));


   return beta_Lambda7;
}

/**
 * Calculates the 3-loop beta function of Lambda7.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

/**
 * Calculates the 4-loop beta function of Lambda7.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

/**
 * Calculates the 5-loop beta function of Lambda7.
 *
 * @return 5-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return beta_Lambda7;
}

} // namespace flexiblesusy
