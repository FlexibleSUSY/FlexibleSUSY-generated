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

// File generated at Tue 22 Jan 2019 16:22:11

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda6.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(0.1*oneOver16PiSqr*(240*Lambda1*Lambda6 + 60*Lambda3*
      Lambda6 + 80*Lambda4*Lambda6 + 60*Lambda3*Lambda7 + 40*Lambda4*Lambda7 +
      90*Lambda6*traceYdAdjYd + 30*Lambda6*traceYeAdjYe + 30*Lambda6*
      traceYuAdjYu + 100*Conj(Lambda5)*Conj(Lambda6) + 20*Conj(Lambda5)*Conj(
      Lambda7) - 18*Lambda6*Sqr(g1) + 45*Lambda6*Sqr(g1d) + 15*Lambda6*Sqr(g1dp
      ) - 90*Lambda6*Sqr(g2) + 15*Lambda6*Sqr(g2u) + 5*Lambda6*Sqr(g2up)));


   return beta_Lambda6;
}

/**
 * Calculates the 2-loop beta function of Lambda6.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda6;

   const double beta_Lambda6_1 = Re(0.0025*twoLoop*(3200*g1d*g1dp*g2u*g2up*
      Lambda6 - 50400*Lambda1*Lambda3*Lambda6 - 14400*Lambda2*Lambda3*Lambda6 -
      53600*Lambda1*Lambda4*Lambda6 - 11200*Lambda2*Lambda4*Lambda6 - 26000*
      Lambda3*Lambda4*Lambda6 + 2400*Lambda1*Lambda2*Lambda7 - 13200*Lambda1*
      Lambda3*Lambda7 - 13200*Lambda2*Lambda3*Lambda7 - 10000*Lambda1*Lambda4*
      Lambda7 - 10000*Lambda2*Lambda4*Lambda7 - 21200*Lambda3*Lambda4*Lambda7 -
      57600*Lambda1*Lambda6*traceYdAdjYd - 7200*Lambda3*Lambda6*traceYdAdjYd -
      9600*Lambda4*Lambda6*traceYdAdjYd - 900*Lambda6*traceYdAdjYdYdAdjYd -
      8400*Lambda6*traceYdAdjYuYuAdjYd - 19200*Lambda1*Lambda6*traceYeAdjYe -
      2400*Lambda3*Lambda6*traceYeAdjYe - 3200*Lambda4*Lambda6*traceYeAdjYe -
      300*Lambda6*traceYeAdjYeYeAdjYe - 7200*Lambda3*Lambda6*traceYuAdjYu -
      9600*Lambda4*Lambda6*traceYuAdjYu - 14400*Lambda3*Lambda7*traceYuAdjYu -
      9600*Lambda4*Lambda7*traceYuAdjYu - 2700*Lambda6*traceYuAdjYuYuAdjYu -
      13800*Lambda6*AbsSqr(Lambda5) - 16200*Lambda7*AbsSqr(Lambda5) - 33600*
      Lambda7*AbsSqr(Lambda6) - 56800*Lambda1*Conj(Lambda5)*Conj(Lambda6) -
      8000*Lambda2*Conj(Lambda5)*Conj(Lambda6) - 27600*Lambda3*Conj(Lambda5)*
      Conj(Lambda6) + 3606*Lambda6*Quad(g1) + 540*Lambda7*Quad(g1) - 2375*
      Lambda6*Quad(g1d) - 475*Lambda6*Quad(g1dp) - 1050*Lambda6*Quad(g2) + 4500
      *Lambda7*Quad(g2) - 1125*Lambda6*Quad(g2u) - 225*Lambda6*Quad(g2up) +
      8640*Lambda1*Lambda6*Sqr(g1) + 1440*Lambda3*Lambda6*Sqr(g1) + 2400*
      Lambda4*Lambda6*Sqr(g1) + 2880*Lambda3*Lambda7*Sqr(g1) + 1920*Lambda4*
      Lambda7*Sqr(g1) + 750*Lambda6*traceYdAdjYd*Sqr(g1) + 2250*Lambda6*
      traceYeAdjYe*Sqr(g1) + 850*Lambda6*traceYuAdjYu*Sqr(g1) + 4800*Conj(
      Lambda5)*Conj(Lambda6)*Sqr(g1) - 28800*Lambda1*Lambda6*Sqr(g1d) - 3600*
      Lambda3*Lambda6*Sqr(g1d) - 4800*Lambda4*Lambda6*Sqr(g1d) - 6000*Conj(
      Lambda5)*Conj(Lambda6)*Sqr(g1d) + 675*Lambda6*Sqr(g1)*Sqr(g1d) - 9600*
      Lambda1*Lambda6*Sqr(g1dp) - 1200*Lambda3*Lambda6*Sqr(g1dp) - 1600*Lambda4
      *Lambda6*Sqr(g1dp) - 2000*Conj(Lambda5)*Conj(Lambda6)*Sqr(g1dp) + 225*
      Lambda6*Sqr(g1)*Sqr(g1dp) - 950*Lambda6*Sqr(g1d)*Sqr(g1dp) + 43200*
      Lambda1*Lambda6*Sqr(g2) + 7200*Lambda3*Lambda6*Sqr(g2) + 14400*Lambda4*
      Lambda6*Sqr(g2) + 14400*Lambda3*Lambda7*Sqr(g2) + 7200*Lambda4*Lambda7*
      Sqr(g2) + 6750*Lambda6*traceYdAdjYd*Sqr(g2) + 2250*Lambda6*traceYeAdjYe*
      Sqr(g2) + 2250*Lambda6*traceYuAdjYu*Sqr(g2) + 21600*Conj(Lambda5)*Conj(
      Lambda6)*Sqr(g2) + 1740*Lambda6*Sqr(g1)*Sqr(g2) + 600*Lambda7*Sqr(g1)*Sqr
      (g2) + 12375*Lambda6*Sqr(g1d)*Sqr(g2) + 1125*Lambda6*Sqr(g1dp)*Sqr(g2) -
      3600*Lambda3*Lambda6*Sqr(g2u) - 4800*Lambda4*Lambda6*Sqr(g2u) - 7200*
      Lambda3*Lambda7*Sqr(g2u) - 4800*Lambda4*Lambda7*Sqr(g2u) - 6000*Conj(
      Lambda5)*Conj(Lambda6)*Sqr(g2u) + 225*Lambda6*Sqr(g1)*Sqr(g2u) - 2600*
      Lambda6*Sqr(g1d)*Sqr(g2u) + 4125*Lambda6*Sqr(g2)*Sqr(g2u) - 1200*Lambda3*
      Lambda6*Sqr(g2up) - 1600*Lambda4*Lambda6*Sqr(g2up) - 2400*Lambda3*Lambda7
      *Sqr(g2up) - 1600*Lambda4*Lambda7*Sqr(g2up) - 2000*Conj(Lambda5)*Conj(
      Lambda6)*Sqr(g2up) + 75*Lambda6*Sqr(g1)*Sqr(g2up) + 200*Lambda6*Sqr(g1dp)
      *Sqr(g2up) + 375*Lambda6*Sqr(g2)*Sqr(g2up) - 450*Lambda6*Sqr(g2u)*Sqr(
      g2up) + 24000*Lambda6*traceYdAdjYd*Sqr(g3) + 8000*Lambda6*traceYuAdjYu*
      Sqr(g3) - 124800*Lambda6*Sqr(Lambda1) + 2400*Lambda6*Sqr(Lambda2) - 12200
      *Lambda6*Sqr(Lambda3) - 13800*Lambda7*Sqr(Lambda3) - 13000*Lambda6*Sqr(
      Lambda4) - 13000*Lambda7*Sqr(Lambda4) - 44400*Conj(Lambda6)*Sqr(Lambda6)
      - 8800*Conj(Lambda6)*Sqr(Lambda7)));
   const double beta_Lambda6_2 = Re(-0.2*twoLoop*(55*Lambda6*AbsSqr(Lambda7) +
      365*Lambda4*Conj(Lambda5)*Conj(Lambda6) + 150*traceYdAdjYd*Conj(Lambda5)*
      Conj(Lambda6) + 50*traceYeAdjYe*Conj(Lambda5)*Conj(Lambda6) + 150*
      traceYuAdjYu*Conj(Lambda5)*Conj(Lambda6) + 85*Lambda1*Conj(Lambda5)*Conj(
      Lambda7) + 85*Lambda2*Conj(Lambda5)*Conj(Lambda7) + 185*Lambda3*Conj(
      Lambda5)*Conj(Lambda7) + 205*Lambda4*Conj(Lambda5)*Conj(Lambda7) + 60*
      traceYuAdjYu*Conj(Lambda5)*Conj(Lambda7) + 6*Conj(Lambda5)*Conj(Lambda7)*
      Sqr(g1) + 30*Conj(Lambda5)*Conj(Lambda7)*Sqr(g2u) + 10*Conj(Lambda5)*Conj
      (Lambda7)*Sqr(g2up) + 210*Conj(Lambda7)*Sqr(Lambda6) + 210*Conj(Lambda7)*
      Sqr(Lambda7)));

   beta_Lambda6 = beta_Lambda6_1 + beta_Lambda6_2;


   return beta_Lambda6;
}

/**
 * Calculates the 3-loop beta function of Lambda6.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

/**
 * Calculates the 4-loop beta function of Lambda6.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

/**
 * Calculates the 5-loop beta function of Lambda6.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda6_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

} // namespace flexiblesusy
