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


#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda7.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda7;

   beta_Lambda7 = Re(0.1*(60*Lambda3*Lambda6 + 40*Lambda4*Lambda6 + 240*Lambda2
      *Lambda7 + 60*Lambda3*Lambda7 + 80*Lambda4*Lambda7 + 30*Lambda7*
      traceYdAdjYd + 10*Lambda7*traceYeAdjYe + 90*Lambda7*traceYuAdjYu + 20*
      Conj(Lambda5)*Conj(Lambda6) + 100*Conj(Lambda5)*Conj(Lambda7) - 18*
      Lambda7*Sqr(g1) + 15*Lambda7*Sqr(g1d) + 5*Lambda7*Sqr(g1dp) - 90*Lambda7*
      Sqr(g2) + 45*Lambda7*Sqr(g2u) + 15*Lambda7*Sqr(g2up)));


   return oneLoop * beta_Lambda7;
}

/**
 * Calculates the 2-loop beta function of Lambda7.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda7;

   const double beta_Lambda7_1 = Re(0.0025*(2400*Lambda1*Lambda2*Lambda6 -
      13200*Lambda1*Lambda3*Lambda6 - 13200*Lambda2*Lambda3*Lambda6 - 10000*
      Lambda1*Lambda4*Lambda6 - 10000*Lambda2*Lambda4*Lambda6 - 21200*Lambda3*
      Lambda4*Lambda6 + 3200*g1d*g1dp*g2u*g2up*Lambda7 - 14400*Lambda1*Lambda3*
      Lambda7 - 50400*Lambda2*Lambda3*Lambda7 - 11200*Lambda1*Lambda4*Lambda7 -
      53600*Lambda2*Lambda4*Lambda7 - 26000*Lambda3*Lambda4*Lambda7 - 14400*
      Lambda3*Lambda6*traceYdAdjYd - 9600*Lambda4*Lambda6*traceYdAdjYd - 7200*
      Lambda3*Lambda7*traceYdAdjYd - 9600*Lambda4*Lambda7*traceYdAdjYd - 2700*
      Lambda7*traceYdAdjYdYdAdjYd - 8400*Lambda7*traceYdAdjYuYuAdjYd - 4800*
      Lambda3*Lambda6*traceYeAdjYe - 3200*Lambda4*Lambda6*traceYeAdjYe - 2400*
      Lambda3*Lambda7*traceYeAdjYe - 3200*Lambda4*Lambda7*traceYeAdjYe - 900*
      Lambda7*traceYeAdjYeYeAdjYe - 57600*Lambda2*Lambda7*traceYuAdjYu - 7200*
      Lambda3*Lambda7*traceYuAdjYu - 9600*Lambda4*Lambda7*traceYuAdjYu - 900*
      Lambda7*traceYuAdjYuYuAdjYu - 16200*Lambda6*AbsSqr(Lambda5) - 13800*
      Lambda7*AbsSqr(Lambda5) - 4400*Lambda7*AbsSqr(Lambda6) - 6800*Lambda1*
      Conj(Lambda5)*Conj(Lambda6) - 6800*Lambda2*Conj(Lambda5)*Conj(Lambda6) -
      14800*Lambda3*Conj(Lambda5)*Conj(Lambda6) - 16400*Lambda4*Conj(Lambda5)*
      Conj(Lambda6) - 4800*traceYdAdjYd*Conj(Lambda5)*Conj(Lambda6) + 540*
      Lambda6*Quad(g1) + 3606*Lambda7*Quad(g1) - 1125*Lambda7*Quad(g1d) - 225*
      Lambda7*Quad(g1dp) + 4500*Lambda6*Quad(g2) - 1050*Lambda7*Quad(g2) - 2375
      *Lambda7*Quad(g2u) - 475*Lambda7*Quad(g2up) + 2880*Lambda3*Lambda6*Sqr(g1
      ) + 1920*Lambda4*Lambda6*Sqr(g1) + 8640*Lambda2*Lambda7*Sqr(g1) + 1440*
      Lambda3*Lambda7*Sqr(g1) + 2400*Lambda4*Lambda7*Sqr(g1) + 250*Lambda7*
      traceYdAdjYd*Sqr(g1) + 750*Lambda7*traceYeAdjYe*Sqr(g1) + 2550*Lambda7*
      traceYuAdjYu*Sqr(g1) - 480*Conj(Lambda5)*Conj(Lambda6)*Sqr(g1) - 7200*
      Lambda3*Lambda6*Sqr(g1d) - 4800*Lambda4*Lambda6*Sqr(g1d) - 3600*Lambda3*
      Lambda7*Sqr(g1d) - 4800*Lambda4*Lambda7*Sqr(g1d) - 2400*Conj(Lambda5)*
      Conj(Lambda6)*Sqr(g1d) + 225*Lambda7*Sqr(g1)*Sqr(g1d) - 2400*Lambda3*
      Lambda6*Sqr(g1dp) - 1600*Lambda4*Lambda6*Sqr(g1dp) - 1200*Lambda3*Lambda7
      *Sqr(g1dp) - 1600*Lambda4*Lambda7*Sqr(g1dp) - 800*Conj(Lambda5)*Conj(
      Lambda6)*Sqr(g1dp) + 75*Lambda7*Sqr(g1)*Sqr(g1dp) - 450*Lambda7*Sqr(g1d)*
      Sqr(g1dp) + 14400*Lambda3*Lambda6*Sqr(g2) + 7200*Lambda4*Lambda6*Sqr(g2)
      + 43200*Lambda2*Lambda7*Sqr(g2) + 7200*Lambda3*Lambda7*Sqr(g2) + 14400*
      Lambda4*Lambda7*Sqr(g2) + 2250*Lambda7*traceYdAdjYd*Sqr(g2) + 750*Lambda7
      *traceYeAdjYe*Sqr(g2) + 6750*Lambda7*traceYuAdjYu*Sqr(g2) + 600*Lambda6*
      Sqr(g1)*Sqr(g2) + 1740*Lambda7*Sqr(g1)*Sqr(g2) + 4125*Lambda7*Sqr(g1d)*
      Sqr(g2) + 375*Lambda7*Sqr(g1dp)*Sqr(g2) - 28800*Lambda2*Lambda7*Sqr(g2u)
      - 3600*Lambda3*Lambda7*Sqr(g2u) - 4800*Lambda4*Lambda7*Sqr(g2u) + 675*
      Lambda7*Sqr(g1)*Sqr(g2u) - 2600*Lambda7*Sqr(g1d)*Sqr(g2u) + 12375*Lambda7
      *Sqr(g2)*Sqr(g2u) - 9600*Lambda2*Lambda7*Sqr(g2up) - 1200*Lambda3*Lambda7
      *Sqr(g2up) - 1600*Lambda4*Lambda7*Sqr(g2up) + 225*Lambda7*Sqr(g1)*Sqr(
      g2up) + 200*Lambda7*Sqr(g1dp)*Sqr(g2up) + 1125*Lambda7*Sqr(g2)*Sqr(g2up)
      - 950*Lambda7*Sqr(g2u)*Sqr(g2up) + 8000*Lambda7*traceYdAdjYd*Sqr(g3) +
      24000*Lambda7*traceYuAdjYu*Sqr(g3) + 2400*Lambda7*Sqr(Lambda1) - 124800*
      Lambda7*Sqr(Lambda2) - 13800*Lambda6*Sqr(Lambda3) - 12200*Lambda7*Sqr(
      Lambda3) - 13000*Lambda6*Sqr(Lambda4) - 13000*Lambda7*Sqr(Lambda4) -
      16800*Conj(Lambda6)*Sqr(Lambda6) - 16800*Conj(Lambda6)*Sqr(Lambda7)));
   const double beta_Lambda7_2 = Re(-84*Lambda6*AbsSqr(Lambda7) - 4*
      traceYeAdjYe*Conj(Lambda5)*Conj(Lambda6) - 20*Lambda1*Conj(Lambda5)*Conj(
      Lambda7) - 142*Lambda2*Conj(Lambda5)*Conj(Lambda7) - 69*Lambda3*Conj(
      Lambda5)*Conj(Lambda7) - 73*Lambda4*Conj(Lambda5)*Conj(Lambda7) - 30*
      traceYdAdjYd*Conj(Lambda5)*Conj(Lambda7) - 10*traceYeAdjYe*Conj(Lambda5)*
      Conj(Lambda7) - 30*traceYuAdjYu*Conj(Lambda5)*Conj(Lambda7) + 12*Conj(
      Lambda5)*Conj(Lambda7)*Sqr(g1) - 15*Conj(Lambda5)*Conj(Lambda7)*Sqr(g1d)
      - 5*Conj(Lambda5)*Conj(Lambda7)*Sqr(g1dp) + 54*Conj(Lambda5)*Conj(Lambda7
      )*Sqr(g2) - 15*Conj(Lambda5)*Conj(Lambda7)*Sqr(g2u) - 5*Conj(Lambda5)*
      Conj(Lambda7)*Sqr(g2up) - 22*Conj(Lambda7)*Sqr(Lambda6) - 111*Conj(
      Lambda7)*Sqr(Lambda7));

   beta_Lambda7 = beta_Lambda7_1 + beta_Lambda7_2;


   return twoLoop * beta_Lambda7;
}

/**
 * Calculates the 3-loop beta function of Lambda7.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return threeLoop * beta_Lambda7;
}

/**
 * Calculates the 4-loop beta function of Lambda7.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return fourLoop * beta_Lambda7;
}

/**
 * Calculates the 5-loop beta function of Lambda7.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda7_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda7;

   beta_Lambda7 = 0;


   return fiveLoop * beta_Lambda7;
}

} // namespace flexiblesusy
