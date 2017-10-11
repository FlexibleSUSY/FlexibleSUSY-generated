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

// File generated at Tue 10 Oct 2017 20:53:46

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

   beta_Lambda3 = Re(oneOver16PiSqr*(2*g1d*g1dp*g2u*g2up + 12*Lambda1*
      Lambda3 + 12*Lambda2*Lambda3 + 4*Lambda1*Lambda4 + 4*Lambda2*Lambda4 + 16
      *Lambda6*Lambda7 + 6*Lambda3*traceYdAdjYd - 12*traceYdAdjYuYuAdjYd + 2*
      Lambda3*traceYeAdjYe + 6*Lambda3*traceYuAdjYu + 0.27*Quad(g1) + 2.25*Quad
      (g2) + 3*Lambda3*Sqr(g1d) + Lambda3*Sqr(g1dp) - 9*Lambda3*Sqr(g2) - 0.9*
      Sqr(g1)*(2*Lambda3 + Sqr(g2)) + 3*Lambda3*Sqr(g2u) - 5*Sqr(g1d)*Sqr(g2u)
      + Lambda3*Sqr(g2up) - Sqr(g1dp)*Sqr(g2up) + 4*Sqr(Lambda3) + 2*Sqr(
      Lambda4) + 2*Sqr(Lambda5) + 4*Sqr(Lambda6) + 4*Sqr(Lambda7)));


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
   const double traceYdAdjYdYdAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd;


   double beta_Lambda3;

   const double beta_Lambda3_1 = Re(-0.005*twoLoop*(765*Power6(g1) + 3*
      Quad(g1)*(-180*Lambda1 - 180*Lambda2 - 631*Lambda3 - 60*Lambda4 + 9*Sqr(
      g1d) + 3*Sqr(g1dp) - 335*Sqr(g2) + 9*Sqr(g2u) + 3*Sqr(g2up)) - 15*Sqr(g1)
      *(192*Lambda1*Lambda3 + 192*Lambda2*Lambda3 + 64*Lambda1*Lambda4 + 64*
      Lambda2*Lambda4 + 43*Quad(g2) + Sqr(g1d)*(15*Lambda3 - 42*Sqr(g2)) - 40*
      Lambda1*Sqr(g2) - 40*Lambda2*Sqr(g2) + 22*Lambda3*Sqr(g2) - 24*Lambda4*
      Sqr(g2) + Sqr(g1dp)*(5*Lambda3 + 2*Sqr(g2)) + 15*Lambda3*Sqr(g2u) - 42*
      Sqr(g2)*Sqr(g2u) + 5*Lambda3*Sqr(g2up) + 2*Sqr(g2)*Sqr(g2up) + 16*Sqr(
      Lambda3) - 16*Sqr(Lambda4)) - 25*(-256*Lambda1*Lambda3*Lambda4 - 256*
      Lambda2*Lambda3*Lambda4 - 12*g1dp*g2u*g2up*Cube(g1d) - 96*Cube(Lambda3) +
      195*Power6(g2) + 180*Lambda1*Quad(g2) + 180*Lambda2*Quad(g2) + 9*Lambda3
      *Quad(g2) + 60*Lambda4*Quad(g2) - 45*Lambda3*Quad(g2u) - 9*Lambda3*Quad(
      g2up) + 576*Lambda1*Lambda3*Sqr(g2) + 576*Lambda2*Lambda3*Sqr(g2) + 288*
      Lambda1*Lambda4*Sqr(g2) + 288*Lambda2*Lambda4*Sqr(g2) - 96*Lambda3*
      Lambda4*Sqr(g2) - 288*Lambda2*Lambda3*Sqr(g2u) - 96*Lambda2*Lambda4*Sqr(
      g2u) - 153*Quad(g2)*Sqr(g2u) + 165*Lambda3*Sqr(g2)*Sqr(g2u) + 3*Quad(g1d)
      *(-15*Lambda3 + 38*Sqr(g2u)) - 96*Lambda2*Lambda3*Sqr(g2up) - 32*Lambda2*
      Lambda4*Sqr(g2up) - 3*Quad(g2)*Sqr(g2up) + 15*Lambda3*Sqr(g2)*Sqr(g2up) -
      18*Lambda3*Sqr(g2u)*Sqr(g2up) - 4*g1d*g1dp*g2u*g2up*(8*Lambda3 + 5*Sqr(
      g1dp) - 8*Sqr(g2) + 3*Sqr(g2u) + 5*Sqr(g2up)) + Quad(g1dp)*(-9*Lambda3 +
      14*Sqr(g2up)) - 480*Lambda3*Sqr(Lambda1) - 128*Lambda4*Sqr(Lambda1) - 480
      *Lambda3*Sqr(Lambda2) - 128*Lambda4*Sqr(Lambda2) - 576*Lambda1*Sqr(
      Lambda3) - 576*Lambda2*Sqr(Lambda3) - 32*Lambda4*Sqr(Lambda3) + 48*Sqr(g2
      )*Sqr(Lambda3) - 48*Sqr(g2u)*Sqr(Lambda3) - 16*Sqr(g2up)*Sqr(Lambda3) +
      Sqr(g1dp)*(-3*Quad(g2) + 15*Lambda3*Sqr(g2) + 2*(-48*Lambda1*Lambda3 - 16
      *Lambda1*Lambda4 + 7*Quad(g2up) + 2*Lambda3*Sqr(g2up) + 9*Sqr(g2u)*Sqr(
      g2up) - 8*Sqr(Lambda3) - 4*Sqr(Lambda4))) - 224*Lambda1*Sqr(Lambda4) -
      224*Lambda2*Sqr(Lambda4) + 48*Sqr(g2)*Sqr(Lambda4) - 24*Sqr(g2u)*Sqr(
      Lambda4) - 8*Sqr(g2up)*Sqr(Lambda4) + Sqr(g1d)*(-153*Quad(g2) - 5*Sqr(g2)
      *(-33*Lambda3 + 32*Sqr(g2u)) + 2*(57*Quad(g2u) + Sqr(g2u)*(22*Lambda3 + 7
      *Sqr(g2up)) + Sqr(g1dp)*(-9*Lambda3 + 7*Sqr(g2u) + 9*Sqr(g2up)) - 12*(4*
      Lambda1*(3*Lambda3 + Lambda4) + 2*Sqr(Lambda3) + Sqr(Lambda4)))))));
   const double beta_Lambda3_2 = Re(-0.01*twoLoop*(8800*Lambda1*Lambda6*
      Lambda7 + 8800*Lambda2*Lambda6*Lambda7 + 7200*Lambda5*Lambda6*Lambda7 +
      4800*Lambda6*Lambda7*traceYdAdjYd - 1200*traceYdAdjYdYdAdjYuYuAdjYd -
      2400*traceYdAdjYuYuAdjYdYdAdjYd - 3600*traceYdAdjYuYuAdjYuYuAdjYd + 1600*
      Lambda6*Lambda7*traceYeAdjYe + 4800*Lambda6*Lambda7*traceYuAdjYu + 1200*
      Cube(Lambda4) - 45*traceYdAdjYd*Quad(g1) + 225*traceYeAdjYe*Quad(g1) +
      171*traceYuAdjYu*Quad(g1) + 225*traceYdAdjYd*Quad(g2) + 75*traceYeAdjYe*
      Quad(g2) + 225*traceYuAdjYu*Quad(g2) - 1920*Lambda6*Lambda7*Sqr(g1) + 80*
      traceYdAdjYuYuAdjYd*Sqr(g1) + 2400*Lambda6*Lambda7*Sqr(g1d) + 800*Lambda6
      *Lambda7*Sqr(g1dp) - 10800*Lambda6*Lambda7*Sqr(g2) + 270*traceYdAdjYd*Sqr
      (g1)*Sqr(g2) + 330*traceYeAdjYe*Sqr(g1)*Sqr(g2) + 630*traceYuAdjYu*Sqr(g1
      )*Sqr(g2) + 2400*Lambda6*Lambda7*Sqr(g2u) + 800*Lambda6*Lambda7*Sqr(g2up)
      + 6400*traceYdAdjYuYuAdjYd*Sqr(g3) + 400*(3*traceYdAdjYd + traceYeAdjYe
      + 3*traceYuAdjYu)*Sqr(Lambda3) + 200*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu)*Sqr(Lambda4) + 3600*Lambda1*Sqr(Lambda5) + 3600*Lambda2*Sqr
      (Lambda5) + 600*traceYdAdjYd*Sqr(Lambda5) + 200*traceYeAdjYe*Sqr(Lambda5)
      + 600*traceYuAdjYu*Sqr(Lambda5) - 240*Sqr(g1)*Sqr(Lambda5) + 300*Sqr(g1d
      )*Sqr(Lambda5) + 100*Sqr(g1dp)*Sqr(Lambda5) + 300*Sqr(g2u)*Sqr(Lambda5) +
      100*Sqr(g2up)*Sqr(Lambda5) + 12400*Lambda1*Sqr(Lambda6) + 4400*Lambda2*
      Sqr(Lambda6) + 6800*Lambda5*Sqr(Lambda6) + 2400*traceYdAdjYd*Sqr(Lambda6)
      + 800*traceYeAdjYe*Sqr(Lambda6) - 120*Sqr(g1)*Sqr(Lambda6) + 1200*Sqr(
      g1d)*Sqr(Lambda6) + 400*Sqr(g1dp)*Sqr(Lambda6) + 4400*Lambda1*Sqr(Lambda7
      ) + 12400*Lambda2*Sqr(Lambda7) + 6800*Lambda5*Sqr(Lambda7) + 2400*
      traceYuAdjYu*Sqr(Lambda7) - 120*Sqr(g1)*Sqr(Lambda7) + 1200*Sqr(g2u)*Sqr(
      Lambda7) + 400*Sqr(g2up)*Sqr(Lambda7) + 400*Lambda4*(22*Lambda6*Lambda7 +
      6*Lambda1*traceYdAdjYd + 2*Lambda1*traceYeAdjYe + 6*Lambda2*traceYuAdjYu
      + 11*Sqr(Lambda5) + 17*Sqr(Lambda6) + 17*Sqr(Lambda7)) + 25*Lambda3*(704
      *Lambda6*Lambda7 + 288*Lambda1*traceYdAdjYd + 54*traceYdAdjYdYdAdjYd - 60
      *traceYdAdjYuYuAdjYd + 96*Lambda1*traceYeAdjYe + 18*traceYeAdjYeYeAdjYe +
      288*Lambda2*traceYuAdjYu + 54*traceYuAdjYuYuAdjYu - 5*traceYdAdjYd*Sqr(
      g1) - 15*traceYeAdjYe*Sqr(g1) - 17*traceYuAdjYu*Sqr(g1) - 45*traceYdAdjYd
      *Sqr(g2) - 15*traceYeAdjYe*Sqr(g2) - 45*traceYuAdjYu*Sqr(g2) - 160*
      traceYdAdjYd*Sqr(g3) - 160*traceYuAdjYu*Sqr(g3) + 64*Sqr(Lambda4) + 72*
      Sqr(Lambda5) + 240*Sqr(Lambda6) + 240*Sqr(Lambda7))));

   beta_Lambda3 = beta_Lambda3_1 + beta_Lambda3_2;


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

} // namespace flexiblesusy