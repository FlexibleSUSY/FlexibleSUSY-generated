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

// File generated at Sun 26 Aug 2018 14:06:25

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda4.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda4;

   beta_Lambda4 = Re(oneOver16PiSqr*(-4*g1d*g1dp*g2u*g2up + 4*Lambda1*Lambda4 +
      4*Lambda2*Lambda4 + 8*Lambda3*Lambda4 + 6*Lambda4*traceYdAdjYd + 12*
      traceYdAdjYuYuAdjYd + 2*Lambda4*traceYeAdjYe + 6*Lambda4*traceYuAdjYu + 8
      *AbsSqr(Lambda5) + 10*AbsSqr(Lambda7) + 2*(5*Lambda6 + Lambda7)*Conj(
      Lambda6) + 2*Lambda6*Conj(Lambda7) - 1.8*Lambda4*Sqr(g1) + 3*Lambda4*Sqr(
      g1d) + Lambda4*Sqr(g1dp) - 9*Lambda4*Sqr(g2) + 1.8*Sqr(g1)*Sqr(g2) + 3*
      Lambda4*Sqr(g2u) + 4*Sqr(g1d)*Sqr(g2u) + Lambda4*Sqr(g2up) + 4*Sqr(
      Lambda4)));


   return beta_Lambda4;
}

/**
 * Calculates the 2-loop beta function of Lambda4.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_2_loop(const Susy_traces& susy_traces) const
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


   double beta_Lambda4;

   const double beta_Lambda4_1 = Re(-0.005*twoLoop*(3*Quad(g1)*(-511*Lambda4 +
      940*Sqr(g2)) - 5*Sqr(g1)*(192*Lambda1*Lambda4 + 192*Lambda2*Lambda4 + 96*
      Lambda3*Lambda4 + 50*Lambda4*traceYdAdjYd + 32*traceYdAdjYuYuAdjYd + 150*
      Lambda4*traceYeAdjYe + 170*Lambda4*traceYuAdjYu - 528*Quad(g2) + 240*
      Lambda1*Sqr(g2) + 240*Lambda2*Sqr(g2) + 48*Lambda3*Sqr(g2) + 306*Lambda4*
      Sqr(g2) + 216*traceYdAdjYd*Sqr(g2) + 264*traceYeAdjYe*Sqr(g2) + 504*
      traceYuAdjYu*Sqr(g2) - 3*Sqr(g1dp)*(-5*Lambda4 + 4*Sqr(g2)) + 9*Sqr(g1d)*
      (5*Lambda4 + 28*Sqr(g2)) + 45*Lambda4*Sqr(g2u) + 252*Sqr(g2)*Sqr(g2u) +
      15*Lambda4*Sqr(g2up) - 12*Sqr(g2)*Sqr(g2up) + 192*Sqr(Lambda4)) + 25*(640
      *Lambda1*Lambda3*Lambda4 + 640*Lambda2*Lambda3*Lambda4 + 384*Lambda5*
      Lambda6*Lambda7 + 192*Lambda1*Lambda4*traceYdAdjYd + 192*Lambda3*Lambda4*
      traceYdAdjYd + 108*Lambda4*traceYdAdjYdYdAdjYd + 96*
      traceYdAdjYdYdAdjYuYuAdjYd + 192*Lambda3*traceYdAdjYuYuAdjYd + 264*
      Lambda4*traceYdAdjYuYuAdjYd + 96*traceYdAdjYuYuAdjYdYdAdjYd + 192*
      traceYdAdjYuYuAdjYuYuAdjYd + 64*Lambda1*Lambda4*traceYeAdjYe + 64*Lambda3
      *Lambda4*traceYeAdjYe + 36*Lambda4*traceYeAdjYeYeAdjYe + 192*Lambda2*
      Lambda4*traceYuAdjYu + 192*Lambda3*Lambda4*traceYuAdjYu + 108*Lambda4*
      traceYuAdjYuYuAdjYu - 24*g1dp*g2u*g2up*Cube(g1d) + 9*Lambda4*Quad(g1dp) +
      111*Lambda4*Quad(g2) + 45*Lambda4*Quad(g2u) + 9*Lambda4*Quad(g2up) - 288*
      Lambda3*Lambda4*Sqr(g2) - 90*Lambda4*traceYdAdjYd*Sqr(g2) - 30*Lambda4*
      traceYeAdjYe*Sqr(g2) - 90*Lambda4*traceYuAdjYu*Sqr(g2) + 96*Lambda2*
      Lambda4*Sqr(g2u) + 96*Lambda3*Lambda4*Sqr(g2u) - 165*Lambda4*Sqr(g2)*Sqr(
      g2u) + 9*Quad(g1d)*(5*Lambda4 + 8*Sqr(g2u)) + Lambda4*Sqr(g1dp)*(-15*Sqr(
      g2) + 4*(8*Lambda1 + 8*Lambda3 + 4*Lambda4 - Sqr(g2up))) + 32*Lambda2*
      Lambda4*Sqr(g2up) + 32*Lambda3*Lambda4*Sqr(g2up) - 15*Lambda4*Sqr(g2)*Sqr
      (g2up) + 18*Lambda4*Sqr(g2u)*Sqr(g2up) - 8*g1d*g1dp*g2u*g2up*(8*Lambda3 +
      12*Lambda4 + 5*Sqr(g1dp) - 8*Sqr(g2) + 3*Sqr(g2u) + 5*Sqr(g2up)) - 320*
      Lambda4*traceYdAdjYd*Sqr(g3) - 512*traceYdAdjYuYuAdjYd*Sqr(g3) - 320*
      Lambda4*traceYuAdjYu*Sqr(g3) + 224*Lambda4*Sqr(Lambda1) + 224*Lambda4*Sqr
      (Lambda2) + 224*Lambda4*Sqr(Lambda3) + 320*Lambda1*Sqr(Lambda4) + 320*
      Lambda2*Sqr(Lambda4) + 224*Lambda3*Sqr(Lambda4) + 96*traceYdAdjYd*Sqr(
      Lambda4) + 32*traceYeAdjYe*Sqr(Lambda4) + 96*traceYuAdjYu*Sqr(Lambda4) -
      144*Sqr(g2)*Sqr(Lambda4) + 48*Sqr(g2u)*Sqr(Lambda4) + 16*Sqr(g2up)*Sqr(
      Lambda4) + Sqr(g1d)*(96*Lambda1*Lambda4 + 96*Lambda3*Lambda4 + 72*Quad(
      g2u) + Sqr(g1dp)*(18*Lambda4 - 8*Sqr(g2u)) + 64*Lambda3*Sqr(g2u) + 84*
      Lambda4*Sqr(g2u) - Sqr(g2)*(165*Lambda4 + 128*Sqr(g2u)) - 8*Sqr(g2u)*Sqr(
      g2up) + 48*Sqr(Lambda4)) + 320*Lambda5*Sqr(Lambda6) + 320*Lambda5*Sqr(
      Lambda7))));
   const double beta_Lambda4_2 = Re(0.2*twoLoop*(Conj(Lambda6)*(6*(7*Lambda6 +
      2*Lambda7)*Sqr(g1) - 5*(148*Lambda1*Lambda6 + 20*Lambda2*Lambda6 + 72*
      Lambda3*Lambda6 + 68*Lambda4*Lambda6 + 20*Lambda1*Lambda7 + 20*Lambda2*
      Lambda7 + 40*Lambda3*Lambda7 + 80*Lambda4*Lambda7 + 60*Lambda6*
      traceYdAdjYd + 6*Lambda7*traceYdAdjYd + 20*Lambda6*traceYeAdjYe + 2*
      Lambda7*traceYeAdjYe + 6*Lambda7*traceYuAdjYu + 3*(10*Lambda6 + Lambda7)*
      Sqr(g1d) + (10*Lambda6 + Lambda7)*Sqr(g1dp) - 54*Lambda6*Sqr(g2) + 3*
      Lambda7*Sqr(g2u) + Lambda7*Sqr(g2up))) + Conj(Lambda7)*(6*(2*Lambda6 + 7*
      Lambda7)*Sqr(g1) - 5*(20*Lambda1*Lambda6 + 20*Lambda2*Lambda6 + 40*
      Lambda3*Lambda6 + 80*Lambda4*Lambda6 + 20*Lambda1*Lambda7 + 148*Lambda2*
      Lambda7 + 72*Lambda3*Lambda7 + 68*Lambda4*Lambda7 + 6*Lambda6*
      traceYdAdjYd + 2*Lambda6*traceYeAdjYe + 6*Lambda6*traceYuAdjYu + 60*
      Lambda7*traceYuAdjYu + 3*Lambda6*Sqr(g1d) + Lambda6*Sqr(g1dp) - 54*
      Lambda7*Sqr(g2) + 3*Lambda6*Sqr(g2u) + 30*Lambda7*Sqr(g2u) + Lambda6*Sqr(
      g2up) + 10*Lambda7*Sqr(g2up))) + 2*Conj(Lambda5)*(-120*Conj(Lambda6)*Conj
      (Lambda7) + Lambda5*(24*Sqr(g1) - 5*(24*Lambda1 + 24*Lambda2 + 24*Lambda3
       + 13*Lambda4 + 12*traceYdAdjYd + 4*traceYeAdjYe + 12*traceYuAdjYu + 6*
      Sqr(g1d) + 2*Sqr(g1dp) - 27*Sqr(g2) + 6*Sqr(g2u) + 2*Sqr(g2up))) - 100*
      Sqr(Conj(Lambda6)) - 100*Sqr(Conj(Lambda7)))));

   beta_Lambda4 = beta_Lambda4_1 + beta_Lambda4_2;


   return beta_Lambda4;
}

/**
 * Calculates the 3-loop beta function of Lambda4.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return beta_Lambda4;
}

/**
 * Calculates the 4-loop beta function of Lambda4.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return beta_Lambda4;
}

} // namespace flexiblesusy
