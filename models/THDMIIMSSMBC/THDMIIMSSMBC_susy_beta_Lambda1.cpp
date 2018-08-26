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

// File generated at Sun 26 Aug 2018 14:08:40

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

   beta_Lambda1 = Re(oneOver16PiSqr*(2*Lambda3*Lambda4 + 12*Lambda1*
      traceYdAdjYd - 6*traceYdAdjYdYdAdjYd + 4*Lambda1*traceYeAdjYe - 2*
      traceYeAdjYeYeAdjYe + AbsSqr(Lambda5) + 12*AbsSqr(Lambda6) + 0.135*Quad(
      g1) + 1.125*Quad(g2) - 1.8*Lambda1*Sqr(g1) - 9*Lambda1*Sqr(g2) + 0.45*Sqr
      (g1)*Sqr(g2) + 24*Sqr(Lambda1) + 2*Sqr(Lambda3) + Sqr(Lambda4)));


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

   beta_Lambda1 = Re(twoLoop*(-20*Lambda1*Lambda3*Lambda4 - 10*Lambda5*Lambda6*
      Lambda7 - 3*Lambda1*traceYdAdjYdYdAdjYd + 30*traceYdAdjYdYdAdjYdYdAdjYd -
      9*Lambda1*traceYdAdjYuYuAdjYd + 6*traceYdAdjYuYuAdjYdYdAdjYd - Lambda1*
      traceYeAdjYeYeAdjYe + 10*traceYeAdjYeYeAdjYeYeAdjYe - 12*Lambda3*Lambda4*
      traceYuAdjYu + 6*Lambda1*AbsSqr(Lambda7) - 18*Lambda3*AbsSqr(Lambda7) -
      14*Lambda4*AbsSqr(Lambda7) - 18*Lambda3*Lambda6*Conj(Lambda7) - 14*
      Lambda4*Lambda6*Conj(Lambda7) - 312*Cube(Lambda1) - 8*Cube(Lambda3) - 6*
      Cube(Lambda4) - 1.7685*Power6(g1) + 18.1875*Power6(g2) + 9.765*Lambda1*
      Quad(g1) + 0.9*Lambda3*Quad(g1) + 0.45*Lambda4*Quad(g1) + 0.45*
      traceYdAdjYd*Quad(g1) - 2.25*traceYeAdjYe*Quad(g1) - 6.375*Lambda1*Quad(
      g2) + 7.5*Lambda3*Quad(g2) + 3.75*Lambda4*Quad(g2) - 2.25*traceYdAdjYd*
      Quad(g2) - 0.75*traceYeAdjYe*Quad(g2) + 2.4*Lambda3*Lambda4*Sqr(g1) + 2.5
      *Lambda1*traceYdAdjYd*Sqr(g1) + 0.8*traceYdAdjYdYdAdjYd*Sqr(g1) + 7.5*
      Lambda1*traceYeAdjYe*Sqr(g1) - 2.4*traceYeAdjYeYeAdjYe*Sqr(g1) - 3.7875*
      Quad(g2)*Sqr(g1) + 12*Lambda3*Lambda4*Sqr(g2) + 22.5*Lambda1*traceYdAdjYd
      *Sqr(g2) + 7.5*Lambda1*traceYeAdjYe*Sqr(g2) - 4.2975*Quad(g1)*Sqr(g2) +
      5.85*Lambda1*Sqr(g1)*Sqr(g2) + 1.5*Lambda4*Sqr(g1)*Sqr(g2) + 2.7*
      traceYdAdjYd*Sqr(g1)*Sqr(g2) + 3.3*traceYeAdjYe*Sqr(g1)*Sqr(g2) + Conj(
      Lambda6)*(10.8*Lambda6*Sqr(g1) + 2*(-159*Lambda1*Lambda6 - 33*Lambda3*
      Lambda6 - 35*Lambda4*Lambda6 - 9*Lambda3*Lambda7 - 7*Lambda4*Lambda7 - 18
      *Lambda6*traceYdAdjYd - 6*Lambda6*traceYeAdjYe - 18*Lambda6*traceYuAdjYu
      + 27*Lambda6*Sqr(g2))) + 80*Lambda1*traceYdAdjYd*Sqr(g3) - 32*
      traceYdAdjYdYdAdjYd*Sqr(g3) - 144*traceYdAdjYd*Sqr(Lambda1) - 48*
      traceYeAdjYe*Sqr(Lambda1) + 21.6*Sqr(g1)*Sqr(Lambda1) + 108*Sqr(g2)*Sqr(
      Lambda1) - 20*Lambda1*Sqr(Lambda3) - 12*Lambda4*Sqr(Lambda3) - 12*
      traceYuAdjYu*Sqr(Lambda3) + 2.4*Sqr(g1)*Sqr(Lambda3) + 12*Sqr(g2)*Sqr(
      Lambda3) - 12*Lambda1*Sqr(Lambda4) - 16*Lambda3*Sqr(Lambda4) - 6*
      traceYuAdjYu*Sqr(Lambda4) + 1.2*Sqr(g1)*Sqr(Lambda4) + 3*Sqr(g2)*Sqr(
      Lambda4) - 37*Lambda5*Sqr(Lambda6) - 5*Lambda5*Sqr(Lambda7) - 0.2*Conj(
      Lambda5)*(50*Conj(Lambda6)*Conj(Lambda7) + Lambda5*(70*Lambda1 + 100*
      Lambda3 + 110*Lambda4 + 30*traceYuAdjYu + 3*Sqr(g1)) + 185*Sqr(Conj(
      Lambda6)) + 25*Sqr(Conj(Lambda7)))));


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

} // namespace flexiblesusy
