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

// File generated at Sun 26 Aug 2018 14:07:22

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda2.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(oneOver16PiSqr*(2*Lambda3*Lambda4 + 12*Lambda2*
      traceYuAdjYu - 6*traceYuAdjYuYuAdjYu + AbsSqr(Lambda5) + 12*AbsSqr(
      Lambda7) + 0.135*Quad(g1) + 1.125*Quad(g2) - 1.8*Lambda2*Sqr(g1) - 9*
      Lambda2*Sqr(g2) + 0.45*Sqr(g1)*Sqr(g2) + 24*Sqr(Lambda2) + 2*Sqr(Lambda3)
      + Sqr(Lambda4)));


   return beta_Lambda2;
}

/**
 * Calculates the 2-loop beta function of Lambda2.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_2_loop(const Susy_traces& susy_traces) const
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

   beta_Lambda2 = Re(twoLoop*(-20*Lambda2*Lambda3*Lambda4 - 10*Lambda5*Lambda6*
      Lambda7 - 12*Lambda3*Lambda4*traceYdAdjYd - 9*Lambda2*traceYdAdjYuYuAdjYd
       + 6*traceYdAdjYuYuAdjYuYuAdjYd - 4*Lambda3*Lambda4*traceYeAdjYe - 3*
      Lambda2*traceYuAdjYuYuAdjYu + 30*traceYuAdjYuYuAdjYuYuAdjYu - 318*Lambda2
      *AbsSqr(Lambda7) - 66*Lambda3*AbsSqr(Lambda7) - 70*Lambda4*AbsSqr(Lambda7
      ) - 36*traceYdAdjYd*AbsSqr(Lambda7) - 12*traceYeAdjYe*AbsSqr(Lambda7) -
      36*traceYuAdjYu*AbsSqr(Lambda7) + 2*(3*Lambda2*Lambda6 - (9*Lambda3 + 7*
      Lambda4)*(Lambda6 + Lambda7))*Conj(Lambda6) - 18*Lambda3*Lambda6*Conj(
      Lambda7) - 14*Lambda4*Lambda6*Conj(Lambda7) - 312*Cube(Lambda2) - 8*Cube(
      Lambda3) - 6*Cube(Lambda4) - 1.9125*Power6(g1) + 16.1875*Power6(g2) +
      10.365*Lambda2*Quad(g1) + 0.9*Lambda3*Quad(g1) + 0.45*Lambda4*Quad(g1) -
      1.71*traceYuAdjYu*Quad(g1) - 1.375*Lambda2*Quad(g2) + 7.5*Lambda3*Quad(g2
      ) + 3.75*Lambda4*Quad(g2) - 2.25*traceYuAdjYu*Quad(g2) + 2.4*Lambda3*
      Lambda4*Sqr(g1) + 8.5*Lambda2*traceYuAdjYu*Sqr(g1) - 1.6*
      traceYuAdjYuYuAdjYu*Sqr(g1) + 10.8*AbsSqr(Lambda7)*Sqr(g1) - 4.1875*Quad(
      g2)*Sqr(g1) + 12*Lambda3*Lambda4*Sqr(g2) + 22.5*Lambda2*traceYuAdjYu*Sqr(
      g2) + 54*AbsSqr(Lambda7)*Sqr(g2) - 4.5375*Quad(g1)*Sqr(g2) + 5.85*Lambda2
      *Sqr(g1)*Sqr(g2) + 1.5*Lambda4*Sqr(g1)*Sqr(g2) + 6.3*traceYuAdjYu*Sqr(g1)
      *Sqr(g2) + 80*Lambda2*traceYuAdjYu*Sqr(g3) - 32*traceYuAdjYuYuAdjYu*Sqr(
      g3) - 144*traceYuAdjYu*Sqr(Lambda2) + 21.6*Sqr(g1)*Sqr(Lambda2) + 108*Sqr
      (g2)*Sqr(Lambda2) - 20*Lambda2*Sqr(Lambda3) - 12*Lambda4*Sqr(Lambda3) -
      12*traceYdAdjYd*Sqr(Lambda3) - 4*traceYeAdjYe*Sqr(Lambda3) + 2.4*Sqr(g1)*
      Sqr(Lambda3) + 12*Sqr(g2)*Sqr(Lambda3) - 12*Lambda2*Sqr(Lambda4) - 16*
      Lambda3*Sqr(Lambda4) - 6*traceYdAdjYd*Sqr(Lambda4) - 2*traceYeAdjYe*Sqr(
      Lambda4) + 1.2*Sqr(g1)*Sqr(Lambda4) + 3*Sqr(g2)*Sqr(Lambda4) - 5*Lambda5*
      Sqr(Lambda6) - 37*Lambda5*Sqr(Lambda7) - 0.2*Conj(Lambda5)*(50*Conj(
      Lambda6)*Conj(Lambda7) + Lambda5*(10*(7*Lambda2 + 10*Lambda3 + 11*Lambda4
       + 3*traceYdAdjYd + traceYeAdjYe) + 3*Sqr(g1)) + 25*Sqr(Conj(Lambda6)) +
      185*Sqr(Conj(Lambda7)))));


   return beta_Lambda2;
}

/**
 * Calculates the 3-loop beta function of Lambda2.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_3_loop(const Susy_traces& susy_traces) const
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
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

} // namespace flexiblesusy
