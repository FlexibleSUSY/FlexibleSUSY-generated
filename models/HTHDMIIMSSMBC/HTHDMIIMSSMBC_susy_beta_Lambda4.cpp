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

// File generated at Sun 26 Aug 2018 14:07:18

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda4.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda4;

   beta_Lambda4 = Re(oneOver16PiSqr*(4*Lambda1*Lambda4 + 4*Lambda2*Lambda4 + 8*
      Lambda3*Lambda4 + 6*Lambda4*traceYdAdjYd + 12*traceYdAdjYuYuAdjYd + 2*
      Lambda4*traceYeAdjYe + 6*Lambda4*traceYuAdjYu + 8*AbsSqr(Lambda5) + 10*
      AbsSqr(Lambda7) + 2*(5*Lambda6 + Lambda7)*Conj(Lambda6) + 2*Lambda6*Conj(
      Lambda7) - 1.8*Lambda4*Sqr(g1) - 9*Lambda4*Sqr(g2) + 1.8*Sqr(g1)*Sqr(g2)
      + 4*Sqr(Lambda4)));


   return beta_Lambda4;
}

/**
 * Calculates the 2-loop beta function of Lambda4.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_2_loop(const Susy_traces& susy_traces) const
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

   beta_Lambda4 = Re(twoLoop*(-80*Lambda1*Lambda3*Lambda4 - 80*Lambda2*Lambda3*
      Lambda4 - 48*Lambda5*Lambda6*Lambda7 - 24*Lambda1*Lambda4*traceYdAdjYd -
      24*Lambda3*Lambda4*traceYdAdjYd - 13.5*Lambda4*traceYdAdjYdYdAdjYd - 12*
      traceYdAdjYdYdAdjYuYuAdjYd - 24*Lambda3*traceYdAdjYuYuAdjYd - 33*Lambda4*
      traceYdAdjYuYuAdjYd - 12*traceYdAdjYuYuAdjYdYdAdjYd - 24*
      traceYdAdjYuYuAdjYuYuAdjYd - 8*Lambda1*Lambda4*traceYeAdjYe - 8*Lambda3*
      Lambda4*traceYeAdjYe - 4.5*Lambda4*traceYeAdjYeYeAdjYe - 24*Lambda2*
      Lambda4*traceYuAdjYu - 24*Lambda3*Lambda4*traceYuAdjYu - 13.5*Lambda4*
      traceYuAdjYuYuAdjYu - 20*Lambda1*AbsSqr(Lambda7) - 148*Lambda2*AbsSqr(
      Lambda7) - 72*Lambda3*AbsSqr(Lambda7) - 68*Lambda4*AbsSqr(Lambda7) - 60*
      traceYuAdjYu*AbsSqr(Lambda7) - 20*Lambda1*Lambda6*Conj(Lambda7) - 20*
      Lambda2*Lambda6*Conj(Lambda7) - 40*Lambda3*Lambda6*Conj(Lambda7) - 80*
      Lambda4*Lambda6*Conj(Lambda7) - 6*Lambda6*traceYdAdjYd*Conj(Lambda7) - 2*
      Lambda6*traceYeAdjYe*Conj(Lambda7) - 6*Lambda6*traceYuAdjYu*Conj(Lambda7)
      + 7.665*Lambda4*Quad(g1) - 23.875*Lambda4*Quad(g2) + 4.8*Lambda1*Lambda4*
      Sqr(g1) + 4.8*Lambda2*Lambda4*Sqr(g1) + 2.4*Lambda3*Lambda4*Sqr(g1) +
      1.25*Lambda4*traceYdAdjYd*Sqr(g1) + 0.8*traceYdAdjYuYuAdjYd*Sqr(g1) +
      3.75*Lambda4*traceYeAdjYe*Sqr(g1) + 4.25*Lambda4*traceYuAdjYu*Sqr(g1) +
      8.4*AbsSqr(Lambda7)*Sqr(g1) + 2.4*Lambda6*Conj(Lambda7)*Sqr(g1) - 10*Quad
      (g2)*Sqr(g1) + 36*Lambda3*Lambda4*Sqr(g2) + 11.25*Lambda4*traceYdAdjYd*
      Sqr(g2) + 3.75*Lambda4*traceYeAdjYe*Sqr(g2) + 11.25*Lambda4*traceYuAdjYu*
      Sqr(g2) + 54*AbsSqr(Lambda7)*Sqr(g2) - 14.1*Quad(g1)*Sqr(g2) + 6*Lambda1*
      Sqr(g1)*Sqr(g2) + 6*Lambda2*Sqr(g1)*Sqr(g2) + 1.2*Lambda3*Sqr(g1)*Sqr(g2)
      + 7.65*Lambda4*Sqr(g1)*Sqr(g2) + 5.4*traceYdAdjYd*Sqr(g1)*Sqr(g2) + 6.6*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) + 12.6*traceYuAdjYu*Sqr(g1)*Sqr(g2) + Conj(
      Lambda6)*(-2*(36*Lambda3*Lambda6 + 34*Lambda4*Lambda6 + 20*Lambda3*
      Lambda7 + 40*Lambda4*Lambda7 + 10*Lambda2*(Lambda6 + Lambda7) + 2*Lambda1
      *(37*Lambda6 + 5*Lambda7) + 30*Lambda6*traceYdAdjYd + 3*Lambda7*
      traceYdAdjYd + 10*Lambda6*traceYeAdjYe + Lambda7*traceYeAdjYe + 3*Lambda7
      *traceYuAdjYu) + 1.2*(7*Lambda6 + 2*Lambda7)*Sqr(g1) + 54*Lambda6*Sqr(g2)
      ) + 40*Lambda4*traceYdAdjYd*Sqr(g3) + 64*traceYdAdjYuYuAdjYd*Sqr(g3) + 40
      *Lambda4*traceYuAdjYu*Sqr(g3) - 28*Lambda4*Sqr(Lambda1) - 28*Lambda4*Sqr(
      Lambda2) - 28*Lambda4*Sqr(Lambda3) - 40*Lambda1*Sqr(Lambda4) - 40*Lambda2
      *Sqr(Lambda4) - 28*Lambda3*Sqr(Lambda4) - 12*traceYdAdjYd*Sqr(Lambda4) -
      4*traceYeAdjYe*Sqr(Lambda4) - 12*traceYuAdjYu*Sqr(Lambda4) + 4.8*Sqr(g1)*
      Sqr(Lambda4) + 18*Sqr(g2)*Sqr(Lambda4) - 40*Lambda5*Sqr(Lambda6) - 40*
      Lambda5*Sqr(Lambda7) - 0.4*Conj(Lambda5)*(120*Conj(Lambda6)*Conj(Lambda7)
      + Lambda5*(-24*Sqr(g1) + 5*(24*Lambda1 + 24*Lambda2 + 24*Lambda3 + 13*
      Lambda4 + 12*traceYdAdjYd + 4*traceYeAdjYe + 12*traceYuAdjYu - 27*Sqr(g2)
      )) + 100*Sqr(Conj(Lambda6)) + 100*Sqr(Conj(Lambda7)))));


   return beta_Lambda4;
}

/**
 * Calculates the 3-loop beta function of Lambda4.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_3_loop(const Susy_traces& susy_traces) const
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
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return beta_Lambda4;
}

} // namespace flexiblesusy
