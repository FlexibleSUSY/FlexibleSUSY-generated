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

// File generated at Sun 26 Aug 2018 14:08:37

#include "THDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda5.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda5;

   beta_Lambda5 = Re(oneOver16PiSqr*(4*Lambda1*Lambda5 + 4*Lambda2*Lambda5 + 8*
      Lambda3*Lambda5 + 12*Lambda4*Lambda5 + 6*Lambda5*traceYdAdjYd + 2*Lambda5
      *traceYeAdjYe + 6*Lambda5*traceYuAdjYu + 4*Conj(Lambda6)*Conj(Lambda7) -
      1.8*Lambda5*Sqr(g1) - 9*Lambda5*Sqr(g2) + 10*Sqr(Conj(Lambda6)) + 10*Sqr(
      Conj(Lambda7))));


   return beta_Lambda5;
}

/**
 * Calculates the 2-loop beta function of Lambda5.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda5;

   beta_Lambda5 = Re(twoLoop*(-80*Lambda1*Lambda3*Lambda5 - 80*Lambda2*Lambda3*
      Lambda5 - 88*Lambda1*Lambda4*Lambda5 - 88*Lambda2*Lambda4*Lambda5 - 76*
      Lambda3*Lambda4*Lambda5 - 24*Lambda1*Lambda5*traceYdAdjYd - 24*Lambda3*
      Lambda5*traceYdAdjYd - 36*Lambda4*Lambda5*traceYdAdjYd - 7.5*Lambda5*
      traceYdAdjYdYdAdjYd - 33*Lambda5*traceYdAdjYuYuAdjYd - 8*Lambda1*Lambda5*
      traceYeAdjYe - 8*Lambda3*Lambda5*traceYeAdjYe - 12*Lambda4*Lambda5*
      traceYeAdjYe - 2.5*Lambda5*traceYeAdjYeYeAdjYe - 24*Lambda2*Lambda5*
      traceYuAdjYu - 24*Lambda3*Lambda5*traceYuAdjYu - 36*Lambda4*Lambda5*
      traceYuAdjYu - 7.5*Lambda5*traceYuAdjYuYuAdjYu - 72*Lambda5*AbsSqr(
      Lambda7) - 84*Lambda5*Lambda6*Conj(Lambda7) + 7.065*Lambda5*Quad(g1) +
      0.45*traceYdAdjYd*Quad(g1) - 2.25*traceYeAdjYe*Quad(g1) - 28.875*Lambda5*
      Quad(g2) - 2.25*traceYdAdjYd*Quad(g2) - 0.75*traceYeAdjYe*Quad(g2) - 2.4*
      Lambda1*Lambda5*Sqr(g1) - 2.4*Lambda2*Lambda5*Sqr(g1) + 9.6*Lambda3*
      Lambda5*Sqr(g1) + 14.4*Lambda4*Lambda5*Sqr(g1) + 1.25*Lambda5*
      traceYdAdjYd*Sqr(g1) + 3.75*Lambda5*traceYeAdjYe*Sqr(g1) + 4.25*Lambda5*
      traceYuAdjYu*Sqr(g1) - 0.8*Conj(Lambda6)*(15*Lambda5*(6*Lambda6 + 7*
      Lambda7) + Conj(Lambda7)*(5*(10*Lambda1 + 10*Lambda2 + 20*Lambda3 + 22*
      Lambda4 + 3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu) + 3*Sqr(g1))) +
      36*Lambda3*Lambda5*Sqr(g2) + 72*Lambda4*Lambda5*Sqr(g2) + 11.25*Lambda5*
      traceYdAdjYd*Sqr(g2) + 3.75*Lambda5*traceYeAdjYe*Sqr(g2) + 11.25*Lambda5*
      traceYuAdjYu*Sqr(g2) + 2.85*Lambda5*Sqr(g1)*Sqr(g2) + 2.7*traceYdAdjYd*
      Sqr(g1)*Sqr(g2) + 3.3*traceYeAdjYe*Sqr(g1)*Sqr(g2) + 40*Lambda5*
      traceYdAdjYd*Sqr(g3) + 40*Lambda5*traceYuAdjYu*Sqr(g3) - 28*Lambda5*Sqr(
      Lambda1) - 28*Lambda5*Sqr(Lambda2) - 28*Lambda5*Sqr(Lambda3) - 32*Lambda5
      *Sqr(Lambda4) + 6*Conj(Lambda5)*Sqr(Lambda5) + 2*(-2*(37*Lambda1 + 5*
      Lambda2 + 18*Lambda3 + 19*Lambda4 + 15*traceYdAdjYd + 5*traceYeAdjYe) + 6
      *Sqr(g1) + 27*Sqr(g2))*Sqr(Conj(Lambda6)) - 20*Lambda1*Sqr(Conj(Lambda7))
      - 148*Lambda2*Sqr(Conj(Lambda7)) - 72*Lambda3*Sqr(Conj(Lambda7)) - 76*
      Lambda4*Sqr(Conj(Lambda7)) - 60*traceYuAdjYu*Sqr(Conj(Lambda7)) + 12*Sqr(
      g1)*Sqr(Conj(Lambda7)) + 54*Sqr(g2)*Sqr(Conj(Lambda7))));


   return beta_Lambda5;
}

/**
 * Calculates the 3-loop beta function of Lambda5.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

/**
 * Calculates the 4-loop beta function of Lambda5.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

} // namespace flexiblesusy
