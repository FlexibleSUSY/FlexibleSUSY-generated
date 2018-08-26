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

// File generated at Sun 26 Aug 2018 14:09:04

#include "THDMII_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda6.
 *
 * @return 1-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda6_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(oneOver16PiSqr*(24*Lambda1*Lambda6 + 6*Lambda3*Lambda6 + 8
      *Lambda4*Lambda6 + 6*Lambda3*Lambda7 + 4*Lambda4*Lambda7 + 9*Lambda6*
      traceYdAdjYd + 3*Lambda6*traceYeAdjYe + 3*Lambda6*traceYuAdjYu + 2*Conj(
      Lambda5)*(5*Conj(Lambda6) + Conj(Lambda7)) - 1.8*Lambda6*Sqr(g1) - 9*
      Lambda6*Sqr(g2)));


   return beta_Lambda6;
}

/**
 * Calculates the 2-loop beta function of Lambda6.
 *
 * @return 2-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda6_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda6;

   beta_Lambda6 = Re(twoLoop*(-132*Lambda1*Lambda3*Lambda6 - 36*Lambda2*Lambda3
      *Lambda6 - 140*Lambda1*Lambda4*Lambda6 - 28*Lambda2*Lambda4*Lambda6 - 68*
      Lambda3*Lambda4*Lambda6 - 36*Lambda1*Lambda3*Lambda7 - 36*Lambda2*Lambda3
      *Lambda7 - 28*Lambda1*Lambda4*Lambda7 - 28*Lambda2*Lambda4*Lambda7 - 56*
      Lambda3*Lambda4*Lambda7 - 144*Lambda1*Lambda6*traceYdAdjYd - 18*Lambda3*
      Lambda6*traceYdAdjYd - 24*Lambda4*Lambda6*traceYdAdjYd - 2.25*Lambda6*
      traceYdAdjYdYdAdjYd - 21*Lambda6*traceYdAdjYuYuAdjYd - 48*Lambda1*Lambda6
      *traceYeAdjYe - 6*Lambda3*Lambda6*traceYeAdjYe - 8*Lambda4*Lambda6*
      traceYeAdjYe - 0.75*Lambda6*traceYeAdjYeYeAdjYe - 18*Lambda3*Lambda6*
      traceYuAdjYu - 24*Lambda4*Lambda6*traceYuAdjYu - 36*Lambda3*Lambda7*
      traceYuAdjYu - 24*Lambda4*Lambda7*traceYuAdjYu - 6.75*Lambda6*
      traceYuAdjYuYuAdjYu - 11*Lambda6*AbsSqr(Lambda7) + 8.415*Lambda6*Quad(g1)
      + 1.35*Lambda7*Quad(g1) - 17.625*Lambda6*Quad(g2) + 11.25*Lambda7*Quad(g2
      ) + 21.6*Lambda1*Lambda6*Sqr(g1) + 3.6*Lambda3*Lambda6*Sqr(g1) + 6*
      Lambda4*Lambda6*Sqr(g1) + 7.2*Lambda3*Lambda7*Sqr(g1) + 4.8*Lambda4*
      Lambda7*Sqr(g1) + 1.875*Lambda6*traceYdAdjYd*Sqr(g1) + 5.625*Lambda6*
      traceYeAdjYe*Sqr(g1) + 2.125*Lambda6*traceYuAdjYu*Sqr(g1) - 0.4*Conj(
      Lambda5)*(15*Lambda5*(6*Lambda6 + 7*Lambda7) + Conj(Lambda7)*(10*(5*
      Lambda1 + 5*Lambda2 + 10*Lambda3 + 11*Lambda4 + 3*traceYuAdjYu) + 3*Sqr(
      g1)) + 5*Conj(Lambda6)*(74*Lambda1 + 10*Lambda2 + 36*Lambda3 + 38*Lambda4
       + 15*traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu - 6*Sqr(g1) - 27*
      Sqr(g2))) + 108*Lambda1*Lambda6*Sqr(g2) + 18*Lambda3*Lambda6*Sqr(g2) + 36
      *Lambda4*Lambda6*Sqr(g2) + 36*Lambda3*Lambda7*Sqr(g2) + 18*Lambda4*
      Lambda7*Sqr(g2) + 16.875*Lambda6*traceYdAdjYd*Sqr(g2) + 5.625*Lambda6*
      traceYeAdjYe*Sqr(g2) + 5.625*Lambda6*traceYuAdjYu*Sqr(g2) + 4.35*Lambda6*
      Sqr(g1)*Sqr(g2) + 1.5*Lambda7*Sqr(g1)*Sqr(g2) + 60*Lambda6*traceYdAdjYd*
      Sqr(g3) + 20*Lambda6*traceYuAdjYu*Sqr(g3) - 318*Lambda6*Sqr(Lambda1) + 6*
      Lambda6*Sqr(Lambda2) - 32*Lambda6*Sqr(Lambda3) - 36*Lambda7*Sqr(Lambda3)
      - 34*Lambda6*Sqr(Lambda4) - 34*Lambda7*Sqr(Lambda4) - 42*Conj(Lambda7)*
      Sqr(Lambda6) - 42*Conj(Lambda7)*Sqr(Lambda7) - Conj(Lambda6)*(84*Lambda6*
      Lambda7 + 111*Sqr(Lambda6) + 22*Sqr(Lambda7))));


   return beta_Lambda6;
}

/**
 * Calculates the 3-loop beta function of Lambda6.
 *
 * @return 3-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda6_3_loop(const Susy_traces& susy_traces) const
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
double THDMII_susy_parameters::calc_beta_Lambda6_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda6;

   beta_Lambda6 = 0;


   return beta_Lambda6;
}

} // namespace flexiblesusy
