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

// File generated at Sun 26 Aug 2018 14:08:45

#include "THDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda3.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda3;

   beta_Lambda3 = Re(oneOver16PiSqr*(12*Lambda1*Lambda3 + 12*Lambda2*Lambda3 +
      4*Lambda1*Lambda4 + 4*Lambda2*Lambda4 + 6*Lambda3*traceYdAdjYd - 12*
      traceYdAdjYuYuAdjYd + 2*Lambda3*traceYeAdjYe + 6*Lambda3*traceYuAdjYu + 2
      *AbsSqr(Lambda5) + 4*AbsSqr(Lambda7) + 4*(Lambda6 + 2*Lambda7)*Conj(
      Lambda6) + 8*Lambda6*Conj(Lambda7) + 0.27*Quad(g1) + 2.25*Quad(g2) - 1.8*
      Lambda3*Sqr(g1) - 9*Lambda3*Sqr(g2) - 0.9*Sqr(g1)*Sqr(g2) + 4*Sqr(Lambda3
      ) + 2*Sqr(Lambda4)));


   return beta_Lambda3;
}

/**
 * Calculates the 2-loop beta function of Lambda3.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_2_loop(const Susy_traces& susy_traces) const
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


   double beta_Lambda3;

   const double beta_Lambda3_1 = Re(twoLoop*(-32*Lambda1*Lambda3*Lambda4 - 32*
      Lambda2*Lambda3*Lambda4 - 36*Lambda5*Lambda6*Lambda7 - 72*Lambda1*Lambda3
      *traceYdAdjYd - 24*Lambda1*Lambda4*traceYdAdjYd - 13.5*Lambda3*
      traceYdAdjYdYdAdjYd + 12*traceYdAdjYdYdAdjYuYuAdjYd + 15*Lambda3*
      traceYdAdjYuYuAdjYd + 24*traceYdAdjYuYuAdjYdYdAdjYd + 36*
      traceYdAdjYuYuAdjYuYuAdjYd - 24*Lambda1*Lambda3*traceYeAdjYe - 8*Lambda1*
      Lambda4*traceYeAdjYe - 4.5*Lambda3*traceYeAdjYeYeAdjYe - 72*Lambda2*
      Lambda3*traceYuAdjYu - 24*Lambda2*Lambda4*traceYuAdjYu - 13.5*Lambda3*
      traceYuAdjYuYuAdjYu - 12*Cube(Lambda3) - 12*Cube(Lambda4) - 3.537*Power6(
      g1) + 36.375*Power6(g2) + 2.7*Lambda1*Quad(g1) + 2.7*Lambda2*Quad(g1) +
      8.865*Lambda3*Quad(g1) + 0.9*Lambda4*Quad(g1) + 0.45*traceYdAdjYd*Quad(g1
      ) - 2.25*traceYeAdjYe*Quad(g1) - 1.71*traceYuAdjYu*Quad(g1) + 22.5*
      Lambda1*Quad(g2) + 22.5*Lambda2*Quad(g2) - 13.875*Lambda3*Quad(g2) + 7.5*
      Lambda4*Quad(g2) - 2.25*traceYdAdjYd*Quad(g2) - 0.75*traceYeAdjYe*Quad(g2
      ) - 2.25*traceYuAdjYu*Quad(g2) + 14.4*Lambda1*Lambda3*Sqr(g1) + 14.4*
      Lambda2*Lambda3*Sqr(g1) + 4.8*Lambda1*Lambda4*Sqr(g1) + 4.8*Lambda2*
      Lambda4*Sqr(g1) + 1.25*Lambda3*traceYdAdjYd*Sqr(g1) - 0.8*
      traceYdAdjYuYuAdjYd*Sqr(g1) + 3.75*Lambda3*traceYeAdjYe*Sqr(g1) + 4.25*
      Lambda3*traceYuAdjYu*Sqr(g1) + 0.825*Quad(g2)*Sqr(g1) + 0.4*AbsSqr(
      Lambda5)*(-5*(18*Lambda1 + 18*Lambda2 + 9*Lambda3 + 22*Lambda4 + 3*
      traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu) + 6*Sqr(g1)) + 72*Lambda1*
      Lambda3*Sqr(g2) + 72*Lambda2*Lambda3*Sqr(g2) + 36*Lambda1*Lambda4*Sqr(g2)
      + 36*Lambda2*Lambda4*Sqr(g2) - 12*Lambda3*Lambda4*Sqr(g2) + 11.25*Lambda3
      *traceYdAdjYd*Sqr(g2) + 3.75*Lambda3*traceYeAdjYe*Sqr(g2) + 11.25*Lambda3
      *traceYuAdjYu*Sqr(g2) + 4.545*Quad(g1)*Sqr(g2) - 3*Lambda1*Sqr(g1)*Sqr(g2
      ) - 3*Lambda2*Sqr(g1)*Sqr(g2) + 1.65*Lambda3*Sqr(g1)*Sqr(g2) - 1.8*
      Lambda4*Sqr(g1)*Sqr(g2) - 2.7*traceYdAdjYd*Sqr(g1)*Sqr(g2) - 3.3*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) - 6.3*traceYuAdjYu*Sqr(g1)*Sqr(g2) + 0.4*
      Conj(Lambda6)*(3*(Lambda6 + 8*Lambda7)*Sqr(g1) - 5*(62*Lambda1*Lambda6 +
      22*Lambda2*Lambda6 + 30*Lambda3*Lambda6 + 34*Lambda4*Lambda6 + 22*Lambda1
      *Lambda7 - 27*Lambda7*Sqr(g2))) + 40*Lambda3*traceYdAdjYd*Sqr(g3) - 64*
      traceYdAdjYuYuAdjYd*Sqr(g3) + 40*Lambda3*traceYuAdjYu*Sqr(g3) - 60*
      Lambda3*Sqr(Lambda1) - 16*Lambda4*Sqr(Lambda1) - 60*Lambda3*Sqr(Lambda2)
      - 16*Lambda4*Sqr(Lambda2) - 72*Lambda1*Sqr(Lambda3) - 72*Lambda2*Sqr(
      Lambda3) - 4*Lambda4*Sqr(Lambda3) - 12*traceYdAdjYd*Sqr(Lambda3) - 4*
      traceYeAdjYe*Sqr(Lambda3) - 12*traceYuAdjYu*Sqr(Lambda3) + 1.2*Sqr(g1)*
      Sqr(Lambda3) + 6*Sqr(g2)*Sqr(Lambda3) - 28*Lambda1*Sqr(Lambda4) - 28*
      Lambda2*Sqr(Lambda4) - 16*Lambda3*Sqr(Lambda4) - 6*traceYdAdjYd*Sqr(
      Lambda4) - 2*traceYeAdjYe*Sqr(Lambda4) - 6*traceYuAdjYu*Sqr(Lambda4) -
      1.2*Sqr(g1)*Sqr(Lambda4) + 6*Sqr(g2)*Sqr(Lambda4) - 34*Lambda5*Sqr(
      Lambda6) - 34*Lambda5*Sqr(Lambda7)));
   const double beta_Lambda3_2 = Re(-0.4*twoLoop*(10*Conj(Lambda6)*(11*Lambda2*
      Lambda7 + 22*Lambda3*Lambda7 + 11*Lambda4*Lambda7 + 6*Lambda6*
      traceYdAdjYd + 6*Lambda7*traceYdAdjYd + 2*Lambda6*traceYeAdjYe + 2*
      Lambda7*traceYeAdjYe + 6*Lambda7*traceYuAdjYu + 9*Conj(Lambda5)*Conj(
      Lambda7)) + Conj(Lambda7)*(85*Conj(Lambda5)*Conj(Lambda7) - 3*(8*Lambda6
      + Lambda7)*Sqr(g1) + 5*(2*(11*Lambda2*Lambda6 + 22*Lambda3*Lambda6 + 11*
      Lambda4*Lambda6 + 31*Lambda2*Lambda7 + 15*Lambda3*Lambda7 + 17*Lambda4*
      Lambda7 + 11*Lambda1*(Lambda6 + Lambda7) + 6*Lambda6*traceYdAdjYd + 2*
      Lambda6*traceYeAdjYe + 6*Lambda6*traceYuAdjYu + 6*Lambda7*traceYuAdjYu) -
      27*Lambda6*Sqr(g2))) + 85*Conj(Lambda5)*Sqr(Conj(Lambda6))));

   beta_Lambda3 = beta_Lambda3_1 + beta_Lambda3_2;


   return beta_Lambda3;
}

/**
 * Calculates the 3-loop beta function of Lambda3.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

/**
 * Calculates the 4-loop beta function of Lambda3.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda3_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

} // namespace flexiblesusy
