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

// File generated at Wed 16 Oct 2019 21:21:44

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda5.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda5;

   beta_Lambda5 = Re(0.2*oneOver16PiSqr*(20*Lambda1*Lambda5 + 20*Lambda2*
      Lambda5 + 40*Lambda3*Lambda5 + 60*Lambda4*Lambda5 + 30*Lambda5*
      traceYdAdjYd + 10*Lambda5*traceYeAdjYe + 30*Lambda5*traceYuAdjYu + 20*
      Conj(Lambda6)*Conj(Lambda7) - 9*Lambda5*Sqr(g1) - 45*Lambda5*Sqr(g2) + 50
      *Sqr(Conj(Lambda6)) + 50*Sqr(Conj(Lambda7))));


   return beta_Lambda5;
}

/**
 * Calculates the 2-loop beta function of Lambda5.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda5;

   beta_Lambda5 = Re(0.005*twoLoop*(-16000*Lambda1*Lambda3*Lambda5 - 16000*
      Lambda2*Lambda3*Lambda5 - 17600*Lambda1*Lambda4*Lambda5 - 17600*Lambda2*
      Lambda4*Lambda5 - 15200*Lambda3*Lambda4*Lambda5 - 4800*Lambda1*Lambda5*
      traceYdAdjYd - 4800*Lambda3*Lambda5*traceYdAdjYd - 7200*Lambda4*Lambda5*
      traceYdAdjYd - 1500*Lambda5*traceYdAdjYdYdAdjYd - 6600*Lambda5*
      traceYdAdjYuYuAdjYd - 1600*Lambda1*Lambda5*traceYeAdjYe - 1600*Lambda3*
      Lambda5*traceYeAdjYe - 2400*Lambda4*Lambda5*traceYeAdjYe - 500*Lambda5*
      traceYeAdjYeYeAdjYe - 4800*Lambda2*Lambda5*traceYuAdjYu - 4800*Lambda3*
      Lambda5*traceYuAdjYu - 7200*Lambda4*Lambda5*traceYuAdjYu - 1500*Lambda5*
      traceYuAdjYuYuAdjYu - 13800*Lambda5*AbsSqr(Lambda6) - 13800*Lambda5*
      AbsSqr(Lambda7) - 16200*Lambda5*Lambda7*Conj(Lambda6) - 16200*Lambda5*
      Lambda6*Conj(Lambda7) - 6800*Lambda1*Conj(Lambda6)*Conj(Lambda7) - 6800*
      Lambda2*Conj(Lambda6)*Conj(Lambda7) - 14800*Lambda3*Conj(Lambda6)*Conj(
      Lambda7) - 16400*Lambda4*Conj(Lambda6)*Conj(Lambda7) - 2400*traceYdAdjYd*
      Conj(Lambda6)*Conj(Lambda7) - 800*traceYeAdjYe*Conj(Lambda6)*Conj(Lambda7
      ) - 2400*traceYuAdjYu*Conj(Lambda6)*Conj(Lambda7) + 1533*Lambda5*Quad(g1)
      + 90*traceYdAdjYd*Quad(g1) - 450*traceYeAdjYe*Quad(g1) - 4775*Lambda5*
      Quad(g2) - 450*traceYdAdjYd*Quad(g2) - 150*traceYeAdjYe*Quad(g2) - 480*
      Lambda1*Lambda5*Sqr(g1) - 480*Lambda2*Lambda5*Sqr(g1) + 1920*Lambda3*
      Lambda5*Sqr(g1) + 2880*Lambda4*Lambda5*Sqr(g1) + 250*Lambda5*traceYdAdjYd
      *Sqr(g1) + 750*Lambda5*traceYeAdjYe*Sqr(g1) + 850*Lambda5*traceYuAdjYu*
      Sqr(g1) - 480*Conj(Lambda6)*Conj(Lambda7)*Sqr(g1) + 7200*Lambda3*Lambda5*
      Sqr(g2) + 14400*Lambda4*Lambda5*Sqr(g2) + 2250*Lambda5*traceYdAdjYd*Sqr(
      g2) + 750*Lambda5*traceYeAdjYe*Sqr(g2) + 2250*Lambda5*traceYuAdjYu*Sqr(g2
      ) + 570*Lambda5*Sqr(g1)*Sqr(g2) + 540*traceYdAdjYd*Sqr(g1)*Sqr(g2) + 660*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) + 8000*Lambda5*traceYdAdjYd*Sqr(g3) + 8000*
      Lambda5*traceYuAdjYu*Sqr(g3) - 5600*Lambda5*Sqr(Lambda1) - 5600*Lambda5*
      Sqr(Lambda2) - 5600*Lambda5*Sqr(Lambda3) - 6400*Lambda5*Sqr(Lambda4) +
      1200*Conj(Lambda5)*Sqr(Lambda5) - 28400*Lambda1*Sqr(Conj(Lambda6)) - 4000
      *Lambda2*Sqr(Conj(Lambda6)) - 13800*Lambda3*Sqr(Conj(Lambda6)) - 14600*
      Lambda4*Sqr(Conj(Lambda6)) - 12000*traceYdAdjYd*Sqr(Conj(Lambda6)) - 4000
      *traceYeAdjYe*Sqr(Conj(Lambda6)) + 2400*Sqr(g1)*Sqr(Conj(Lambda6)) +
      10800*Sqr(g2)*Sqr(Conj(Lambda6)) - 4000*Lambda1*Sqr(Conj(Lambda7)) -
      28400*Lambda2*Sqr(Conj(Lambda7)) - 13800*Lambda3*Sqr(Conj(Lambda7)) -
      14600*Lambda4*Sqr(Conj(Lambda7)) - 12000*traceYuAdjYu*Sqr(Conj(Lambda7))
      + 2400*Sqr(g1)*Sqr(Conj(Lambda7)) + 10800*Sqr(g2)*Sqr(Conj(Lambda7))));


   return beta_Lambda5;
}

/**
 * Calculates the 3-loop beta function of Lambda5.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_3_loop(const Susy_traces& susy_traces) const
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
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

/**
 * Calculates the 5-loop beta function of Lambda5.
 *
 * @return 5-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda5_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda5;

   beta_Lambda5 = 0;


   return beta_Lambda5;
}

} // namespace flexiblesusy
