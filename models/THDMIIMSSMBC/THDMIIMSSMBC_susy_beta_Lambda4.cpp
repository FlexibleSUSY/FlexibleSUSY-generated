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


#include "THDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda4.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda4;

   beta_Lambda4 = Re(0.2*(20*Lambda1*Lambda4 + 20*Lambda2*Lambda4 + 40*Lambda3*
      Lambda4 + 30*Lambda4*traceYdAdjYd + 60*traceYdAdjYuYuAdjYd + 10*Lambda4*
      traceYeAdjYe + 30*Lambda4*traceYuAdjYu + 40*AbsSqr(Lambda5) + 50*AbsSqr(
      Lambda6) + 50*AbsSqr(Lambda7) + 10*Lambda7*Conj(Lambda6) + 10*Lambda6*
      Conj(Lambda7) - 9*Lambda4*Sqr(g1) - 45*Lambda4*Sqr(g2) + 9*Sqr(g1)*Sqr(g2
      ) + 20*Sqr(Lambda4)));


   return oneLoop * beta_Lambda4;
}

/**
 * Calculates the 2-loop beta function of Lambda4.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_2_loop(const Susy_traces& susy_traces) const
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

   beta_Lambda4 = Re(0.005*(-16000*Lambda1*Lambda3*Lambda4 - 16000*Lambda2*
      Lambda3*Lambda4 - 9000*Lambda5*Lambda6*Lambda7 - 4800*Lambda1*Lambda4*
      traceYdAdjYd - 4800*Lambda3*Lambda4*traceYdAdjYd - 2700*Lambda4*
      traceYdAdjYdYdAdjYd - 2400*traceYdAdjYdYdAdjYuYuAdjYd - 4800*Lambda3*
      traceYdAdjYuYuAdjYd - 6600*Lambda4*traceYdAdjYuYuAdjYd - 2400*
      traceYdAdjYuYuAdjYdYdAdjYd - 4800*traceYdAdjYuYuAdjYuYuAdjYd - 1600*
      Lambda1*Lambda4*traceYeAdjYe - 1600*Lambda3*Lambda4*traceYeAdjYe - 900*
      Lambda4*traceYeAdjYeYeAdjYe - 4800*Lambda2*Lambda4*traceYuAdjYu - 4800*
      Lambda3*Lambda4*traceYuAdjYu - 2700*Lambda4*traceYuAdjYuYuAdjYu - 9600*
      Lambda1*AbsSqr(Lambda5) - 9600*Lambda2*AbsSqr(Lambda5) - 9600*Lambda3*
      AbsSqr(Lambda5) - 5200*Lambda4*AbsSqr(Lambda5) - 4800*traceYdAdjYd*AbsSqr
      (Lambda5) - 1600*traceYeAdjYe*AbsSqr(Lambda5) - 4800*traceYuAdjYu*AbsSqr(
      Lambda5) - 28400*Lambda1*AbsSqr(Lambda6) - 4000*Lambda2*AbsSqr(Lambda6) -
      13800*Lambda3*AbsSqr(Lambda6) - 13000*Lambda4*AbsSqr(Lambda6) - 12000*
      traceYdAdjYd*AbsSqr(Lambda6) - 4000*traceYeAdjYe*AbsSqr(Lambda6) - 4000*
      Lambda1*AbsSqr(Lambda7) - 28400*Lambda2*AbsSqr(Lambda7) - 13800*Lambda3*
      AbsSqr(Lambda7) - 13000*Lambda4*AbsSqr(Lambda7) - 12000*traceYuAdjYu*
      AbsSqr(Lambda7) - 3400*Lambda1*Lambda7*Conj(Lambda6) - 3400*Lambda2*
      Lambda7*Conj(Lambda6) - 7400*Lambda3*Lambda7*Conj(Lambda6) - 15400*
      Lambda4*Lambda7*Conj(Lambda6) - 1200*Lambda7*traceYdAdjYd*Conj(Lambda6) -
      400*Lambda7*traceYeAdjYe*Conj(Lambda6) - 1200*Lambda7*traceYuAdjYu*Conj(
      Lambda6) - 3400*Lambda1*Lambda6*Conj(Lambda7) - 3400*Lambda2*Lambda6*Conj
      (Lambda7) - 7400*Lambda3*Lambda6*Conj(Lambda7) - 15400*Lambda4*Lambda6*
      Conj(Lambda7) - 1200*Lambda6*traceYdAdjYd*Conj(Lambda7) - 400*Lambda6*
      traceYeAdjYe*Conj(Lambda7) - 1200*Lambda6*traceYuAdjYu*Conj(Lambda7) -
      9000*Conj(Lambda5)*Conj(Lambda6)*Conj(Lambda7) + 1413*Lambda4*Quad(g1) -
      5775*Lambda4*Quad(g2) + 960*Lambda1*Lambda4*Sqr(g1) + 960*Lambda2*Lambda4
      *Sqr(g1) + 480*Lambda3*Lambda4*Sqr(g1) + 250*Lambda4*traceYdAdjYd*Sqr(g1)
      + 160*traceYdAdjYuYuAdjYd*Sqr(g1) + 750*Lambda4*traceYeAdjYe*Sqr(g1) +
      850*Lambda4*traceYuAdjYu*Sqr(g1) + 1920*AbsSqr(Lambda5)*Sqr(g1) + 1680*
      AbsSqr(Lambda6)*Sqr(g1) + 1680*AbsSqr(Lambda7)*Sqr(g1) + 480*Lambda7*Conj
      (Lambda6)*Sqr(g1) + 480*Lambda6*Conj(Lambda7)*Sqr(g1) - 1680*Quad(g2)*Sqr
      (g1) + 7200*Lambda3*Lambda4*Sqr(g2) + 2250*Lambda4*traceYdAdjYd*Sqr(g2) +
      750*Lambda4*traceYeAdjYe*Sqr(g2) + 2250*Lambda4*traceYuAdjYu*Sqr(g2) +
      10800*AbsSqr(Lambda5)*Sqr(g2) + 10800*AbsSqr(Lambda6)*Sqr(g2) + 10800*
      AbsSqr(Lambda7)*Sqr(g2) - 2628*Quad(g1)*Sqr(g2) + 1200*Lambda1*Sqr(g1)*
      Sqr(g2) + 1200*Lambda2*Sqr(g1)*Sqr(g2) + 240*Lambda3*Sqr(g1)*Sqr(g2) +
      1530*Lambda4*Sqr(g1)*Sqr(g2) + 1080*traceYdAdjYd*Sqr(g1)*Sqr(g2) + 1320*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) + 2520*traceYuAdjYu*Sqr(g1)*Sqr(g2) + 8000*
      Lambda4*traceYdAdjYd*Sqr(g3) + 12800*traceYdAdjYuYuAdjYd*Sqr(g3) + 8000*
      Lambda4*traceYuAdjYu*Sqr(g3) - 5600*Lambda4*Sqr(Lambda1) - 5600*Lambda4*
      Sqr(Lambda2) - 5600*Lambda4*Sqr(Lambda3) - 8000*Lambda1*Sqr(Lambda4) -
      8000*Lambda2*Sqr(Lambda4) - 5600*Lambda3*Sqr(Lambda4) - 2400*traceYdAdjYd
      *Sqr(Lambda4) - 800*traceYeAdjYe*Sqr(Lambda4) - 2400*traceYuAdjYu*Sqr(
      Lambda4) + 960*Sqr(g1)*Sqr(Lambda4) + 3600*Sqr(g2)*Sqr(Lambda4) - 7700*
      Lambda5*Sqr(Lambda6) - 7700*Lambda5*Sqr(Lambda7) - 7700*Conj(Lambda5)*Sqr
      (Conj(Lambda6)) - 7700*Conj(Lambda5)*Sqr(Conj(Lambda7))));


   return twoLoop * beta_Lambda4;
}

/**
 * Calculates the 3-loop beta function of Lambda4.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return threeLoop * beta_Lambda4;
}

/**
 * Calculates the 4-loop beta function of Lambda4.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return fourLoop * beta_Lambda4;
}

/**
 * Calculates the 5-loop beta function of Lambda4.
 *
 * @return 5-loop beta function
 */
double THDMIIMSSMBC_susy_parameters::calc_beta_Lambda4_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda4;

   beta_Lambda4 = 0;


   return fiveLoop * beta_Lambda4;
}

} // namespace flexiblesusy
