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


#include "THDMII_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda3.
 *
 * @return 1-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda3_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;


   double beta_Lambda3;

   beta_Lambda3 = Re(0.01*(1200*Lambda1*Lambda3 + 1200*Lambda2*Lambda3 + 400*
      Lambda1*Lambda4 + 400*Lambda2*Lambda4 + 600*Lambda3*traceYdAdjYd - 1200*
      traceYdAdjYuYuAdjYd + 200*Lambda3*traceYeAdjYe + 600*Lambda3*traceYuAdjYu
       + 200*AbsSqr(Lambda5) + 400*AbsSqr(Lambda6) + 400*AbsSqr(Lambda7) + 800*
      Lambda7*Conj(Lambda6) + 800*Lambda6*Conj(Lambda7) + 27*Quad(g1) + 225*
      Quad(g2) - 180*Lambda3*Sqr(g1) - 900*Lambda3*Sqr(g2) - 90*Sqr(g1)*Sqr(g2)
      + 400*Sqr(Lambda3) + 200*Sqr(Lambda4)));


   return oneLoop * beta_Lambda3;
}

/**
 * Calculates the 2-loop beta function of Lambda3.
 *
 * @return 2-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda3_2_loop(const Susy_traces& susy_traces) const
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

   const double beta_Lambda3_1 = Re(0.001*(-32000*Lambda1*Lambda3*Lambda4 -
      32000*Lambda2*Lambda3*Lambda4 - 33000*Lambda5*Lambda6*Lambda7 - 72000*
      Lambda1*Lambda3*traceYdAdjYd - 24000*Lambda1*Lambda4*traceYdAdjYd - 13500
      *Lambda3*traceYdAdjYdYdAdjYd + 12000*traceYdAdjYdYdAdjYuYuAdjYd + 15000*
      Lambda3*traceYdAdjYuYuAdjYd + 24000*traceYdAdjYuYuAdjYdYdAdjYd + 36000*
      traceYdAdjYuYuAdjYuYuAdjYd - 24000*Lambda1*Lambda3*traceYeAdjYe - 8000*
      Lambda1*Lambda4*traceYeAdjYe - 4500*Lambda3*traceYeAdjYeYeAdjYe - 72000*
      Lambda2*Lambda3*traceYuAdjYu - 24000*Lambda2*Lambda4*traceYuAdjYu - 13500
      *Lambda3*traceYuAdjYuYuAdjYu - 36000*Lambda1*AbsSqr(Lambda5) - 36000*
      Lambda2*AbsSqr(Lambda5) - 18000*Lambda3*AbsSqr(Lambda5) - 44000*Lambda4*
      AbsSqr(Lambda5) - 6000*traceYdAdjYd*AbsSqr(Lambda5) - 2000*traceYeAdjYe*
      AbsSqr(Lambda5) - 6000*traceYuAdjYu*AbsSqr(Lambda5) - 118000*Lambda1*
      AbsSqr(Lambda6) - 44000*Lambda2*AbsSqr(Lambda6) - 57000*Lambda3*AbsSqr(
      Lambda6) - 65000*Lambda4*AbsSqr(Lambda6) - 41000*Lambda1*Lambda7*Conj(
      Lambda6) - 12000*Cube(Lambda3) - 12000*Cube(Lambda4) - 3537*Power6(g1) +
      36375*Power6(g2) + 2700*Lambda1*Quad(g1) + 2700*Lambda2*Quad(g1) + 8865*
      Lambda3*Quad(g1) + 900*Lambda4*Quad(g1) + 450*traceYdAdjYd*Quad(g1) -
      2250*traceYeAdjYe*Quad(g1) - 1710*traceYuAdjYu*Quad(g1) + 22500*Lambda1*
      Quad(g2) + 22500*Lambda2*Quad(g2) - 13875*Lambda3*Quad(g2) + 7500*Lambda4
      *Quad(g2) - 2250*traceYdAdjYd*Quad(g2) - 750*traceYeAdjYe*Quad(g2) - 2250
      *traceYuAdjYu*Quad(g2) + 14400*Lambda1*Lambda3*Sqr(g1) + 14400*Lambda2*
      Lambda3*Sqr(g1) + 4800*Lambda1*Lambda4*Sqr(g1) + 4800*Lambda2*Lambda4*Sqr
      (g1) + 1250*Lambda3*traceYdAdjYd*Sqr(g1) - 800*traceYdAdjYuYuAdjYd*Sqr(g1
      ) + 3750*Lambda3*traceYeAdjYe*Sqr(g1) + 4250*Lambda3*traceYuAdjYu*Sqr(g1)
      + 2400*AbsSqr(Lambda5)*Sqr(g1) + 1200*AbsSqr(Lambda6)*Sqr(g1) + 9600*
      Lambda7*Conj(Lambda6)*Sqr(g1) + 825*Quad(g2)*Sqr(g1) + 72000*Lambda1*
      Lambda3*Sqr(g2) + 72000*Lambda2*Lambda3*Sqr(g2) + 36000*Lambda1*Lambda4*
      Sqr(g2) + 36000*Lambda2*Lambda4*Sqr(g2) - 12000*Lambda3*Lambda4*Sqr(g2) +
      11250*Lambda3*traceYdAdjYd*Sqr(g2) + 3750*Lambda3*traceYeAdjYe*Sqr(g2) +
      11250*Lambda3*traceYuAdjYu*Sqr(g2) + 54000*Lambda7*Conj(Lambda6)*Sqr(g2)
      + 4545*Quad(g1)*Sqr(g2) - 3000*Lambda1*Sqr(g1)*Sqr(g2) - 3000*Lambda2*Sqr
      (g1)*Sqr(g2) + 1650*Lambda3*Sqr(g1)*Sqr(g2) - 1800*Lambda4*Sqr(g1)*Sqr(g2
      ) - 2700*traceYdAdjYd*Sqr(g1)*Sqr(g2) - 3300*traceYeAdjYe*Sqr(g1)*Sqr(g2)
      - 6300*traceYuAdjYu*Sqr(g1)*Sqr(g2) + 40000*Lambda3*traceYdAdjYd*Sqr(g3)
      - 64000*traceYdAdjYuYuAdjYd*Sqr(g3) + 40000*Lambda3*traceYuAdjYu*Sqr(g3)
      - 60000*Lambda3*Sqr(Lambda1) - 16000*Lambda4*Sqr(Lambda1) - 60000*Lambda3
      *Sqr(Lambda2) - 16000*Lambda4*Sqr(Lambda2) - 72000*Lambda1*Sqr(Lambda3) -
      72000*Lambda2*Sqr(Lambda3) - 4000*Lambda4*Sqr(Lambda3) - 12000*
      traceYdAdjYd*Sqr(Lambda3) - 4000*traceYeAdjYe*Sqr(Lambda3) - 12000*
      traceYuAdjYu*Sqr(Lambda3) + 1200*Sqr(g1)*Sqr(Lambda3) + 6000*Sqr(g2)*Sqr(
      Lambda3) - 28000*Lambda1*Sqr(Lambda4) - 28000*Lambda2*Sqr(Lambda4) -
      16000*Lambda3*Sqr(Lambda4) - 6000*traceYdAdjYd*Sqr(Lambda4) - 2000*
      traceYeAdjYe*Sqr(Lambda4) - 6000*traceYuAdjYu*Sqr(Lambda4) - 1200*Sqr(g1)
      *Sqr(Lambda4) + 6000*Sqr(g2)*Sqr(Lambda4) - 32500*Lambda5*Sqr(Lambda6) -
      32500*Lambda5*Sqr(Lambda7)));
   const double beta_Lambda3_2 = Re(0.1*(-240*traceYdAdjYd*AbsSqr(Lambda6) - 80
      *traceYeAdjYe*AbsSqr(Lambda6) - 440*Lambda1*AbsSqr(Lambda7) - 1180*
      Lambda2*AbsSqr(Lambda7) - 570*Lambda3*AbsSqr(Lambda7) - 650*Lambda4*
      AbsSqr(Lambda7) - 240*traceYuAdjYu*AbsSqr(Lambda7) - 410*Lambda2*Lambda7*
      Conj(Lambda6) - 850*Lambda3*Lambda7*Conj(Lambda6) - 410*Lambda4*Lambda7*
      Conj(Lambda6) - 240*Lambda7*traceYdAdjYd*Conj(Lambda6) - 80*Lambda7*
      traceYeAdjYe*Conj(Lambda6) - 240*Lambda7*traceYuAdjYu*Conj(Lambda6) - 410
      *Lambda1*Lambda6*Conj(Lambda7) - 410*Lambda2*Lambda6*Conj(Lambda7) - 850*
      Lambda3*Lambda6*Conj(Lambda7) - 410*Lambda4*Lambda6*Conj(Lambda7) - 240*
      Lambda6*traceYdAdjYd*Conj(Lambda7) - 80*Lambda6*traceYeAdjYe*Conj(Lambda7
      ) - 240*Lambda6*traceYuAdjYu*Conj(Lambda7) - 330*Conj(Lambda5)*Conj(
      Lambda6)*Conj(Lambda7) + 12*AbsSqr(Lambda7)*Sqr(g1) + 96*Lambda6*Conj(
      Lambda7)*Sqr(g1) + 540*Lambda6*Conj(Lambda7)*Sqr(g2) - 325*Conj(Lambda5)*
      Sqr(Conj(Lambda6)) - 325*Conj(Lambda5)*Sqr(Conj(Lambda7))));

   beta_Lambda3 = beta_Lambda3_1 + beta_Lambda3_2;


   return twoLoop * beta_Lambda3;
}

/**
 * Calculates the 3-loop beta function of Lambda3.
 *
 * @return 3-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return threeLoop * beta_Lambda3;
}

/**
 * Calculates the 4-loop beta function of Lambda3.
 *
 * @return 4-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda3_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return fourLoop * beta_Lambda3;
}

/**
 * Calculates the 5-loop beta function of Lambda3.
 *
 * @return 5-loop beta function
 */
double THDMII_susy_parameters::calc_beta_Lambda3_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return fiveLoop * beta_Lambda3;
}

} // namespace flexiblesusy
