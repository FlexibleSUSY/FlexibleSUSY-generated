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

// File generated at Tue 10 Oct 2017 21:35:04

#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vd.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = Re(0.1*oneOver16PiSqr*vd*(-30*traceYdAdjYd - 10*traceYeAdjYe
      - 10*AbsSqr(LamSD) - 15*AbsSqr(LamTD) + 3*Sqr(g1) + 15*Sqr(g2)));


   return beta_vd;
}

/**
 * Calculates the 2-loop beta function of vd.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_vd;

   beta_vd = Re(-0.025*twoLoop*vd*(-360*traceYdAdjYdYdAdjYd - 120*
      traceYdAdjYuYuAdjYd - 120*traceYeAdjYeYeAdjYe + 45*Quad(g1) + 145*Quad(g2
      ) + 20*traceYdAdjYd*Sqr(g1) + 60*traceYeAdjYe*Sqr(g1) + 180*traceYdAdjYd*
      Sqr(g2) + 60*traceYeAdjYe*Sqr(g2) + 18*Sqr(g1)*Sqr(g2) + 6*AbsSqr(LamTD)*
      (-10*AbsSqr(LamTU) + 3*Sqr(g1) + 55*Sqr(g2)) + 4*AbsSqr(LamSD)*(-20*
      AbsSqr(LamSU) + 3*(-10*AbsSqr(LamTD) + Sqr(g1) + 5*Sqr(g2))) + 640*
      traceYdAdjYd*Sqr(g3) - 120*Sqr(LamSD)*Sqr(Conj(LamSD)) - 150*Sqr(LamTD)*
      Sqr(Conj(LamTD))));


   return beta_vd;
}

/**
 * Calculates the 3-loop beta function of vd.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vd;

   beta_vd = 0;


   return beta_vd;
}

} // namespace flexiblesusy
