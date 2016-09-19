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

// File generated at Mon 19 Sep 2016 09:37:07

#include "MRSSMtower_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vu.
 *
 * @return one-loop beta function
 */
double MRSSMtower_susy_parameters::calc_beta_vu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vu;

   beta_vu = Re(0.1*oneOver16PiSqr*vu*(-10*AbsSqr(LamSU) + 3*(-10*
      traceYuAdjYu - 5*AbsSqr(LamTU) + Sqr(g1) + 5*Sqr(g2))));


   return beta_vu;
}

/**
 * Calculates the two-loop beta function of vu.
 *
 * @return two-loop beta function
 */
double MRSSMtower_susy_parameters::calc_beta_vu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = Re(-0.025*twoLoop*vu*(45*Power(g1,4) + 145*Power(g2,4) - 120
      *traceYdAdjYuYuAdjYd - 360*traceYuAdjYuYuAdjYu + 68*traceYuAdjYu*Sqr(g1)
      + 180*traceYuAdjYu*Sqr(g2) + 18*Sqr(g1)*Sqr(g2) + 6*AbsSqr(LamTU)*(-10*
      AbsSqr(LamTD) + 3*Sqr(g1) + 55*Sqr(g2)) + 4*AbsSqr(LamSU)*(-20*AbsSqr(
      LamSD) + 3*(-10*AbsSqr(LamTU) + Sqr(g1) + 5*Sqr(g2))) + 640*traceYuAdjYu*
      Sqr(g3) - 120*Sqr(LamSU)*Sqr(Conj(LamSU)) - 150*Sqr(LamTU)*Sqr(Conj(LamTU
      ))));


   return beta_vu;
}

/**
 * Calculates the three-loop beta function of vu.
 *
 * @return three-loop beta function
 */
double MRSSMtower_susy_parameters::calc_beta_vu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return beta_vu;
}

} // namespace flexiblesusy