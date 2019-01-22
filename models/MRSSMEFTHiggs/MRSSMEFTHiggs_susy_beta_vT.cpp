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

// File generated at Tue 22 Jan 2019 13:53:08

#include "MRSSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vT.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vT_1_loop(const Susy_traces& susy_traces) const
{


   double beta_vT;

   beta_vT = Re(oneOver16PiSqr*vT*(-AbsSqr(LamTD) - AbsSqr(LamTU) + 4*Sqr(g2)))
      ;


   return beta_vT;
}

/**
 * Calculates the 2-loop beta function of vT.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vT_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vT;

   beta_vT = Re(-0.06666666666666667*twoLoop*vT*(-45*traceYdAdjYd*AbsSqr(LamTD)
      - 15*traceYeAdjYe*AbsSqr(LamTD) - 30*AbsSqr(LamSD)*AbsSqr(LamTD) - 45*
      traceYuAdjYu*AbsSqr(LamTU) - 30*AbsSqr(LamSU)*AbsSqr(LamTU) + 220*Quad(g2
      ) + 9*AbsSqr(LamTD)*Sqr(g1) + 9*AbsSqr(LamTU)*Sqr(g1) + 45*AbsSqr(LamTD)*
      Sqr(g2) + 45*AbsSqr(LamTU)*Sqr(g2) - 45*Sqr(LamTD)*Sqr(Conj(LamTD)) - 45*
      Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_vT;
}

/**
 * Calculates the 3-loop beta function of vT.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vT_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vT;

   beta_vT = 0;


   return beta_vT;
}

/**
 * Calculates the 4-loop beta function of vT.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vT_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vT;

   beta_vT = 0;


   return beta_vT;
}

/**
 * Calculates the 5-loop beta function of vT.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vT_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vT;

   beta_vT = 0;


   return beta_vT;
}

} // namespace flexiblesusy
