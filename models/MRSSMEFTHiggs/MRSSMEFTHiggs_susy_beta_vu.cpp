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


#include "MRSSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vu.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vu;

   beta_vu = Re(0.1*vu*(-30*traceYuAdjYu - 10*AbsSqr(LamSU) - 15*AbsSqr(LamTU)
      + 3*Sqr(g1) + 15*Sqr(g2)));


   return oneLoop * beta_vu;
}

/**
 * Calculates the 2-loop beta function of vu.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = Re(-0.025*vu*(-120*traceYdAdjYuYuAdjYd - 360*traceYuAdjYuYuAdjYu -
      80*AbsSqr(LamSD)*AbsSqr(LamSU) - 120*AbsSqr(LamSU)*AbsSqr(LamTU) - 60*
      AbsSqr(LamTD)*AbsSqr(LamTU) + 45*Quad(g1) + 145*Quad(g2) + 68*
      traceYuAdjYu*Sqr(g1) + 12*AbsSqr(LamSU)*Sqr(g1) + 18*AbsSqr(LamTU)*Sqr(g1
      ) + 180*traceYuAdjYu*Sqr(g2) + 60*AbsSqr(LamSU)*Sqr(g2) + 330*AbsSqr(
      LamTU)*Sqr(g2) + 18*Sqr(g1)*Sqr(g2) + 640*traceYuAdjYu*Sqr(g3) - 120*Sqr(
      LamSU)*Sqr(Conj(LamSU)) - 150*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return twoLoop * beta_vu;
}

/**
 * Calculates the 3-loop beta function of vu.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return threeLoop * beta_vu;
}

/**
 * Calculates the 4-loop beta function of vu.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return fourLoop * beta_vu;
}

/**
 * Calculates the 5-loop beta function of vu.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_vu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return fiveLoop * beta_vu;
}

} // namespace flexiblesusy
