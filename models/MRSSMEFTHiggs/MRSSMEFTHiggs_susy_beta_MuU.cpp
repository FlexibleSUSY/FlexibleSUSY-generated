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
 * Calculates the 1-loop beta function of MuU.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_MuU_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MuU;

   beta_MuU = Re(-0.2*MuU*(-15*traceYuAdjYu - 10*AbsSqr(LamSU) - 15*AbsSqr(
      LamTU) + 3*Sqr(g1) + 15*Sqr(g2)));


   return oneLoop * beta_MuU;
}

/**
 * Calculates the 2-loop beta function of MuU.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_MuU_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_MuU;

   beta_MuU = Re(0.1*MuU*(-30*traceYdAdjYuYuAdjYd - 90*traceYuAdjYuYuAdjYu - 30
      *traceYuAdjYu*AbsSqr(LamSU) - 40*AbsSqr(LamSD)*AbsSqr(LamSU) - 45*
      traceYuAdjYu*AbsSqr(LamTU) - 60*AbsSqr(LamSU)*AbsSqr(LamTU) - 30*AbsSqr(
      LamTD)*AbsSqr(LamTU) + 45*Quad(g1) + 165*Quad(g2) + 8*traceYuAdjYu*Sqr(g1
      ) + 120*AbsSqr(LamTU)*Sqr(g2) + 18*Sqr(g1)*Sqr(g2) + 160*traceYuAdjYu*Sqr
      (g3) - 60*Sqr(LamSU)*Sqr(Conj(LamSU)) - 75*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return twoLoop * beta_MuU;
}

/**
 * Calculates the 3-loop beta function of MuU.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_MuU_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuU;

   beta_MuU = 0;


   return threeLoop * beta_MuU;
}

/**
 * Calculates the 4-loop beta function of MuU.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_MuU_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuU;

   beta_MuU = 0;


   return fourLoop * beta_MuU;
}

/**
 * Calculates the 5-loop beta function of MuU.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_MuU_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuU;

   beta_MuU = 0;


   return fiveLoop * beta_MuU;
}

} // namespace flexiblesusy
