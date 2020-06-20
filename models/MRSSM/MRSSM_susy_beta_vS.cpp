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


#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vS.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vS_1_loop(const Susy_traces& susy_traces) const
{


   double beta_vS;

   beta_vS = Re(-2*vS*(AbsSqr(LamSD) + AbsSqr(LamSU)));


   return oneLoop * beta_vS;
}

/**
 * Calculates the 2-loop beta function of vS.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vS_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vS;

   beta_vS = Re(-0.4*vS*(-15*traceYdAdjYd*AbsSqr(LamSD) - 5*traceYeAdjYe*AbsSqr
      (LamSD) - 15*traceYuAdjYu*AbsSqr(LamSU) - 15*AbsSqr(LamSD)*AbsSqr(LamTD)
      - 15*AbsSqr(LamSU)*AbsSqr(LamTU) + 3*AbsSqr(LamSD)*Sqr(g1) + 3*AbsSqr(
      LamSU)*Sqr(g1) + 15*AbsSqr(LamSD)*Sqr(g2) + 15*AbsSqr(LamSU)*Sqr(g2) - 10
      *Sqr(LamSD)*Sqr(Conj(LamSD)) - 10*Sqr(LamSU)*Sqr(Conj(LamSU))));


   return twoLoop * beta_vS;
}

/**
 * Calculates the 3-loop beta function of vS.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vS_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return threeLoop * beta_vS;
}

/**
 * Calculates the 4-loop beta function of vS.
 *
 * @return 4-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vS_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return fourLoop * beta_vS;
}

/**
 * Calculates the 5-loop beta function of vS.
 *
 * @return 5-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_vS_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vS;

   beta_vS = 0;


   return fiveLoop * beta_vS;
}

} // namespace flexiblesusy
