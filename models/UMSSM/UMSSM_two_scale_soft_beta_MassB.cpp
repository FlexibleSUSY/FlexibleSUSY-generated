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

// File generated at Mon 8 Jun 2015 17:49:20

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassB.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassB;

   beta_MassB = Re(13.2*MassB*oneOver16PiSqr*Sqr(g1));


   return beta_MassB;
}

/**
 * Calculates the two-loop beta function of MassB.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassB_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_MassB;

   beta_MassB = Re(0.08*twoLoop*Sqr(g1)*(70*traceAdjYdTYd + 90*
      traceAdjYeTYe + 130*traceAdjYuTYu - 70*MassB*traceYdAdjYd - 90*MassB*
      traceYeAdjYe - 130*MassB*traceYuAdjYu + 398*MassB*Sqr(g1) + 135*MassB*Sqr
      (g2) + 135*MassWB*Sqr(g2) + 440*MassB*Sqr(g3) + 440*MassG*Sqr(g3) + 60*
      MassB*Sqr(gp)*Sqr(Qd) + 60*MassU*Sqr(gp)*Sqr(Qd) + 180*MassB*Sqr(gp)*Sqr(
      Qe) + 180*MassU*Sqr(gp)*Sqr(Qe) + 30*MassB*Sqr(gp)*Sqr(QHd) + 30*MassU*
      Sqr(gp)*Sqr(QHd) + 30*MassB*Sqr(gp)*Sqr(QHu) + 30*MassU*Sqr(gp)*Sqr(QHu)
      + 90*MassB*Sqr(gp)*Sqr(Ql) + 90*MassU*Sqr(gp)*Sqr(Ql) + 30*MassB*Sqr(gp)*
      Sqr(Qq) + 30*MassU*Sqr(gp)*Sqr(Qq) + 240*MassB*Sqr(gp)*Sqr(Qu) + 240*
      MassU*Sqr(gp)*Sqr(Qu) - 30*Conj(Lambdax)*(MassB*Lambdax - TLambdax)));


   return beta_MassB;
}

} // namespace flexiblesusy
