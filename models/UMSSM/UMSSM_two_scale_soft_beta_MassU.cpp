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

// File generated at Fri 26 Jun 2015 19:03:50

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of MassU.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassU_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);


   double beta_MassU;

   beta_MassU = Re(2*MassU*oneOver16PiSqr*Sqr(gp)*(9*Sqr(Qd) + 3*Sqr(Qe)
      + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu))
      );


   return beta_MassU;
}

/**
 * Calculates the two-loop beta function of MassU.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassU_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_MassU;

   beta_MassU = Re(0.8*twoLoop*Sqr(gp)*(180*MassU*Power(Qd,4)*Sqr(gp) +
      60*MassU*Power(Qe,4)*Sqr(gp) + 40*MassU*Power(QHd,4)*Sqr(gp) + 40*MassU*
      Power(QHu,4)*Sqr(gp) + 120*MassU*Power(Ql,4)*Sqr(gp) + 360*MassU*Power(Qq
      ,4)*Sqr(gp) + 20*MassU*Power(Qs,4)*Sqr(gp) + 180*MassU*Power(Qu,4)*Sqr(gp
      ) + 30*traceAdjYdTYd*Sqr(Qd) + 6*MassB*Sqr(g1)*Sqr(Qd) + 6*MassU*Sqr(g1)*
      Sqr(Qd) + 120*MassG*Sqr(g3)*Sqr(Qd) + 120*MassU*Sqr(g3)*Sqr(Qd) + 10*
      traceAdjYeTYe*Sqr(Qe) - 10*MassU*traceYeAdjYe*Sqr(Qe) + 18*MassB*Sqr(g1)*
      Sqr(Qe) + 18*MassU*Sqr(g1)*Sqr(Qe) + 30*traceAdjYdTYd*Sqr(QHd) + 10*
      traceAdjYeTYe*Sqr(QHd) - 10*MassU*traceYeAdjYe*Sqr(QHd) + 3*MassB*Sqr(g1)
      *Sqr(QHd) + 3*MassU*Sqr(g1)*Sqr(QHd) + 15*MassU*Sqr(g2)*Sqr(QHd) + 15*
      MassWB*Sqr(g2)*Sqr(QHd) + 30*traceAdjYuTYu*Sqr(QHu) - 30*MassU*
      traceYuAdjYu*Sqr(QHu) + 3*MassB*Sqr(g1)*Sqr(QHu) + 3*MassU*Sqr(g1)*Sqr(
      QHu) + 15*MassU*Sqr(g2)*Sqr(QHu) + 15*MassWB*Sqr(g2)*Sqr(QHu) + 10*
      traceAdjYeTYe*Sqr(Ql) - 10*MassU*traceYeAdjYe*Sqr(Ql) + 9*MassB*Sqr(g1)*
      Sqr(Ql) + 9*MassU*Sqr(g1)*Sqr(Ql) + 45*MassU*Sqr(g2)*Sqr(Ql) + 45*MassWB*
      Sqr(g2)*Sqr(Ql) + 30*traceAdjYdTYd*Sqr(Qq) + 30*traceAdjYuTYu*Sqr(Qq) -
      30*MassU*traceYuAdjYu*Sqr(Qq) + 3*MassB*Sqr(g1)*Sqr(Qq) + 3*MassU*Sqr(g1)
      *Sqr(Qq) + 135*MassU*Sqr(g2)*Sqr(Qq) + 135*MassWB*Sqr(g2)*Sqr(Qq) + 240*
      MassG*Sqr(g3)*Sqr(Qq) + 240*MassU*Sqr(g3)*Sqr(Qq) - 30*MassU*traceYdAdjYd
      *(Sqr(Qd) + Sqr(QHd) + Sqr(Qq)) + 30*traceAdjYuTYu*Sqr(Qu) - 30*MassU*
      traceYuAdjYu*Sqr(Qu) + 24*MassB*Sqr(g1)*Sqr(Qu) + 24*MassU*Sqr(g1)*Sqr(Qu
      ) + 120*MassG*Sqr(g3)*Sqr(Qu) + 120*MassU*Sqr(g3)*Sqr(Qu) - 10*Conj(
      Lambdax)*(Sqr(QHd) + Sqr(QHu) + Sqr(Qs))*(MassU*Lambdax - TLambdax)));


   return beta_MassU;
}

/**
 * Calculates the three-loop beta function of MassU.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassU_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassU;

   beta_MassU = 0;


   return beta_MassU;
}

} // namespace flexiblesusy
