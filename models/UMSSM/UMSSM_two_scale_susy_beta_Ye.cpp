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

// File generated at Wed 12 Apr 2017 12:25:08

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(-0.2*Ye*(-5*AbsSqr(Lambdax) + 9*Sqr(g1) + 5
      *(-3*traceYdAdjYd - traceYeAdjYe + 3*Sqr(g2) + 2*Sqr(gp)*(Sqr(Qe) + Sqr(
      QHd) + Sqr(Ql)))) + 3*(Ye*Ye.adjoint()*Ye) + Ye*Yv.conjugate()*
      Yv.transpose())).real();


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.1*Ye*(135*Power(g1,4) + 2*Sqr(g1)*(9*Sqr(g2) + 2
      *(-traceYdAdjYd + 3*traceYeAdjYe + 3*Sqr(gp)*(-(QHd*QHu) + 4*QHd*Ql - QHu
      *Ql + Qd*(6*Qe - 3*(QHd + Ql)) - 3*QHd*Qq - 3*Ql*Qq + Qe*(-5*QHd + 2*QHu
      - 9*Ql + 6*Qq - 12*Qu) + 6*QHd*Qu + 6*Ql*Qu + 10*Sqr(Qe) + 2*Sqr(QHd) + 4
      *Sqr(Ql)))) - 10*AbsSqr(Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 2*Sqr(
      gp)*(Sqr(QHd) - Sqr(QHu) - Sqr(Qs))) + 5*(15*Power(g2,4) + 12*Sqr(g2)*Sqr
      (gp)*(Sqr(QHd) + Sqr(Ql)) + 2*(-9*traceYdAdjYdYdAdjYd - 3*
      traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 16
      *traceYdAdjYd*Sqr(g3) + 2*Sqr(gp)*(3*traceYdAdjYd*Sqr(Qd) + traceYeAdjYe*
      Sqr(Qe) - (3*traceYdAdjYd + traceYeAdjYe)*Sqr(QHd) + traceYeAdjYe*Sqr(Ql)
      + 3*traceYdAdjYd*Sqr(Qq)) + 2*Power(gp,4)*(5*Power(Qe,4) + 4*Power(QHd,4
      ) + 8*Power(Ql,4) + 2*Sqr(QHd)*Sqr(QHu) + 8*Sqr(QHd)*Sqr(Ql) + 2*Sqr(QHu)
      *Sqr(Ql) + 9*Sqr(Qd)*(Sqr(Qe) + Sqr(QHd) + Sqr(Ql)) + 18*Sqr(QHd)*Sqr(Qq)
      + 18*Sqr(Ql)*Sqr(Qq) + Sqr(QHd)*Sqr(Qs) + Sqr(Ql)*Sqr(Qs) + 9*Sqr(QHd)*
      Sqr(Qu) + 9*Sqr(Ql)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv) + 3*Sqr(Ql)*Sqr(Qv) +
      Sqr(Qe)*(5*Sqr(QHd) + 2*Sqr(QHu) + 9*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*
      Sqr(Qu) + 3*Sqr(Qv))))) - 30*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-9*
      traceYdAdjYd - 3*traceYeAdjYe - 3*AbsSqr(Lambdax) + 6*Sqr(g2) - 2*Sqr(gp)
      *Sqr(Qe) + 6*Sqr(gp)*Sqr(QHd) + 2*Sqr(gp)*Sqr(Ql))*(Ye*Ye.adjoint()*Ye) +
      (-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax) + 2*Sqr(gp)*(Sqr(QHu)
      - Sqr(Ql) + Sqr(Qv)))*(Ye*Yv.conjugate()*Yv.transpose()) - 4*(Ye*
      Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 2*(Ye*Yv.conjugate()*Yv.transpose()*
      Ye.adjoint()*Ye) - 2*(Ye*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*
      Yv.transpose()))).real();


   return beta_Ye;
}

/**
 * Calculates the three-loop beta function of Ye.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
