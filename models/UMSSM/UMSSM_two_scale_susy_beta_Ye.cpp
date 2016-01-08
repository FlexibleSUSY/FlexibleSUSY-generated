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

// File generated at Fri 8 Jan 2016 12:29:21

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

   beta_Ye = (oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 1.8*Sqr(g1) - 3*Sqr(g2) - 2*Sqr(gp)*Sqr(Qe) - 2*Sqr(gp)*Sqr(
      QHd) - 2*Sqr(gp)*Sqr(Ql)) + 3*(Ye*Ye.adjoint()*Ye) + Ye*Yv.conjugate()*
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

   beta_Ye = (twoLoop*(Ye*(16.74*Power(g1,4) + 7.5*Power(g2,4) + 10*Power
      (gp,4)*Power(Qe,4) + 8*Power(gp,4)*Power(QHd,4) + 16*Power(gp,4)*Power(Ql
      ,4) - 9*traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 1.2*traceYeAdjYe*Sqr(g1) +
      1.8*Sqr(g1)*Sqr(g2) + 7.2*Qd*Qe*Sqr(g1)*Sqr(gp) - 3.6*Qd*QHd*Sqr(g1)*Sqr
      (gp) - 6*Qe*QHd*Sqr(g1)*Sqr(gp) + 2.4*Qe*QHu*Sqr(g1)*Sqr(gp) - 1.2*QHd*
      QHu*Sqr(g1)*Sqr(gp) - 3.6*Qd*Ql*Sqr(g1)*Sqr(gp) - 10.8*Qe*Ql*Sqr(g1)*Sqr(
      gp) + 4.8*QHd*Ql*Sqr(g1)*Sqr(gp) - 1.2*QHu*Ql*Sqr(g1)*Sqr(gp) + 7.2*Qe*Qq
      *Sqr(g1)*Sqr(gp) - 3.6*QHd*Qq*Sqr(g1)*Sqr(gp) - 3.6*Ql*Qq*Sqr(g1)*Sqr(gp)
      - 14.4*Qe*Qu*Sqr(g1)*Sqr(gp) + 7.2*QHd*Qu*Sqr(g1)*Sqr(gp) + 7.2*Ql*Qu*
      Sqr(g1)*Sqr(gp) + 7.2*Qe*Qv*Sqr(g1)*Sqr(gp) - 3.6*QHd*Qv*Sqr(g1)*Sqr(gp)
      - 3.6*Ql*Qv*Sqr(g1)*Sqr(gp) + 2*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 12*Sqr(g1)
      *Sqr(gp)*Sqr(Qe) + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qe) - 2*traceYeAdjYe*Sqr(gp
      )*Sqr(QHd) + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd) +
      18*Power(gp,4)*Sqr(Qd)*Sqr(QHd) + 10*Power(gp,4)*Sqr(Qe)*Sqr(QHd) + 4*
      Power(gp,4)*Sqr(Qe)*Sqr(QHu) + 4*Power(gp,4)*Sqr(QHd)*Sqr(QHu) + 2*
      traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 4.8*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 6*Sqr(g2)*
      Sqr(gp)*Sqr(Ql) + 18*Power(gp,4)*Sqr(Qd)*Sqr(Ql) + 18*Power(gp,4)*Sqr(Qe)
      *Sqr(Ql) + 16*Power(gp,4)*Sqr(QHd)*Sqr(Ql) + 4*Power(gp,4)*Sqr(QHu)*Sqr(
      Ql) + 36*Power(gp,4)*Sqr(Qe)*Sqr(Qq) + 36*Power(gp,4)*Sqr(QHd)*Sqr(Qq) +
      36*Power(gp,4)*Sqr(Ql)*Sqr(Qq) - 0.4*traceYdAdjYd*(Sqr(g1) - 5*(8*Sqr(g3)
      + 3*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) + Sqr(Qq)))) - AbsSqr(Lambdax)*(3*
      traceYuAdjYu + traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) - Sqr(Qs)))
      + 2*Power(gp,4)*Sqr(Qe)*Sqr(Qs) + 2*Power(gp,4)*Sqr(QHd)*Sqr(Qs) + 2*
      Power(gp,4)*Sqr(Ql)*Sqr(Qs) + 18*Power(gp,4)*Sqr(Qe)*Sqr(Qu) + 18*Power(
      gp,4)*Sqr(QHd)*Sqr(Qu) + 18*Power(gp,4)*Sqr(Ql)*Sqr(Qu) + 6*Power(gp,4)*
      Sqr(Qe)*Sqr(Qv) + 6*Power(gp,4)*Sqr(QHd)*Sqr(Qv) + 6*Power(gp,4)*Sqr(Ql)*
      Sqr(Qv) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYdAdjYd - 3*
      traceYeAdjYe - 3*AbsSqr(Lambdax) + 6*Sqr(g2) - 2*Sqr(gp)*Sqr(Qe) + 6*Sqr(
      gp)*Sqr(QHd) + 2*Sqr(gp)*Sqr(Ql))*(Ye*Ye.adjoint()*Ye) - 3*traceYuAdjYu*(
      Ye*Yv.conjugate()*Yv.transpose()) - traceYvAdjYv*(Ye*Yv.conjugate()*
      Yv.transpose()) - AbsSqr(Lambdax)*(Ye*Yv.conjugate()*Yv.transpose()) +
      1.2*Sqr(g1)*(Ye*Yv.conjugate()*Yv.transpose()) + 2*Sqr(gp)*Sqr(QHu)*(Ye*
      Yv.conjugate()*Yv.transpose()) - 2*Sqr(gp)*Sqr(Ql)*(Ye*Yv.conjugate()*
      Yv.transpose()) + 2*Sqr(gp)*Sqr(Qv)*(Ye*Yv.conjugate()*Yv.transpose()) -
      4*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 2*(Ye*Yv.conjugate()*
      Yv.transpose()*Ye.adjoint()*Ye) - 2*(Ye*Yv.conjugate()*Yv.transpose()*
      Yv.conjugate()*Yv.transpose()))).real();


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
