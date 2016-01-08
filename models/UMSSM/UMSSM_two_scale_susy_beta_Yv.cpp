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

// File generated at Fri 8 Jan 2016 15:11:08

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yv.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_one_loop(const Susy_traces& susy_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qv = INPUT(Qv);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (oneOver16PiSqr*(Yv*(3*traceYuAdjYu + traceYvAdjYv + AbsSqr(
      Lambdax) - 1.8*Sqr(g1) - 3*Sqr(g2) - 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(
      Ql) - 2*Sqr(gp)*Sqr(Qv)) + 3*(Yv*Yv.adjoint()*Yv) + Ye.transpose()*
      Ye.conjugate()*Yv)).real();


   return beta_Yv;
}

/**
 * Calculates the two-loop beta function of Yv.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (twoLoop*(Yv*(16.74*Power(g1,4) + 7.5*Power(g2,4) + 8*Power(
      gp,4)*Power(QHu,4) + 16*Power(gp,4)*Power(Ql,4) + 10*Power(gp,4)*Power(Qv
      ,4) - 3*traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu -
      traceYvAdjYvTpYeconjYe - 3*traceYvAdjYvYvAdjYv + 1.2*traceYvAdjYv*Sqr(g1)
      + 1.8*Sqr(g1)*Sqr(g2) + 3.6*Qd*QHu*Sqr(g1)*Sqr(gp) + 3.6*Qe*QHu*Sqr(g1)*
      Sqr(gp) - 1.2*QHd*QHu*Sqr(g1)*Sqr(gp) - 3.6*Qd*Ql*Sqr(g1)*Sqr(gp) - 3.6*
      Qe*Ql*Sqr(g1)*Sqr(gp) + 1.2*QHd*Ql*Sqr(g1)*Sqr(gp) - 4.8*QHu*Ql*Sqr(g1)*
      Sqr(gp) + 3.6*QHu*Qq*Sqr(g1)*Sqr(gp) - 3.6*Ql*Qq*Sqr(g1)*Sqr(gp) - 7.2*
      QHu*Qu*Sqr(g1)*Sqr(gp) + 7.2*Ql*Qu*Sqr(g1)*Sqr(gp) + 7.2*Qd*Qv*Sqr(g1)*
      Sqr(gp) + 7.2*Qe*Qv*Sqr(g1)*Sqr(gp) - 2.4*QHd*Qv*Sqr(g1)*Sqr(gp) + 6*QHu*
      Qv*Sqr(g1)*Sqr(gp) - 10.8*Ql*Qv*Sqr(g1)*Sqr(gp) + 7.2*Qq*Qv*Sqr(g1)*Sqr(
      gp) - 14.4*Qu*Qv*Sqr(g1)*Sqr(gp) - 2*traceYvAdjYv*Sqr(gp)*Sqr(QHu) + 2.4*
      Sqr(g1)*Sqr(gp)*Sqr(QHu) + 6*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 18*Power(gp,4)*
      Sqr(Qd)*Sqr(QHu) + 6*Power(gp,4)*Sqr(Qe)*Sqr(QHu) + 4*Power(gp,4)*Sqr(QHd
      )*Sqr(QHu) + 2*traceYvAdjYv*Sqr(gp)*Sqr(Ql) + 4.8*Sqr(g1)*Sqr(gp)*Sqr(Ql)
      + 6*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 18*Power(gp,4)*Sqr(Qd)*Sqr(Ql) + 6*Power(
      gp,4)*Sqr(Qe)*Sqr(Ql) + 4*Power(gp,4)*Sqr(QHd)*Sqr(Ql) + 16*Power(gp,4)*
      Sqr(QHu)*Sqr(Ql) + 36*Power(gp,4)*Sqr(QHu)*Sqr(Qq) + 36*Power(gp,4)*Sqr(
      Ql)*Sqr(Qq) + 2*Power(gp,4)*Sqr(QHu)*Sqr(Qs) + 2*Power(gp,4)*Sqr(Ql)*Sqr(
      Qs) + AbsSqr(Lambdax)*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(
      QHd) - Sqr(QHu) + Sqr(Qs))) + 18*Power(gp,4)*Sqr(QHu)*Sqr(Qu) + 18*Power(
      gp,4)*Sqr(Ql)*Sqr(Qu) + 0.4*traceYuAdjYu*(2*Sqr(g1) + 5*(8*Sqr(g3) + 3*
      Sqr(gp)*(-Sqr(QHu) + Sqr(Qq) + Sqr(Qu)))) + 2*traceYvAdjYv*Sqr(gp)*Sqr(Qv
      ) + 12*Sqr(g1)*Sqr(gp)*Sqr(Qv) + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qv) + 6*Power
      (gp,4)*Sqr(Qe)*Sqr(Qv) + 4*Power(gp,4)*Sqr(QHd)*Sqr(Qv) + 10*Power(gp,4)*
      Sqr(QHu)*Sqr(Qv) + 18*Power(gp,4)*Sqr(Ql)*Sqr(Qv) + 36*Power(gp,4)*Sqr(Qq
      )*Sqr(Qv) + 2*Power(gp,4)*Sqr(Qs)*Sqr(Qv) + 18*Power(gp,4)*Sqr(Qu)*Sqr(Qv
      ) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYuAdjYu - 3*
      traceYvAdjYv - 3*AbsSqr(Lambdax) + 6*Sqr(g2) + 6*Sqr(gp)*Sqr(QHu) + 2*Sqr
      (gp)*Sqr(Ql) - 2*Sqr(gp)*Sqr(Qv))*(Yv*Yv.adjoint()*Yv) - 3*traceYdAdjYd*(
      Ye.transpose()*Ye.conjugate()*Yv) - traceYeAdjYe*(Ye.transpose()*
      Ye.conjugate()*Yv) - AbsSqr(Lambdax)*(Ye.transpose()*Ye.conjugate()*Yv) +
      1.2*Sqr(g1)*(Ye.transpose()*Ye.conjugate()*Yv) + 2*Sqr(gp)*Sqr(Qe)*(
      Ye.transpose()*Ye.conjugate()*Yv) + 2*Sqr(gp)*Sqr(QHd)*(Ye.transpose()*
      Ye.conjugate()*Yv) - 2*Sqr(gp)*Sqr(Ql)*(Ye.transpose()*Ye.conjugate()*Yv)
      - 4*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*Yv) - 2*(Yv*Yv.adjoint()*
      Ye.transpose()*Ye.conjugate()*Yv) - 2*(Ye.transpose()*Ye.conjugate()*
      Ye.transpose()*Ye.conjugate()*Yv))).real();


   return beta_Yv;
}

/**
 * Calculates the three-loop beta function of Yv.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return beta_Yv;
}

} // namespace flexiblesusy
