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

// File generated at Fri 20 Oct 2017 08:51:40

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_1_loop(const Susy_traces& susy_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qv = INPUT(Qv);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (oneOver16PiSqr*(-0.2*Yv*(-5*AbsSqr(Lambdax) + 3*Sqr(g1) + 5
      *(-3*traceYuAdjYu - traceYvAdjYv + 3*Sqr(g2) + 2*Sqr(gp)*(Sqr(QHu) + Sqr(
      Ql) + Sqr(Qv)))) + 3*(Yv*Yv.adjoint()*Yv) + Ye.transpose()*Ye.conjugate()
      *Yv)).real();


   return beta_Yv;
}

/**
 * Calculates the 2-loop beta function of Yv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_2_loop(const Susy_traces& susy_traces) const
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

   beta_Yv = (twoLoop*(0.02*Yv*(207*Quad(g1) + 10*Sqr(g1)*(9*Sqr(g2) + 2*
      (2*traceYuAdjYu + 3*Sqr(gp)*(3*Qd*QHu + 3*Qe*QHu - QHd*QHu - 3*Qd*Ql - 3*
      Qe*Ql + QHd*Ql - 4*QHu*Ql + 3*QHu*Qq - 3*Ql*Qq - 6*QHu*Qu + 6*Ql*Qu + 2*
      Sqr(QHu) + 4*Sqr(Ql)))) + 50*AbsSqr(Lambdax)*(-3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs))) + 25*(15*Quad(
      g2) + 12*Sqr(g2)*Sqr(gp)*(Sqr(QHu) + Sqr(Ql)) + 2*(-3*traceYdAdjYuYuAdjYd
      - 9*traceYuAdjYuYuAdjYu - traceYvAdjYvTpYeconjYe - 3*traceYvAdjYvYvAdjYv
      + 16*traceYuAdjYu*Sqr(g3) + 2*Sqr(gp)*(-((3*traceYuAdjYu + traceYvAdjYv)
      *Sqr(QHu)) + traceYvAdjYv*Sqr(Ql) + 3*traceYuAdjYu*Sqr(Qq) + 3*
      traceYuAdjYu*Sqr(Qu) + traceYvAdjYv*Sqr(Qv)) + 2*Quad(gp)*(4*Quad(QHu) +
      8*Quad(Ql) + 5*Quad(Qv) + 2*Sqr(QHd)*Sqr(QHu) + 2*Sqr(QHd)*Sqr(Ql) + 8*
      Sqr(QHu)*Sqr(Ql) + 18*Sqr(QHu)*Sqr(Qq) + 18*Sqr(Ql)*Sqr(Qq) + Sqr(QHu)*
      Sqr(Qs) + Sqr(Ql)*Sqr(Qs) + 9*Sqr(QHu)*Sqr(Qu) + 9*Sqr(Ql)*Sqr(Qu) + 2*
      Sqr(QHd)*Sqr(Qv) + 5*Sqr(QHu)*Sqr(Qv) + 9*Sqr(Ql)*Sqr(Qv) + 18*Sqr(Qq)*
      Sqr(Qv) + Sqr(Qs)*Sqr(Qv) + 9*Sqr(Qu)*Sqr(Qv) + 9*Sqr(Qd)*(Sqr(QHu) + Sqr
      (Ql) + Sqr(Qv)) + 3*Sqr(Qe)*(Sqr(QHu) + Sqr(Ql) + Sqr(Qv))))) - 150*Sqr(
      Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYuAdjYu - 3*traceYvAdjYv - 3*
      AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*Sqr(g2) + 6*Sqr(gp)*Sqr(QHu) + 2*Sqr(gp
      )*Sqr(Ql) - 2*Sqr(gp)*Sqr(Qv))*(Yv*Yv.adjoint()*Yv) + (-3*traceYdAdjYd -
      traceYeAdjYe - AbsSqr(Lambdax) + 1.2*Sqr(g1) + 2*Sqr(gp)*Sqr(Qe) + 2*Sqr(
      gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(Ql))*(Ye.transpose()*Ye.conjugate()*Yv) - 4*
      (Yv*Yv.adjoint()*Yv*Yv.adjoint()*Yv) - 2*(Yv*Yv.adjoint()*Ye.transpose()*
      Ye.conjugate()*Yv) - 2*(Ye.transpose()*Ye.conjugate()*Ye.transpose()*
      Ye.conjugate()*Yv))).real();


   return beta_Yv;
}

/**
 * Calculates the 3-loop beta function of Yv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return beta_Yv;
}

} // namespace flexiblesusy
