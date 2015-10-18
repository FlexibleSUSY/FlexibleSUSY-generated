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

// File generated at Sun 18 Oct 2015 12:18:37

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(3*traceYdAdjYd*TYe + traceYeAdjYe*TYe +
      AbsSqr(Lambdax)*TYe - 1.8*Sqr(g1)*TYe - 3*Sqr(g2)*TYe - 2*Sqr(gp)*Sqr(Qe)
      *TYe - 2*Sqr(gp)*Sqr(QHd)*TYe - 2*Sqr(gp)*Sqr(Ql)*TYe + Ye*(6*
      traceAdjYdTYd + 2*traceAdjYeTYe + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      4*MassU*Sqr(gp)*Sqr(Qe) + 4*MassU*Sqr(gp)*Sqr(QHd) + 4*MassU*Sqr(gp)*Sqr(
      Ql) + 2*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.adjoint()*TYe) + 5*(TYe*
      Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(13.5*Power(g1,4)*TYe + 7.5*Power(g2,4)*TYe + 10*
      Power(gp,4)*Power(Qe,4)*TYe + 8*Power(gp,4)*Power(QHd,4)*TYe + 16*Power(
      gp,4)*Power(Ql,4)*TYe - 9*traceYdAdjYdYdAdjYd*TYe - 3*traceYdAdjYuYuAdjYd
      *TYe - 3*traceYeAdjYeYeAdjYe*TYe - 3*traceYuAdjYu*AbsSqr(Lambdax)*TYe -
      0.4*traceYdAdjYd*Sqr(g1)*TYe + 1.2*traceYeAdjYe*Sqr(g1)*TYe + 1.8*Sqr(g1)
      *Sqr(g2)*TYe + 16*traceYdAdjYd*Sqr(g3)*TYe + 7.2*Qd*Qe*Sqr(g1)*Sqr(gp)*
      TYe - 3.6*Qd*QHd*Sqr(g1)*Sqr(gp)*TYe - 6*Qe*QHd*Sqr(g1)*Sqr(gp)*TYe + 2.4
      *Qe*QHu*Sqr(g1)*Sqr(gp)*TYe - 1.2*QHd*QHu*Sqr(g1)*Sqr(gp)*TYe - 3.6*Qd*Ql
      *Sqr(g1)*Sqr(gp)*TYe - 10.8*Qe*Ql*Sqr(g1)*Sqr(gp)*TYe + 4.8*QHd*Ql*Sqr(g1
      )*Sqr(gp)*TYe - 1.2*QHu*Ql*Sqr(g1)*Sqr(gp)*TYe + 7.2*Qe*Qq*Sqr(g1)*Sqr(gp
      )*TYe - 3.6*QHd*Qq*Sqr(g1)*Sqr(gp)*TYe - 3.6*Ql*Qq*Sqr(g1)*Sqr(gp)*TYe -
      14.4*Qe*Qu*Sqr(g1)*Sqr(gp)*TYe + 7.2*QHd*Qu*Sqr(g1)*Sqr(gp)*TYe + 7.2*Ql*
      Qu*Sqr(g1)*Sqr(gp)*TYe + 6*traceYdAdjYd*Sqr(gp)*Sqr(Qd)*TYe + 2*
      traceYeAdjYe*Sqr(gp)*Sqr(Qe)*TYe + 12*Sqr(g1)*Sqr(gp)*Sqr(Qe)*TYe + 18*
      Power(gp,4)*Sqr(Qd)*Sqr(Qe)*TYe - 6*traceYdAdjYd*Sqr(gp)*Sqr(QHd)*TYe - 2
      *traceYeAdjYe*Sqr(gp)*Sqr(QHd)*TYe - 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd)*
      TYe + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHd)*TYe + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd)*TYe +
      18*Power(gp,4)*Sqr(Qd)*Sqr(QHd)*TYe + 10*Power(gp,4)*Sqr(Qe)*Sqr(QHd)*
      TYe + 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu)*TYe + 4*Power(gp,4)*Sqr(Qe)*Sqr(
      QHu)*TYe + 4*Power(gp,4)*Sqr(QHd)*Sqr(QHu)*TYe + 2*traceYeAdjYe*Sqr(gp)*
      Sqr(Ql)*TYe + 4.8*Sqr(g1)*Sqr(gp)*Sqr(Ql)*TYe + 6*Sqr(g2)*Sqr(gp)*Sqr(Ql)
      *TYe + 18*Power(gp,4)*Sqr(Qd)*Sqr(Ql)*TYe + 18*Power(gp,4)*Sqr(Qe)*Sqr(Ql
      )*TYe + 16*Power(gp,4)*Sqr(QHd)*Sqr(Ql)*TYe + 4*Power(gp,4)*Sqr(QHu)*Sqr(
      Ql)*TYe + 6*traceYdAdjYd*Sqr(gp)*Sqr(Qq)*TYe + 36*Power(gp,4)*Sqr(Qe)*Sqr
      (Qq)*TYe + 36*Power(gp,4)*Sqr(QHd)*Sqr(Qq)*TYe + 36*Power(gp,4)*Sqr(Ql)*
      Sqr(Qq)*TYe + 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)*TYe + 2*Power(gp,4)*Sqr(
      Qe)*Sqr(Qs)*TYe + 2*Power(gp,4)*Sqr(QHd)*Sqr(Qs)*TYe + 2*Power(gp,4)*Sqr(
      Ql)*Sqr(Qs)*TYe + 18*Power(gp,4)*Sqr(Qe)*Sqr(Qu)*TYe + 18*Power(gp,4)*Sqr
      (QHd)*Sqr(Qu)*TYe + 18*Power(gp,4)*Sqr(Ql)*Sqr(Qu)*TYe - 3*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)*TYe - 0.4*Ye*(135*Power(g1,4)*MassB + 75*Power(g2,
      4)*MassWB + 100*Power(gp,4)*MassU*Power(Qe,4) + 80*Power(gp,4)*MassU*
      Power(QHd,4) + 160*Power(gp,4)*MassU*Power(Ql,4) + 90*
      traceYdAdjYdTYdAdjYd + 15*traceYdAdjYuTYuAdjYd + 30*traceYeAdjYeTYeAdjYe
      + 15*traceYuAdjYdTYdAdjYu + 2*traceAdjYdTYd*Sqr(g1) - 6*traceAdjYeTYe*Sqr
      (g1) + 6*MassB*traceYeAdjYe*Sqr(g1) + 9*MassB*Sqr(g1)*Sqr(g2) + 9*MassWB*
      Sqr(g1)*Sqr(g2) - 80*traceAdjYdTYd*Sqr(g3) + 36*MassB*Qd*Qe*Sqr(g1)*Sqr(
      gp) + 36*MassU*Qd*Qe*Sqr(g1)*Sqr(gp) - 18*MassB*Qd*QHd*Sqr(g1)*Sqr(gp) -
      18*MassU*Qd*QHd*Sqr(g1)*Sqr(gp) - 30*MassB*Qe*QHd*Sqr(g1)*Sqr(gp) - 30*
      MassU*Qe*QHd*Sqr(g1)*Sqr(gp) + 12*MassB*Qe*QHu*Sqr(g1)*Sqr(gp) + 12*MassU
      *Qe*QHu*Sqr(g1)*Sqr(gp) - 6*MassB*QHd*QHu*Sqr(g1)*Sqr(gp) - 6*MassU*QHd*
      QHu*Sqr(g1)*Sqr(gp) - 18*MassB*Qd*Ql*Sqr(g1)*Sqr(gp) - 18*MassU*Qd*Ql*Sqr
      (g1)*Sqr(gp) - 54*MassB*Qe*Ql*Sqr(g1)*Sqr(gp) - 54*MassU*Qe*Ql*Sqr(g1)*
      Sqr(gp) + 24*MassB*QHd*Ql*Sqr(g1)*Sqr(gp) + 24*MassU*QHd*Ql*Sqr(g1)*Sqr(
      gp) - 6*MassB*QHu*Ql*Sqr(g1)*Sqr(gp) - 6*MassU*QHu*Ql*Sqr(g1)*Sqr(gp) +
      36*MassB*Qe*Qq*Sqr(g1)*Sqr(gp) + 36*MassU*Qe*Qq*Sqr(g1)*Sqr(gp) - 18*
      MassB*QHd*Qq*Sqr(g1)*Sqr(gp) - 18*MassU*QHd*Qq*Sqr(g1)*Sqr(gp) - 18*MassB
      *Ql*Qq*Sqr(g1)*Sqr(gp) - 18*MassU*Ql*Qq*Sqr(g1)*Sqr(gp) - 72*MassB*Qe*Qu*
      Sqr(g1)*Sqr(gp) - 72*MassU*Qe*Qu*Sqr(g1)*Sqr(gp) + 36*MassB*QHd*Qu*Sqr(g1
      )*Sqr(gp) + 36*MassU*QHd*Qu*Sqr(g1)*Sqr(gp) + 36*MassB*Ql*Qu*Sqr(g1)*Sqr(
      gp) + 36*MassU*Ql*Qu*Sqr(g1)*Sqr(gp) - 30*traceAdjYdTYd*Sqr(gp)*Sqr(Qd) -
      10*traceAdjYeTYe*Sqr(gp)*Sqr(Qe) + 10*MassU*traceYeAdjYe*Sqr(gp)*Sqr(Qe)
      + 60*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qe) + 60*MassU*Sqr(g1)*Sqr(gp)*Sqr(Qe) +
      180*Power(gp,4)*MassU*Sqr(Qd)*Sqr(Qe) + 30*traceAdjYdTYd*Sqr(gp)*Sqr(QHd)
      + 10*traceAdjYeTYe*Sqr(gp)*Sqr(QHd) - 10*MassU*traceYeAdjYe*Sqr(gp)*Sqr(
      QHd) + 12*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 12*MassU*Sqr(g1)*Sqr(gp)*Sqr(
      QHd) + 30*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 30*MassWB*Sqr(g2)*Sqr(gp)*Sqr(
      QHd) + 180*Power(gp,4)*MassU*Sqr(Qd)*Sqr(QHd) + 100*Power(gp,4)*MassU*Sqr
      (Qe)*Sqr(QHd) + 40*Power(gp,4)*MassU*Sqr(Qe)*Sqr(QHu) + 40*Power(gp,4)*
      MassU*Sqr(QHd)*Sqr(QHu) - 10*traceAdjYeTYe*Sqr(gp)*Sqr(Ql) + 10*MassU*
      traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 24*MassB*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 24*
      MassU*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 30*MassU*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 30*
      MassWB*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 180*Power(gp,4)*MassU*Sqr(Qd)*Sqr(Ql) +
      180*Power(gp,4)*MassU*Sqr(Qe)*Sqr(Ql) + 160*Power(gp,4)*MassU*Sqr(QHd)*
      Sqr(Ql) + 40*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Ql) - 30*traceAdjYdTYd*Sqr(gp
      )*Sqr(Qq) + 360*Power(gp,4)*MassU*Sqr(Qe)*Sqr(Qq) + 360*Power(gp,4)*MassU
      *Sqr(QHd)*Sqr(Qq) + 360*Power(gp,4)*MassU*Sqr(Ql)*Sqr(Qq) + traceYdAdjYd*
      (-2*MassB*Sqr(g1) + 10*(8*MassG*Sqr(g3) + 3*MassU*Sqr(gp)*(Sqr(Qd) - Sqr(
      QHd) + Sqr(Qq)))) + 20*Power(gp,4)*MassU*Sqr(Qe)*Sqr(Qs) + 20*Power(gp,4)
      *MassU*Sqr(QHd)*Sqr(Qs) + 20*Power(gp,4)*MassU*Sqr(Ql)*Sqr(Qs) + 180*
      Power(gp,4)*MassU*Sqr(Qe)*Sqr(Qu) + 180*Power(gp,4)*MassU*Sqr(QHd)*Sqr(Qu
      ) + 180*Power(gp,4)*MassU*Sqr(Ql)*Sqr(Qu) + 30*Lambdax*Sqr(Conj(Lambdax))
      *TLambdax - 5*Conj(Lambdax)*(Lambdax*(-3*traceAdjYuTYu + 2*MassU*Sqr(gp)*
      (Sqr(QHd) - Sqr(QHu) - Sqr(Qs))) + (-3*traceYuAdjYu + 2*Sqr(gp)*(-Sqr(QHd
      ) + Sqr(QHu) + Sqr(Qs)))*TLambdax)) - 2*(9*traceAdjYdTYd + 3*
      traceAdjYeTYe + 6*MassWB*Sqr(g2) - 2*MassU*Sqr(gp)*Sqr(Qe) + 6*MassU*Sqr(
      gp)*Sqr(QHd) + 2*MassU*Sqr(gp)*Sqr(Ql) + 3*Conj(Lambdax)*TLambdax)*(Ye*
      Ye.adjoint()*Ye) - 12*traceYdAdjYd*(Ye*Ye.adjoint()*TYe) - 4*traceYeAdjYe
      *(Ye*Ye.adjoint()*TYe) - 4*AbsSqr(Lambdax)*(Ye*Ye.adjoint()*TYe) + 1.2*
      Sqr(g1)*(Ye*Ye.adjoint()*TYe) + 6*Sqr(g2)*(Ye*Ye.adjoint()*TYe) + 8*Sqr(
      gp)*Sqr(QHd)*(Ye*Ye.adjoint()*TYe) - 15*traceYdAdjYd*(TYe*Ye.adjoint()*Ye
      ) - 5*traceYeAdjYe*(TYe*Ye.adjoint()*Ye) - 5*AbsSqr(Lambdax)*(TYe*
      Ye.adjoint()*Ye) - 1.2*Sqr(g1)*(TYe*Ye.adjoint()*Ye) + 12*Sqr(g2)*(TYe*
      Ye.adjoint()*Ye) - 6*Sqr(gp)*Sqr(Qe)*(TYe*Ye.adjoint()*Ye) + 10*Sqr(gp)*
      Sqr(QHd)*(TYe*Ye.adjoint()*Ye) + 6*Sqr(gp)*Sqr(Ql)*(TYe*Ye.adjoint()*Ye)
      - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*
      Ye.adjoint()*Ye) - 6*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the three-loop beta function of TYe.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
