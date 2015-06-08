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

// File generated at Mon 8 Jun 2015 17:48:58

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qq = INPUT(Qq);
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*(3*traceYdAdjYd*TYd + traceYeAdjYe*TYd +
      AbsSqr(Lambdax)*TYd - 0.4666666666666667*Sqr(g1)*TYd - 3*Sqr(g2)*TYd -
      5.333333333333333*Sqr(g3)*TYd - 2*Sqr(gp)*Sqr(Qd)*TYd - 2*Sqr(gp)*Sqr(QHd
      )*TYd - 2*Sqr(gp)*Sqr(Qq)*TYd + Yd*(6*traceAdjYdTYd + 2*traceAdjYeTYe +
      0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3) + 4*MassU*Sqr(gp)*Sqr(Qd) + 4*MassU*Sqr(gp)*Sqr(QHd) + 4*
      MassU*Sqr(gp)*Sqr(Qq) + 2*Conj(Lambdax)*TLambdax) + 4*(Yd*Yd.adjoint()*
      TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint
      ()*Yu)).real();


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(3.188888888888889*Power(g1,4)*TYd + 7.5*Power(g2,
      4)*TYd - 1.7777777777777777*Power(g3,4)*TYd + 22*Power(gp,4)*Power(Qd,4)*
      TYd + 8*Power(gp,4)*Power(QHd,4)*TYd + 40*Power(gp,4)*Power(Qq,4)*TYd - 9
      *traceYdAdjYdYdAdjYd*TYd - 3*traceYdAdjYuYuAdjYd*TYd - 3*
      traceYeAdjYeYeAdjYe*TYd - 3*traceYuAdjYu*AbsSqr(Lambdax)*TYd - 0.4*
      traceYdAdjYd*Sqr(g1)*TYd + 1.2*traceYeAdjYe*Sqr(g1)*TYd + Sqr(g1)*Sqr(g2)
      *TYd + 16*traceYdAdjYd*Sqr(g3)*TYd + 0.8888888888888888*Sqr(g1)*Sqr(g3)*
      TYd + 8*Sqr(g2)*Sqr(g3)*TYd + 6*traceYdAdjYd*Sqr(gp)*Sqr(Qd)*TYd +
      0.5333333333333333*Sqr(g1)*Sqr(gp)*Sqr(Qd)*TYd + 10.666666666666666*Sqr(
      g3)*Sqr(gp)*Sqr(Qd)*TYd + 2*traceYeAdjYe*Sqr(gp)*Sqr(Qe)*TYd + 6*Power(gp
      ,4)*Sqr(Qd)*Sqr(Qe)*TYd - 6*traceYdAdjYd*Sqr(gp)*Sqr(QHd)*TYd - 2*
      traceYeAdjYe*Sqr(gp)*Sqr(QHd)*TYd - 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd)*
      TYd + 1.2*Sqr(g1)*Sqr(gp)*Sqr(QHd)*TYd + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd)*TYd +
      22*Power(gp,4)*Sqr(Qd)*Sqr(QHd)*TYd + 6*Power(gp,4)*Sqr(Qe)*Sqr(QHd)*TYd
      + 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu)*TYd + 4*Power(gp,4)*Sqr(Qd)*Sqr(QHu
      )*TYd + 4*Power(gp,4)*Sqr(QHd)*Sqr(QHu)*TYd + 2*traceYeAdjYe*Sqr(gp)*Sqr(
      Ql)*TYd + 12*Power(gp,4)*Sqr(Qd)*Sqr(Ql)*TYd + 12*Power(gp,4)*Sqr(QHd)*
      Sqr(Ql)*TYd + 6*traceYdAdjYd*Sqr(gp)*Sqr(Qq)*TYd + 0.13333333333333333*
      Sqr(g1)*Sqr(gp)*Sqr(Qq)*TYd + 6*Sqr(g2)*Sqr(gp)*Sqr(Qq)*TYd +
      10.666666666666666*Sqr(g3)*Sqr(gp)*Sqr(Qq)*TYd + 54*Power(gp,4)*Sqr(Qd)*
      Sqr(Qq)*TYd + 6*Power(gp,4)*Sqr(Qe)*Sqr(Qq)*TYd + 40*Power(gp,4)*Sqr(QHd)
      *Sqr(Qq)*TYd + 4*Power(gp,4)*Sqr(QHu)*Sqr(Qq)*TYd + 12*Power(gp,4)*Sqr(Ql
      )*Sqr(Qq)*TYd + 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)*TYd + 2*Power(gp,4)*Sqr
      (Qd)*Sqr(Qs)*TYd + 2*Power(gp,4)*Sqr(QHd)*Sqr(Qs)*TYd + 2*Power(gp,4)*Sqr
      (Qq)*Sqr(Qs)*TYd + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qu)*TYd + 18*Power(gp,4)*
      Sqr(QHd)*Sqr(Qu)*TYd + 18*Power(gp,4)*Sqr(Qq)*Sqr(Qu)*TYd - 3*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)*TYd - 0.044444444444444446*Yd*(287*Power(g1,4)*
      MassB - 160*Power(g3,4)*MassG + 675*Power(g2,4)*MassWB + 1980*Power(gp,4)
      *MassU*Power(Qd,4) + 720*Power(gp,4)*MassU*Power(QHd,4) + 3600*Power(gp,4
      )*MassU*Power(Qq,4) + 810*traceYdAdjYdTYdAdjYd + 135*traceYdAdjYuTYuAdjYd
      + 270*traceYeAdjYeTYeAdjYe + 135*traceYuAdjYdTYdAdjYu + 18*traceAdjYdTYd
      *Sqr(g1) - 54*traceAdjYeTYe*Sqr(g1) + 54*MassB*traceYeAdjYe*Sqr(g1) + 45*
      MassB*Sqr(g1)*Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(g2) - 720*traceAdjYdTYd*Sqr
      (g3) + 40*MassB*Sqr(g1)*Sqr(g3) + 40*MassG*Sqr(g1)*Sqr(g3) + 360*MassG*
      Sqr(g2)*Sqr(g3) + 360*MassWB*Sqr(g2)*Sqr(g3) - 270*traceAdjYdTYd*Sqr(gp)*
      Sqr(Qd) + 24*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qd) + 24*MassU*Sqr(g1)*Sqr(gp)*Sqr
      (Qd) + 480*MassG*Sqr(g3)*Sqr(gp)*Sqr(Qd) + 480*MassU*Sqr(g3)*Sqr(gp)*Sqr(
      Qd) - 90*traceAdjYeTYe*Sqr(gp)*Sqr(Qe) + 90*MassU*traceYeAdjYe*Sqr(gp)*
      Sqr(Qe) + 540*Power(gp,4)*MassU*Sqr(Qd)*Sqr(Qe) + 270*traceAdjYdTYd*Sqr(
      gp)*Sqr(QHd) + 90*traceAdjYeTYe*Sqr(gp)*Sqr(QHd) - 90*MassU*traceYeAdjYe*
      Sqr(gp)*Sqr(QHd) + 54*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 54*MassU*Sqr(g1)*
      Sqr(gp)*Sqr(QHd) + 270*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 270*MassWB*Sqr(g2
      )*Sqr(gp)*Sqr(QHd) + 1980*Power(gp,4)*MassU*Sqr(Qd)*Sqr(QHd) + 540*Power(
      gp,4)*MassU*Sqr(Qe)*Sqr(QHd) + 360*Power(gp,4)*MassU*Sqr(Qd)*Sqr(QHu) +
      360*Power(gp,4)*MassU*Sqr(QHd)*Sqr(QHu) - 90*traceAdjYeTYe*Sqr(gp)*Sqr(Ql
      ) + 90*MassU*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 1080*Power(gp,4)*MassU*Sqr(Qd
      )*Sqr(Ql) + 1080*Power(gp,4)*MassU*Sqr(QHd)*Sqr(Ql) - 270*traceAdjYdTYd*
      Sqr(gp)*Sqr(Qq) + 6*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 6*MassU*Sqr(g1)*Sqr(
      gp)*Sqr(Qq) + 270*MassU*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 270*MassWB*Sqr(g2)*Sqr(
      gp)*Sqr(Qq) + 480*MassG*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 480*MassU*Sqr(g3)*Sqr(
      gp)*Sqr(Qq) + 4860*Power(gp,4)*MassU*Sqr(Qd)*Sqr(Qq) + 540*Power(gp,4)*
      MassU*Sqr(Qe)*Sqr(Qq) + 3600*Power(gp,4)*MassU*Sqr(QHd)*Sqr(Qq) + 360*
      Power(gp,4)*MassU*Sqr(QHu)*Sqr(Qq) + 1080*Power(gp,4)*MassU*Sqr(Ql)*Sqr(
      Qq) - 18*traceYdAdjYd*(MassB*Sqr(g1) - 5*(8*MassG*Sqr(g3) + 3*MassU*Sqr(
      gp)*(Sqr(Qd) - Sqr(QHd) + Sqr(Qq)))) + 180*Power(gp,4)*MassU*Sqr(Qd)*Sqr(
      Qs) + 180*Power(gp,4)*MassU*Sqr(QHd)*Sqr(Qs) + 180*Power(gp,4)*MassU*Sqr(
      Qq)*Sqr(Qs) + 1620*Power(gp,4)*MassU*Sqr(Qd)*Sqr(Qu) + 1620*Power(gp,4)*
      MassU*Sqr(QHd)*Sqr(Qu) + 1620*Power(gp,4)*MassU*Sqr(Qq)*Sqr(Qu) + 270*
      Lambdax*Sqr(Conj(Lambdax))*TLambdax - 45*Conj(Lambdax)*(Lambdax*(-3*
      traceAdjYuTYu + 2*MassU*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) - Sqr(Qs))) + (-3*
      traceYuAdjYu + 2*Sqr(gp)*(-Sqr(QHd) + Sqr(QHu) + Sqr(Qs)))*TLambdax)) -
      0.4*(45*traceAdjYdTYd + 15*traceAdjYeTYe + 4*MassB*Sqr(g1) + 30*MassWB*
      Sqr(g2) - 10*MassU*Sqr(gp)*Sqr(Qd) + 30*MassU*Sqr(gp)*Sqr(QHd) + 10*MassU
      *Sqr(gp)*Sqr(Qq) + 15*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) - 12*
      traceYdAdjYd*(Yd*Yd.adjoint()*TYd) - 4*traceYeAdjYe*(Yd*Yd.adjoint()*TYd)
      - 4*AbsSqr(Lambdax)*(Yd*Yd.adjoint()*TYd) + 1.2*Sqr(g1)*(Yd*Yd.adjoint()
      *TYd) + 6*Sqr(g2)*(Yd*Yd.adjoint()*TYd) + 8*Sqr(gp)*Sqr(QHd)*(Yd*
      Yd.adjoint()*TYd) - 6*traceAdjYuTYu*(Yd*Yu.adjoint()*Yu) - 1.6*MassB*Sqr(
      g1)*(Yd*Yu.adjoint()*Yu) - 4*MassU*Sqr(gp)*Sqr(QHu)*(Yd*Yu.adjoint()*Yu)
      + 4*MassU*Sqr(gp)*Sqr(Qq)*(Yd*Yu.adjoint()*Yu) - 4*MassU*Sqr(gp)*Sqr(Qu)*
      (Yd*Yu.adjoint()*Yu) - 2*Conj(Lambdax)*TLambdax*(Yd*Yu.adjoint()*Yu) - 6*
      traceYuAdjYu*(Yd*Yu.adjoint()*TYu) - 2*AbsSqr(Lambdax)*(Yd*Yu.adjoint()*
      TYu) + 1.6*Sqr(g1)*(Yd*Yu.adjoint()*TYu) + 4*Sqr(gp)*Sqr(QHu)*(Yd*
      Yu.adjoint()*TYu) - 4*Sqr(gp)*Sqr(Qq)*(Yd*Yu.adjoint()*TYu) + 4*Sqr(gp)*
      Sqr(Qu)*(Yd*Yu.adjoint()*TYu) - 15*traceYdAdjYd*(TYd*Yd.adjoint()*Yd) - 5
      *traceYeAdjYe*(TYd*Yd.adjoint()*Yd) - 5*AbsSqr(Lambdax)*(TYd*Yd.adjoint()
      *Yd) + 1.2*Sqr(g1)*(TYd*Yd.adjoint()*Yd) + 12*Sqr(g2)*(TYd*Yd.adjoint()*
      Yd) - 6*Sqr(gp)*Sqr(Qd)*(TYd*Yd.adjoint()*Yd) + 10*Sqr(gp)*Sqr(QHd)*(TYd*
      Yd.adjoint()*Yd) + 6*Sqr(gp)*Sqr(Qq)*(TYd*Yd.adjoint()*Yd) - 3*
      traceYuAdjYu*(TYd*Yu.adjoint()*Yu) - AbsSqr(Lambdax)*(TYd*Yu.adjoint()*Yu
      ) + 0.8*Sqr(g1)*(TYd*Yu.adjoint()*Yu) + 2*Sqr(gp)*Sqr(QHu)*(TYd*
      Yu.adjoint()*Yu) - 2*Sqr(gp)*Sqr(Qq)*(TYd*Yu.adjoint()*Yu) + 2*Sqr(gp)*
      Sqr(Qu)*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) -
      8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*
      Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*
      Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYd;
}

} // namespace flexiblesusy
