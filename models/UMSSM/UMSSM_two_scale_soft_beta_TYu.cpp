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

// File generated at Fri 26 Jun 2015 19:03:34

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of TYu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*(3*traceYuAdjYu*TYu + AbsSqr(Lambdax)*TYu -
      0.8666666666666667*Sqr(g1)*TYu - 3*Sqr(g2)*TYu - 5.333333333333333*Sqr(
      g3)*TYu - 2*Sqr(gp)*Sqr(QHu)*TYu - 2*Sqr(gp)*Sqr(Qq)*TYu - 2*Sqr(gp)*Sqr(
      Qu)*TYu + Yu*(6*traceAdjYuTYu + 1.7333333333333334*MassB*Sqr(g1) + 6*
      MassWB*Sqr(g2) + 10.666666666666666*MassG*Sqr(g3) + 4*MassU*Sqr(gp)*Sqr(
      QHu) + 4*MassU*Sqr(gp)*Sqr(Qq) + 4*MassU*Sqr(gp)*Sqr(Qu) + 2*Conj(Lambdax
      )*TLambdax) + 2*(Yu*Yd.adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) + TYu*
      Yd.adjoint()*Yd + 5*(TYu*Yu.adjoint()*Yu))).real();


   return beta_TYu;
}

/**
 * Calculates the two-loop beta function of TYu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (twoLoop*(6.095555555555555*Power(g1,4)*TYu + 7.5*Power(g2,
      4)*TYu - 1.7777777777777777*Power(g3,4)*TYu + 8*Power(gp,4)*Power(QHu,4)*
      TYu + 40*Power(gp,4)*Power(Qq,4)*TYu + 22*Power(gp,4)*Power(Qu,4)*TYu - 3
      *traceYdAdjYuYuAdjYd*TYu - 9*traceYuAdjYuYuAdjYu*TYu - 3*traceYdAdjYd*
      AbsSqr(Lambdax)*TYu - traceYeAdjYe*AbsSqr(Lambdax)*TYu + 0.8*traceYuAdjYu
      *Sqr(g1)*TYu + Sqr(g1)*Sqr(g2)*TYu + 16*traceYuAdjYu*Sqr(g3)*TYu +
      3.022222222222222*Sqr(g1)*Sqr(g3)*TYu + 8*Sqr(g2)*Sqr(g3)*TYu + 3.6*Qd*
      QHu*Sqr(g1)*Sqr(gp)*TYu + 3.6*Qe*QHu*Sqr(g1)*Sqr(gp)*TYu - 1.2*QHd*QHu*
      Sqr(g1)*Sqr(gp)*TYu - 3.6*QHu*Ql*Sqr(g1)*Sqr(gp)*TYu + 1.2*Qd*Qq*Sqr(g1)*
      Sqr(gp)*TYu + 1.2*Qe*Qq*Sqr(g1)*Sqr(gp)*TYu - 0.4*QHd*Qq*Sqr(g1)*Sqr(gp)*
      TYu + 4*QHu*Qq*Sqr(g1)*Sqr(gp)*TYu - 1.2*Ql*Qq*Sqr(g1)*Sqr(gp)*TYu - 4.8*
      Qd*Qu*Sqr(g1)*Sqr(gp)*TYu - 4.8*Qe*Qu*Sqr(g1)*Sqr(gp)*TYu + 1.6*QHd*Qu*
      Sqr(g1)*Sqr(gp)*TYu - 8.8*QHu*Qu*Sqr(g1)*Sqr(gp)*TYu + 4.8*Ql*Qu*Sqr(g1)*
      Sqr(gp)*TYu - 7.2*Qq*Qu*Sqr(g1)*Sqr(gp)*TYu + 2*AbsSqr(Lambdax)*Sqr(gp)*
      Sqr(QHd)*TYu - 6*traceYuAdjYu*Sqr(gp)*Sqr(QHu)*TYu - 2*AbsSqr(Lambdax)*
      Sqr(gp)*Sqr(QHu)*TYu + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHu)*TYu + 6*Sqr(g2)*Sqr(
      gp)*Sqr(QHu)*TYu + 18*Power(gp,4)*Sqr(Qd)*Sqr(QHu)*TYu + 6*Power(gp,4)*
      Sqr(Qe)*Sqr(QHu)*TYu + 4*Power(gp,4)*Sqr(QHd)*Sqr(QHu)*TYu + 12*Power(gp,
      4)*Sqr(QHu)*Sqr(Ql)*TYu + 6*traceYuAdjYu*Sqr(gp)*Sqr(Qq)*TYu +
      1.3333333333333333*Sqr(g1)*Sqr(gp)*Sqr(Qq)*TYu + 6*Sqr(g2)*Sqr(gp)*Sqr(Qq
      )*TYu + 10.666666666666666*Sqr(g3)*Sqr(gp)*Sqr(Qq)*TYu + 18*Power(gp,4)*
      Sqr(Qd)*Sqr(Qq)*TYu + 6*Power(gp,4)*Sqr(Qe)*Sqr(Qq)*TYu + 4*Power(gp,4)*
      Sqr(QHd)*Sqr(Qq)*TYu + 40*Power(gp,4)*Sqr(QHu)*Sqr(Qq)*TYu + 12*Power(gp,
      4)*Sqr(Ql)*Sqr(Qq)*TYu + 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)*TYu + 2*Power(
      gp,4)*Sqr(QHu)*Sqr(Qs)*TYu + 2*Power(gp,4)*Sqr(Qq)*Sqr(Qs)*TYu + 6*
      traceYuAdjYu*Sqr(gp)*Sqr(Qu)*TYu + 11.733333333333333*Sqr(g1)*Sqr(gp)*Sqr
      (Qu)*TYu + 10.666666666666666*Sqr(g3)*Sqr(gp)*Sqr(Qu)*TYu + 18*Power(gp,4
      )*Sqr(Qd)*Sqr(Qu)*TYu + 6*Power(gp,4)*Sqr(Qe)*Sqr(Qu)*TYu + 4*Power(gp,4)
      *Sqr(QHd)*Sqr(Qu)*TYu + 22*Power(gp,4)*Sqr(QHu)*Sqr(Qu)*TYu + 12*Power(gp
      ,4)*Sqr(Ql)*Sqr(Qu)*TYu + 54*Power(gp,4)*Sqr(Qq)*Sqr(Qu)*TYu + 2*Power(gp
      ,4)*Sqr(Qs)*Sqr(Qu)*TYu - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYu -
      0.008888888888888889*Yu*(2743*Power(g1,4)*MassB - 800*Power(g3,4)*MassG +
      3375*Power(g2,4)*MassWB + 3600*Power(gp,4)*MassU*Power(QHu,4) + 18000*
      Power(gp,4)*MassU*Power(Qq,4) + 9900*Power(gp,4)*MassU*Power(Qu,4) + 675*
      traceYdAdjYuTYuAdjYd + 675*traceYuAdjYdTYdAdjYu + 4050*
      traceYuAdjYuTYuAdjYu - 180*traceAdjYuTYu*Sqr(g1) + 225*MassB*Sqr(g1)*Sqr(
      g2) + 225*MassWB*Sqr(g1)*Sqr(g2) - 3600*traceAdjYuTYu*Sqr(g3) + 680*MassB
      *Sqr(g1)*Sqr(g3) + 680*MassG*Sqr(g1)*Sqr(g3) + 1800*MassG*Sqr(g2)*Sqr(g3)
      + 1800*MassWB*Sqr(g2)*Sqr(g3) + 810*MassB*Qd*QHu*Sqr(g1)*Sqr(gp) + 810*
      MassU*Qd*QHu*Sqr(g1)*Sqr(gp) + 810*MassB*Qe*QHu*Sqr(g1)*Sqr(gp) + 810*
      MassU*Qe*QHu*Sqr(g1)*Sqr(gp) - 270*MassB*QHd*QHu*Sqr(g1)*Sqr(gp) - 270*
      MassU*QHd*QHu*Sqr(g1)*Sqr(gp) - 810*MassB*QHu*Ql*Sqr(g1)*Sqr(gp) - 810*
      MassU*QHu*Ql*Sqr(g1)*Sqr(gp) + 270*MassB*Qd*Qq*Sqr(g1)*Sqr(gp) + 270*
      MassU*Qd*Qq*Sqr(g1)*Sqr(gp) + 270*MassB*Qe*Qq*Sqr(g1)*Sqr(gp) + 270*MassU
      *Qe*Qq*Sqr(g1)*Sqr(gp) - 90*MassB*QHd*Qq*Sqr(g1)*Sqr(gp) - 90*MassU*QHd*
      Qq*Sqr(g1)*Sqr(gp) + 900*MassB*QHu*Qq*Sqr(g1)*Sqr(gp) + 900*MassU*QHu*Qq*
      Sqr(g1)*Sqr(gp) - 270*MassB*Ql*Qq*Sqr(g1)*Sqr(gp) - 270*MassU*Ql*Qq*Sqr(
      g1)*Sqr(gp) - 1080*MassB*Qd*Qu*Sqr(g1)*Sqr(gp) - 1080*MassU*Qd*Qu*Sqr(g1)
      *Sqr(gp) - 1080*MassB*Qe*Qu*Sqr(g1)*Sqr(gp) - 1080*MassU*Qe*Qu*Sqr(g1)*
      Sqr(gp) + 360*MassB*QHd*Qu*Sqr(g1)*Sqr(gp) + 360*MassU*QHd*Qu*Sqr(g1)*Sqr
      (gp) - 1980*MassB*QHu*Qu*Sqr(g1)*Sqr(gp) - 1980*MassU*QHu*Qu*Sqr(g1)*Sqr(
      gp) + 1080*MassB*Ql*Qu*Sqr(g1)*Sqr(gp) + 1080*MassU*Ql*Qu*Sqr(g1)*Sqr(gp)
      - 1620*MassB*Qq*Qu*Sqr(g1)*Sqr(gp) - 1620*MassU*Qq*Qu*Sqr(g1)*Sqr(gp) +
      1350*traceAdjYuTYu*Sqr(gp)*Sqr(QHu) + 540*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHu)
      + 540*MassU*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 1350*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHu
      ) + 1350*MassWB*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 8100*Power(gp,4)*MassU*Sqr(Qd)
      *Sqr(QHu) + 2700*Power(gp,4)*MassU*Sqr(Qe)*Sqr(QHu) + 1800*Power(gp,4)*
      MassU*Sqr(QHd)*Sqr(QHu) + 5400*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Ql) - 1350*
      traceAdjYuTYu*Sqr(gp)*Sqr(Qq) + 300*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 300*
      MassU*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 1350*MassU*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 1350
      *MassWB*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 2400*MassG*Sqr(g3)*Sqr(gp)*Sqr(Qq) +
      2400*MassU*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 8100*Power(gp,4)*MassU*Sqr(Qd)*Sqr(
      Qq) + 2700*Power(gp,4)*MassU*Sqr(Qe)*Sqr(Qq) + 1800*Power(gp,4)*MassU*Sqr
      (QHd)*Sqr(Qq) + 18000*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Qq) + 5400*Power(gp,
      4)*MassU*Sqr(Ql)*Sqr(Qq) + 900*Power(gp,4)*MassU*Sqr(QHu)*Sqr(Qs) + 900*
      Power(gp,4)*MassU*Sqr(Qq)*Sqr(Qs) - 1350*traceAdjYuTYu*Sqr(gp)*Sqr(Qu) +
      2640*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qu) + 2640*MassU*Sqr(g1)*Sqr(gp)*Sqr(Qu) +
      2400*MassG*Sqr(g3)*Sqr(gp)*Sqr(Qu) + 2400*MassU*Sqr(g3)*Sqr(gp)*Sqr(Qu)
      + 8100*Power(gp,4)*MassU*Sqr(Qd)*Sqr(Qu) + 2700*Power(gp,4)*MassU*Sqr(Qe)
      *Sqr(Qu) + 1800*Power(gp,4)*MassU*Sqr(QHd)*Sqr(Qu) + 9900*Power(gp,4)*
      MassU*Sqr(QHu)*Sqr(Qu) + 5400*Power(gp,4)*MassU*Sqr(Ql)*Sqr(Qu) + 24300*
      Power(gp,4)*MassU*Sqr(Qq)*Sqr(Qu) + 900*Power(gp,4)*MassU*Sqr(Qs)*Sqr(Qu)
      + 90*traceYuAdjYu*(2*MassB*Sqr(g1) + 5*(8*MassG*Sqr(g3) + 3*MassU*Sqr(gp
      )*(-Sqr(QHu) + Sqr(Qq) + Sqr(Qu)))) + 1350*Lambdax*Sqr(Conj(Lambdax))*
      TLambdax + 225*Conj(Lambdax)*(Lambdax*(3*traceAdjYdTYd + traceAdjYeTYe +
      2*MassU*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs))) + (3*traceYdAdjYd +
      traceYeAdjYe - 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs)))*TLambdax)) + (
      -6*traceAdjYdTYd - 2*traceAdjYeTYe - 0.8*MassB*Sqr(g1) - 4*MassU*Sqr(gp)*
      Sqr(Qd) - 4*MassU*Sqr(gp)*Sqr(QHd) + 4*MassU*Sqr(gp)*Sqr(Qq) - 2*Conj(
      Lambdax)*TLambdax)*(Yu*Yd.adjoint()*Yd) - 6*traceYdAdjYd*(Yu*Yd.adjoint()
      *TYd) - 2*traceYeAdjYe*(Yu*Yd.adjoint()*TYd) - 2*AbsSqr(Lambdax)*(Yu*
      Yd.adjoint()*TYd) + 0.8*Sqr(g1)*(Yu*Yd.adjoint()*TYd) + 4*Sqr(gp)*Sqr(Qd)
      *(Yu*Yd.adjoint()*TYd) + 4*Sqr(gp)*Sqr(QHd)*(Yu*Yd.adjoint()*TYd) - 4*Sqr
      (gp)*Sqr(Qq)*(Yu*Yd.adjoint()*TYd) - 18*traceAdjYuTYu*(Yu*Yu.adjoint()*Yu
      ) - 0.8*MassB*Sqr(g1)*(Yu*Yu.adjoint()*Yu) - 12*MassWB*Sqr(g2)*(Yu*
      Yu.adjoint()*Yu) - 12*MassU*Sqr(gp)*Sqr(QHu)*(Yu*Yu.adjoint()*Yu) - 4*
      MassU*Sqr(gp)*Sqr(Qq)*(Yu*Yu.adjoint()*Yu) + 4*MassU*Sqr(gp)*Sqr(Qu)*(Yu*
      Yu.adjoint()*Yu) - 6*Conj(Lambdax)*TLambdax*(Yu*Yu.adjoint()*Yu) - 12*
      traceYuAdjYu*(Yu*Yu.adjoint()*TYu) - 4*AbsSqr(Lambdax)*(Yu*Yu.adjoint()*
      TYu) + 1.2*Sqr(g1)*(Yu*Yu.adjoint()*TYu) + 6*Sqr(g2)*(Yu*Yu.adjoint()*TYu
      ) + 8*Sqr(gp)*Sqr(QHu)*(Yu*Yu.adjoint()*TYu) - 3*traceYdAdjYd*(TYu*
      Yd.adjoint()*Yd) - traceYeAdjYe*(TYu*Yd.adjoint()*Yd) - AbsSqr(Lambdax)*(
      TYu*Yd.adjoint()*Yd) + 0.4*Sqr(g1)*(TYu*Yd.adjoint()*Yd) + 2*Sqr(gp)*Sqr(
      Qd)*(TYu*Yd.adjoint()*Yd) + 2*Sqr(gp)*Sqr(QHd)*(TYu*Yd.adjoint()*Yd) - 2*
      Sqr(gp)*Sqr(Qq)*(TYu*Yd.adjoint()*Yd) - 15*traceYuAdjYu*(TYu*Yu.adjoint()
      *Yu) - 5*AbsSqr(Lambdax)*(TYu*Yu.adjoint()*Yu) + 12*Sqr(g2)*(TYu*
      Yu.adjoint()*Yu) + 10*Sqr(gp)*Sqr(QHu)*(TYu*Yu.adjoint()*Yu) + 6*Sqr(gp)*
      Sqr(Qq)*(TYu*Yu.adjoint()*Yu) - 6*Sqr(gp)*Sqr(Qu)*(TYu*Yu.adjoint()*Yu) -
      4*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 2*(Yu*Yd.adjoint()*Yd*
      Yu.adjoint()*TYu) - 4*(Yu*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 4*(Yu*
      Yd.adjoint()*TYd*Yu.adjoint()*Yu) - 6*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*
      TYu) - 8*(Yu*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 2*(TYu*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) - 4*(TYu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 6*(TYu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYu;
}

/**
 * Calculates the three-loop beta function of TYu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

} // namespace flexiblesusy
