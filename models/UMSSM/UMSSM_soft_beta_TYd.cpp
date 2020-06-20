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


#include "UMSSM_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the 1-loop beta function of TYd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qq = INPUT(Qq);
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (0.06666666666666667*(90*traceAdjYdTYd*Yd + 30*traceAdjYeTYe*Yd +
      14*MassB*Yd*Sqr(g1) + 90*MassWB*Yd*Sqr(g2) + 160*MassG*Yd*Sqr(g3) + 60*
      MassU*Yd*Sqr(gp)*Sqr(Qd) + 60*MassU*Yd*Sqr(gp)*Sqr(QHd) + 60*MassU*Yd*Sqr
      (gp)*Sqr(Qq) + 45*traceYdAdjYd*TYd + 15*traceYeAdjYe*TYd + 15*AbsSqr(
      Lambdax)*TYd - 7*Sqr(g1)*TYd - 45*Sqr(g2)*TYd - 80*Sqr(g3)*TYd - 30*Sqr(
      gp)*Sqr(Qd)*TYd - 30*Sqr(gp)*Sqr(QHd)*TYd - 30*Sqr(gp)*Sqr(Qq)*TYd + 30*
      Yd*Conj(Lambdax)*TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()
      *TYu) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu).real();


   return oneLoop * beta_TYd;
}

/**
 * Calculates the 2-loop beta function of TYd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_2_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceAdjYeTYeconjYvTpYv = TRACE_STRUCT.traceAdjYeTYeconjYvTpYv;
   const double traceAdjYvTpYeconjYeTYv = TRACE_STRUCT.traceAdjYvTpYeconjYeTYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   const Eigen::Matrix<double,3,3> beta_TYd_1 = (-0.044444444444444446*Yd*(45*
      traceAdjYeTYeconjYvTpYv + 45*traceAdjYvTpYeconjYeTYv + 810*
      traceYdAdjYdTYdAdjYd + 135*traceYdAdjYuTYuAdjYd + 270*
      traceYeAdjYeTYeAdjYe + 287*MassB*Quad(g1) + 675*MassWB*Quad(g2) - 160*
      MassG*Quad(g3) + 1980*MassU*Quad(gp)*Quad(Qd) + 720*MassU*Quad(gp)*Quad(
      QHd) + 3600*MassU*Quad(gp)*Quad(Qq) + 18*traceAdjYdTYd*Sqr(g1) - 54*
      traceAdjYeTYe*Sqr(g1) - 18*MassB*traceYdAdjYd*Sqr(g1) + 54*MassB*
      traceYeAdjYe*Sqr(g1) + 45*MassB*Sqr(g1)*Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(
      g2) - 720*traceAdjYdTYd*Sqr(g3) + 720*MassG*traceYdAdjYd*Sqr(g3) + 40*
      MassB*Sqr(g1)*Sqr(g3) + 40*MassG*Sqr(g1)*Sqr(g3) + 360*MassG*Sqr(g2)*Sqr(
      g3) + 360*MassWB*Sqr(g2)*Sqr(g3) + 108*MassB*Qd*Qe*Sqr(g1)*Sqr(gp) + 108*
      MassU*Qd*Qe*Sqr(g1)*Sqr(gp) - 198*MassB*Qd*QHd*Sqr(g1)*Sqr(gp) - 198*
      MassU*Qd*QHd*Sqr(g1)*Sqr(gp) - 162*MassB*Qe*QHd*Sqr(g1)*Sqr(gp) - 162*
      MassU*Qe*QHd*Sqr(g1)*Sqr(gp) + 36*MassB*Qd*QHu*Sqr(g1)*Sqr(gp) + 36*MassU
      *Qd*QHu*Sqr(g1)*Sqr(gp) - 54*MassB*QHd*QHu*Sqr(g1)*Sqr(gp) - 54*MassU*QHd
      *QHu*Sqr(g1)*Sqr(gp) - 108*MassB*Qd*Ql*Sqr(g1)*Sqr(gp) - 108*MassU*Qd*Ql*
      Sqr(g1)*Sqr(gp) + 162*MassB*QHd*Ql*Sqr(g1)*Sqr(gp) + 162*MassU*QHd*Ql*Sqr
      (g1)*Sqr(gp) + 162*MassB*Qd*Qq*Sqr(g1)*Sqr(gp) + 162*MassU*Qd*Qq*Sqr(g1)*
      Sqr(gp) + 54*MassB*Qe*Qq*Sqr(g1)*Sqr(gp) + 54*MassU*Qe*Qq*Sqr(g1)*Sqr(gp)
      - 180*MassB*QHd*Qq*Sqr(g1)*Sqr(gp) - 180*MassU*QHd*Qq*Sqr(g1)*Sqr(gp) +
      18*MassB*QHu*Qq*Sqr(g1)*Sqr(gp) + 18*MassU*QHu*Qq*Sqr(g1)*Sqr(gp) - 54*
      MassB*Ql*Qq*Sqr(g1)*Sqr(gp) - 54*MassU*Ql*Qq*Sqr(g1)*Sqr(gp) - 216*MassB*
      Qd*Qu*Sqr(g1)*Sqr(gp) - 216*MassU*Qd*Qu*Sqr(g1)*Sqr(gp) + 324*MassB*QHd*
      Qu*Sqr(g1)*Sqr(gp) + 324*MassU*QHd*Qu*Sqr(g1)*Sqr(gp) - 108*MassB*Qq*Qu*
      Sqr(g1)*Sqr(gp) - 108*MassU*Qq*Qu*Sqr(g1)*Sqr(gp) - 270*traceAdjYdTYd*Sqr
      (gp)*Sqr(Qd) + 270*MassU*traceYdAdjYd*Sqr(gp)*Sqr(Qd) + 132*MassB*Sqr(g1)
      *Sqr(gp)*Sqr(Qd) + 132*MassU*Sqr(g1)*Sqr(gp)*Sqr(Qd) + 480*MassG*Sqr(g3)*
      Sqr(gp)*Sqr(Qd) + 480*MassU*Sqr(g3)*Sqr(gp)*Sqr(Qd) - 90*traceAdjYeTYe*
      Sqr(gp)*Sqr(Qe) + 90*MassU*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 540*MassU*Quad(
      gp)*Sqr(Qd)*Sqr(Qe) + 270*traceAdjYdTYd*Sqr(gp)*Sqr(QHd) + 90*
      traceAdjYeTYe*Sqr(gp)*Sqr(QHd) - 270*MassU*traceYdAdjYd*Sqr(gp)*Sqr(QHd)
      - 90*MassU*traceYeAdjYe*Sqr(gp)*Sqr(QHd) + 108*MassB*Sqr(g1)*Sqr(gp)*Sqr(
      QHd) + 108*MassU*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 270*MassU*Sqr(g2)*Sqr(gp)*Sqr
      (QHd) + 270*MassWB*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 1980*MassU*Quad(gp)*Sqr(Qd)
      *Sqr(QHd) + 540*MassU*Quad(gp)*Sqr(Qe)*Sqr(QHd) + 360*MassU*Quad(gp)*Sqr(
      Qd)*Sqr(QHu) + 360*MassU*Quad(gp)*Sqr(QHd)*Sqr(QHu) - 90*traceAdjYeTYe*
      Sqr(gp)*Sqr(Ql) + 90*MassU*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 1080*MassU*Quad
      (gp)*Sqr(Qd)*Sqr(Ql) + 1080*MassU*Quad(gp)*Sqr(QHd)*Sqr(Ql) - 270*
      traceAdjYdTYd*Sqr(gp)*Sqr(Qq) + 270*MassU*traceYdAdjYd*Sqr(gp)*Sqr(Qq) +
      60*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 60*MassU*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 270
      *MassU*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 270*MassWB*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 480
      *MassG*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 480*MassU*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 4860
      *MassU*Quad(gp)*Sqr(Qd)*Sqr(Qq) + 540*MassU*Quad(gp)*Sqr(Qe)*Sqr(Qq) +
      3600*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qq) + 360*MassU*Quad(gp)*Sqr(QHu)*Sqr(Qq
      ) + 1080*MassU*Quad(gp)*Sqr(Ql)*Sqr(Qq) + 180*MassU*Quad(gp)*Sqr(Qd)*Sqr(
      Qs) + 180*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qs) + 180*MassU*Quad(gp)*Sqr(Qq)*
      Sqr(Qs) + 1620*MassU*Quad(gp)*Sqr(Qd)*Sqr(Qu) + 1620*MassU*Quad(gp)*Sqr(
      QHd)*Sqr(Qu) + 1620*MassU*Quad(gp)*Sqr(Qq)*Sqr(Qu) + 540*MassU*Quad(gp)*
      Sqr(Qd)*Sqr(Qv) + 540*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qv) + 540*MassU*Quad(gp
      )*Sqr(Qq)*Sqr(Qv))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_2 = ((0.011111111111111112*(-540*
      traceYuAdjYdTYdAdjYu*Yd - 540*traceAdjYuTYu*Yd*AbsSqr(Lambdax) - 180*
      traceAdjYvTYv*Yd*AbsSqr(Lambdax) + 360*MassU*Yd*AbsSqr(Lambdax)*Sqr(gp)*
      Sqr(QHd) - 360*MassU*Yd*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 360*MassU*Yd*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 287*Quad(g1)*TYd + 675*Quad(g2)*TYd -
      160*Quad(g3)*TYd + 1980*Quad(gp)*Quad(Qd)*TYd + 720*Quad(gp)*Quad(QHd)*
      TYd + 3600*Quad(gp)*Quad(Qq)*TYd + 90*Sqr(g1)*Sqr(g2)*TYd + 80*Sqr(g1)*
      Sqr(g3)*TYd + 720*Sqr(g2)*Sqr(g3)*TYd + 216*Qd*Qe*Sqr(g1)*Sqr(gp)*TYd -
      396*Qd*QHd*Sqr(g1)*Sqr(gp)*TYd - 324*Qe*QHd*Sqr(g1)*Sqr(gp)*TYd + 72*Qd*
      QHu*Sqr(g1)*Sqr(gp)*TYd - 108*QHd*QHu*Sqr(g1)*Sqr(gp)*TYd - 216*Qd*Ql*Sqr
      (g1)*Sqr(gp)*TYd + 324*QHd*Ql*Sqr(g1)*Sqr(gp)*TYd + 324*Qd*Qq*Sqr(g1)*Sqr
      (gp)*TYd + 108*Qe*Qq*Sqr(g1)*Sqr(gp)*TYd - 360*QHd*Qq*Sqr(g1)*Sqr(gp)*TYd
       + 36*QHu*Qq*Sqr(g1)*Sqr(gp)*TYd - 108*Ql*Qq*Sqr(g1)*Sqr(gp)*TYd - 432*Qd
      *Qu*Sqr(g1)*Sqr(gp)*TYd + 264*Sqr(g1)*Sqr(gp)*Sqr(Qd)*TYd + 960*Sqr(g3)*
      Sqr(gp)*Sqr(Qd)*TYd + 540*Quad(gp)*Sqr(Qd)*Sqr(Qe)*TYd + 216*Sqr(g1)*Sqr(
      gp)*Sqr(QHd)*TYd + 540*Sqr(g2)*Sqr(gp)*Sqr(QHd)*TYd + 1980*Quad(gp)*Sqr(
      Qd)*Sqr(QHd)*TYd + 540*Quad(gp)*Sqr(Qe)*Sqr(QHd)*TYd + 360*Quad(gp)*Sqr(
      Qd)*Sqr(QHu)*TYd + 360*Quad(gp)*Sqr(QHd)*Sqr(QHu)*TYd + 1080*Quad(gp)*Sqr
      (Qd)*Sqr(Ql)*TYd + 1080*Quad(gp)*Sqr(QHd)*Sqr(Ql)*TYd + 120*Sqr(g1)*Sqr(
      gp)*Sqr(Qq)*TYd + 540*Sqr(g2)*Sqr(gp)*Sqr(Qq)*TYd + 960*Sqr(g3)*Sqr(gp)*
      Sqr(Qq)*TYd + 4860*Quad(gp)*Sqr(Qd)*Sqr(Qq)*TYd + 540*Quad(gp)*Sqr(Qe)*
      Sqr(Qq)*TYd + 3600*Quad(gp)*Sqr(QHd)*Sqr(Qq)*TYd + 360*Quad(gp)*Sqr(QHu)*
      Sqr(Qq)*TYd + 1080*Quad(gp)*Sqr(Ql)*Sqr(Qq)*TYd + 180*Quad(gp)*Sqr(Qd)*
      Sqr(Qs)*TYd + 180*Quad(gp)*Sqr(QHd)*Sqr(Qs)*TYd + 180*Quad(gp)*Sqr(Qq)*
      Sqr(Qs)*TYd) - 0.4*(45*traceAdjYdTYd + 15*traceAdjYeTYe + 4*MassB*Sqr(g1)
      + 30*MassWB*Sqr(g2) - 10*MassU*Sqr(gp)*Sqr(Qd) + 30*MassU*Sqr(gp)*Sqr(QHd
      ) + 10*MassU*Sqr(gp)*Sqr(Qq))*(Yd*Yd.adjoint()*Yd) + 0.4*(-30*
      traceYdAdjYd - 10*traceYeAdjYe - 10*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(
      g2) + 20*Sqr(gp)*Sqr(QHd))*(Yd*Yd.adjoint()*TYd) - 0.4*(15*traceAdjYuTYu
      + 5*traceAdjYvTYv + 4*MassB*Sqr(g1) + 10*MassU*Sqr(gp)*Sqr(QHu) - 10*
      MassU*Sqr(gp)*Sqr(Qq) + 10*MassU*Sqr(gp)*Sqr(Qu))*(Yd*Yu.adjoint()*Yu) +
      0.4*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 4*Sqr(g1) +
      10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq) + 10*Sqr(gp)*Sqr(Qu))*(Yd*Yu.
      adjoint()*TYu) + 0.2*(-75*traceYdAdjYd - 25*traceYeAdjYe - 25*AbsSqr(
      Lambdax) + 6*Sqr(g1) + 60*Sqr(g2) - 30*Sqr(gp)*Sqr(Qd) + 50*Sqr(gp)*Sqr(
      QHd) + 30*Sqr(gp)*Sqr(Qq))*(TYd*Yd.adjoint()*Yd) + 0.2*(-15*traceYuAdjYu
      - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + 10*Sqr(gp)*Sqr(QHu) -
      10*Sqr(gp)*Sqr(Qq) + 10*Sqr(gp)*Sqr(Qu))*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd
      .adjoint()*Yd*Yd.adjoint()*TYd) - 8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd)
      - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.
      adjoint()*TYu) - 4*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.
      adjoint()*TYu*Yu.adjoint()*Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd)
      - 4*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.
      adjoint()*Yu))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_3 = ((0.2*(-45*traceYdAdjYdYdAdjYd*
      TYd - 15*traceYdAdjYuYuAdjYd*TYd - 15*traceYeAdjYeYeAdjYe*TYd - 5*
      traceYvAdjYvTpYeconjYe*TYd - 15*traceYuAdjYu*AbsSqr(Lambdax)*TYd - 5*
      traceYvAdjYv*AbsSqr(Lambdax)*TYd - 2*traceYdAdjYd*Sqr(g1)*TYd + 6*
      traceYeAdjYe*Sqr(g1)*TYd + 80*traceYdAdjYd*Sqr(g3)*TYd + 36*QHd*Qu*Sqr(g1
      )*Sqr(gp)*TYd - 12*Qq*Qu*Sqr(g1)*Sqr(gp)*TYd + 30*traceYdAdjYd*Sqr(gp)*
      Sqr(Qd)*TYd + 10*traceYeAdjYe*Sqr(gp)*Sqr(Qe)*TYd - 30*traceYdAdjYd*Sqr(
      gp)*Sqr(QHd)*TYd - 10*traceYeAdjYe*Sqr(gp)*Sqr(QHd)*TYd - 10*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(QHd)*TYd + 10*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu)*TYd +
      10*traceYeAdjYe*Sqr(gp)*Sqr(Ql)*TYd + 30*traceYdAdjYd*Sqr(gp)*Sqr(Qq)*TYd
       + 10*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)*TYd + 90*Quad(gp)*Sqr(Qd)*Sqr(Qu)*
      TYd + 90*Quad(gp)*Sqr(QHd)*Sqr(Qu)*TYd + 90*Quad(gp)*Sqr(Qq)*Sqr(Qu)*TYd
      + 30*Quad(gp)*Sqr(Qd)*Sqr(Qv)*TYd + 30*Quad(gp)*Sqr(QHd)*Sqr(Qv)*TYd + 30
      *Quad(gp)*Sqr(Qq)*Sqr(Qv)*TYd - 15*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYd -
      30*traceYuAdjYu*Yd*Conj(Lambdax)*TLambdax - 10*traceYvAdjYv*Yd*Conj(
      Lambdax)*TLambdax - 20*Yd*Conj(Lambdax)*Sqr(gp)*Sqr(QHd)*TLambdax + 20*Yd
      *Conj(Lambdax)*Sqr(gp)*Sqr(QHu)*TLambdax + 20*Yd*Conj(Lambdax)*Sqr(gp)*
      Sqr(Qs)*TLambdax - 60*Yd*Lambdax*Sqr(Conj(Lambdax))*TLambdax) - 6*Conj(
      Lambdax)*TLambdax*(Yd*Yd.adjoint()*Yd) - 2*Conj(Lambdax)*TLambdax*(Yd*Yu.
      adjoint()*Yu))*UNITMATRIX(3)).real();

   beta_TYd = beta_TYd_1 + beta_TYd_2 + beta_TYd_3;


   return twoLoop * beta_TYd;
}

/**
 * Calculates the 3-loop beta function of TYd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return threeLoop * beta_TYd;
}

/**
 * Calculates the 4-loop beta function of TYd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return fourLoop * beta_TYd;
}

/**
 * Calculates the 5-loop beta function of TYd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return fiveLoop * beta_TYd;
}

} // namespace flexiblesusy
