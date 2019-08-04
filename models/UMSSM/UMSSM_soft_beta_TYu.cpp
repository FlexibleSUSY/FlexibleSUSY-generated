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

// File generated at Sun 4 Aug 2019 19:34:26

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
 * Calculates the 1-loop beta function of TYu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_1_loop(const Soft_traces& soft_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*(0.06666666666666667*(90*traceAdjYuTYu*Yu + 30*
      traceAdjYvTYv*Yu + 26*MassB*Yu*Sqr(g1) + 90*MassWB*Yu*Sqr(g2) + 160*MassG
      *Yu*Sqr(g3) + 60*MassU*Yu*Sqr(gp)*Sqr(QHu) + 60*MassU*Yu*Sqr(gp)*Sqr(Qq)
      + 60*MassU*Yu*Sqr(gp)*Sqr(Qu) + 45*traceYuAdjYu*TYu + 15*traceYvAdjYv*TYu
       + 15*AbsSqr(Lambdax)*TYu - 13*Sqr(g1)*TYu - 45*Sqr(g2)*TYu - 80*Sqr(g3)*
      TYu - 30*Sqr(gp)*Sqr(QHu)*TYu - 30*Sqr(gp)*Sqr(Qq)*TYu - 30*Sqr(gp)*Sqr(
      Qu)*TYu + 30*Yu*Conj(Lambdax)*TLambdax) + 2*(Yu*Yd.adjoint()*TYd) + 4*(Yu
      *Yu.adjoint()*TYu) + TYu*Yd.adjoint()*Yd + 5*(TYu*Yu.adjoint()*Yu))).real
      ();


   return beta_TYu;
}

/**
 * Calculates the 2-loop beta function of TYu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_2_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYvAdjYvTYvAdjYv = TRACE_STRUCT.traceYvAdjYvTYvAdjYv;
   const double traceAdjYeTYeconjYvTpYv = TRACE_STRUCT.traceAdjYeTYeconjYvTpYv;
   const double traceAdjYvTpYeconjYeTYv = TRACE_STRUCT.traceAdjYvTpYeconjYeTYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_TYu;

   const Eigen::Matrix<double,3,3> beta_TYu_1 = (-0.008888888888888889*twoLoop*
      Yu*(225*traceAdjYeTYeconjYvTpYv + 225*traceAdjYvTpYeconjYeTYv + 675*
      traceYdAdjYuTYuAdjYd + 675*traceYuAdjYdTYdAdjYu + 4050*
      traceYuAdjYuTYuAdjYu + 1350*traceYvAdjYvTYvAdjYv + 2743*MassB*Quad(g1) +
      3375*MassWB*Quad(g2) - 800*MassG*Quad(g3) + 3600*MassU*Quad(gp)*Quad(QHu)
      + 18000*MassU*Quad(gp)*Quad(Qq) + 9900*MassU*Quad(gp)*Quad(Qu) - 180*
      traceAdjYuTYu*Sqr(g1) + 180*MassB*traceYuAdjYu*Sqr(g1) + 225*MassB*Sqr(g1
      )*Sqr(g2) + 225*MassWB*Sqr(g1)*Sqr(g2) - 3600*traceAdjYuTYu*Sqr(g3) +
      3600*MassG*traceYuAdjYu*Sqr(g3) + 680*MassB*Sqr(g1)*Sqr(g3) + 680*MassG*
      Sqr(g1)*Sqr(g3) + 1800*MassG*Sqr(g2)*Sqr(g3) + 1800*MassWB*Sqr(g2)*Sqr(g3
      ) + 810*MassB*Qd*QHu*Sqr(g1)*Sqr(gp) + 810*MassU*Qd*QHu*Sqr(g1)*Sqr(gp) +
      810*MassB*Qe*QHu*Sqr(g1)*Sqr(gp) + 810*MassU*Qe*QHu*Sqr(g1)*Sqr(gp) - 270
      *MassB*QHd*QHu*Sqr(g1)*Sqr(gp) - 270*MassU*QHd*QHu*Sqr(g1)*Sqr(gp) - 810*
      MassB*QHu*Ql*Sqr(g1)*Sqr(gp) - 810*MassU*QHu*Ql*Sqr(g1)*Sqr(gp) + 270*
      MassB*Qd*Qq*Sqr(g1)*Sqr(gp) + 270*MassU*Qd*Qq*Sqr(g1)*Sqr(gp) + 270*MassB
      *Qe*Qq*Sqr(g1)*Sqr(gp) + 270*MassU*Qe*Qq*Sqr(g1)*Sqr(gp) - 90*MassB*QHd*
      Qq*Sqr(g1)*Sqr(gp) - 90*MassU*QHd*Qq*Sqr(g1)*Sqr(gp) + 900*MassB*QHu*Qq*
      Sqr(g1)*Sqr(gp) + 900*MassU*QHu*Qq*Sqr(g1)*Sqr(gp) - 270*MassB*Ql*Qq*Sqr(
      g1)*Sqr(gp) - 270*MassU*Ql*Qq*Sqr(g1)*Sqr(gp) - 1080*MassB*Qd*Qu*Sqr(g1)*
      Sqr(gp) - 1080*MassU*Qd*Qu*Sqr(g1)*Sqr(gp) - 1080*MassB*Qe*Qu*Sqr(g1)*Sqr
      (gp) - 1080*MassU*Qe*Qu*Sqr(g1)*Sqr(gp) + 360*MassB*QHd*Qu*Sqr(g1)*Sqr(gp
      ) + 360*MassU*QHd*Qu*Sqr(g1)*Sqr(gp) - 1980*MassB*QHu*Qu*Sqr(g1)*Sqr(gp)
      - 1980*MassU*QHu*Qu*Sqr(g1)*Sqr(gp) + 1080*MassB*Ql*Qu*Sqr(g1)*Sqr(gp) +
      1080*MassU*Ql*Qu*Sqr(g1)*Sqr(gp) - 1620*MassB*Qq*Qu*Sqr(g1)*Sqr(gp) -
      1620*MassU*Qq*Qu*Sqr(g1)*Sqr(gp) + 450*MassU*AbsSqr(Lambdax)*Sqr(gp)*Sqr(
      QHd) + 1350*traceAdjYuTYu*Sqr(gp)*Sqr(QHu) + 450*traceAdjYvTYv*Sqr(gp)*
      Sqr(QHu) - 1350*MassU*traceYuAdjYu*Sqr(gp)*Sqr(QHu) - 450*MassU*
      traceYvAdjYv*Sqr(gp)*Sqr(QHu) + 540*MassB*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 540*
      MassU*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 1350*MassU*Sqr(g2)*Sqr(gp)*Sqr(QHu) +
      1350*MassWB*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 8100*MassU*Quad(gp)*Sqr(Qd)*Sqr(
      QHu) + 2700*MassU*Quad(gp)*Sqr(Qe)*Sqr(QHu) + 1800*MassU*Quad(gp)*Sqr(QHd
      )*Sqr(QHu) - 450*traceAdjYvTYv*Sqr(gp)*Sqr(Ql) + 450*MassU*traceYvAdjYv*
      Sqr(gp)*Sqr(Ql) + 5400*MassU*Quad(gp)*Sqr(QHu)*Sqr(Ql) - 1350*
      traceAdjYuTYu*Sqr(gp)*Sqr(Qq) + 1350*MassU*traceYuAdjYu*Sqr(gp)*Sqr(Qq) +
      300*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 300*MassU*Sqr(g1)*Sqr(gp)*Sqr(Qq) +
      1350*MassU*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 1350*MassWB*Sqr(g2)*Sqr(gp)*Sqr(Qq)
      + 2400*MassG*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 2400*MassU*Sqr(g3)*Sqr(gp)*Sqr(Qq)
      + 8100*MassU*Quad(gp)*Sqr(Qd)*Sqr(Qq) + 2700*MassU*Quad(gp)*Sqr(Qe)*Sqr(
      Qq) + 1800*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qq) + 18000*MassU*Quad(gp)*Sqr(QHu
      )*Sqr(Qq) + 5400*MassU*Quad(gp)*Sqr(Ql)*Sqr(Qq) + 900*MassU*Quad(gp)*Sqr(
      QHu)*Sqr(Qs) + 900*MassU*Quad(gp)*Sqr(Qq)*Sqr(Qs) - 1350*traceAdjYuTYu*
      Sqr(gp)*Sqr(Qu) + 1350*MassU*traceYuAdjYu*Sqr(gp)*Sqr(Qu) + 2640*MassB*
      Sqr(g1)*Sqr(gp)*Sqr(Qu) + 2640*MassU*Sqr(g1)*Sqr(gp)*Sqr(Qu) + 2400*MassG
      *Sqr(g3)*Sqr(gp)*Sqr(Qu) + 2400*MassU*Sqr(g3)*Sqr(gp)*Sqr(Qu) + 8100*
      MassU*Quad(gp)*Sqr(Qd)*Sqr(Qu) + 2700*MassU*Quad(gp)*Sqr(Qe)*Sqr(Qu) +
      1800*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qu) + 9900*MassU*Quad(gp)*Sqr(QHu)*Sqr(
      Qu) + 5400*MassU*Quad(gp)*Sqr(Ql)*Sqr(Qu) + 24300*MassU*Quad(gp)*Sqr(Qq)*
      Sqr(Qu) + 900*MassU*Quad(gp)*Sqr(Qs)*Sqr(Qu) - 450*traceAdjYvTYv*Sqr(gp)*
      Sqr(Qv) + 450*MassU*traceYvAdjYv*Sqr(gp)*Sqr(Qv) + 2700*MassU*Quad(gp)*
      Sqr(QHu)*Sqr(Qv) + 2700*MassU*Quad(gp)*Sqr(Qq)*Sqr(Qv) + 2700*MassU*Quad(
      gp)*Sqr(Qu)*Sqr(Qv))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYu_2 = ((0.0022222222222222222*twoLoop
      *(-2700*traceAdjYdTYd*Yu*AbsSqr(Lambdax) - 900*traceAdjYeTYe*Yu*AbsSqr(
      Lambdax) + 1800*MassU*Yu*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 1800*MassU*Yu
      *AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 2743*Quad(g1)*TYu + 3375*Quad(g2)*TYu
      - 800*Quad(g3)*TYu + 3600*Quad(gp)*Quad(QHu)*TYu + 18000*Quad(gp)*Quad(Qq
      )*TYu + 450*Sqr(g1)*Sqr(g2)*TYu + 1360*Sqr(g1)*Sqr(g3)*TYu + 3600*Sqr(g2)
      *Sqr(g3)*TYu + 1620*Qd*QHu*Sqr(g1)*Sqr(gp)*TYu + 1620*Qe*QHu*Sqr(g1)*Sqr(
      gp)*TYu - 540*QHd*QHu*Sqr(g1)*Sqr(gp)*TYu - 1620*QHu*Ql*Sqr(g1)*Sqr(gp)*
      TYu + 540*Qd*Qq*Sqr(g1)*Sqr(gp)*TYu + 540*Qe*Qq*Sqr(g1)*Sqr(gp)*TYu - 180
      *QHd*Qq*Sqr(g1)*Sqr(gp)*TYu + 1800*QHu*Qq*Sqr(g1)*Sqr(gp)*TYu - 540*Ql*Qq
      *Sqr(g1)*Sqr(gp)*TYu - 2160*Qd*Qu*Sqr(g1)*Sqr(gp)*TYu - 2160*Qe*Qu*Sqr(g1
      )*Sqr(gp)*TYu + 720*QHd*Qu*Sqr(g1)*Sqr(gp)*TYu - 3960*QHu*Qu*Sqr(g1)*Sqr(
      gp)*TYu + 2160*Ql*Qu*Sqr(g1)*Sqr(gp)*TYu - 3240*Qq*Qu*Sqr(g1)*Sqr(gp)*TYu
       + 1080*Sqr(g1)*Sqr(gp)*Sqr(QHu)*TYu + 2700*Sqr(g2)*Sqr(gp)*Sqr(QHu)*TYu
      + 8100*Quad(gp)*Sqr(Qd)*Sqr(QHu)*TYu + 2700*Quad(gp)*Sqr(Qe)*Sqr(QHu)*TYu
       + 1800*Quad(gp)*Sqr(QHd)*Sqr(QHu)*TYu + 5400*Quad(gp)*Sqr(QHu)*Sqr(Ql)*
      TYu + 600*Sqr(g1)*Sqr(gp)*Sqr(Qq)*TYu + 2700*Sqr(g2)*Sqr(gp)*Sqr(Qq)*TYu
      + 4800*Sqr(g3)*Sqr(gp)*Sqr(Qq)*TYu + 8100*Quad(gp)*Sqr(Qd)*Sqr(Qq)*TYu +
      2700*Quad(gp)*Sqr(Qe)*Sqr(Qq)*TYu + 1800*Quad(gp)*Sqr(QHd)*Sqr(Qq)*TYu +
      18000*Quad(gp)*Sqr(QHu)*Sqr(Qq)*TYu + 5400*Quad(gp)*Sqr(Ql)*Sqr(Qq)*TYu +
      900*Quad(gp)*Sqr(QHu)*Sqr(Qs)*TYu + 900*Quad(gp)*Sqr(Qq)*Sqr(Qs)*TYu +
      5280*Sqr(g1)*Sqr(gp)*Sqr(Qu)*TYu + 4800*Sqr(g3)*Sqr(gp)*Sqr(Qu)*TYu +
      8100*Quad(gp)*Sqr(Qd)*Sqr(Qu)*TYu + 2700*Quad(gp)*Sqr(Qe)*Sqr(Qu)*TYu +
      1800*Quad(gp)*Sqr(QHd)*Sqr(Qu)*TYu + 9900*Quad(gp)*Sqr(QHu)*Sqr(Qu)*TYu +
      5400*Quad(gp)*Sqr(Ql)*Sqr(Qu)*TYu + 24300*Quad(gp)*Sqr(Qq)*Sqr(Qu)*TYu) -
      0.4*twoLoop*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 2*MassB*Sqr(g1) + 10*
      MassU*Sqr(gp)*Sqr(Qd) + 10*MassU*Sqr(gp)*Sqr(QHd) - 10*MassU*Sqr(gp)*Sqr(
      Qq))*(Yu*Yd.adjoint()*Yd) + 0.4*twoLoop*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1) + 10*Sqr(gp)*Sqr(Qd) + 10*
      Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Qq))*(Yu*Yd.adjoint()*TYd) - 0.4*
      twoLoop*(45*traceAdjYuTYu + 15*traceAdjYvTYv + 2*MassB*Sqr(g1) + 30*
      MassWB*Sqr(g2) + 30*MassU*Sqr(gp)*Sqr(QHu) + 10*MassU*Sqr(gp)*Sqr(Qq) -
      10*MassU*Sqr(gp)*Sqr(Qu))*(Yu*Yu.adjoint()*Yu) + 0.4*twoLoop*(-30*
      traceYuAdjYu - 10*traceYvAdjYv - 10*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(
      g2) + 20*Sqr(gp)*Sqr(QHu))*(Yu*Yu.adjoint()*TYu) + 0.2*twoLoop*(-15*
      traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1) + 10*Sqr(gp
      )*Sqr(Qd) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Qq))*(TYu*Yd.adjoint()*
      Yd) + twoLoop*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 12
      *Sqr(g2) + 10*Sqr(gp)*Sqr(QHu) + 6*Sqr(gp)*Sqr(Qq) - 6*Sqr(gp)*Sqr(Qu))*(
      TYu*Yu.adjoint()*Yu) - 4*twoLoop*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd) -
      2*twoLoop*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*twoLoop*(Yu*Yd.
      adjoint()*TYd*Yd.adjoint()*Yd) - 4*twoLoop*(Yu*Yd.adjoint()*TYd*Yu.
      adjoint()*Yu) - 6*twoLoop*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 8*
      twoLoop*(Yu*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 2*twoLoop*(TYu*Yd.adjoint
      ()*Yd*Yd.adjoint()*Yd) - 4*twoLoop*(TYu*Yd.adjoint()*Yd*Yu.adjoint()*Yu)
      - 6*twoLoop*(TYu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYu_3 = ((0.2*twoLoop*(-15*
      traceYdAdjYuYuAdjYd*TYu - 45*traceYuAdjYuYuAdjYu*TYu - 5*
      traceYvAdjYvTpYeconjYe*TYu - 15*traceYvAdjYvYvAdjYv*TYu - 15*traceYdAdjYd
      *AbsSqr(Lambdax)*TYu - 5*traceYeAdjYe*AbsSqr(Lambdax)*TYu + 110*Quad(gp)*
      Quad(Qu)*TYu + 4*traceYuAdjYu*Sqr(g1)*TYu + 80*traceYuAdjYu*Sqr(g3)*TYu +
      10*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd)*TYu - 30*traceYuAdjYu*Sqr(gp)*Sqr(QHu
      )*TYu - 10*traceYvAdjYv*Sqr(gp)*Sqr(QHu)*TYu - 10*AbsSqr(Lambdax)*Sqr(gp)
      *Sqr(QHu)*TYu + 10*traceYvAdjYv*Sqr(gp)*Sqr(Ql)*TYu + 30*traceYuAdjYu*Sqr
      (gp)*Sqr(Qq)*TYu + 10*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)*TYu + 30*
      traceYuAdjYu*Sqr(gp)*Sqr(Qu)*TYu + 10*Quad(gp)*Sqr(Qs)*Sqr(Qu)*TYu + 10*
      traceYvAdjYv*Sqr(gp)*Sqr(Qv)*TYu + 30*Quad(gp)*Sqr(QHu)*Sqr(Qv)*TYu + 30*
      Quad(gp)*Sqr(Qq)*Sqr(Qv)*TYu + 30*Quad(gp)*Sqr(Qu)*Sqr(Qv)*TYu - 15*Sqr(
      Conj(Lambdax))*Sqr(Lambdax)*TYu - 30*traceYdAdjYd*Yu*Conj(Lambdax)*
      TLambdax - 10*traceYeAdjYe*Yu*Conj(Lambdax)*TLambdax + 20*Yu*Conj(Lambdax
      )*Sqr(gp)*Sqr(QHd)*TLambdax - 20*Yu*Conj(Lambdax)*Sqr(gp)*Sqr(QHu)*
      TLambdax + 20*Yu*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*TLambdax - 60*Yu*Lambdax*
      Sqr(Conj(Lambdax))*TLambdax) - 2*twoLoop*Conj(Lambdax)*TLambdax*(Yu*Yd.
      adjoint()*Yd) - 6*twoLoop*Conj(Lambdax)*TLambdax*(Yu*Yu.adjoint()*Yu))*
      UNITMATRIX(3)).real();

   beta_TYu = beta_TYu_1 + beta_TYu_2 + beta_TYu_3;


   return beta_TYu;
}

/**
 * Calculates the 3-loop beta function of TYu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

/**
 * Calculates the 4-loop beta function of TYu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

/**
 * Calculates the 5-loop beta function of TYu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYu_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

} // namespace flexiblesusy
