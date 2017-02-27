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

// File generated at Mon 27 Feb 2017 13:35:26

#include "UMSSM_two_scale_soft_parameters.hpp"
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
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*((3*traceYuAdjYu + traceYvAdjYv + AbsSqr(
      Lambdax) - 0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr
      (g3) - 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Qq) - 2*Sqr(gp)*Sqr(Qu))*TYu +
      0.13333333333333333*Yu*(13*MassB*Sqr(g1) + 5*(9*traceAdjYuTYu + 3*
      traceAdjYvTYv + 9*MassWB*Sqr(g2) + 16*MassG*Sqr(g3) + 6*MassU*Sqr(gp)*Sqr
      (QHu) + 6*MassU*Sqr(gp)*Sqr(Qq) + 6*MassU*Sqr(gp)*Sqr(Qu)) + 15*Conj(
      Lambdax)*TLambdax) + 2*(Yu*Yd.adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) +
      TYu*Yd.adjoint()*Yd + 5*(TYu*Yu.adjoint()*Yu))).real();


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
   const auto Qv = INPUT(Qv);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYvAdjYvTYvAdjYv = TRACE_STRUCT.traceYvAdjYvTYvAdjYv;
   const double traceAdjYeTYeconjYvTpYv =
      TRACE_STRUCT.traceAdjYeTYeconjYvTpYv;
   const double traceAdjYvTpYeconjYeTYv =
      TRACE_STRUCT.traceAdjYvTpYeconjYeTYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_TYu;

   const Eigen::Matrix<double,3,3> beta_TYu_1 = (-0.008888888888888889*
      twoLoop*Yu*(2743*Power(g1,4)*MassB + 450*MassU*AbsSqr(Lambdax)*Sqr(gp)*
      Sqr(QHd) + 5*Sqr(g1)*(45*(MassB + MassWB)*Sqr(g2) + 2*(68*(MassB + MassG)
      *Sqr(g3) + 3*(-6*traceAdjYuTYu + 6*MassB*traceYuAdjYu + (MassB + MassU)*
      Sqr(gp)*(-9*QHd*QHu - 27*QHu*Ql - 3*QHd*Qq + 30*QHu*Qq - 9*Ql*Qq + 9*Qd*(
      3*QHu + Qq - 4*Qu) + 9*Qe*(3*QHu + Qq - 4*Qu) + 12*QHd*Qu - 66*QHu*Qu +
      36*Ql*Qu - 54*Qq*Qu + 18*Sqr(QHu) + 10*Sqr(Qq) + 88*Sqr(Qu))))) + 25*(-32
      *Power(g3,4)*MassG + 135*Power(g2,4)*MassWB + 18*Sqr(g2)*(4*(MassG +
      MassWB)*Sqr(g3) + 3*(MassU + MassWB)*Sqr(gp)*(Sqr(QHu) + Sqr(Qq))) + 48*
      Sqr(g3)*(-3*traceAdjYuTYu + 3*MassG*traceYuAdjYu + 2*(MassG + MassU)*Sqr(
      gp)*(Sqr(Qq) + Sqr(Qu))) + 9*(traceAdjYeTYeconjYvTpYv +
      traceAdjYvTpYeconjYeTYv + 3*traceYdAdjYuTYuAdjYd + 3*traceYuAdjYdTYdAdjYu
      + 18*traceYuAdjYuTYuAdjYu + 6*traceYvAdjYvTYvAdjYv + 2*Sqr(gp)*((3*
      traceAdjYuTYu + traceAdjYvTYv - MassU*(3*traceYuAdjYu + traceYvAdjYv))*
      Sqr(QHu) - traceAdjYvTYv*Sqr(Ql) + MassU*traceYvAdjYv*Sqr(Ql) - 3*(
      traceAdjYuTYu - MassU*traceYuAdjYu)*Sqr(Qq) - 3*traceAdjYuTYu*Sqr(Qu) + 3
      *MassU*traceYuAdjYu*Sqr(Qu) - traceAdjYvTYv*Sqr(Qv) + MassU*traceYvAdjYv*
      Sqr(Qv)) + 4*Power(gp,4)*MassU*(4*Power(QHu,4) + 20*Power(Qq,4) + 11*
      Power(Qu,4) + 2*Sqr(QHd)*Sqr(QHu) + 6*Sqr(QHu)*Sqr(Ql) + 2*Sqr(QHd)*Sqr(
      Qq) + 20*Sqr(QHu)*Sqr(Qq) + 6*Sqr(Ql)*Sqr(Qq) + Sqr(QHu)*Sqr(Qs) + Sqr(Qq
      )*Sqr(Qs) + 2*Sqr(QHd)*Sqr(Qu) + 11*Sqr(QHu)*Sqr(Qu) + 6*Sqr(Ql)*Sqr(Qu)
      + 27*Sqr(Qq)*Sqr(Qu) + Sqr(Qs)*Sqr(Qu) + 9*Sqr(Qd)*(Sqr(QHu) + Sqr(Qq) +
      Sqr(Qu)) + 3*Sqr(Qe)*(Sqr(QHu) + Sqr(Qq) + Sqr(Qu)) + 3*Sqr(QHu)*Sqr(Qv)
      + 3*Sqr(Qq)*Sqr(Qv) + 3*Sqr(Qu)*Sqr(Qv)))))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYu_2 = ((0.0022222222222222222*
      twoLoop*(900*Yu*AbsSqr(Lambdax)*(-3*traceAdjYdTYd - traceAdjYeTYe + 2*
      MassU*Sqr(gp)*(Sqr(QHu) - Sqr(Qs))) + (2743*Power(g1,4) + 10*Sqr(g1)*(45*
      Sqr(g2) + 2*(68*Sqr(g3) + 3*Sqr(gp)*(-9*QHd*QHu - 27*QHu*Ql - 3*QHd*Qq +
      30*QHu*Qq - 9*Ql*Qq + 9*Qd*(3*QHu + Qq - 4*Qu) + 9*Qe*(3*QHu + Qq - 4*Qu)
      + 12*QHd*Qu - 66*QHu*Qu + 36*Ql*Qu - 54*Qq*Qu + 18*Sqr(QHu) + 10*Sqr(Qq)
      + 88*Sqr(Qu)))) + 25*(135*Power(g2,4) + 36*Sqr(g2)*(4*Sqr(g3) + 3*Sqr(gp
      )*(Sqr(QHu) + Sqr(Qq))) + 4*(-8*Power(g3,4) + 48*Sqr(g3)*Sqr(gp)*(Sqr(Qq)
      + Sqr(Qu)) + 9*Power(gp,4)*(4*Power(QHu,4) + 20*Power(Qq,4) + 2*Sqr(QHd)
      *Sqr(QHu) + 6*Sqr(QHu)*Sqr(Ql) + 2*Sqr(QHd)*Sqr(Qq) + 20*Sqr(QHu)*Sqr(Qq)
      + 6*Sqr(Ql)*Sqr(Qq) + Sqr(QHu)*Sqr(Qs) + Sqr(Qq)*Sqr(Qs) + 2*Sqr(QHd)*
      Sqr(Qu) + 11*Sqr(QHu)*Sqr(Qu) + 6*Sqr(Ql)*Sqr(Qu) + 27*Sqr(Qq)*Sqr(Qu) +
      9*Sqr(Qd)*(Sqr(QHu) + Sqr(Qq) + Sqr(Qu)) + 3*Sqr(Qe)*(Sqr(QHu) + Sqr(Qq)
      + Sqr(Qu))))))*TYu) - 0.4*twoLoop*(2*MassB*Sqr(g1) + 5*(3*traceAdjYdTYd +
      traceAdjYeTYe + 2*MassU*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq))))*(Yu*
      Yd.adjoint()*Yd) + 0.4*twoLoop*(-5*AbsSqr(Lambdax) + 2*Sqr(g1) + 5*(-3*
      traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq))))*
      (Yu*Yd.adjoint()*TYd) - 0.4*twoLoop*(2*MassB*Sqr(g1) + 5*(3*(3*
      traceAdjYuTYu + traceAdjYvTYv) + 6*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)*(3*
      Sqr(QHu) + Sqr(Qq) - Sqr(Qu))))*(Yu*Yu.adjoint()*Yu) + twoLoop*(-12*
      traceYuAdjYu - 4*traceYvAdjYv - 4*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*Sqr(
      g2) + 8*Sqr(gp)*Sqr(QHu))*(Yu*Yu.adjoint()*TYu) + 0.2*twoLoop*(-5*AbsSqr(
      Lambdax) + 2*Sqr(g1) + 5*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr
      (Qd) + Sqr(QHd) - Sqr(Qq))))*(TYu*Yd.adjoint()*Yd) + twoLoop*(-15*
      traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 12*Sqr(g2) + 10*Sqr(
      gp)*Sqr(QHu) + 6*Sqr(gp)*Sqr(Qq) - 6*Sqr(gp)*Sqr(Qu))*(TYu*Yu.adjoint()*
      Yu) - 4*twoLoop*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 2*twoLoop*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*twoLoop*(Yu*Yd.adjoint()*TYd*
      Yd.adjoint()*Yd) - 4*twoLoop*(Yu*Yd.adjoint()*TYd*Yu.adjoint()*Yu) - 6*
      twoLoop*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 8*twoLoop*(Yu*Yu.adjoint(
      )*TYu*Yu.adjoint()*Yu) - 2*twoLoop*(TYu*Yd.adjoint()*Yd*Yd.adjoint()*Yd)
      - 4*twoLoop*(TYu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 6*twoLoop*(TYu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYu_3 = ((0.2*twoLoop*((-15*
      traceYdAdjYuYuAdjYd - 45*traceYuAdjYuYuAdjYu - 5*traceYvAdjYvTpYeconjYe -
      15*traceYvAdjYvYvAdjYv + 4*traceYuAdjYu*Sqr(g1) + 80*traceYuAdjYu*Sqr(g3
      ) + 5*AbsSqr(Lambdax)*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(
      QHd) - Sqr(QHu) + Sqr(Qs))) + 10*Sqr(gp)*(-((3*traceYuAdjYu +
      traceYvAdjYv)*Sqr(QHu)) + traceYvAdjYv*Sqr(Ql) + 3*traceYuAdjYu*Sqr(Qq) +
      3*traceYuAdjYu*Sqr(Qu) + traceYvAdjYv*Sqr(Qv)) + 10*Power(gp,4)*(11*
      Power(Qu,4) + Sqr(Qs)*Sqr(Qu) + 3*(Sqr(QHu) + Sqr(Qq))*Sqr(Qv) + 3*Sqr(Qu
      )*Sqr(Qv)) - 15*Sqr(Conj(Lambdax))*Sqr(Lambdax))*TYu - 10*Yu*Conj(Lambdax
      )*(3*traceYdAdjYd + traceYeAdjYe + 6*AbsSqr(Lambdax) - 2*Sqr(gp)*(Sqr(QHd
      ) - Sqr(QHu) + Sqr(Qs)))*TLambdax) - 2*twoLoop*Conj(Lambdax)*TLambdax*(Yu
      *Yd.adjoint()*Yd) - 6*twoLoop*Conj(Lambdax)*TLambdax*(Yu*Yu.adjoint()*Yu)
      )*UNITMATRIX(3)).real();

   beta_TYu = beta_TYu_1 + beta_TYu_2 + beta_TYu_3;


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
