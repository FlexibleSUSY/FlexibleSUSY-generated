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

// File generated at Tue 12 Jul 2016 11:16:31

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
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceAdjYeTYeconjYvTpYv =
      TRACE_STRUCT.traceAdjYeTYeconjYvTpYv;
   const double traceAdjYvTpYeconjYeTYv =
      TRACE_STRUCT.traceAdjYvTpYeconjYeTYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   const Eigen::Matrix<double,3,3> beta_TYd_1 = (-0.044444444444444446*
      twoLoop*Yd*(287*Power(g1,4)*MassB + Sqr(g1)*(45*(MassB + MassWB)*Sqr(g2)
      + 2*(20*(MassB + MassG)*Sqr(g3) + 3*(3*(traceAdjYdTYd - 3*traceAdjYeTYe -
      MassB*traceYdAdjYd + 3*MassB*traceYeAdjYe) + (MassB + MassU)*Sqr(gp)*(-9
      *QHd*QHu + 27*QHd*Ql - 30*QHd*Qq + 3*QHu*Qq - 9*Ql*Qq + 9*Qe*(-3*QHd + Qq
      ) + 3*Qd*(6*Qe - 11*QHd + 2*QHu - 6*Ql + 9*Qq - 12*Qu) + 54*QHd*Qu - 18*
      Qq*Qu + 22*Sqr(Qd) + 18*Sqr(QHd) + 10*Sqr(Qq))))) + 5*(-32*Power(g3,4)*
      MassG + 135*Power(g2,4)*MassWB + 48*Sqr(g3)*(-3*traceAdjYdTYd + 3*MassG*
      traceYdAdjYd + 2*(MassG + MassU)*Sqr(gp)*(Sqr(Qd) + Sqr(Qq))) + 18*Sqr(g2
      )*(4*(MassG + MassWB)*Sqr(g3) + 3*(MassU + MassWB)*Sqr(gp)*(Sqr(QHd) +
      Sqr(Qq))) + 9*(traceAdjYeTYeconjYvTpYv + traceAdjYvTpYeconjYeTYv + 18*
      traceYdAdjYdTYdAdjYd + 3*traceYdAdjYuTYuAdjYd + 6*traceYeAdjYeTYeAdjYe -
      2*Sqr(gp)*(3*(traceAdjYdTYd - MassU*traceYdAdjYd)*Sqr(Qd) + traceAdjYeTYe
      *Sqr(Qe) - MassU*traceYeAdjYe*Sqr(Qe) + (-3*traceAdjYdTYd - traceAdjYeTYe
      + 3*MassU*traceYdAdjYd + MassU*traceYeAdjYe)*Sqr(QHd) + traceAdjYeTYe*
      Sqr(Ql) - MassU*traceYeAdjYe*Sqr(Ql) + 3*traceAdjYdTYd*Sqr(Qq) - 3*MassU*
      traceYdAdjYd*Sqr(Qq)) + 4*Power(gp,4)*MassU*(11*Power(Qd,4) + 4*Power(QHd
      ,4) + 20*Power(Qq,4) + 2*Sqr(QHd)*Sqr(QHu) + 6*Sqr(QHd)*Sqr(Ql) + 20*Sqr(
      QHd)*Sqr(Qq) + 2*Sqr(QHu)*Sqr(Qq) + 6*Sqr(Ql)*Sqr(Qq) + 3*Sqr(Qe)*(Sqr(
      QHd) + Sqr(Qq)) + Sqr(QHd)*Sqr(Qs) + Sqr(Qq)*Sqr(Qs) + 9*Sqr(QHd)*Sqr(Qu)
      + 9*Sqr(Qq)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv) + 3*Sqr(Qq)*Sqr(Qv) + Sqr(Qd)*(
      3*Sqr(Qe) + 11*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 27*Sqr(Qq) + Sqr(Qs) +
      9*Sqr(Qu) + 3*Sqr(Qv))))))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_2 = (twoLoop*(-6*
      traceYuAdjYdTYdAdjYu*Yd + 3.188888888888889*Power(g1,4)*TYd + 7.5*Power(
      g2,4)*TYd - 1.7777777777777777*Power(g3,4)*TYd + 22*Power(gp,4)*Power(Qd,
      4)*TYd + 8*Power(gp,4)*Power(QHd,4)*TYd + 40*Power(gp,4)*Power(Qq,4)*TYd
      + Sqr(g1)*Sqr(g2)*TYd + 0.8888888888888888*Sqr(g1)*Sqr(g3)*TYd + 8*Sqr(g2
      )*Sqr(g3)*TYd + 2.4*Qd*Qe*Sqr(g1)*Sqr(gp)*TYd - 4.4*Qd*QHd*Sqr(g1)*Sqr(gp
      )*TYd - 3.6*Qe*QHd*Sqr(g1)*Sqr(gp)*TYd + 0.8*Qd*QHu*Sqr(g1)*Sqr(gp)*TYd -
      1.2*QHd*QHu*Sqr(g1)*Sqr(gp)*TYd - 2.4*Qd*Ql*Sqr(g1)*Sqr(gp)*TYd + 3.6*
      QHd*Ql*Sqr(g1)*Sqr(gp)*TYd + 3.6*Qd*Qq*Sqr(g1)*Sqr(gp)*TYd + 1.2*Qe*Qq*
      Sqr(g1)*Sqr(gp)*TYd - 4*QHd*Qq*Sqr(g1)*Sqr(gp)*TYd + 0.4*QHu*Qq*Sqr(g1)*
      Sqr(gp)*TYd - 1.2*Ql*Qq*Sqr(g1)*Sqr(gp)*TYd - 4.8*Qd*Qu*Sqr(g1)*Sqr(gp)*
      TYd + 2.933333333333333*Sqr(g1)*Sqr(gp)*Sqr(Qd)*TYd + 10.666666666666666*
      Sqr(g3)*Sqr(gp)*Sqr(Qd)*TYd + 6*Power(gp,4)*Sqr(Qd)*Sqr(Qe)*TYd + 2.4*Sqr
      (g1)*Sqr(gp)*Sqr(QHd)*TYd + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd)*TYd + 22*Power(gp,
      4)*Sqr(Qd)*Sqr(QHd)*TYd + 6*Power(gp,4)*Sqr(Qe)*Sqr(QHd)*TYd + 4*Power(gp
      ,4)*Sqr(Qd)*Sqr(QHu)*TYd + 4*Power(gp,4)*Sqr(QHd)*Sqr(QHu)*TYd + 12*Power
      (gp,4)*Sqr(Qd)*Sqr(Ql)*TYd + 12*Power(gp,4)*Sqr(QHd)*Sqr(Ql)*TYd +
      1.3333333333333333*Sqr(g1)*Sqr(gp)*Sqr(Qq)*TYd + 6*Sqr(g2)*Sqr(gp)*Sqr(Qq
      )*TYd + 10.666666666666666*Sqr(g3)*Sqr(gp)*Sqr(Qq)*TYd + 54*Power(gp,4)*
      Sqr(Qd)*Sqr(Qq)*TYd + 6*Power(gp,4)*Sqr(Qe)*Sqr(Qq)*TYd + 40*Power(gp,4)*
      Sqr(QHd)*Sqr(Qq)*TYd + 4*Power(gp,4)*Sqr(QHu)*Sqr(Qq)*TYd + 12*Power(gp,4
      )*Sqr(Ql)*Sqr(Qq)*TYd + 2*Power(gp,4)*Sqr(Qd)*Sqr(Qs)*TYd + 2*Power(gp,4)
      *Sqr(QHd)*Sqr(Qs)*TYd + 2*Power(gp,4)*Sqr(Qq)*Sqr(Qs)*TYd - 0.4*(4*MassB*
      Sqr(g1) + 5*(3*(3*traceAdjYdTYd + traceAdjYeTYe) + 6*MassWB*Sqr(g2) + 2*
      MassU*Sqr(gp)*(-Sqr(Qd) + 3*Sqr(QHd) + Sqr(Qq))))*(Yd*Yd.adjoint()*Yd) -
      12*traceYdAdjYd*(Yd*Yd.adjoint()*TYd) - 4*traceYeAdjYe*(Yd*Yd.adjoint()*
      TYd) + 1.2*Sqr(g1)*(Yd*Yd.adjoint()*TYd) + 6*Sqr(g2)*(Yd*Yd.adjoint()*TYd
      ) + 8*Sqr(gp)*Sqr(QHd)*(Yd*Yd.adjoint()*TYd) - 6*traceAdjYuTYu*(Yd*
      Yu.adjoint()*Yu) - 2*traceAdjYvTYv*(Yd*Yu.adjoint()*Yu) - 1.6*MassB*Sqr(
      g1)*(Yd*Yu.adjoint()*Yu) - 4*MassU*Sqr(gp)*Sqr(QHu)*(Yd*Yu.adjoint()*Yu)
      + 4*MassU*Sqr(gp)*Sqr(Qq)*(Yd*Yu.adjoint()*Yu) - 4*MassU*Sqr(gp)*Sqr(Qu)*
      (Yd*Yu.adjoint()*Yu) - 6*traceYuAdjYu*(Yd*Yu.adjoint()*TYu) - 2*
      traceYvAdjYv*(Yd*Yu.adjoint()*TYu) + 1.6*Sqr(g1)*(Yd*Yu.adjoint()*TYu) +
      4*Sqr(gp)*Sqr(QHu)*(Yd*Yu.adjoint()*TYu) - 4*Sqr(gp)*Sqr(Qq)*(Yd*
      Yu.adjoint()*TYu) + 4*Sqr(gp)*Sqr(Qu)*(Yd*Yu.adjoint()*TYu) - 15*
      traceYdAdjYd*(TYd*Yd.adjoint()*Yd) - 5*traceYeAdjYe*(TYd*Yd.adjoint()*Yd)
      + 1.2*Sqr(g1)*(TYd*Yd.adjoint()*Yd) + 12*Sqr(g2)*(TYd*Yd.adjoint()*Yd) -
      6*Sqr(gp)*Sqr(Qd)*(TYd*Yd.adjoint()*Yd) + 10*Sqr(gp)*Sqr(QHd)*(TYd*
      Yd.adjoint()*Yd) + 6*Sqr(gp)*Sqr(Qq)*(TYd*Yd.adjoint()*Yd) + AbsSqr(
      Lambdax)*(-6*traceAdjYuTYu*Yd - 2*traceAdjYvTYv*Yd + 4*MassU*Yd*Sqr(gp)*
      Sqr(QHd) - 4*MassU*Yd*Sqr(gp)*Sqr(QHu) - 4*MassU*Yd*Sqr(gp)*Sqr(Qs) - 4*(
      Yd*Yd.adjoint()*TYd) - 2*(Yd*Yu.adjoint()*TYu) - 5*(TYd*Yd.adjoint()*Yd)
      - TYd*Yu.adjoint()*Yu) - 3*traceYuAdjYu*(TYd*Yu.adjoint()*Yu) -
      traceYvAdjYv*(TYd*Yu.adjoint()*Yu) + 0.8*Sqr(g1)*(TYd*Yu.adjoint()*Yu) +
      2*Sqr(gp)*Sqr(QHu)*(TYd*Yu.adjoint()*Yu) - 2*Sqr(gp)*Sqr(Qq)*(TYd*
      Yu.adjoint()*Yu) + 2*Sqr(gp)*Sqr(Qu)*(TYd*Yu.adjoint()*Yu) - 6*(Yd*
      Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*
      Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*
      Yu.adjoint()*TYu) - 4*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*
      Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*
      Yd) - 4*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_3 = (0.2*twoLoop*((2*Sqr(g1)*
      (-traceYdAdjYd + 3*traceYeAdjYe + 6*(3*QHd - Qq)*Qu*Sqr(gp)) - 5*AbsSqr(
      Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu)
      - Sqr(Qs))) + 5*(-9*traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 16*traceYdAdjYd*Sqr(g3) +
      2*Sqr(gp)*(3*traceYdAdjYd*Sqr(Qd) + traceYeAdjYe*Sqr(Qe) - (3*
      traceYdAdjYd + traceYeAdjYe)*Sqr(QHd) + traceYeAdjYe*Sqr(Ql) + 3*
      traceYdAdjYd*Sqr(Qq)) + 6*Power(gp,4)*(Sqr(Qd) + Sqr(QHd) + Sqr(Qq))*(3*
      Sqr(Qu) + Sqr(Qv))) - 15*Sqr(Conj(Lambdax))*Sqr(Lambdax))*TYd - 10*Conj(
      Lambdax)*TLambdax*(3*traceYuAdjYu*Yd + traceYvAdjYv*Yd + 6*Yd*AbsSqr(
      Lambdax) + 2*Yd*Sqr(gp)*Sqr(QHd) - 2*Yd*Sqr(gp)*Sqr(QHu) - 2*Yd*Sqr(gp)*
      Sqr(Qs) + 3*(Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu))*UNITMATRIX(3))
      .real();

   beta_TYd = beta_TYd_1 + beta_TYd_2 + beta_TYd_3;


   return beta_TYd;
}

/**
 * Calculates the three-loop beta function of TYd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

} // namespace flexiblesusy
