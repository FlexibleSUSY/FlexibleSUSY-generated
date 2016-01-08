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

// File generated at Fri 8 Jan 2016 15:11:40

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
 * Calculates the one-loop beta function of TYv.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qv = INPUT(Qv);
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = (oneOver16PiSqr*(3*traceYuAdjYu*TYv + traceYvAdjYv*TYv +
      AbsSqr(Lambdax)*TYv - 1.8*Sqr(g1)*TYv - 3*Sqr(g2)*TYv - 2*Sqr(gp)*Sqr(QHu
      )*TYv - 2*Sqr(gp)*Sqr(Ql)*TYv - 2*Sqr(gp)*Sqr(Qv)*TYv + Yv*(6*
      traceAdjYuTYu + 2*traceAdjYvTYv + 3.6*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      4*MassU*Sqr(gp)*Sqr(QHu) + 4*MassU*Sqr(gp)*Sqr(Ql) + 4*MassU*Sqr(gp)*Sqr(
      Qv) + 2*Conj(Lambdax)*TLambdax) + 5*(Yv*Yv.adjoint()*TYv) + 4*(TYv*
      Yv.adjoint()*Yv) + Ye.transpose()*Ye.conjugate()*TYv + 2*((TYe).transpose
      ()*Ye.conjugate()*Yv))).real();


   return beta_TYv;
}

/**
 * Calculates the two-loop beta function of TYv.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYv;

   const Eigen::Matrix<double,3,3> beta_TYv_1 = (-0.08*twoLoop*Yv*(837*
      Power(g1,4)*MassB + 50*MassU*AbsSqr(Lambdax)*Sqr(gp)*(Sqr(QHd) - Sqr(QHu)
      ) + 25*(15*Power(g2,4)*MassWB + traceAdjYeTYeconjYvTpYv +
      traceAdjYvTpYeconjYeTYv + 3*traceYdAdjYuTYuAdjYd + 3*traceYuAdjYdTYdAdjYu
      + 18*traceYuAdjYuTYuAdjYu + 6*traceYvAdjYvTYvAdjYv - 16*traceAdjYuTYu*
      Sqr(g3) + 16*MassG*traceYuAdjYu*Sqr(g3) + 6*(MassU + MassWB)*Sqr(g2)*Sqr(
      gp)*(Sqr(QHu) + Sqr(Ql)) + 2*Sqr(gp)*((3*traceAdjYuTYu + traceAdjYvTYv -
      MassU*(3*traceYuAdjYu + traceYvAdjYv))*Sqr(QHu) - traceAdjYvTYv*Sqr(Ql) +
      MassU*traceYvAdjYv*Sqr(Ql) - 3*(traceAdjYuTYu - MassU*traceYuAdjYu)*Sqr(
      Qq) - 3*traceAdjYuTYu*Sqr(Qu) + 3*MassU*traceYuAdjYu*Sqr(Qu) -
      traceAdjYvTYv*Sqr(Qv) + MassU*traceYvAdjYv*Sqr(Qv)) + 4*Power(gp,4)*MassU
      *(4*Power(QHu,4) + 8*Power(Ql,4) + 5*Power(Qv,4) + 2*Sqr(QHd)*Sqr(QHu) +
      2*Sqr(QHd)*Sqr(Ql) + 8*Sqr(QHu)*Sqr(Ql) + 18*Sqr(QHu)*Sqr(Qq) + 18*Sqr(Ql
      )*Sqr(Qq) + Sqr(QHu)*Sqr(Qs) + Sqr(Ql)*Sqr(Qs) + 9*Sqr(QHu)*Sqr(Qu) + 9*
      Sqr(Ql)*Sqr(Qu) + 2*Sqr(QHd)*Sqr(Qv) + 5*Sqr(QHu)*Sqr(Qv) + 9*Sqr(Ql)*Sqr
      (Qv) + 18*Sqr(Qq)*Sqr(Qv) + Sqr(Qs)*Sqr(Qv) + 9*Sqr(Qu)*Sqr(Qv) + 9*Sqr(
      Qd)*(Sqr(QHu) + Sqr(Ql) + Sqr(Qv)) + 3*Sqr(Qe)*(Sqr(QHu) + Sqr(Ql) + Sqr(
      Qv)))) + 5*Sqr(g1)*(9*(MassB + MassWB)*Sqr(g2) + 2*(-2*traceAdjYuTYu - 3*
      traceAdjYvTYv + 2*MassB*traceYuAdjYu + 3*MassB*traceYvAdjYv + 3*(MassB +
      MassU)*Sqr(gp)*(-(QHd*QHu) + QHd*Ql - 4*QHu*Ql + 3*QHu*Qq - 3*Ql*Qq - 6*
      QHu*Qu + 6*Ql*Qu - 2*QHd*Qv + 5*QHu*Qv - 9*Ql*Qv + 6*Qq*Qv - 12*Qu*Qv + 3
      *Qd*(QHu - Ql + 2*Qv) + 3*Qe*(QHu - Ql + 2*Qv) + 2*Sqr(QHu) + 4*Sqr(Ql) +
      10*Sqr(Qv)))))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYv_2 = (twoLoop*(16.74*Power(g1,
      4)*TYv + 7.5*Power(g2,4)*TYv + 8*Power(gp,4)*Power(QHu,4)*TYv + 16*Power(
      gp,4)*Power(Ql,4)*TYv + 1.8*Sqr(g1)*Sqr(g2)*TYv + 3.6*Qd*QHu*Sqr(g1)*Sqr(
      gp)*TYv + 3.6*Qe*QHu*Sqr(g1)*Sqr(gp)*TYv - 1.2*QHd*QHu*Sqr(g1)*Sqr(gp)*
      TYv - 3.6*Qd*Ql*Sqr(g1)*Sqr(gp)*TYv - 3.6*Qe*Ql*Sqr(g1)*Sqr(gp)*TYv + 1.2
      *QHd*Ql*Sqr(g1)*Sqr(gp)*TYv - 4.8*QHu*Ql*Sqr(g1)*Sqr(gp)*TYv + 3.6*QHu*Qq
      *Sqr(g1)*Sqr(gp)*TYv - 3.6*Ql*Qq*Sqr(g1)*Sqr(gp)*TYv - 7.2*QHu*Qu*Sqr(g1)
      *Sqr(gp)*TYv + 7.2*Ql*Qu*Sqr(g1)*Sqr(gp)*TYv + 7.2*Qd*Qv*Sqr(g1)*Sqr(gp)*
      TYv + 7.2*Qe*Qv*Sqr(g1)*Sqr(gp)*TYv - 2.4*QHd*Qv*Sqr(g1)*Sqr(gp)*TYv + 6*
      QHu*Qv*Sqr(g1)*Sqr(gp)*TYv - 10.8*Ql*Qv*Sqr(g1)*Sqr(gp)*TYv + 7.2*Qq*Qv*
      Sqr(g1)*Sqr(gp)*TYv - 14.4*Qu*Qv*Sqr(g1)*Sqr(gp)*TYv + 2.4*Sqr(g1)*Sqr(gp
      )*Sqr(QHu)*TYv + 6*Sqr(g2)*Sqr(gp)*Sqr(QHu)*TYv + 18*Power(gp,4)*Sqr(Qd)*
      Sqr(QHu)*TYv + 6*Power(gp,4)*Sqr(Qe)*Sqr(QHu)*TYv + 4*Power(gp,4)*Sqr(QHd
      )*Sqr(QHu)*TYv + 4.8*Sqr(g1)*Sqr(gp)*Sqr(Ql)*TYv + 6*Sqr(g2)*Sqr(gp)*Sqr(
      Ql)*TYv + 18*Power(gp,4)*Sqr(Qd)*Sqr(Ql)*TYv + 6*Power(gp,4)*Sqr(Qe)*Sqr(
      Ql)*TYv + 4*Power(gp,4)*Sqr(QHd)*Sqr(Ql)*TYv + 16*Power(gp,4)*Sqr(QHu)*
      Sqr(Ql)*TYv + 36*Power(gp,4)*Sqr(QHu)*Sqr(Qq)*TYv + 36*Power(gp,4)*Sqr(Ql
      )*Sqr(Qq)*TYv + 2*Power(gp,4)*Sqr(QHu)*Sqr(Qs)*TYv + 2*Power(gp,4)*Sqr(Ql
      )*Sqr(Qs)*TYv + 18*Power(gp,4)*Sqr(QHu)*Sqr(Qu)*TYv + 18*Power(gp,4)*Sqr(
      Ql)*Sqr(Qu)*TYv + 12*Sqr(g1)*Sqr(gp)*Sqr(Qv)*TYv + 18*Power(gp,4)*Sqr(Qd)
      *Sqr(Qv)*TYv + 6*Power(gp,4)*Sqr(Qe)*Sqr(Qv)*TYv + 4*Power(gp,4)*Sqr(QHd)
      *Sqr(Qv)*TYv + 10*Power(gp,4)*Sqr(QHu)*Sqr(Qv)*TYv + 18*Power(gp,4)*Sqr(
      Ql)*Sqr(Qv)*TYv + 36*Power(gp,4)*Sqr(Qq)*Sqr(Qv)*TYv + 2*Power(gp,4)*Sqr(
      Qs)*Sqr(Qv)*TYv - 2*(3*(3*traceAdjYuTYu + traceAdjYvTYv) + 6*MassWB*Sqr(
      g2) + 2*MassU*Sqr(gp)*(3*Sqr(QHu) + Sqr(Ql) - Sqr(Qv)))*(Yv*Yv.adjoint()*
      Yv) - 15*traceYuAdjYu*(Yv*Yv.adjoint()*TYv) - 5*traceYvAdjYv*(Yv*
      Yv.adjoint()*TYv) - 1.2*Sqr(g1)*(Yv*Yv.adjoint()*TYv) + 12*Sqr(g2)*(Yv*
      Yv.adjoint()*TYv) + 10*Sqr(gp)*Sqr(QHu)*(Yv*Yv.adjoint()*TYv) + 6*Sqr(gp)
      *Sqr(Ql)*(Yv*Yv.adjoint()*TYv) - 6*Sqr(gp)*Sqr(Qv)*(Yv*Yv.adjoint()*TYv)
      - 12*traceYuAdjYu*(TYv*Yv.adjoint()*Yv) - 4*traceYvAdjYv*(TYv*Yv.adjoint(
      )*Yv) + 1.2*Sqr(g1)*(TYv*Yv.adjoint()*Yv) + 6*Sqr(g2)*(TYv*Yv.adjoint()*
      Yv) + 8*Sqr(gp)*Sqr(QHu)*(TYv*Yv.adjoint()*Yv) - 6*traceAdjYdTYd*(
      Ye.transpose()*Ye.conjugate()*Yv) - 2*traceAdjYeTYe*(Ye.transpose()*
      Ye.conjugate()*Yv) - 2.4*MassB*Sqr(g1)*(Ye.transpose()*Ye.conjugate()*Yv)
      - 4*MassU*Sqr(gp)*Sqr(Qe)*(Ye.transpose()*Ye.conjugate()*Yv) - 4*MassU*
      Sqr(gp)*Sqr(QHd)*(Ye.transpose()*Ye.conjugate()*Yv) + 4*MassU*Sqr(gp)*Sqr
      (Ql)*(Ye.transpose()*Ye.conjugate()*Yv) - 3*traceYdAdjYd*(Ye.transpose()*
      Ye.conjugate()*TYv) - traceYeAdjYe*(Ye.transpose()*Ye.conjugate()*TYv) +
      1.2*Sqr(g1)*(Ye.transpose()*Ye.conjugate()*TYv) + 2*Sqr(gp)*Sqr(Qe)*(
      Ye.transpose()*Ye.conjugate()*TYv) + 2*Sqr(gp)*Sqr(QHd)*(Ye.transpose()*
      Ye.conjugate()*TYv) - 2*Sqr(gp)*Sqr(Ql)*(Ye.transpose()*Ye.conjugate()*
      TYv) - 6*traceYdAdjYd*((TYe).transpose()*Ye.conjugate()*Yv) - 2*
      traceYeAdjYe*((TYe).transpose()*Ye.conjugate()*Yv) + 2.4*Sqr(g1)*((TYe)
      .transpose()*Ye.conjugate()*Yv) + 4*Sqr(gp)*Sqr(Qe)*((TYe).transpose()*
      Ye.conjugate()*Yv) + 4*Sqr(gp)*Sqr(QHd)*((TYe).transpose()*Ye.conjugate()
      *Yv) - 4*Sqr(gp)*Sqr(Ql)*((TYe).transpose()*Ye.conjugate()*Yv) - AbsSqr(
      Lambdax)*(6*traceAdjYdTYd*Yv + 2*traceAdjYeTYe*Yv + 4*MassU*Yv*Sqr(gp)*
      Sqr(Qs) + 5*(Yv*Yv.adjoint()*TYv) + 4*(TYv*Yv.adjoint()*Yv) +
      Ye.transpose()*Ye.conjugate()*TYv + 2*((TYe).transpose()*Ye.conjugate()*
      Yv)) - 6*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*TYv) - 8*(Yv*Yv.adjoint()*TYv*
      Yv.adjoint()*Yv) - 4*(Yv*Yv.adjoint()*Ye.transpose()*Ye.conjugate()*TYv)
      - 4*(Yv*Yv.adjoint()*(TYe).transpose()*Ye.conjugate()*Yv) - 6*(TYv*
      Yv.adjoint()*Yv*Yv.adjoint()*Yv) - 2*(TYv*Yv.adjoint()*Ye.transpose()*
      Ye.conjugate()*Yv) - 2*(Ye.transpose()*Ye.conjugate()*Ye.transpose()*
      Ye.conjugate()*TYv) - 4*(Ye.transpose()*Ye.conjugate()*(TYe).transpose()*
      Ye.conjugate()*Yv) - 4*((TYe).transpose()*Ye.conjugate()*Ye.transpose()*
      Ye.conjugate()*Yv))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYv_3 = (0.2*twoLoop*((-15*
      traceYdAdjYuYuAdjYd - 45*traceYuAdjYuYuAdjYu - 5*traceYvAdjYvTpYeconjYe -
      15*traceYvAdjYvYvAdjYv + 4*traceYuAdjYu*Sqr(g1) + 6*traceYvAdjYv*Sqr(g1)
      + 80*traceYuAdjYu*Sqr(g3) + 5*AbsSqr(Lambdax)*(-3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs))) + 10*Sqr(gp)*(-
      ((3*traceYuAdjYu + traceYvAdjYv)*Sqr(QHu)) + traceYvAdjYv*Sqr(Ql) + 3*
      traceYuAdjYu*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qu) + traceYvAdjYv*Sqr(Qv)) +
      10*Power(gp,4)*(5*Power(Qv,4) + 9*Sqr(Qu)*Sqr(Qv)) - 15*Sqr(Conj(Lambdax)
      )*Sqr(Lambdax))*TYv - 10*Conj(Lambdax)*TLambdax*(3*traceYdAdjYd*Yv +
      traceYeAdjYe*Yv + 6*Yv*AbsSqr(Lambdax) - 2*Yv*Sqr(gp)*Sqr(QHd) + 2*Yv*Sqr
      (gp)*Sqr(QHu) - 2*Yv*Sqr(gp)*Sqr(Qs) + 3*(Yv*Yv.adjoint()*Yv) +
      Ye.transpose()*Ye.conjugate()*Yv))*UNITMATRIX(3)).real();

   beta_TYv = beta_TYv_1 + beta_TYv_2 + beta_TYv_3;


   return beta_TYv;
}

/**
 * Calculates the three-loop beta function of TYv.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

} // namespace flexiblesusy
