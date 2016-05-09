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

// File generated at Mon 9 May 2016 12:42:43

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
      Ql) + 2*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.adjoint()*TYe) + 2*(Ye*
      Yv.conjugate()*(TYv).transpose()) + 5*(TYe*Ye.adjoint()*Ye) + TYe*
      Yv.conjugate()*Yv.transpose())).real();


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


   Eigen::Matrix<double,3,3> beta_TYe;

   const Eigen::Matrix<double,3,3> beta_TYe_1 = (-0.4*twoLoop*(5*Ye*
      AbsSqr(Lambdax)*(3*traceAdjYuTYu + traceAdjYvTYv + 2*MassU*Sqr(gp)*(-Sqr(
      QHd) + Sqr(QHu) + Sqr(Qs))) + Ye*(135*Power(g1,4)*MassB + Sqr(g1)*(9*(
      MassB + MassWB)*Sqr(g2) + 2*(traceAdjYdTYd - 3*traceAdjYeTYe - MassB*
      traceYdAdjYd + 3*MassB*traceYeAdjYe + 3*(MassB + MassU)*Sqr(gp)*(-(QHd*
      QHu) + 4*QHd*Ql - QHu*Ql + Qd*(6*Qe - 3*(QHd + Ql)) - 3*QHd*Qq - 3*Ql*Qq
      + Qe*(-5*QHd + 2*QHu - 9*Ql + 6*Qq - 12*Qu) + 6*QHd*Qu + 6*Ql*Qu + 10*Sqr
      (Qe) + 2*Sqr(QHd) + 4*Sqr(Ql)))) + 5*(15*Power(g2,4)*MassWB +
      traceAdjYeTYeconjYvTpYv + traceAdjYvTpYeconjYeTYv + 18*
      traceYdAdjYdTYdAdjYd + 3*traceYdAdjYuTYuAdjYd + 6*traceYeAdjYeTYeAdjYe +
      3*traceYuAdjYdTYdAdjYu - 16*traceAdjYdTYd*Sqr(g3) + 16*MassG*traceYdAdjYd
      *Sqr(g3) + 6*(MassU + MassWB)*Sqr(g2)*Sqr(gp)*(Sqr(QHd) + Sqr(Ql)) - 2*
      Sqr(gp)*(3*(traceAdjYdTYd - MassU*traceYdAdjYd)*Sqr(Qd) + traceAdjYeTYe*
      Sqr(Qe) - MassU*traceYeAdjYe*Sqr(Qe) + (-3*traceAdjYdTYd - traceAdjYeTYe
      + 3*MassU*traceYdAdjYd + MassU*traceYeAdjYe)*Sqr(QHd) + traceAdjYeTYe*Sqr
      (Ql) - MassU*traceYeAdjYe*Sqr(Ql) + 3*traceAdjYdTYd*Sqr(Qq) - 3*MassU*
      traceYdAdjYd*Sqr(Qq)) + 4*Power(gp,4)*MassU*(5*Power(Qe,4) + 4*Power(QHd,
      4) + 8*Power(Ql,4) + 2*Sqr(QHd)*Sqr(QHu) + 8*Sqr(QHd)*Sqr(Ql) + 2*Sqr(QHu
      )*Sqr(Ql) + 9*Sqr(Qd)*(Sqr(Qe) + Sqr(QHd) + Sqr(Ql)) + 18*Sqr(QHd)*Sqr(Qq
      ) + 18*Sqr(Ql)*Sqr(Qq) + Sqr(QHd)*Sqr(Qs) + Sqr(Ql)*Sqr(Qs) + 9*Sqr(QHd)*
      Sqr(Qu) + 9*Sqr(Ql)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv) + 3*Sqr(Ql)*Sqr(Qv) +
      Sqr(Qe)*(5*Sqr(QHd) + 2*Sqr(QHu) + 9*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*
      Sqr(Qu) + 3*Sqr(Qv))))) + 10*(3*MassWB*Sqr(g2) - MassU*Sqr(gp)*(Sqr(Qe) -
      3*Sqr(QHd)))*(Ye*Ye.adjoint()*Ye))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYe_2 = (twoLoop*(13.5*Power(g1,4
      )*TYe + 7.5*Power(g2,4)*TYe + 10*Power(gp,4)*Power(Qe,4)*TYe + 8*Power(gp
      ,4)*Power(QHd,4)*TYe + 16*Power(gp,4)*Power(Ql,4)*TYe - 9*
      traceYdAdjYdYdAdjYd*TYe - 3*traceYdAdjYuYuAdjYd*TYe - 0.4*traceYdAdjYd*
      Sqr(g1)*TYe + 1.2*traceYeAdjYe*Sqr(g1)*TYe + 1.8*Sqr(g1)*Sqr(g2)*TYe + 16
      *traceYdAdjYd*Sqr(g3)*TYe + 7.2*Qd*Qe*Sqr(g1)*Sqr(gp)*TYe - 3.6*Qd*QHd*
      Sqr(g1)*Sqr(gp)*TYe - 6*Qe*QHd*Sqr(g1)*Sqr(gp)*TYe + 2.4*Qe*QHu*Sqr(g1)*
      Sqr(gp)*TYe - 1.2*QHd*QHu*Sqr(g1)*Sqr(gp)*TYe - 3.6*Qd*Ql*Sqr(g1)*Sqr(gp)
      *TYe - 10.8*Qe*Ql*Sqr(g1)*Sqr(gp)*TYe + 4.8*QHd*Ql*Sqr(g1)*Sqr(gp)*TYe -
      1.2*QHu*Ql*Sqr(g1)*Sqr(gp)*TYe + 7.2*Qe*Qq*Sqr(g1)*Sqr(gp)*TYe - 3.6*QHd*
      Qq*Sqr(g1)*Sqr(gp)*TYe - 3.6*Ql*Qq*Sqr(g1)*Sqr(gp)*TYe - 14.4*Qe*Qu*Sqr(
      g1)*Sqr(gp)*TYe + 7.2*QHd*Qu*Sqr(g1)*Sqr(gp)*TYe + 7.2*Ql*Qu*Sqr(g1)*Sqr(
      gp)*TYe + 6*traceYdAdjYd*Sqr(gp)*Sqr(Qd)*TYe + 2*traceYeAdjYe*Sqr(gp)*Sqr
      (Qe)*TYe + 12*Sqr(g1)*Sqr(gp)*Sqr(Qe)*TYe + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qe
      )*TYe - 6*traceYdAdjYd*Sqr(gp)*Sqr(QHd)*TYe - 2*traceYeAdjYe*Sqr(gp)*Sqr(
      QHd)*TYe + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHd)*TYe + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd)*
      TYe + 18*Power(gp,4)*Sqr(Qd)*Sqr(QHd)*TYe + 10*Power(gp,4)*Sqr(Qe)*Sqr(
      QHd)*TYe + 4*Power(gp,4)*Sqr(Qe)*Sqr(QHu)*TYe + 4*Power(gp,4)*Sqr(QHd)*
      Sqr(QHu)*TYe + 4.8*Sqr(g1)*Sqr(gp)*Sqr(Ql)*TYe + 6*Sqr(g2)*Sqr(gp)*Sqr(Ql
      )*TYe + 18*Power(gp,4)*Sqr(Qd)*Sqr(Ql)*TYe + 18*Power(gp,4)*Sqr(Qe)*Sqr(
      Ql)*TYe + 16*Power(gp,4)*Sqr(QHd)*Sqr(Ql)*TYe + 4*Power(gp,4)*Sqr(QHu)*
      Sqr(Ql)*TYe + 6*traceYdAdjYd*Sqr(gp)*Sqr(Qq)*TYe + 36*Power(gp,4)*Sqr(Qe)
      *Sqr(Qq)*TYe + 36*Power(gp,4)*Sqr(QHd)*Sqr(Qq)*TYe + 36*Power(gp,4)*Sqr(
      Ql)*Sqr(Qq)*TYe + 2*Power(gp,4)*Sqr(Qe)*Sqr(Qs)*TYe + 2*Power(gp,4)*Sqr(
      QHd)*Sqr(Qs)*TYe + 2*Power(gp,4)*Sqr(Ql)*Sqr(Qs)*TYe + 18*Power(gp,4)*Sqr
      (Qe)*Sqr(Qu)*TYe + 18*Power(gp,4)*Sqr(QHd)*Sqr(Qu)*TYe + 18*Power(gp,4)*
      Sqr(Ql)*Sqr(Qu)*TYe + 6*Power(gp,4)*Sqr(Qe)*Sqr(Qv)*TYe + 6*Power(gp,4)*
      Sqr(QHd)*Sqr(Qv)*TYe + 6*Power(gp,4)*Sqr(Ql)*Sqr(Qv)*TYe - 2*(9*
      traceAdjYdTYd + 3*traceAdjYeTYe + 2*MassU*Sqr(gp)*Sqr(Ql))*(Ye*Ye.adjoint
      ()*Ye) + (-12*traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(Lambdax) + 1.2*Sqr
      (g1) + 6*Sqr(g2) + 8*Sqr(gp)*Sqr(QHd))*(Ye*Ye.adjoint()*TYe) - 6*
      traceAdjYuTYu*(Ye*Yv.conjugate()*Yv.transpose()) - 2*traceAdjYvTYv*(Ye*
      Yv.conjugate()*Yv.transpose()) - 4*MassU*Sqr(gp)*Sqr(QHu)*(Ye*
      Yv.conjugate()*Yv.transpose()) + 4*MassU*Sqr(gp)*Sqr(Ql)*(Ye*Yv.conjugate
      ()*Yv.transpose()) - 4*MassU*Sqr(gp)*Sqr(Qv)*(Ye*Yv.conjugate()*
      Yv.transpose()) - 6*traceYuAdjYu*(Ye*Yv.conjugate()*(TYv).transpose()) -
      2*traceYvAdjYv*(Ye*Yv.conjugate()*(TYv).transpose()) - 2*AbsSqr(Lambdax)*
      (Ye*Yv.conjugate()*(TYv).transpose()) + 4*Sqr(gp)*Sqr(QHu)*(Ye*
      Yv.conjugate()*(TYv).transpose()) - 4*Sqr(gp)*Sqr(Ql)*(Ye*Yv.conjugate()*
      (TYv).transpose()) + 4*Sqr(gp)*Sqr(Qv)*(Ye*Yv.conjugate()*(TYv).transpose
      ()) - 15*traceYdAdjYd*(TYe*Ye.adjoint()*Ye) - 5*traceYeAdjYe*(TYe*
      Ye.adjoint()*Ye) - 5*AbsSqr(Lambdax)*(TYe*Ye.adjoint()*Ye) - 1.2*Sqr(g1)*
      (TYe*Ye.adjoint()*Ye) + 12*Sqr(g2)*(TYe*Ye.adjoint()*Ye) - 6*Sqr(gp)*Sqr(
      Qe)*(TYe*Ye.adjoint()*Ye) + 10*Sqr(gp)*Sqr(QHd)*(TYe*Ye.adjoint()*Ye) + 6
      *Sqr(gp)*Sqr(Ql)*(TYe*Ye.adjoint()*Ye) - 3*traceYuAdjYu*(TYe*Yv.conjugate
      ()*Yv.transpose()) - traceYvAdjYv*(TYe*Yv.conjugate()*Yv.transpose()) -
      AbsSqr(Lambdax)*(TYe*Yv.conjugate()*Yv.transpose()) + 2*Sqr(gp)*Sqr(QHu)*
      (TYe*Yv.conjugate()*Yv.transpose()) - 2*Sqr(gp)*Sqr(Ql)*(TYe*Yv.conjugate
      ()*Yv.transpose()) + 2*Sqr(gp)*Sqr(Qv)*(TYe*Yv.conjugate()*Yv.transpose()
      ) - 6*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*
      Ye.adjoint()*Ye) - 2*(Ye*Yv.conjugate()*Yv.transpose()*Ye.adjoint()*TYe)
      - 4*(Ye*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*(TYv).transpose()) -
      4*(Ye*Yv.conjugate()*(TYv).transpose()*Ye.adjoint()*Ye) - 4*(Ye*
      Yv.conjugate()*(TYv).transpose()*Yv.conjugate()*Yv.transpose()) - 6*(TYe*
      Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 4*(TYe*Yv.conjugate()*Yv.transpose()*
      Ye.adjoint()*Ye) - 2*(TYe*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*
      Yv.transpose()))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYe_3 = (twoLoop*((-3*
      traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 2*traceYeAdjYe*Sqr(gp)*Sqr
      (Ql) - AbsSqr(Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 2*Sqr(gp)*(Sqr(
      QHd) - Sqr(QHu) - Sqr(Qs))) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax))*TYe - 2*
      Conj(Lambdax)*TLambdax*(3*traceYuAdjYu*Ye + traceYvAdjYv*Ye + 6*Ye*AbsSqr
      (Lambdax) + 2*Ye*Sqr(gp)*Sqr(QHd) - 2*Ye*Sqr(gp)*Sqr(QHu) - 2*Ye*Sqr(gp)*
      Sqr(Qs) + 3*(Ye*Ye.adjoint()*Ye) + Ye*Yv.conjugate()*Yv.transpose()))*
      UNITMATRIX(3)).real();

   beta_TYe = beta_TYe_1 + beta_TYe_2 + beta_TYe_3;


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
