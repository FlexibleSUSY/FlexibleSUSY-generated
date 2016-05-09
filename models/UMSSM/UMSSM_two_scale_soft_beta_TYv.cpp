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

// File generated at Mon 9 May 2016 12:42:52

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
      AbsSqr(Lambdax)*TYv - 0.6*Sqr(g1)*TYv - 3*Sqr(g2)*TYv - 2*Sqr(gp)*Sqr(QHu
      )*TYv - 2*Sqr(gp)*Sqr(Ql)*TYv - 2*Sqr(gp)*Sqr(Qv)*TYv + Yv*(6*
      traceAdjYuTYu + 2*traceAdjYvTYv + 1.2*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
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

   const Eigen::Matrix<double,3,3> beta_TYv_1 = (-0.04*twoLoop*(414*Power
      (g1,4)*MassB*Yv + 750*Power(g2,4)*MassWB*Yv + 800*Power(gp,4)*MassU*Power
      (QHu,4)*Yv + 1600*Power(gp,4)*MassU*Power(Ql,4)*Yv + 1000*Power(gp,4)*
      MassU*Power(Qv,4)*Yv + 50*traceAdjYeTYeconjYvTpYv*Yv + 50*
      traceAdjYvTpYeconjYeTYv*Yv + 150*traceYdAdjYuTYuAdjYd*Yv + 150*
      traceYuAdjYdTYdAdjYu*Yv + 900*traceYuAdjYuTYuAdjYu*Yv + 300*
      traceYvAdjYvTYvAdjYv*Yv - 40*traceAdjYuTYu*Yv*Sqr(g1) + 40*MassB*
      traceYuAdjYu*Yv*Sqr(g1) + 90*MassB*Yv*Sqr(g1)*Sqr(g2) + 90*MassWB*Yv*Sqr(
      g1)*Sqr(g2) - 800*traceAdjYuTYu*Yv*Sqr(g3) + 800*MassG*traceYuAdjYu*Yv*
      Sqr(g3) + 180*MassB*Qd*QHu*Yv*Sqr(g1)*Sqr(gp) + 180*MassU*Qd*QHu*Yv*Sqr(
      g1)*Sqr(gp) + 180*MassB*Qe*QHu*Yv*Sqr(g1)*Sqr(gp) + 180*MassU*Qe*QHu*Yv*
      Sqr(g1)*Sqr(gp) - 60*MassB*QHd*QHu*Yv*Sqr(g1)*Sqr(gp) - 60*MassU*QHd*QHu*
      Yv*Sqr(g1)*Sqr(gp) - 180*MassB*Qd*Ql*Yv*Sqr(g1)*Sqr(gp) - 180*MassU*Qd*Ql
      *Yv*Sqr(g1)*Sqr(gp) - 180*MassB*Qe*Ql*Yv*Sqr(g1)*Sqr(gp) - 180*MassU*Qe*
      Ql*Yv*Sqr(g1)*Sqr(gp) + 60*MassB*QHd*Ql*Yv*Sqr(g1)*Sqr(gp) + 60*MassU*QHd
      *Ql*Yv*Sqr(g1)*Sqr(gp) - 240*MassB*QHu*Ql*Yv*Sqr(g1)*Sqr(gp) - 240*MassU*
      QHu*Ql*Yv*Sqr(g1)*Sqr(gp) + 180*MassB*QHu*Qq*Yv*Sqr(g1)*Sqr(gp) + 180*
      MassU*QHu*Qq*Yv*Sqr(g1)*Sqr(gp) - 180*MassB*Ql*Qq*Yv*Sqr(g1)*Sqr(gp) -
      180*MassU*Ql*Qq*Yv*Sqr(g1)*Sqr(gp) - 360*MassB*QHu*Qu*Yv*Sqr(g1)*Sqr(gp)
      - 360*MassU*QHu*Qu*Yv*Sqr(g1)*Sqr(gp) + 360*MassB*Ql*Qu*Yv*Sqr(g1)*Sqr(gp
      ) + 360*MassU*Ql*Qu*Yv*Sqr(g1)*Sqr(gp) + 300*traceAdjYuTYu*Yv*Sqr(gp)*Sqr
      (QHu) + 100*traceAdjYvTYv*Yv*Sqr(gp)*Sqr(QHu) - 300*MassU*traceYuAdjYu*Yv
      *Sqr(gp)*Sqr(QHu) - 100*MassU*traceYvAdjYv*Yv*Sqr(gp)*Sqr(QHu) + 120*
      MassB*Yv*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 120*MassU*Yv*Sqr(g1)*Sqr(gp)*Sqr(QHu)
      + 300*MassU*Yv*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 300*MassWB*Yv*Sqr(g2)*Sqr(gp)*
      Sqr(QHu) + 1800*Power(gp,4)*MassU*Yv*Sqr(Qd)*Sqr(QHu) + 600*Power(gp,4)*
      MassU*Yv*Sqr(Qe)*Sqr(QHu) + 400*Power(gp,4)*MassU*Yv*Sqr(QHd)*Sqr(QHu) -
      100*traceAdjYvTYv*Yv*Sqr(gp)*Sqr(Ql) + 100*MassU*traceYvAdjYv*Yv*Sqr(gp)*
      Sqr(Ql) + 240*MassB*Yv*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 240*MassU*Yv*Sqr(g1)*Sqr
      (gp)*Sqr(Ql) + 300*MassU*Yv*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 300*MassWB*Yv*Sqr(
      g2)*Sqr(gp)*Sqr(Ql) + 1800*Power(gp,4)*MassU*Yv*Sqr(Qd)*Sqr(Ql) + 600*
      Power(gp,4)*MassU*Yv*Sqr(Qe)*Sqr(Ql) + 400*Power(gp,4)*MassU*Yv*Sqr(QHd)*
      Sqr(Ql) + 1600*Power(gp,4)*MassU*Yv*Sqr(QHu)*Sqr(Ql) - 300*traceAdjYuTYu*
      Yv*Sqr(gp)*Sqr(Qq) + 300*MassU*traceYuAdjYu*Yv*Sqr(gp)*Sqr(Qq) + 3600*
      Power(gp,4)*MassU*Yv*Sqr(QHu)*Sqr(Qq) + 3600*Power(gp,4)*MassU*Yv*Sqr(Ql)
      *Sqr(Qq) + 200*Power(gp,4)*MassU*Yv*Sqr(QHu)*Sqr(Qs) + 200*Power(gp,4)*
      MassU*Yv*Sqr(Ql)*Sqr(Qs) - 300*traceAdjYuTYu*Yv*Sqr(gp)*Sqr(Qu) + 300*
      MassU*traceYuAdjYu*Yv*Sqr(gp)*Sqr(Qu) + 1800*Power(gp,4)*MassU*Yv*Sqr(QHu
      )*Sqr(Qu) + 1800*Power(gp,4)*MassU*Yv*Sqr(Ql)*Sqr(Qu) - 100*traceAdjYvTYv
      *Yv*Sqr(gp)*Sqr(Qv) + 100*MassU*traceYvAdjYv*Yv*Sqr(gp)*Sqr(Qv) + 1800*
      Power(gp,4)*MassU*Yv*Sqr(Qd)*Sqr(Qv) + 600*Power(gp,4)*MassU*Yv*Sqr(Qe)*
      Sqr(Qv) + 400*Power(gp,4)*MassU*Yv*Sqr(QHd)*Sqr(Qv) + 1000*Power(gp,4)*
      MassU*Yv*Sqr(QHu)*Sqr(Qv) + 1800*Power(gp,4)*MassU*Yv*Sqr(Ql)*Sqr(Qv) +
      3600*Power(gp,4)*MassU*Yv*Sqr(Qq)*Sqr(Qv) + 200*Power(gp,4)*MassU*Yv*Sqr(
      Qs)*Sqr(Qv) + 1800*Power(gp,4)*MassU*Yv*Sqr(Qu)*Sqr(Qv) + 10*(6*MassB*Sqr
      (g1) + 5*(3*(3*traceAdjYuTYu + traceAdjYvTYv) + 6*MassWB*Sqr(g2) + 2*
      MassU*Sqr(gp)*(3*Sqr(QHu) + Sqr(Ql) - Sqr(Qv))))*(Yv*Yv.adjoint()*Yv) +
      375*traceYuAdjYu*(Yv*Yv.adjoint()*TYv) + 125*traceYvAdjYv*(Yv*Yv.adjoint(
      )*TYv) - 60*Sqr(g1)*(Yv*Yv.adjoint()*TYv) - 300*Sqr(g2)*(Yv*Yv.adjoint()*
      TYv) - 250*Sqr(gp)*Sqr(QHu)*(Yv*Yv.adjoint()*TYv) - 150*Sqr(gp)*Sqr(Ql)*(
      Yv*Yv.adjoint()*TYv) + 150*Sqr(gp)*Sqr(Qv)*(Yv*Yv.adjoint()*TYv) + 25*
      AbsSqr(Lambdax)*(2*Yv*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassU*Sqr(gp)*
      (Sqr(QHd) - Sqr(QHu) + Sqr(Qs))) + 5*(Yv*Yv.adjoint()*TYv)))*UNITMATRIX(3
      )).real();
   const Eigen::Matrix<double,3,3> beta_TYv_2 = (twoLoop*(4.14*Power(g1,4
      )*TYv + 7.5*Power(g2,4)*TYv + 8*Power(gp,4)*Power(QHu,4)*TYv + 16*Power(
      gp,4)*Power(Ql,4)*TYv + 10*Power(gp,4)*Power(Qv,4)*TYv - 3*
      traceYdAdjYuYuAdjYd*TYv - 9*traceYuAdjYuYuAdjYu*TYv -
      traceYvAdjYvTpYeconjYe*TYv - 3*traceYvAdjYvYvAdjYv*TYv - 3*traceYdAdjYd*
      AbsSqr(Lambdax)*TYv - traceYeAdjYe*AbsSqr(Lambdax)*TYv + 0.8*traceYuAdjYu
      *Sqr(g1)*TYv + 1.8*Sqr(g1)*Sqr(g2)*TYv + 16*traceYuAdjYu*Sqr(g3)*TYv +
      3.6*Qd*QHu*Sqr(g1)*Sqr(gp)*TYv + 3.6*Qe*QHu*Sqr(g1)*Sqr(gp)*TYv - 1.2*QHd
      *QHu*Sqr(g1)*Sqr(gp)*TYv - 3.6*Qd*Ql*Sqr(g1)*Sqr(gp)*TYv - 3.6*Qe*Ql*Sqr(
      g1)*Sqr(gp)*TYv + 1.2*QHd*Ql*Sqr(g1)*Sqr(gp)*TYv - 4.8*QHu*Ql*Sqr(g1)*Sqr
      (gp)*TYv + 3.6*QHu*Qq*Sqr(g1)*Sqr(gp)*TYv - 3.6*Ql*Qq*Sqr(g1)*Sqr(gp)*TYv
      - 7.2*QHu*Qu*Sqr(g1)*Sqr(gp)*TYv + 7.2*Ql*Qu*Sqr(g1)*Sqr(gp)*TYv + 2*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd)*TYv - 6*traceYuAdjYu*Sqr(gp)*Sqr(QHu)*
      TYv - 2*traceYvAdjYv*Sqr(gp)*Sqr(QHu)*TYv - 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr
      (QHu)*TYv + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHu)*TYv + 6*Sqr(g2)*Sqr(gp)*Sqr(QHu)
      *TYv + 18*Power(gp,4)*Sqr(Qd)*Sqr(QHu)*TYv + 6*Power(gp,4)*Sqr(Qe)*Sqr(
      QHu)*TYv + 4*Power(gp,4)*Sqr(QHd)*Sqr(QHu)*TYv + 2*traceYvAdjYv*Sqr(gp)*
      Sqr(Ql)*TYv + 4.8*Sqr(g1)*Sqr(gp)*Sqr(Ql)*TYv + 6*Sqr(g2)*Sqr(gp)*Sqr(Ql)
      *TYv + 18*Power(gp,4)*Sqr(Qd)*Sqr(Ql)*TYv + 6*Power(gp,4)*Sqr(Qe)*Sqr(Ql)
      *TYv + 4*Power(gp,4)*Sqr(QHd)*Sqr(Ql)*TYv + 16*Power(gp,4)*Sqr(QHu)*Sqr(
      Ql)*TYv + 6*traceYuAdjYu*Sqr(gp)*Sqr(Qq)*TYv + 36*Power(gp,4)*Sqr(QHu)*
      Sqr(Qq)*TYv + 36*Power(gp,4)*Sqr(Ql)*Sqr(Qq)*TYv + 2*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(Qs)*TYv + 2*Power(gp,4)*Sqr(QHu)*Sqr(Qs)*TYv + 2*Power(gp,4)*Sqr(
      Ql)*Sqr(Qs)*TYv + 6*traceYuAdjYu*Sqr(gp)*Sqr(Qu)*TYv + 18*Power(gp,4)*Sqr
      (QHu)*Sqr(Qu)*TYv + 18*Power(gp,4)*Sqr(Ql)*Sqr(Qu)*TYv + 2*traceYvAdjYv*
      Sqr(gp)*Sqr(Qv)*TYv + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qv)*TYv + 6*Power(gp,4)*
      Sqr(Qe)*Sqr(Qv)*TYv + 4*Power(gp,4)*Sqr(QHd)*Sqr(Qv)*TYv + 10*Power(gp,4)
      *Sqr(QHu)*Sqr(Qv)*TYv + 18*Power(gp,4)*Sqr(Ql)*Sqr(Qv)*TYv + 36*Power(gp,
      4)*Sqr(Qq)*Sqr(Qv)*TYv + 2*Power(gp,4)*Sqr(Qs)*Sqr(Qv)*TYv + 18*Power(gp,
      4)*Sqr(Qu)*Sqr(Qv)*TYv - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYv - 6*
      traceYdAdjYd*Yv*Conj(Lambdax)*TLambdax - 2*traceYeAdjYe*Yv*Conj(Lambdax)*
      TLambdax + 4*Yv*Conj(Lambdax)*Sqr(gp)*Sqr(QHd)*TLambdax - 4*Yv*Conj(
      Lambdax)*Sqr(gp)*Sqr(QHu)*TLambdax + 4*Yv*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*
      TLambdax + (-12*traceYuAdjYu - 4*traceYvAdjYv - 4*AbsSqr(Lambdax) + 1.2*
      Sqr(g1) + 6*Sqr(g2) + 8*Sqr(gp)*Sqr(QHu))*(TYv*Yv.adjoint()*Yv) - 0.4*(6*
      MassB*Sqr(g1) + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassU*Sqr(gp)*(Sqr
      (Qe) + Sqr(QHd) - Sqr(Ql))))*(Ye.transpose()*Ye.conjugate()*Yv) - 3*
      traceYdAdjYd*(Ye.transpose()*Ye.conjugate()*TYv) - traceYeAdjYe*(
      Ye.transpose()*Ye.conjugate()*TYv) - AbsSqr(Lambdax)*(Ye.transpose()*
      Ye.conjugate()*TYv) + 1.2*Sqr(g1)*(Ye.transpose()*Ye.conjugate()*TYv) + 2
      *Sqr(gp)*Sqr(Qe)*(Ye.transpose()*Ye.conjugate()*TYv) + 2*Sqr(gp)*Sqr(QHd)
      *(Ye.transpose()*Ye.conjugate()*TYv) - 2*Sqr(gp)*Sqr(Ql)*(Ye.transpose()*
      Ye.conjugate()*TYv) - 6*traceYdAdjYd*((TYe).transpose()*Ye.conjugate()*Yv
      ) - 2*traceYeAdjYe*((TYe).transpose()*Ye.conjugate()*Yv) - 2*AbsSqr(
      Lambdax)*((TYe).transpose()*Ye.conjugate()*Yv) + 2.4*Sqr(g1)*((TYe)
      .transpose()*Ye.conjugate()*Yv) + 4*Sqr(gp)*Sqr(Qe)*((TYe).transpose()*
      Ye.conjugate()*Yv) + 4*Sqr(gp)*Sqr(QHd)*((TYe).transpose()*Ye.conjugate()
      *Yv) - 4*Sqr(gp)*Sqr(Ql)*((TYe).transpose()*Ye.conjugate()*Yv) - 6*(Yv*
      Yv.adjoint()*Yv*Yv.adjoint()*TYv) - 8*(Yv*Yv.adjoint()*TYv*Yv.adjoint()*
      Yv) - 4*(Yv*Yv.adjoint()*Ye.transpose()*Ye.conjugate()*TYv) - 4*(Yv*
      Yv.adjoint()*(TYe).transpose()*Ye.conjugate()*Yv) - 6*(TYv*Yv.adjoint()*
      Yv*Yv.adjoint()*Yv) - 2*(TYv*Yv.adjoint()*Ye.transpose()*Ye.conjugate()*
      Yv) - 2*(Ye.transpose()*Ye.conjugate()*Ye.transpose()*Ye.conjugate()*TYv)
      - 4*(Ye.transpose()*Ye.conjugate()*(TYe).transpose()*Ye.conjugate()*Yv)
      - 4*((TYe).transpose()*Ye.conjugate()*Ye.transpose()*Ye.conjugate()*Yv))*
      UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYv_3 = (-2*twoLoop*Conj(Lambdax)
      *TLambdax*(6*Yv*AbsSqr(Lambdax) + 3*(Yv*Yv.adjoint()*Yv) + Ye.transpose()
      *Ye.conjugate()*Yv)*UNITMATRIX(3)).real();

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
