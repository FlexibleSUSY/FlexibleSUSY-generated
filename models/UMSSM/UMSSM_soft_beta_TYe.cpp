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

// File generated at Sun 26 Aug 2018 14:21:39

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
 * Calculates the 1-loop beta function of TYe.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(3*traceYdAdjYd*TYe + traceYeAdjYe*TYe + AbsSqr(
      Lambdax)*TYe - 1.8*Sqr(g1)*TYe - 3*Sqr(g2)*TYe - 2*Sqr(gp)*Sqr(Qe)*TYe -
      2*Sqr(gp)*Sqr(QHd)*TYe - 2*Sqr(gp)*Sqr(Ql)*TYe + 0.4*Ye*(9*MassB*Sqr(g1)
      + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 3*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)
      *(Sqr(Qe) + Sqr(QHd) + Sqr(Ql))) + 5*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.
      adjoint()*TYe) + 2*(Ye*Yv.conjugate()*(TYv).transpose()) + 5*(TYe*Ye.
      adjoint()*Ye) + TYe*Yv.conjugate()*Yv.transpose())).real();


   return beta_TYe;
}

/**
 * Calculates the 2-loop beta function of TYe.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_2_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYeTYeconjYvTpYv = TRACE_STRUCT.traceAdjYeTYeconjYvTpYv;
   const double traceAdjYvTpYeconjYeTYv = TRACE_STRUCT.traceAdjYvTpYeconjYeTYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   const Eigen::Matrix<double,3,3> beta_TYe_1 = ((-0.4*twoLoop*Ye*(135*MassB*
      Quad(g1) + Sqr(g1)*(9*(MassB + MassWB)*Sqr(g2) + 2*(traceAdjYdTYd - 3*
      traceAdjYeTYe - MassB*traceYdAdjYd + 3*MassB*traceYeAdjYe + 3*(MassB +
      MassU)*Sqr(gp)*(-(QHd*QHu) + 4*QHd*Ql - QHu*Ql + Qd*(6*Qe - 3*(QHd + Ql))
      - 3*QHd*Qq - 3*Ql*Qq + Qe*(-5*QHd + 2*QHu - 9*Ql + 6*Qq - 12*Qu) + 6*QHd*
      Qu + 6*Ql*Qu + 10*Sqr(Qe) + 2*Sqr(QHd) + 4*Sqr(Ql)))) + 5*AbsSqr(Lambdax)
      *(3*traceAdjYuTYu + traceAdjYvTYv + 2*MassU*Sqr(gp)*(-Sqr(QHd) + Sqr(QHu)
      + Sqr(Qs))) + 5*(traceAdjYeTYeconjYvTpYv + traceAdjYvTpYeconjYeTYv + 18*
      traceYdAdjYdTYdAdjYd + 3*traceYdAdjYuTYuAdjYd + 6*traceYeAdjYeTYeAdjYe +
      3*traceYuAdjYdTYdAdjYu + 15*MassWB*Quad(g2) - 16*traceAdjYdTYd*Sqr(g3) +
      16*MassG*traceYdAdjYd*Sqr(g3) + 6*(MassU + MassWB)*Sqr(g2)*Sqr(gp)*(Sqr(
      QHd) + Sqr(Ql)) - 2*Sqr(gp)*(3*(traceAdjYdTYd - MassU*traceYdAdjYd)*Sqr(
      Qd) + traceAdjYeTYe*Sqr(Qe) - MassU*traceYeAdjYe*Sqr(Qe) + (-3*
      traceAdjYdTYd - traceAdjYeTYe + 3*MassU*traceYdAdjYd + MassU*traceYeAdjYe
      )*Sqr(QHd) + traceAdjYeTYe*Sqr(Ql) - MassU*traceYeAdjYe*Sqr(Ql) + 3*
      traceAdjYdTYd*Sqr(Qq) - 3*MassU*traceYdAdjYd*Sqr(Qq)) + 4*MassU*Quad(gp)*
      (5*Quad(Qe) + 4*Quad(QHd) + 8*Quad(Ql) + 2*Sqr(QHd)*Sqr(QHu) + 8*Sqr(QHd)
      *Sqr(Ql) + 2*Sqr(QHu)*Sqr(Ql) + 9*Sqr(Qd)*(Sqr(Qe) + Sqr(QHd) + Sqr(Ql))
      + 18*Sqr(QHd)*Sqr(Qq) + 18*Sqr(Ql)*Sqr(Qq) + Sqr(QHd)*Sqr(Qs) + Sqr(Ql)*
      Sqr(Qs) + 9*Sqr(QHd)*Sqr(Qu) + 9*Sqr(Ql)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv) + 3
      *Sqr(Ql)*Sqr(Qv) + Sqr(Qe)*(5*Sqr(QHd) + 2*Sqr(QHu) + 9*Sqr(Ql) + 18*Sqr(
      Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv))))) - 4*twoLoop*(3*MassWB*Sqr(g2) -
      MassU*Sqr(gp)*(Sqr(Qe) - 3*Sqr(QHd)))*(Ye*Ye.adjoint()*Ye))*UNITMATRIX(3)
      ).real();
   const Eigen::Matrix<double,3,3> beta_TYe_2 = ((0.1*twoLoop*(135*Quad(g1) + 2
      *Sqr(g1)*(9*Sqr(g2) + 2*(-traceYdAdjYd + 3*traceYeAdjYe + 3*Sqr(gp)*(-(
      QHd*QHu) + 4*QHd*Ql - QHu*Ql + Qd*(6*Qe - 3*(QHd + Ql)) - 3*QHd*Qq - 3*Ql
      *Qq + Qe*(-5*QHd + 2*QHu - 9*Ql + 6*Qq - 12*Qu) + 6*QHd*Qu + 6*Ql*Qu + 10
      *Sqr(Qe) + 2*Sqr(QHd) + 4*Sqr(Ql)))) + 5*(15*Quad(g2) + 12*Sqr(g2)*Sqr(gp
      )*(Sqr(QHd) + Sqr(Ql)) + 2*(-9*traceYdAdjYdYdAdjYd - 3*
      traceYdAdjYuYuAdjYd + 16*traceYdAdjYd*Sqr(g3) + 2*Sqr(gp)*(3*traceYdAdjYd
      *Sqr(Qd) + traceYeAdjYe*Sqr(Qe) - (3*traceYdAdjYd + traceYeAdjYe)*Sqr(QHd
      ) + 3*traceYdAdjYd*Sqr(Qq)) + 2*Quad(gp)*(5*Quad(Qe) + 4*Quad(QHd) + 8*
      Quad(Ql) + 2*Sqr(QHd)*Sqr(QHu) + 8*Sqr(QHd)*Sqr(Ql) + 2*Sqr(QHu)*Sqr(Ql)
      + 9*Sqr(Qd)*(Sqr(Qe) + Sqr(QHd) + Sqr(Ql)) + 18*Sqr(QHd)*Sqr(Qq) + 18*Sqr
      (Ql)*Sqr(Qq) + Sqr(QHd)*Sqr(Qs) + Sqr(Ql)*Sqr(Qs) + 9*Sqr(QHd)*Sqr(Qu) +
      9*Sqr(Ql)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv) + 3*Sqr(Ql)*Sqr(Qv) + Sqr(Qe)*(5*
      Sqr(QHd) + 2*Sqr(QHu) + 9*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*
      Sqr(Qv))))))*TYe - 2*twoLoop*(9*traceAdjYdTYd + 3*traceAdjYeTYe + 2*MassU
      *Sqr(gp)*Sqr(Ql))*(Ye*Ye.adjoint()*Ye) + twoLoop*(-12*traceYdAdjYd - 4*
      traceYeAdjYe - 4*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*Sqr(g2) + 8*Sqr(gp)*
      Sqr(QHd))*(Ye*Ye.adjoint()*TYe) - 2*twoLoop*(3*traceAdjYuTYu +
      traceAdjYvTYv + 2*MassU*Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(Ye*Yv.
      conjugate()*Yv.transpose()) + twoLoop*(-2*(3*traceYuAdjYu + traceYvAdjYv)
      - 2*AbsSqr(Lambdax) + 4*Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(Ye*Yv.
      conjugate()*(TYv).transpose()) + twoLoop*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) - 1.2*Sqr(g1) + 12*Sqr(g2) - 6*Sqr(gp)*
      Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) + 6*Sqr(gp)*Sqr(Ql))*(TYe*Ye.adjoint()*Ye)
      + twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax) + 2*Sqr(gp)*(
      Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(TYe*Yv.conjugate()*Yv.transpose()) - 6*
      twoLoop*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*twoLoop*(Ye*Ye.adjoint(
      )*TYe*Ye.adjoint()*Ye) - 2*twoLoop*(Ye*Yv.conjugate()*Yv.transpose()*Ye.
      adjoint()*TYe) - 4*twoLoop*(Ye*Yv.conjugate()*Yv.transpose()*Yv.conjugate
      ()*(TYv).transpose()) - 4*twoLoop*(Ye*Yv.conjugate()*(TYv).transpose()*Ye
      .adjoint()*Ye) - 4*twoLoop*(Ye*Yv.conjugate()*(TYv).transpose()*Yv.
      conjugate()*Yv.transpose()) - 6*twoLoop*(TYe*Ye.adjoint()*Ye*Ye.adjoint()
      *Ye) - 4*twoLoop*(TYe*Yv.conjugate()*Yv.transpose()*Ye.adjoint()*Ye) - 2*
      twoLoop*(TYe*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose())
      )*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYe_3 = ((twoLoop*((-3*
      traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 2*traceYeAdjYe*Sqr(gp)*Sqr
      (Ql) - AbsSqr(Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 2*Sqr(gp)*(Sqr(
      QHd) - Sqr(QHu) - Sqr(Qs))) - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax))*TYe - 2*
      Ye*Conj(Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 6*AbsSqr(Lambdax) + 2*
      Sqr(gp)*(Sqr(QHd) - Sqr(QHu) - Sqr(Qs)))*TLambdax) - 6*twoLoop*Conj(
      Lambdax)*TLambdax*(Ye*Ye.adjoint()*Ye) - 2*twoLoop*Conj(Lambdax)*TLambdax
      *(Ye*Yv.conjugate()*Yv.transpose()))*UNITMATRIX(3)).real();

   beta_TYe = beta_TYe_1 + beta_TYe_2 + beta_TYe_3;


   return beta_TYe;
}

/**
 * Calculates the 3-loop beta function of TYe.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 4-loop beta function of TYe.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
