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

// File generated at Tue 10 Oct 2017 22:14:39

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
 * Calculates the 1-loop beta function of TYv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_1_loop(const Soft_traces& soft_traces) const
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
      )*TYv - 2*Sqr(gp)*Sqr(Ql)*TYv - 2*Sqr(gp)*Sqr(Qv)*TYv + 0.4*Yv*(3*MassB*
      Sqr(g1) + 5*(3*traceAdjYuTYu + traceAdjYvTYv + 3*MassWB*Sqr(g2) + 2*MassU
      *Sqr(gp)*(Sqr(QHu) + Sqr(Ql) + Sqr(Qv))) + 5*Conj(Lambdax)*TLambdax) + 5*
      (Yv*Yv.adjoint()*TYv) + 4*(TYv*Yv.adjoint()*Yv) + Ye.transpose()*
      Ye.conjugate()*TYv + 2*((TYe).transpose()*Ye.conjugate()*Yv))).real();


   return beta_TYv;
}

/**
 * Calculates the 2-loop beta function of TYv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_2_loop(const Soft_traces& soft_traces) const
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

   const Eigen::Matrix<double,3,3> beta_TYv_1 = ((-0.08*twoLoop*Yv*(207*
      MassB*Quad(g1) + 5*Sqr(g1)*(9*(MassB + MassWB)*Sqr(g2) + 2*(-2*
      traceAdjYuTYu + 2*MassB*traceYuAdjYu + 3*(MassB + MassU)*Sqr(gp)*(-(QHd*
      QHu) + 3*Qd*(QHu - Ql) + 3*Qe*(QHu - Ql) + QHd*Ql - 4*QHu*Ql + 3*QHu*Qq -
      3*Ql*Qq - 6*QHu*Qu + 6*Ql*Qu + 2*Sqr(QHu) + 4*Sqr(Ql)))) + 25*AbsSqr(
      Lambdax)*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassU*Sqr(gp)*(Sqr(QHd) -
      Sqr(QHu) + Sqr(Qs))) + 25*(traceAdjYeTYeconjYvTpYv +
      traceAdjYvTpYeconjYeTYv + 3*traceYdAdjYuTYuAdjYd + 3*traceYuAdjYdTYdAdjYu
      + 18*traceYuAdjYuTYuAdjYu + 6*traceYvAdjYvTYvAdjYv + 15*MassWB*Quad(g2)
      - 16*traceAdjYuTYu*Sqr(g3) + 16*MassG*traceYuAdjYu*Sqr(g3) + 6*(MassU +
      MassWB)*Sqr(g2)*Sqr(gp)*(Sqr(QHu) + Sqr(Ql)) + 2*Sqr(gp)*((3*
      traceAdjYuTYu + traceAdjYvTYv - MassU*(3*traceYuAdjYu + traceYvAdjYv))*
      Sqr(QHu) - traceAdjYvTYv*Sqr(Ql) + MassU*traceYvAdjYv*Sqr(Ql) - 3*(
      traceAdjYuTYu - MassU*traceYuAdjYu)*Sqr(Qq) - 3*traceAdjYuTYu*Sqr(Qu) + 3
      *MassU*traceYuAdjYu*Sqr(Qu) - traceAdjYvTYv*Sqr(Qv) + MassU*traceYvAdjYv*
      Sqr(Qv)) + 4*MassU*Quad(gp)*(4*Quad(QHu) + 8*Quad(Ql) + 5*Quad(Qv) + 2*
      Sqr(QHd)*Sqr(QHu) + 2*Sqr(QHd)*Sqr(Ql) + 8*Sqr(QHu)*Sqr(Ql) + 18*Sqr(QHu)
      *Sqr(Qq) + 18*Sqr(Ql)*Sqr(Qq) + Sqr(QHu)*Sqr(Qs) + Sqr(Ql)*Sqr(Qs) + 9*
      Sqr(QHu)*Sqr(Qu) + 9*Sqr(Ql)*Sqr(Qu) + 2*Sqr(QHd)*Sqr(Qv) + 5*Sqr(QHu)*
      Sqr(Qv) + 9*Sqr(Ql)*Sqr(Qv) + 18*Sqr(Qq)*Sqr(Qv) + Sqr(Qs)*Sqr(Qv) + 9*
      Sqr(Qu)*Sqr(Qv) + 9*Sqr(Qd)*(Sqr(QHu) + Sqr(Ql) + Sqr(Qv)) + 3*Sqr(Qe)*(
      Sqr(QHu) + Sqr(Ql) + Sqr(Qv))))) - 0.4*twoLoop*(6*MassB*Sqr(g1) + 5*(3*(3
      *traceAdjYuTYu + traceAdjYvTYv) + 6*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)*(3*
      Sqr(QHu) + Sqr(Ql) - Sqr(Qv))))*(Yv*Yv.adjoint()*Yv) + 0.2*twoLoop*(-25*
      AbsSqr(Lambdax) + 12*Sqr(g1) + 5*(-15*traceYuAdjYu - 5*traceYvAdjYv + 12*
      Sqr(g2) + 10*Sqr(gp)*Sqr(QHu) + 6*Sqr(gp)*Sqr(Ql) - 6*Sqr(gp)*Sqr(Qv)))*(
      Yv*Yv.adjoint()*TYv))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYv_2 = ((0.02*twoLoop*((207*Quad
      (g1) + 10*Sqr(g1)*(9*Sqr(g2) + 2*(2*traceYuAdjYu + 3*Sqr(gp)*(3*Qd*QHu +
      3*Qe*QHu - QHd*QHu - 3*Qd*Ql - 3*Qe*Ql + QHd*Ql - 4*QHu*Ql + 3*QHu*Qq - 3
      *Ql*Qq - 6*QHu*Qu + 6*Ql*Qu + 2*Sqr(QHu) + 4*Sqr(Ql)))) + 50*AbsSqr(
      Lambdax)*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu)
      + Sqr(Qs))) + 25*(15*Quad(g2) + 12*Sqr(g2)*Sqr(gp)*(Sqr(QHu) + Sqr(Ql))
      + 2*(-3*traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu -
      traceYvAdjYvTpYeconjYe - 3*traceYvAdjYvYvAdjYv + 16*traceYuAdjYu*Sqr(g3)
      + 2*Sqr(gp)*(-((3*traceYuAdjYu + traceYvAdjYv)*Sqr(QHu)) + traceYvAdjYv*
      Sqr(Ql) + 3*traceYuAdjYu*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qu) + traceYvAdjYv*
      Sqr(Qv)) + 2*Quad(gp)*(4*Quad(QHu) + 8*Quad(Ql) + 5*Quad(Qv) + 2*Sqr(QHd)
      *Sqr(QHu) + 2*Sqr(QHd)*Sqr(Ql) + 8*Sqr(QHu)*Sqr(Ql) + 18*Sqr(QHu)*Sqr(Qq)
      + 18*Sqr(Ql)*Sqr(Qq) + Sqr(QHu)*Sqr(Qs) + Sqr(Ql)*Sqr(Qs) + 9*Sqr(QHu)*
      Sqr(Qu) + 9*Sqr(Ql)*Sqr(Qu) + 2*Sqr(QHd)*Sqr(Qv) + 5*Sqr(QHu)*Sqr(Qv) + 9
      *Sqr(Ql)*Sqr(Qv) + 18*Sqr(Qq)*Sqr(Qv) + Sqr(Qs)*Sqr(Qv) + 9*Sqr(Qu)*Sqr(
      Qv) + 9*Sqr(Qd)*(Sqr(QHu) + Sqr(Ql) + Sqr(Qv)) + 3*Sqr(Qe)*(Sqr(QHu) +
      Sqr(Ql) + Sqr(Qv))))) - 150*Sqr(Conj(Lambdax))*Sqr(Lambdax))*TYv + 100*Yv
      *Conj(Lambdax)*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(QHd) -
      Sqr(QHu) + Sqr(Qs)))*TLambdax) + twoLoop*(-12*traceYuAdjYu - 4*
      traceYvAdjYv - 4*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*Sqr(g2) + 8*Sqr(gp)*
      Sqr(QHu))*(TYv*Yv.adjoint()*Yv) - 0.4*twoLoop*(6*MassB*Sqr(g1) + 5*(3*
      traceAdjYdTYd + traceAdjYeTYe + 2*MassU*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr
      (Ql))))*(Ye.transpose()*Ye.conjugate()*Yv) + 0.2*twoLoop*(-5*AbsSqr(
      Lambdax) + 6*Sqr(g1) + 5*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr
      (Qe) + Sqr(QHd) - Sqr(Ql))))*(Ye.transpose()*Ye.conjugate()*TYv) + 0.4*
      twoLoop*(-5*AbsSqr(Lambdax) + 6*Sqr(g1) + 5*(-3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*((TYe)
      .transpose()*Ye.conjugate()*Yv) - 6*twoLoop*(Yv*Yv.adjoint()*Yv*
      Yv.adjoint()*TYv) - 8*twoLoop*(Yv*Yv.adjoint()*TYv*Yv.adjoint()*Yv) - 4*
      twoLoop*(Yv*Yv.adjoint()*Ye.transpose()*Ye.conjugate()*TYv) - 4*twoLoop*(
      Yv*Yv.adjoint()*(TYe).transpose()*Ye.conjugate()*Yv) - 6*twoLoop*(TYv*
      Yv.adjoint()*Yv*Yv.adjoint()*Yv) - 2*twoLoop*(TYv*Yv.adjoint()*
      Ye.transpose()*Ye.conjugate()*Yv) - 2*twoLoop*(Ye.transpose()*
      Ye.conjugate()*Ye.transpose()*Ye.conjugate()*TYv) - 4*twoLoop*(
      Ye.transpose()*Ye.conjugate()*(TYe).transpose()*Ye.conjugate()*Yv) - 4*
      twoLoop*((TYe).transpose()*Ye.conjugate()*Ye.transpose()*Ye.conjugate()*
      Yv))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYv_3 = ((-12*twoLoop*Yv*Lambdax*
      Sqr(Conj(Lambdax))*TLambdax - 6*twoLoop*Conj(Lambdax)*TLambdax*(Yv*
      Yv.adjoint()*Yv) - 2*twoLoop*Conj(Lambdax)*TLambdax*(Ye.transpose()*
      Ye.conjugate()*Yv))*UNITMATRIX(3)).real();

   beta_TYv = beta_TYv_1 + beta_TYv_2 + beta_TYv_3;


   return beta_TYv;
}

/**
 * Calculates the 3-loop beta function of TYv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

} // namespace flexiblesusy
