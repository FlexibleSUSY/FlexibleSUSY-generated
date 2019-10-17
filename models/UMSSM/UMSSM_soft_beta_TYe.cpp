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

// File generated at Wed 16 Oct 2019 22:26:53

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

   beta_TYe = (oneOver16PiSqr*(0.2*(30*traceAdjYdTYd*Ye + 10*traceAdjYeTYe*Ye +
      18*MassB*Ye*Sqr(g1) + 30*MassWB*Ye*Sqr(g2) + 20*MassU*Ye*Sqr(gp)*Sqr(Qe)
      + 20*MassU*Ye*Sqr(gp)*Sqr(QHd) + 20*MassU*Ye*Sqr(gp)*Sqr(Ql) + 15*
      traceYdAdjYd*TYe + 5*traceYeAdjYe*TYe + 5*AbsSqr(Lambdax)*TYe - 9*Sqr(g1)
      *TYe - 15*Sqr(g2)*TYe - 10*Sqr(gp)*Sqr(Qe)*TYe - 10*Sqr(gp)*Sqr(QHd)*TYe
      - 10*Sqr(gp)*Sqr(Ql)*TYe + 10*Ye*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.
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


   Eigen::Matrix<double,3,3> beta_TYe;

   const Eigen::Matrix<double,3,3> beta_TYe_1 = ((-0.4*twoLoop*Ye*(5*
      traceAdjYeTYeconjYvTpYv + 5*traceAdjYvTpYeconjYeTYv + 90*
      traceYdAdjYdTYdAdjYd + 15*traceYdAdjYuTYuAdjYd + 30*traceYeAdjYeTYeAdjYe
      + 15*traceYuAdjYdTYdAdjYu + 15*traceAdjYuTYu*AbsSqr(Lambdax) + 5*
      traceAdjYvTYv*AbsSqr(Lambdax) + 135*MassB*Quad(g1) + 75*MassWB*Quad(g2) +
      100*MassU*Quad(gp)*Quad(Qe) + 80*MassU*Quad(gp)*Quad(QHd) + 160*MassU*
      Quad(gp)*Quad(Ql) + 2*traceAdjYdTYd*Sqr(g1) - 6*traceAdjYeTYe*Sqr(g1) - 2
      *MassB*traceYdAdjYd*Sqr(g1) + 6*MassB*traceYeAdjYe*Sqr(g1) + 9*MassB*Sqr(
      g1)*Sqr(g2) + 9*MassWB*Sqr(g1)*Sqr(g2) - 80*traceAdjYdTYd*Sqr(g3) + 80*
      MassG*traceYdAdjYd*Sqr(g3) + 36*MassB*Qd*Qe*Sqr(g1)*Sqr(gp) + 36*MassU*Qd
      *Qe*Sqr(g1)*Sqr(gp) - 18*MassB*Qd*QHd*Sqr(g1)*Sqr(gp) - 18*MassU*Qd*QHd*
      Sqr(g1)*Sqr(gp) - 30*MassB*Qe*QHd*Sqr(g1)*Sqr(gp) - 30*MassU*Qe*QHd*Sqr(
      g1)*Sqr(gp) + 12*MassB*Qe*QHu*Sqr(g1)*Sqr(gp) + 12*MassU*Qe*QHu*Sqr(g1)*
      Sqr(gp) - 6*MassB*QHd*QHu*Sqr(g1)*Sqr(gp) - 6*MassU*QHd*QHu*Sqr(g1)*Sqr(
      gp) - 18*MassB*Qd*Ql*Sqr(g1)*Sqr(gp) - 18*MassU*Qd*Ql*Sqr(g1)*Sqr(gp) -
      54*MassB*Qe*Ql*Sqr(g1)*Sqr(gp) - 54*MassU*Qe*Ql*Sqr(g1)*Sqr(gp) + 24*
      MassB*QHd*Ql*Sqr(g1)*Sqr(gp) + 24*MassU*QHd*Ql*Sqr(g1)*Sqr(gp) - 6*MassB*
      QHu*Ql*Sqr(g1)*Sqr(gp) - 6*MassU*QHu*Ql*Sqr(g1)*Sqr(gp) + 36*MassB*Qe*Qq*
      Sqr(g1)*Sqr(gp) + 36*MassU*Qe*Qq*Sqr(g1)*Sqr(gp) - 18*MassB*QHd*Qq*Sqr(g1
      )*Sqr(gp) - 18*MassU*QHd*Qq*Sqr(g1)*Sqr(gp) - 18*MassB*Ql*Qq*Sqr(g1)*Sqr(
      gp) - 18*MassU*Ql*Qq*Sqr(g1)*Sqr(gp) - 72*MassB*Qe*Qu*Sqr(g1)*Sqr(gp) -
      72*MassU*Qe*Qu*Sqr(g1)*Sqr(gp) + 36*MassB*QHd*Qu*Sqr(g1)*Sqr(gp) + 36*
      MassU*QHd*Qu*Sqr(g1)*Sqr(gp) + 36*MassB*Ql*Qu*Sqr(g1)*Sqr(gp) + 36*MassU*
      Ql*Qu*Sqr(g1)*Sqr(gp) - 30*traceAdjYdTYd*Sqr(gp)*Sqr(Qd) + 30*MassU*
      traceYdAdjYd*Sqr(gp)*Sqr(Qd) - 10*traceAdjYeTYe*Sqr(gp)*Sqr(Qe) + 10*
      MassU*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 60*MassB*Sqr(g1)*Sqr(gp)*Sqr(Qe) +
      60*MassU*Sqr(g1)*Sqr(gp)*Sqr(Qe) + 180*MassU*Quad(gp)*Sqr(Qd)*Sqr(Qe) +
      30*traceAdjYdTYd*Sqr(gp)*Sqr(QHd) + 10*traceAdjYeTYe*Sqr(gp)*Sqr(QHd) -
      30*MassU*traceYdAdjYd*Sqr(gp)*Sqr(QHd) - 10*MassU*traceYeAdjYe*Sqr(gp)*
      Sqr(QHd) - 10*MassU*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 12*MassB*Sqr(g1)*
      Sqr(gp)*Sqr(QHd) + 12*MassU*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 30*MassU*Sqr(g2)*
      Sqr(gp)*Sqr(QHd) + 30*MassWB*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 180*MassU*Quad(gp
      )*Sqr(Qd)*Sqr(QHd) + 100*MassU*Quad(gp)*Sqr(Qe)*Sqr(QHd) + 10*MassU*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 40*MassU*Quad(gp)*Sqr(Qe)*Sqr(QHu) +
      40*MassU*Quad(gp)*Sqr(QHd)*Sqr(QHu) - 10*traceAdjYeTYe*Sqr(gp)*Sqr(Ql) +
      10*MassU*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 24*MassB*Sqr(g1)*Sqr(gp)*Sqr(Ql)
      + 24*MassU*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 30*MassU*Sqr(g2)*Sqr(gp)*Sqr(Ql) +
      30*MassWB*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 180*MassU*Quad(gp)*Sqr(Qd)*Sqr(Ql) +
      180*MassU*Quad(gp)*Sqr(Qe)*Sqr(Ql) + 160*MassU*Quad(gp)*Sqr(QHd)*Sqr(Ql)
      + 40*MassU*Quad(gp)*Sqr(QHu)*Sqr(Ql) - 30*traceAdjYdTYd*Sqr(gp)*Sqr(Qq) +
      30*MassU*traceYdAdjYd*Sqr(gp)*Sqr(Qq) + 360*MassU*Quad(gp)*Sqr(Qe)*Sqr(Qq
      ) + 360*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qq) + 360*MassU*Quad(gp)*Sqr(Ql)*Sqr(
      Qq) + 10*MassU*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*MassU*Quad(gp)*Sqr(Qe
      )*Sqr(Qs) + 20*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qs) + 20*MassU*Quad(gp)*Sqr(Ql
      )*Sqr(Qs) + 180*MassU*Quad(gp)*Sqr(Qe)*Sqr(Qu) + 180*MassU*Quad(gp)*Sqr(
      QHd)*Sqr(Qu) + 180*MassU*Quad(gp)*Sqr(Ql)*Sqr(Qu) + 60*MassU*Quad(gp)*Sqr
      (Qe)*Sqr(Qv) + 60*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qv) + 60*MassU*Quad(gp)*Sqr
      (Ql)*Sqr(Qv)) - 4*twoLoop*(3*MassWB*Sqr(g2) - MassU*Sqr(gp)*Sqr(Qe) + 3*
      MassU*Sqr(gp)*Sqr(QHd))*(Ye*Ye.adjoint()*Ye))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYe_2 = ((0.1*twoLoop*(-90*
      traceYdAdjYdYdAdjYd - 30*traceYdAdjYuYuAdjYd + 135*Quad(g1) + 75*Quad(g2)
      + 100*Quad(gp)*Quad(Qe) + 80*Quad(gp)*Quad(QHd) + 160*Quad(gp)*Quad(Ql) -
      4*traceYdAdjYd*Sqr(g1) + 12*traceYeAdjYe*Sqr(g1) + 18*Sqr(g1)*Sqr(g2) +
      160*traceYdAdjYd*Sqr(g3) + 72*Qd*Qe*Sqr(g1)*Sqr(gp) - 36*Qd*QHd*Sqr(g1)*
      Sqr(gp) - 60*Qe*QHd*Sqr(g1)*Sqr(gp) + 24*Qe*QHu*Sqr(g1)*Sqr(gp) - 12*QHd*
      QHu*Sqr(g1)*Sqr(gp) - 36*Qd*Ql*Sqr(g1)*Sqr(gp) - 108*Qe*Ql*Sqr(g1)*Sqr(gp
      ) + 48*QHd*Ql*Sqr(g1)*Sqr(gp) - 12*QHu*Ql*Sqr(g1)*Sqr(gp) + 72*Qe*Qq*Sqr(
      g1)*Sqr(gp) - 36*QHd*Qq*Sqr(g1)*Sqr(gp) - 36*Ql*Qq*Sqr(g1)*Sqr(gp) - 144*
      Qe*Qu*Sqr(g1)*Sqr(gp) + 72*QHd*Qu*Sqr(g1)*Sqr(gp) + 72*Ql*Qu*Sqr(g1)*Sqr(
      gp) + 60*traceYdAdjYd*Sqr(gp)*Sqr(Qd) + 20*traceYeAdjYe*Sqr(gp)*Sqr(Qe) +
      120*Sqr(g1)*Sqr(gp)*Sqr(Qe) + 180*Quad(gp)*Sqr(Qd)*Sqr(Qe) - 60*
      traceYdAdjYd*Sqr(gp)*Sqr(QHd) - 20*traceYeAdjYe*Sqr(gp)*Sqr(QHd) + 24*Sqr
      (g1)*Sqr(gp)*Sqr(QHd) + 60*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 180*Quad(gp)*Sqr(Qd
      )*Sqr(QHd) + 100*Quad(gp)*Sqr(Qe)*Sqr(QHd) + 40*Quad(gp)*Sqr(Qe)*Sqr(QHu)
      + 40*Quad(gp)*Sqr(QHd)*Sqr(QHu) + 48*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 60*Sqr(g2)
      *Sqr(gp)*Sqr(Ql) + 180*Quad(gp)*Sqr(Qd)*Sqr(Ql) + 180*Quad(gp)*Sqr(Qe)*
      Sqr(Ql) + 160*Quad(gp)*Sqr(QHd)*Sqr(Ql) + 40*Quad(gp)*Sqr(QHu)*Sqr(Ql) +
      60*traceYdAdjYd*Sqr(gp)*Sqr(Qq) + 360*Quad(gp)*Sqr(Qe)*Sqr(Qq) + 360*Quad
      (gp)*Sqr(QHd)*Sqr(Qq) + 360*Quad(gp)*Sqr(Ql)*Sqr(Qq) + 20*Quad(gp)*Sqr(Qe
      )*Sqr(Qs) + 20*Quad(gp)*Sqr(QHd)*Sqr(Qs) + 20*Quad(gp)*Sqr(Ql)*Sqr(Qs) +
      180*Quad(gp)*Sqr(Qe)*Sqr(Qu) + 180*Quad(gp)*Sqr(QHd)*Sqr(Qu) + 180*Quad(
      gp)*Sqr(Ql)*Sqr(Qu) + 60*Quad(gp)*Sqr(Qe)*Sqr(Qv) + 60*Quad(gp)*Sqr(QHd)*
      Sqr(Qv) + 60*Quad(gp)*Sqr(Ql)*Sqr(Qv))*TYe - 2*twoLoop*(9*traceAdjYdTYd +
      3*traceAdjYeTYe + 2*MassU*Sqr(gp)*Sqr(Ql))*(Ye*Ye.adjoint()*Ye) + 0.4*
      twoLoop*(-30*traceYdAdjYd - 10*traceYeAdjYe - 10*AbsSqr(Lambdax) + 3*Sqr(
      g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*(Ye*Ye.adjoint()*TYe) - 2*twoLoop
      *(3*traceAdjYuTYu + traceAdjYvTYv + 2*MassU*Sqr(gp)*Sqr(QHu) - 2*MassU*
      Sqr(gp)*Sqr(Ql) + 2*MassU*Sqr(gp)*Sqr(Qv))*(Ye*Yv.conjugate()*Yv.
      transpose()) + 2*twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax
      ) + 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Ql) + 2*Sqr(gp)*Sqr(Qv))*(Ye*Yv.
      conjugate()*(TYv).transpose()) - 0.2*twoLoop*(75*traceYdAdjYd + 25*
      traceYeAdjYe + 25*AbsSqr(Lambdax) + 6*Sqr(g1) - 60*Sqr(g2) + 30*Sqr(gp)*
      Sqr(Qe) - 50*Sqr(gp)*Sqr(QHd) - 30*Sqr(gp)*Sqr(Ql))*(TYe*Ye.adjoint()*Ye)
      + twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax) + 2*Sqr(gp)*
      Sqr(QHu) - 2*Sqr(gp)*Sqr(Ql) + 2*Sqr(gp)*Sqr(Qv))*(TYe*Yv.conjugate()*Yv.
      transpose()) - 6*twoLoop*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe) - 8*
      twoLoop*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 2*twoLoop*(Ye*Yv.
      conjugate()*Yv.transpose()*Ye.adjoint()*TYe) - 4*twoLoop*(Ye*Yv.conjugate
      ()*Yv.transpose()*Yv.conjugate()*(TYv).transpose()) - 4*twoLoop*(Ye*Yv.
      conjugate()*(TYv).transpose()*Ye.adjoint()*Ye) - 4*twoLoop*(Ye*Yv.
      conjugate()*(TYv).transpose()*Yv.conjugate()*Yv.transpose()) - 6*twoLoop*
      (TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 4*twoLoop*(TYe*Yv.conjugate()*Yv.
      transpose()*Ye.adjoint()*Ye) - 2*twoLoop*(TYe*Yv.conjugate()*Yv.transpose
      ()*Yv.conjugate()*Yv.transpose()))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYe_3 = ((twoLoop*(-3*
      traceYeAdjYeYeAdjYe*TYe - traceYvAdjYvTpYeconjYe*TYe - 3*traceYuAdjYu*
      AbsSqr(Lambdax)*TYe - traceYvAdjYv*AbsSqr(Lambdax)*TYe - 2*AbsSqr(Lambdax
      )*Sqr(gp)*Sqr(QHd)*TYe + 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu)*TYe + 2*
      traceYeAdjYe*Sqr(gp)*Sqr(Ql)*TYe + 2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)*TYe
      - 3*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYe - 6*traceYuAdjYu*Ye*Conj(Lambdax)
      *TLambdax - 2*traceYvAdjYv*Ye*Conj(Lambdax)*TLambdax - 4*Ye*Conj(Lambdax)
      *Sqr(gp)*Sqr(QHd)*TLambdax + 4*Ye*Conj(Lambdax)*Sqr(gp)*Sqr(QHu)*TLambdax
       + 4*Ye*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*TLambdax - 12*Ye*Lambdax*Sqr(Conj(
      Lambdax))*TLambdax) - 6*twoLoop*Conj(Lambdax)*TLambdax*(Ye*Ye.adjoint()*
      Ye) - 2*twoLoop*Conj(Lambdax)*TLambdax*(Ye*Yv.conjugate()*Yv.transpose())
      )*UNITMATRIX(3)).real();

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

/**
 * Calculates the 5-loop beta function of TYe.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYe_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
