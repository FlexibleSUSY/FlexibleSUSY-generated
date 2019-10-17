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

// File generated at Wed 16 Oct 2019 22:26:59

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

   beta_TYv = (oneOver16PiSqr*(0.2*(30*traceAdjYuTYu*Yv + 10*traceAdjYvTYv*Yv +
      6*MassB*Yv*Sqr(g1) + 30*MassWB*Yv*Sqr(g2) + 20*MassU*Yv*Sqr(gp)*Sqr(QHu)
      + 20*MassU*Yv*Sqr(gp)*Sqr(Ql) + 20*MassU*Yv*Sqr(gp)*Sqr(Qv) + 15*
      traceYuAdjYu*TYv + 5*traceYvAdjYv*TYv + 5*AbsSqr(Lambdax)*TYv - 3*Sqr(g1)
      *TYv - 15*Sqr(g2)*TYv - 10*Sqr(gp)*Sqr(QHu)*TYv - 10*Sqr(gp)*Sqr(Ql)*TYv
      - 10*Sqr(gp)*Sqr(Qv)*TYv + 10*Yv*Conj(Lambdax)*TLambdax) + 5*(Yv*Yv.
      adjoint()*TYv) + 4*(TYv*Yv.adjoint()*Yv) + Ye.transpose()*Ye.conjugate()*
      TYv + 2*((TYe).transpose()*Ye.conjugate()*Yv))).real();


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


   Eigen::Matrix<double,3,3> beta_TYv;

   const Eigen::Matrix<double,3,3> beta_TYv_1 = ((-0.08*twoLoop*Yv*(25*
      traceAdjYeTYeconjYvTpYv + 25*traceAdjYvTpYeconjYeTYv + 75*
      traceYdAdjYuTYuAdjYd + 75*traceYuAdjYdTYdAdjYu + 450*traceYuAdjYuTYuAdjYu
       + 150*traceYvAdjYvTYvAdjYv + 75*traceAdjYdTYd*AbsSqr(Lambdax) + 25*
      traceAdjYeTYe*AbsSqr(Lambdax) + 207*MassB*Quad(g1) + 375*MassWB*Quad(g2)
      + 400*MassU*Quad(gp)*Quad(QHu) + 800*MassU*Quad(gp)*Quad(Ql) + 500*MassU*
      Quad(gp)*Quad(Qv) - 20*traceAdjYuTYu*Sqr(g1) + 20*MassB*traceYuAdjYu*Sqr(
      g1) + 45*MassB*Sqr(g1)*Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(g2) - 400*
      traceAdjYuTYu*Sqr(g3) + 400*MassG*traceYuAdjYu*Sqr(g3) + 90*MassB*Qd*QHu*
      Sqr(g1)*Sqr(gp) + 90*MassU*Qd*QHu*Sqr(g1)*Sqr(gp) + 90*MassB*Qe*QHu*Sqr(
      g1)*Sqr(gp) + 90*MassU*Qe*QHu*Sqr(g1)*Sqr(gp) - 30*MassB*QHd*QHu*Sqr(g1)*
      Sqr(gp) - 30*MassU*QHd*QHu*Sqr(g1)*Sqr(gp) - 90*MassB*Qd*Ql*Sqr(g1)*Sqr(
      gp) - 90*MassU*Qd*Ql*Sqr(g1)*Sqr(gp) - 90*MassB*Qe*Ql*Sqr(g1)*Sqr(gp) -
      90*MassU*Qe*Ql*Sqr(g1)*Sqr(gp) + 30*MassB*QHd*Ql*Sqr(g1)*Sqr(gp) + 30*
      MassU*QHd*Ql*Sqr(g1)*Sqr(gp) - 120*MassB*QHu*Ql*Sqr(g1)*Sqr(gp) - 120*
      MassU*QHu*Ql*Sqr(g1)*Sqr(gp) + 90*MassB*QHu*Qq*Sqr(g1)*Sqr(gp) + 90*MassU
      *QHu*Qq*Sqr(g1)*Sqr(gp) - 90*MassB*Ql*Qq*Sqr(g1)*Sqr(gp) - 90*MassU*Ql*Qq
      *Sqr(g1)*Sqr(gp) - 180*MassB*QHu*Qu*Sqr(g1)*Sqr(gp) - 180*MassU*QHu*Qu*
      Sqr(g1)*Sqr(gp) + 180*MassB*Ql*Qu*Sqr(g1)*Sqr(gp) + 180*MassU*Ql*Qu*Sqr(
      g1)*Sqr(gp) + 50*MassU*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 150*
      traceAdjYuTYu*Sqr(gp)*Sqr(QHu) + 50*traceAdjYvTYv*Sqr(gp)*Sqr(QHu) - 150*
      MassU*traceYuAdjYu*Sqr(gp)*Sqr(QHu) - 50*MassU*traceYvAdjYv*Sqr(gp)*Sqr(
      QHu) - 50*MassU*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 60*MassB*Sqr(g1)*Sqr(
      gp)*Sqr(QHu) + 60*MassU*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 150*MassU*Sqr(g2)*Sqr(
      gp)*Sqr(QHu) + 150*MassWB*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 900*MassU*Quad(gp)*
      Sqr(Qd)*Sqr(QHu) + 300*MassU*Quad(gp)*Sqr(Qe)*Sqr(QHu) + 200*MassU*Quad(
      gp)*Sqr(QHd)*Sqr(QHu) - 50*traceAdjYvTYv*Sqr(gp)*Sqr(Ql) + 50*MassU*
      traceYvAdjYv*Sqr(gp)*Sqr(Ql) + 120*MassB*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 120*
      MassU*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 150*MassU*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 150*
      MassWB*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 900*MassU*Quad(gp)*Sqr(Qd)*Sqr(Ql) + 300
      *MassU*Quad(gp)*Sqr(Qe)*Sqr(Ql) + 200*MassU*Quad(gp)*Sqr(QHd)*Sqr(Ql) +
      800*MassU*Quad(gp)*Sqr(QHu)*Sqr(Ql) - 150*traceAdjYuTYu*Sqr(gp)*Sqr(Qq) +
      150*MassU*traceYuAdjYu*Sqr(gp)*Sqr(Qq) + 1800*MassU*Quad(gp)*Sqr(QHu)*Sqr
      (Qq) + 1800*MassU*Quad(gp)*Sqr(Ql)*Sqr(Qq) + 50*MassU*AbsSqr(Lambdax)*Sqr
      (gp)*Sqr(Qs) + 100*MassU*Quad(gp)*Sqr(QHu)*Sqr(Qs) + 100*MassU*Quad(gp)*
      Sqr(Ql)*Sqr(Qs) - 150*traceAdjYuTYu*Sqr(gp)*Sqr(Qu) + 150*MassU*
      traceYuAdjYu*Sqr(gp)*Sqr(Qu) + 900*MassU*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 900*
      MassU*Quad(gp)*Sqr(Ql)*Sqr(Qu) - 50*traceAdjYvTYv*Sqr(gp)*Sqr(Qv) + 50*
      MassU*traceYvAdjYv*Sqr(gp)*Sqr(Qv) + 900*MassU*Quad(gp)*Sqr(Qd)*Sqr(Qv) +
      300*MassU*Quad(gp)*Sqr(Qe)*Sqr(Qv) + 200*MassU*Quad(gp)*Sqr(QHd)*Sqr(Qv)
      + 500*MassU*Quad(gp)*Sqr(QHu)*Sqr(Qv) + 900*MassU*Quad(gp)*Sqr(Ql)*Sqr(Qv
      ) + 1800*MassU*Quad(gp)*Sqr(Qq)*Sqr(Qv) + 100*MassU*Quad(gp)*Sqr(Qs)*Sqr(
      Qv) + 900*MassU*Quad(gp)*Sqr(Qu)*Sqr(Qv)) - 0.4*twoLoop*(45*traceAdjYuTYu
       + 15*traceAdjYvTYv + 6*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 30*MassU*Sqr(
      gp)*Sqr(QHu) + 10*MassU*Sqr(gp)*Sqr(Ql) - 10*MassU*Sqr(gp)*Sqr(Qv))*(Yv*
      Yv.adjoint()*Yv) + 0.2*twoLoop*(-75*traceYuAdjYu - 25*traceYvAdjYv - 25*
      AbsSqr(Lambdax) + 12*Sqr(g1) + 60*Sqr(g2) + 50*Sqr(gp)*Sqr(QHu) + 30*Sqr(
      gp)*Sqr(Ql) - 30*Sqr(gp)*Sqr(Qv))*(Yv*Yv.adjoint()*TYv))*UNITMATRIX(3)).
      real();
   const Eigen::Matrix<double,3,3> beta_TYv_2 = ((0.02*twoLoop*(-150*
      traceYdAdjYuYuAdjYd*TYv - 450*traceYuAdjYuYuAdjYu*TYv - 50*
      traceYvAdjYvTpYeconjYe*TYv - 150*traceYvAdjYvYvAdjYv*TYv - 150*
      traceYdAdjYd*AbsSqr(Lambdax)*TYv - 50*traceYeAdjYe*AbsSqr(Lambdax)*TYv +
      207*Quad(g1)*TYv + 375*Quad(g2)*TYv + 400*Quad(gp)*Quad(QHu)*TYv + 800*
      Quad(gp)*Quad(Ql)*TYv + 500*Quad(gp)*Quad(Qv)*TYv + 40*traceYuAdjYu*Sqr(
      g1)*TYv + 90*Sqr(g1)*Sqr(g2)*TYv + 800*traceYuAdjYu*Sqr(g3)*TYv + 180*Qd*
      QHu*Sqr(g1)*Sqr(gp)*TYv + 180*Qe*QHu*Sqr(g1)*Sqr(gp)*TYv - 60*QHd*QHu*Sqr
      (g1)*Sqr(gp)*TYv - 180*Qd*Ql*Sqr(g1)*Sqr(gp)*TYv - 180*Qe*Ql*Sqr(g1)*Sqr(
      gp)*TYv + 60*QHd*Ql*Sqr(g1)*Sqr(gp)*TYv - 240*QHu*Ql*Sqr(g1)*Sqr(gp)*TYv
      + 180*QHu*Qq*Sqr(g1)*Sqr(gp)*TYv - 180*Ql*Qq*Sqr(g1)*Sqr(gp)*TYv - 360*
      QHu*Qu*Sqr(g1)*Sqr(gp)*TYv + 360*Ql*Qu*Sqr(g1)*Sqr(gp)*TYv + 100*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(QHd)*TYv - 300*traceYuAdjYu*Sqr(gp)*Sqr(QHu)*TYv -
      100*traceYvAdjYv*Sqr(gp)*Sqr(QHu)*TYv - 100*AbsSqr(Lambdax)*Sqr(gp)*Sqr(
      QHu)*TYv + 120*Sqr(g1)*Sqr(gp)*Sqr(QHu)*TYv + 300*Sqr(g2)*Sqr(gp)*Sqr(QHu
      )*TYv + 900*Quad(gp)*Sqr(Qd)*Sqr(QHu)*TYv + 300*Quad(gp)*Sqr(Qe)*Sqr(QHu)
      *TYv + 200*Quad(gp)*Sqr(QHd)*Sqr(QHu)*TYv + 100*traceYvAdjYv*Sqr(gp)*Sqr(
      Ql)*TYv + 240*Sqr(g1)*Sqr(gp)*Sqr(Ql)*TYv + 300*Sqr(g2)*Sqr(gp)*Sqr(Ql)*
      TYv + 900*Quad(gp)*Sqr(Qd)*Sqr(Ql)*TYv + 300*Quad(gp)*Sqr(Qe)*Sqr(Ql)*TYv
       + 200*Quad(gp)*Sqr(QHd)*Sqr(Ql)*TYv + 800*Quad(gp)*Sqr(QHu)*Sqr(Ql)*TYv
      + 300*traceYuAdjYu*Sqr(gp)*Sqr(Qq)*TYv + 1800*Quad(gp)*Sqr(QHu)*Sqr(Qq)*
      TYv + 1800*Quad(gp)*Sqr(Ql)*Sqr(Qq)*TYv + 100*AbsSqr(Lambdax)*Sqr(gp)*Sqr
      (Qs)*TYv + 100*Quad(gp)*Sqr(QHu)*Sqr(Qs)*TYv + 100*Quad(gp)*Sqr(Ql)*Sqr(
      Qs)*TYv + 300*traceYuAdjYu*Sqr(gp)*Sqr(Qu)*TYv + 900*Quad(gp)*Sqr(QHu)*
      Sqr(Qu)*TYv + 900*Quad(gp)*Sqr(Ql)*Sqr(Qu)*TYv + 100*traceYvAdjYv*Sqr(gp)
      *Sqr(Qv)*TYv + 900*Quad(gp)*Sqr(Qd)*Sqr(Qv)*TYv + 300*Quad(gp)*Sqr(Qe)*
      Sqr(Qv)*TYv + 200*Quad(gp)*Sqr(QHd)*Sqr(Qv)*TYv + 500*Quad(gp)*Sqr(QHu)*
      Sqr(Qv)*TYv + 900*Quad(gp)*Sqr(Ql)*Sqr(Qv)*TYv + 1800*Quad(gp)*Sqr(Qq)*
      Sqr(Qv)*TYv + 100*Quad(gp)*Sqr(Qs)*Sqr(Qv)*TYv + 900*Quad(gp)*Sqr(Qu)*Sqr
      (Qv)*TYv - 150*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYv - 300*traceYdAdjYd*Yv*
      Conj(Lambdax)*TLambdax - 100*traceYeAdjYe*Yv*Conj(Lambdax)*TLambdax + 200
      *Yv*Conj(Lambdax)*Sqr(gp)*Sqr(QHd)*TLambdax - 200*Yv*Conj(Lambdax)*Sqr(gp
      )*Sqr(QHu)*TLambdax + 200*Yv*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*TLambdax) +
      0.4*twoLoop*(-30*traceYuAdjYu - 10*traceYvAdjYv - 10*AbsSqr(Lambdax) + 3*
      Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*(TYv*Yv.adjoint()*Yv) - 0.4*
      twoLoop*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 6*MassB*Sqr(g1) + 10*MassU*
      Sqr(gp)*Sqr(Qe) + 10*MassU*Sqr(gp)*Sqr(QHd) - 10*MassU*Sqr(gp)*Sqr(Ql))*(
      Ye.transpose()*Ye.conjugate()*Yv) + 0.2*twoLoop*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + 6*Sqr(g1) + 10*Sqr(gp)*Sqr(Qe) + 10*
      Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql))*(Ye.transpose()*Ye.conjugate()*TYv
      ) + 0.4*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) +
      6*Sqr(g1) + 10*Sqr(gp)*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql)
      )*((TYe).transpose()*Ye.conjugate()*Yv) - 6*twoLoop*(Yv*Yv.adjoint()*Yv*
      Yv.adjoint()*TYv) - 8*twoLoop*(Yv*Yv.adjoint()*TYv*Yv.adjoint()*Yv) - 4*
      twoLoop*(Yv*Yv.adjoint()*Ye.transpose()*Ye.conjugate()*TYv) - 4*twoLoop*(
      Yv*Yv.adjoint()*(TYe).transpose()*Ye.conjugate()*Yv) - 6*twoLoop*(TYv*Yv.
      adjoint()*Yv*Yv.adjoint()*Yv) - 2*twoLoop*(TYv*Yv.adjoint()*Ye.transpose(
      )*Ye.conjugate()*Yv) - 2*twoLoop*(Ye.transpose()*Ye.conjugate()*Ye.
      transpose()*Ye.conjugate()*TYv) - 4*twoLoop*(Ye.transpose()*Ye.conjugate(
      )*(TYe).transpose()*Ye.conjugate()*Yv) - 4*twoLoop*((TYe).transpose()*Ye.
      conjugate()*Ye.transpose()*Ye.conjugate()*Yv))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYv_3 = ((-12*twoLoop*Yv*Lambdax*Sqr(
      Conj(Lambdax))*TLambdax - 6*twoLoop*Conj(Lambdax)*TLambdax*(Yv*Yv.adjoint
      ()*Yv) - 2*twoLoop*Conj(Lambdax)*TLambdax*(Ye.transpose()*Ye.conjugate()*
      Yv))*UNITMATRIX(3)).real();

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

/**
 * Calculates the 4-loop beta function of TYv.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

/**
 * Calculates the 5-loop beta function of TYv.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYv_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYv;

   beta_TYv = ZEROMATRIX(3,3);


   return beta_TYv;
}

} // namespace flexiblesusy
