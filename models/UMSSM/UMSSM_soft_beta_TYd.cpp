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

// File generated at Fri 20 Oct 2017 08:52:14

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
 * Calculates the 1-loop beta function of TYd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qq = INPUT(Qq);
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*((3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr
      (g3) - 2*Sqr(gp)*Sqr(Qd) - 2*Sqr(gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(Qq))*TYd +
      0.13333333333333333*Yd*(7*MassB*Sqr(g1) + 5*(9*traceAdjYdTYd + 3*
      traceAdjYeTYe + 9*MassWB*Sqr(g2) + 16*MassG*Sqr(g3) + 6*MassU*Sqr(gp)*Sqr
      (Qd) + 6*MassU*Sqr(gp)*Sqr(QHd) + 6*MassU*Sqr(gp)*Sqr(Qq)) + 15*Conj(
      Lambdax)*TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) +
      5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu)).real();


   return beta_TYd;
}

/**
 * Calculates the 2-loop beta function of TYd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_2_loop(const Soft_traces& soft_traces) const
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
      twoLoop*Yd*(287*MassB*Quad(g1) + Sqr(g1)*(45*(MassB + MassWB)*Sqr(g2) + 2
      *(20*(MassB + MassG)*Sqr(g3) + 3*(3*(traceAdjYdTYd - 3*traceAdjYeTYe -
      MassB*traceYdAdjYd + 3*MassB*traceYeAdjYe) + (MassB + MassU)*Sqr(gp)*(-9*
      QHd*QHu + 27*QHd*Ql - 30*QHd*Qq + 3*QHu*Qq - 9*Ql*Qq + 9*Qe*(-3*QHd + Qq)
      + 3*Qd*(6*Qe - 11*QHd + 2*QHu - 6*Ql + 9*Qq - 12*Qu) + 54*QHd*Qu - 18*Qq
      *Qu + 22*Sqr(Qd) + 18*Sqr(QHd) + 10*Sqr(Qq))))) + 5*(135*MassWB*Quad(g2)
      - 32*MassG*Quad(g3) + 48*Sqr(g3)*(-3*traceAdjYdTYd + 3*MassG*traceYdAdjYd
      + 2*(MassG + MassU)*Sqr(gp)*(Sqr(Qd) + Sqr(Qq))) + 18*Sqr(g2)*(4*(MassG
      + MassWB)*Sqr(g3) + 3*(MassU + MassWB)*Sqr(gp)*(Sqr(QHd) + Sqr(Qq))) + 9*
      (traceAdjYeTYeconjYvTpYv + traceAdjYvTpYeconjYeTYv + 18*
      traceYdAdjYdTYdAdjYd + 3*traceYdAdjYuTYuAdjYd + 6*traceYeAdjYeTYeAdjYe -
      2*Sqr(gp)*(3*(traceAdjYdTYd - MassU*traceYdAdjYd)*Sqr(Qd) + traceAdjYeTYe
      *Sqr(Qe) - MassU*traceYeAdjYe*Sqr(Qe) + (-3*traceAdjYdTYd - traceAdjYeTYe
      + 3*MassU*traceYdAdjYd + MassU*traceYeAdjYe)*Sqr(QHd) + traceAdjYeTYe*
      Sqr(Ql) - MassU*traceYeAdjYe*Sqr(Ql) + 3*traceAdjYdTYd*Sqr(Qq) - 3*MassU*
      traceYdAdjYd*Sqr(Qq)) + 4*MassU*Quad(gp)*(11*Quad(Qd) + 4*Quad(QHd) + 20*
      Quad(Qq) + 2*Sqr(QHd)*Sqr(QHu) + 6*Sqr(QHd)*Sqr(Ql) + 20*Sqr(QHd)*Sqr(Qq)
      + 2*Sqr(QHu)*Sqr(Qq) + 6*Sqr(Ql)*Sqr(Qq) + 3*Sqr(Qe)*(Sqr(QHd) + Sqr(Qq)
      ) + Sqr(QHd)*Sqr(Qs) + Sqr(Qq)*Sqr(Qs) + 9*Sqr(QHd)*Sqr(Qu) + 9*Sqr(Qq)*
      Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv) + 3*Sqr(Qq)*Sqr(Qv) + Sqr(Qd)*(3*Sqr(Qe) +
      11*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 27*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) +
      3*Sqr(Qv))))))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_2 = ((-0.011111111111111112*
      twoLoop*(540*traceYuAdjYdTYdAdjYu*Yd - 180*Yd*AbsSqr(Lambdax)*(-3*
      traceAdjYuTYu - traceAdjYvTYv + 2*MassU*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) -
      Sqr(Qs))) - (287*Quad(g1) + 2*Sqr(g1)*(45*Sqr(g2) + 2*(20*Sqr(g3) + 3*Sqr
      (gp)*(-9*QHd*QHu + 27*QHd*Ql - 30*QHd*Qq + 3*QHu*Qq - 9*Ql*Qq + 9*Qe*(-3*
      QHd + Qq) + 3*Qd*(6*Qe - 11*QHd + 2*QHu - 6*Ql + 9*Qq - 12*Qu) + 22*Sqr(
      Qd) + 18*Sqr(QHd) + 10*Sqr(Qq)))) + 5*(135*Quad(g2) + 36*Sqr(g2)*(4*Sqr(
      g3) + 3*Sqr(gp)*(Sqr(QHd) + Sqr(Qq))) + 4*(-8*Quad(g3) + 48*Sqr(g3)*Sqr(
      gp)*(Sqr(Qd) + Sqr(Qq)) + 9*Quad(gp)*(11*Quad(Qd) + 4*Quad(QHd) + 20*Quad
      (Qq) + 2*Sqr(QHd)*Sqr(QHu) + 6*Sqr(QHd)*Sqr(Ql) + 20*Sqr(QHd)*Sqr(Qq) + 2
      *Sqr(QHu)*Sqr(Qq) + 6*Sqr(Ql)*Sqr(Qq) + 3*Sqr(Qe)*(Sqr(QHd) + Sqr(Qq)) +
      Sqr(QHd)*Sqr(Qs) + Sqr(Qq)*Sqr(Qs) + Sqr(Qd)*(3*Sqr(Qe) + 11*Sqr(QHd) + 2
      *Sqr(QHu) + 6*Sqr(Ql) + 27*Sqr(Qq) + Sqr(Qs))))))*TYd) - 0.4*twoLoop*(4*
      MassB*Sqr(g1) + 5*(3*(3*traceAdjYdTYd + traceAdjYeTYe) + 6*MassWB*Sqr(g2)
      + 2*MassU*Sqr(gp)*(-Sqr(Qd) + 3*Sqr(QHd) + Sqr(Qq))))*(Yd*Yd.adjoint()*
      Yd) + twoLoop*(-12*traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(Lambdax) +
      1.2*Sqr(g1) + 6*Sqr(g2) + 8*Sqr(gp)*Sqr(QHd))*(Yd*Yd.adjoint()*TYd) - 0.4
      *twoLoop*(4*MassB*Sqr(g1) + 5*(3*traceAdjYuTYu + traceAdjYvTYv + 2*MassU*
      Sqr(gp)*(Sqr(QHu) - Sqr(Qq) + Sqr(Qu))))*(Yd*Yu.adjoint()*Yu) + 0.4*
      twoLoop*(-5*AbsSqr(Lambdax) + 4*Sqr(g1) + 5*(-3*traceYuAdjYu -
      traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Qq) + Sqr(Qu))))*(Yd*Yu.adjoint(
      )*TYu) + twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) +
      1.2*Sqr(g1) + 12*Sqr(g2) - 6*Sqr(gp)*Sqr(Qd) + 10*Sqr(gp)*Sqr(QHd) + 6*
      Sqr(gp)*Sqr(Qq))*(TYd*Yd.adjoint()*Yd) + 0.2*twoLoop*(-5*AbsSqr(Lambdax)
      + 4*Sqr(g1) + 5*(-3*traceYuAdjYu - traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHu) -
      Sqr(Qq) + Sqr(Qu))))*(TYd*Yu.adjoint()*Yu) - 6*twoLoop*(Yd*Yd.adjoint()*
      Yd*Yd.adjoint()*TYd) - 8*twoLoop*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) -
      2*twoLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd) - 4*twoLoop*(Yd*
      Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*twoLoop*(Yd*Yu.adjoint()*TYu*
      Yd.adjoint()*Yd) - 4*twoLoop*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*
      twoLoop*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*twoLoop*(TYd*Yu.adjoint
      ()*Yu*Yd.adjoint()*Yd) - 2*twoLoop*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))
      *UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_3 = ((0.2*twoLoop*((2*Sqr(g1)
      *(-traceYdAdjYd + 3*traceYeAdjYe + 6*(3*QHd - Qq)*Qu*Sqr(gp)) - 5*AbsSqr(
      Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu)
      - Sqr(Qs))) + 5*(-9*traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd - 3*
      traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 16*traceYdAdjYd*Sqr(g3) +
      2*Sqr(gp)*(3*traceYdAdjYd*Sqr(Qd) + traceYeAdjYe*Sqr(Qe) - (3*
      traceYdAdjYd + traceYeAdjYe)*Sqr(QHd) + traceYeAdjYe*Sqr(Ql) + 3*
      traceYdAdjYd*Sqr(Qq)) + 6*Quad(gp)*(Sqr(Qd) + Sqr(QHd) + Sqr(Qq))*(3*Sqr(
      Qu) + Sqr(Qv))) - 15*Sqr(Conj(Lambdax))*Sqr(Lambdax))*TYd - 10*Yd*Conj(
      Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 6*AbsSqr(Lambdax) + 2*Sqr(gp)*(
      Sqr(QHd) - Sqr(QHu) - Sqr(Qs)))*TLambdax) - 6*twoLoop*Conj(Lambdax)*
      TLambdax*(Yd*Yd.adjoint()*Yd) - 2*twoLoop*Conj(Lambdax)*TLambdax*(Yd*
      Yu.adjoint()*Yu))*UNITMATRIX(3)).real();

   beta_TYd = beta_TYd_1 + beta_TYd_2 + beta_TYd_3;


   return beta_TYd;
}

/**
 * Calculates the 3-loop beta function of TYd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_TYd_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

} // namespace flexiblesusy
