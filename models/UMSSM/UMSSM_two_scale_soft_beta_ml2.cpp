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

// File generated at Mon 27 Feb 2017 13:35:47

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
 * Calculates the one-loop beta function of ml2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Ql = INPUT(Ql);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + 2*((TYe)
      .adjoint()*TYe) + 2*mHu2*(Yv.conjugate()*Yv.transpose()) + 2*(
      TYv.conjugate()*(TYv).transpose()) + ml2*Ye.adjoint()*Ye + ml2*
      Yv.conjugate()*Yv.transpose() + 2*(Ye.adjoint()*me2*Ye) + Ye.adjoint()*Ye
      *ml2 + 2*(Yv.conjugate()*mvR2*Yv.transpose()) + Yv.conjugate()*
      Yv.transpose()*ml2 - 0.2*(3.872983346207417*g1*Tr11 - 10*gp*Ql*Tr14 + 6*
      AbsSqr(MassB)*Sqr(g1) + 30*AbsSqr(MassWB)*Sqr(g2) + 40*AbsSqr(MassU)*Sqr(
      gp)*Sqr(Ql))*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the two-loop beta function of ml2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Ql = INPUT(Ql);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const auto Qs = INPUT(Qs);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_ml2;

   const Eigen::Matrix<double,3,3> beta_ml2_1 = ((0.4*twoLoop*(-15*
      traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd - 5*
      traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*AbsSqr(Lambdax) - 5*mHu2*
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) + 6*mHd2*Sqr(g1) + 12*AbsSqr(
      MassB)*Sqr(g1) + 10*mHd2*Sqr(gp)*Sqr(Qe) + 10*mHd2*Sqr(gp)*Sqr(QHd) + 20*
      AbsSqr(MassU)*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql)) - 10*mHd2*Sqr(gp)*
      Sqr(Ql))*(Ye.adjoint()*Ye) - 0.4*twoLoop*(6*Conj(MassB)*Sqr(g1) + 5*(3*
      traceconjTYdTpYd + traceconjTYeTpYe + Conj(TLambdax)*Lambdax + 2*Conj(
      MassU)*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*(Ye.adjoint()*TYe) - 0.4*
      twoLoop*(6*MassB*Sqr(g1) + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassU*
      Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*((TYe).adjoint()*Ye) + 0.4*
      twoLoop*(-5*AbsSqr(Lambdax) + 6*Sqr(g1) + 5*(-3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*((TYe).adjoint(
      )*TYe) + 2*twoLoop*(-3*traceconjTYuTpTYu - traceconjTYvTpTYv - 3*
      tracemq2AdjYuYu - 3*tracemu2YuAdjYu - 6*mHu2*traceYuAdjYu - 2*mHu2*
      traceYvAdjYv - traceYvAdjYvconjml2 - traceYvconjmvR2AdjYv - (mHd2 + 2*
      mHu2 + ms2)*AbsSqr(Lambdax) + 2*mHu2*Sqr(gp)*Sqr(QHu) - 2*mHu2*Sqr(gp)*
      Sqr(Ql) + 2*mHu2*Sqr(gp)*Sqr(Qv) + 4*AbsSqr(MassU)*Sqr(gp)*(Sqr(QHu) -
      Sqr(Ql) + Sqr(Qv)))*(Yv.conjugate()*Yv.transpose()) - 2*twoLoop*(3*
      traceconjTYuTpYu + traceconjTYvTpYv + Conj(TLambdax)*Lambdax + 2*Conj(
      MassU)*Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(Yv.conjugate()*(TYv)
      .transpose()) - 2*twoLoop*(3*traceAdjYuTYu + traceAdjYvTYv + 2*MassU*Sqr(
      gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(TYv.conjugate()*Yv.transpose()) + 2*
      twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax) + 2*Sqr(gp)*(
      Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(TYv.conjugate()*(TYv).transpose()) + 0.2*
      twoLoop*(-5*AbsSqr(Lambdax) + 6*Sqr(g1) + 5*(-3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*(ml2*Ye.adjoint
      ()*Ye) + twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(Lambdax) + 2*
      Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(ml2*Yv.conjugate()*Yv.transpose(
      )) + 0.4*twoLoop*(-5*AbsSqr(Lambdax) + 6*Sqr(g1) + 5*(-3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*(Ye.adjoint()*
      me2*Ye) + 0.2*twoLoop*(-5*AbsSqr(Lambdax) + 6*Sqr(g1) + 5*(-3*
      traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qe) + Sqr(QHd) - Sqr(Ql))))*
      (Ye.adjoint()*Ye*ml2))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_ml2_2 = (UNITMATRIX(3)*(-2*
      twoLoop*AbsSqr(TLambdax)*(Ye.adjoint()*Ye) - 2*twoLoop*Conj(Lambdax)*
      TLambdax*((TYe).adjoint()*Ye) - 2*twoLoop*AbsSqr(TLambdax)*(Yv.conjugate(
      )*Yv.transpose()) - 2*twoLoop*Conj(Lambdax)*TLambdax*(TYv.conjugate()*
      Yv.transpose()) + 2*twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(
      Lambdax) + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(Yv.conjugate()*mvR2
      *Yv.transpose()) + twoLoop*(-3*traceYuAdjYu - traceYvAdjYv - AbsSqr(
      Lambdax) + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Ql) + Sqr(Qv)))*(Yv.conjugate()*
      Yv.transpose()*ml2) - 8*mHd2*twoLoop*(Ye.adjoint()*Ye*Ye.adjoint()*Ye) -
      4*twoLoop*(Ye.adjoint()*Ye*(TYe).adjoint()*TYe) - 4*twoLoop*(Ye.adjoint()
      *TYe*(TYe).adjoint()*Ye) - 4*twoLoop*((TYe).adjoint()*Ye*Ye.adjoint()*TYe
      ) - 4*twoLoop*((TYe).adjoint()*TYe*Ye.adjoint()*Ye) - 8*mHu2*twoLoop*(
      Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose()) - 4*twoLoop*
      (Yv.conjugate()*Yv.transpose()*TYv.conjugate()*(TYv).transpose()) - 4*
      twoLoop*(Yv.conjugate()*(TYv).transpose()*TYv.conjugate()*Yv.transpose())
      - 4*twoLoop*(TYv.conjugate()*Yv.transpose()*Yv.conjugate()*(TYv)
      .transpose()) - 4*twoLoop*(TYv.conjugate()*(TYv).transpose()*Yv.conjugate
      ()*Yv.transpose()) - 2*twoLoop*(ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 2*
      twoLoop*(ml2*Yv.conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose())
      - 4*twoLoop*(Ye.adjoint()*me2*Ye*Ye.adjoint()*Ye) - 4*twoLoop*(
      Ye.adjoint()*Ye*ml2*Ye.adjoint()*Ye) - 4*twoLoop*(Ye.adjoint()*Ye*
      Ye.adjoint()*me2*Ye) - 2*twoLoop*(Ye.adjoint()*Ye*Ye.adjoint()*Ye*ml2) -
      4*twoLoop*(Yv.conjugate()*mvR2*Yv.transpose()*Yv.conjugate()*Yv.transpose
      ()) - 4*twoLoop*(Yv.conjugate()*Yv.transpose()*ml2*Yv.conjugate()*
      Yv.transpose()) - 4*twoLoop*(Yv.conjugate()*Yv.transpose()*Yv.conjugate()
      *mvR2*Yv.transpose()) - 2*twoLoop*(Yv.conjugate()*Yv.transpose()*
      Yv.conjugate()*Yv.transpose()*ml2) + 0.04*twoLoop*(3*Conj(MassB)*Sqr(g1)*
      (207*MassB*Sqr(g1) + 5*(3*(2*MassB + MassWB)*Sqr(g2) + 4*(2*MassB + MassU
      )*Ql*(-3*Qd - 3*Qe + QHd - QHu + 4*Ql - 3*Qq + 6*Qu)*Sqr(gp))) + 5*(30*
      Power(g2,4)*Tr22 - 15.491933384829668*g1*gp*Ql*Tr2U114 -
      15.491933384829668*g1*gp*Ql*Tr2U141 - 15.491933384829668*g1*Tr31 + 40*gp*
      Ql*Tr34 + 165*Power(g2,4)*AbsSqr(MassWB) + 6*Tr2U111*Sqr(g1) + 18*AbsSqr(
      MassWB)*Sqr(g1)*Sqr(g2) + 9*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 40*
      Tr2U144*Sqr(gp)*Sqr(Ql) + 120*AbsSqr(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 60
      *MassU*Conj(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 12*Ql*Conj(MassU)*Sqr(gp)*(
      -((MassB + 2*MassU)*(3*Qd + 3*Qe - QHd + QHu - 4*Ql + 3*Qq - 6*Qu)*Sqr(g1
      )) + 5*Ql*((2*MassU + MassWB)*Sqr(g2) + 2*MassU*Sqr(gp)*(9*Sqr(Qd) + 3*
      Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 8*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*
      Sqr(Qu) + 3*Sqr(Qv))))))*UNITMATRIX(3))).real();

   beta_ml2 = beta_ml2_1 + beta_ml2_2;


   return beta_ml2;
}

/**
 * Calculates the three-loop beta function of ml2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_ml2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
