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

// File generated at Thu 15 Dec 2016 12:53:55

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
 * Calculates the one-loop beta function of me2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qe = INPUT(Qe);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (oneOver16PiSqr*(4*mHd2*(Ye*Ye.adjoint()) + 4*(TYe*(TYe)
      .adjoint()) + 2*(me2*Ye*Ye.adjoint()) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*
      Ye.adjoint()*me2) + 0.4*(3.872983346207417*g1*Tr11 + 5*gp*Qe*Tr14 - 12*
      AbsSqr(MassB)*Sqr(g1) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qe))*UNITMATRIX(3)))
      .real();


   return beta_me2;
}

/**
 * Calculates the two-loop beta function of me2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_me2;

   const Eigen::Matrix<double,3,3> beta_me2_1 = (UNITMATRIX(3)*(-0.8*
      twoLoop*(15*traceconjTYdTpTYd + 5*traceconjTYeTpTYe + 15*tracemd2YdAdjYd
      + 5*traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*AbsSqr(Lambdax) + 5*mHu2*
      AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*AbsSqr(TLambdax) + 3*mHd2*Sqr
      (g1) + 6*AbsSqr(MassB)*Sqr(g1) - 15*mHd2*Sqr(g2) - 30*AbsSqr(MassWB)*Sqr(
      g2) + 10*mHd2*Sqr(gp)*Sqr(Qe) - 10*mHd2*Sqr(gp)*Sqr(QHd) + 20*AbsSqr(
      MassU)*Sqr(gp)*(Sqr(Qe) - Sqr(QHd) - Sqr(Ql)) - 10*mHd2*Sqr(gp)*Sqr(Ql))*
      (Ye*Ye.adjoint()) + 0.8*twoLoop*(3*MassB*Sqr(g1) - 5*(3*traceAdjYdTYd +
      traceAdjYeTYe + 3*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)*(-Sqr(Qe) + Sqr(QHd) +
      Sqr(Ql))) - 5*Conj(Lambdax)*TLambdax)*(Ye*(TYe).adjoint()) - 0.8*twoLoop
      *(-3*Conj(MassB)*Sqr(g1) + 5*(3*traceconjTYdTpYd + traceconjTYeTpYe +
      Conj(TLambdax)*Lambdax + 3*Conj(MassWB)*Sqr(g2) + 2*Conj(MassU)*Sqr(gp)*(
      -Sqr(Qe) + Sqr(QHd) + Sqr(Ql))))*(TYe*Ye.adjoint()) - 0.8*twoLoop*(5*
      AbsSqr(Lambdax) + 3*Sqr(g1) - 5*(-3*traceYdAdjYd - traceYeAdjYe + 3*Sqr(
      g2) + 2*Sqr(gp)*(-Sqr(Qe) + Sqr(QHd) + Sqr(Ql))))*(TYe*(TYe).adjoint()) -
      0.4*twoLoop*(5*AbsSqr(Lambdax) + 3*Sqr(g1) - 5*(-3*traceYdAdjYd -
      traceYeAdjYe + 3*Sqr(g2) + 2*Sqr(gp)*(-Sqr(Qe) + Sqr(QHd) + Sqr(Ql))))*(
      me2*Ye*Ye.adjoint()) - 0.8*twoLoop*(5*AbsSqr(Lambdax) + 3*Sqr(g1) - 5*(-3
      *traceYdAdjYd - traceYeAdjYe + 3*Sqr(g2) + 2*Sqr(gp)*(-Sqr(Qe) + Sqr(QHd)
      + Sqr(Ql))))*(Ye*ml2*Ye.adjoint()) - 0.4*twoLoop*(5*AbsSqr(Lambdax) + 3*
      Sqr(g1) - 5*(-3*traceYdAdjYd - traceYeAdjYe + 3*Sqr(g2) + 2*Sqr(gp)*(-Sqr
      (Qe) + Sqr(QHd) + Sqr(Ql))))*(Ye*Ye.adjoint()*me2) - 8*mHd2*twoLoop*(Ye*
      Ye.adjoint()*Ye*Ye.adjoint()) - 4*twoLoop*(Ye*Ye.adjoint()*TYe*(TYe)
      .adjoint()) - 4*twoLoop*(Ye*(TYe).adjoint()*TYe*Ye.adjoint()) - 4*(mHd2 +
      mHu2)*twoLoop*(Ye*Yv.conjugate()*Yv.transpose()*Ye.adjoint()) - 4*
      twoLoop*(Ye*Yv.conjugate()*(TYv).transpose()*(TYe).adjoint()) - 4*twoLoop
      *(Ye*TYv.conjugate()*(TYv).transpose()*Ye.adjoint()) - 4*twoLoop*(TYe*
      Ye.adjoint()*Ye*(TYe).adjoint()) - 4*twoLoop*(TYe*(TYe).adjoint()*Ye*
      Ye.adjoint()) - 4*twoLoop*(TYe*Yv.conjugate()*Yv.transpose()*(TYe)
      .adjoint()) - 4*twoLoop*(TYe*TYv.conjugate()*Yv.transpose()*Ye.adjoint())
      - 2*twoLoop*(me2*Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 2*twoLoop*(me2*Ye*
      Yv.conjugate()*Yv.transpose()*Ye.adjoint()) - 4*twoLoop*(Ye*ml2*
      Ye.adjoint()*Ye*Ye.adjoint()) - 4*twoLoop*(Ye*ml2*Yv.conjugate()*
      Yv.transpose()*Ye.adjoint()) - 4*twoLoop*(Ye*Ye.adjoint()*me2*Ye*
      Ye.adjoint()) - 4*twoLoop*(Ye*Ye.adjoint()*Ye*ml2*Ye.adjoint()) - 2*
      twoLoop*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*me2) - 4*twoLoop*(Ye*
      Yv.conjugate()*mvR2*Yv.transpose()*Ye.adjoint()) - 4*twoLoop*(Ye*
      Yv.conjugate()*Yv.transpose()*ml2*Ye.adjoint()) - 2*twoLoop*(Ye*
      Yv.conjugate()*Yv.transpose()*Ye.adjoint()*me2) + 0.32*twoLoop*(5*(
      3.872983346207417*g1*(gp*Qe*(Tr2U114 + Tr2U141) + Tr31) + 5*gp*Qe*(gp*Qe*
      Tr2U144 + Tr34) + 3*Tr2U111*Sqr(g1)) + 9*Conj(MassB)*Sqr(g1)*(39*MassB*
      Sqr(g1) + 5*(2*MassB + MassU)*Qd*Qe*Sqr(gp)))*UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_me2_2 = (4.8*Qe*twoLoop*Sqr(gp)*(
      (2*MassB + MassU)*(5*Qe - QHd + QHu - 3*Ql + 3*Qq - 6*Qu)*Conj(MassB)*Sqr
      (g1) + Conj(MassU)*((MassB + 2*MassU)*(3*Qd + 5*Qe - QHd + QHu - 3*Ql + 3
      *Qq - 6*Qu)*Sqr(g1) + 5*MassU*Qe*Sqr(gp)*(9*Sqr(Qd) + 5*Sqr(Qe) + 2*Sqr(
      QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(
      Qv))))*UNITMATRIX(3)).real();

   beta_me2 = beta_me2_1 + beta_me2_2;


   return beta_me2;
}

/**
 * Calculates the three-loop beta function of me2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_me2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

} // namespace flexiblesusy
