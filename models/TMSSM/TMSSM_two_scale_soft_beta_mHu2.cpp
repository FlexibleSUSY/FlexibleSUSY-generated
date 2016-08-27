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

// File generated at Sat 27 Aug 2016 11:55:17

#include "TMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mHu2.
 *
 * @return one-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mHu2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHu2;

   beta_mHu2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 6*
      traceconjTYuTpTYu + 6*tracemq2AdjYuYu + 6*tracemu2YuAdjYu + 6*mHu2*
      traceYuAdjYu + 3*mHd2*AbsSqr(Lambdax) + 3*mHu2*AbsSqr(Lambdax) + 3*mT2*
      AbsSqr(Lambdax) + 3*AbsSqr(TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*
      AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHu2;
}

/**
 * Calculates the two-loop beta function of mHu2.
 *
 * @return two-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mHu2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu =
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu =
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHu2;

   beta_mHu2 = Re(twoLoop*(6*Power(g2,4)*Tr22 + 3.0983866769659336*g1*
      Tr31 - 6*tracemd2YdAdjYuYuAdjYd - 6*tracemq2AdjYdYdAdjYuYu - 6*
      tracemq2AdjYuYuAdjYdYd - 36*tracemq2AdjYuYuAdjYuYu - 6*
      tracemu2YuAdjYdYdAdjYu - 36*tracemu2YuAdjYuYuAdjYu - 6*
      traceYdAdjTYuTYuAdjYd - 6*traceYdAdjYuTYuAdjTYd - 6*mHd2*
      traceYdAdjYuYuAdjYd - 6*mHu2*traceYdAdjYuYuAdjYd - 6*
      traceYuAdjTYdTYdAdjYu - 36*traceYuAdjTYuTYuAdjYu - 6*
      traceYuAdjYdTYdAdjTYu - 36*traceYuAdjYuTYuAdjTYu - 36*mHu2*
      traceYuAdjYuYuAdjYu - 9*traceconjTYdTpTYd*AbsSqr(Lambdax) - 3*
      traceconjTYeTpTYe*AbsSqr(Lambdax) - 9*tracemd2YdAdjYd*AbsSqr(Lambdax) - 3
      *traceme2YeAdjYe*AbsSqr(Lambdax) - 3*traceml2AdjYeYe*AbsSqr(Lambdax) - 9*
      tracemq2AdjYdYd*AbsSqr(Lambdax) - 18*mHd2*traceYdAdjYd*AbsSqr(Lambdax) -
      9*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 9*mT2*traceYdAdjYd*AbsSqr(Lambdax)
      - 6*mHd2*traceYeAdjYe*AbsSqr(Lambdax) - 3*mHu2*traceYeAdjYe*AbsSqr(
      Lambdax) - 3*mT2*traceYeAdjYe*AbsSqr(Lambdax) - 9*traceYdAdjYd*AbsSqr(
      TLambdax) - 3*traceYeAdjYe*AbsSqr(TLambdax) - 30*AbsSqr(Lambdax)*AbsSqr(
      TLambdax) - 9*traceAdjYdTYd*Conj(TLambdax)*Lambdax - 3*traceAdjYeTYe*Conj
      (TLambdax)*Lambdax + 1.2*Tr2U111*Sqr(g1) + 1.6*traceconjTYuTpTYu*Sqr(g1)
      - 1.6*MassB*traceconjTYuTpYu*Sqr(g1) + 1.6*tracemq2AdjYuYu*Sqr(g1) + 1.6*
      tracemu2YuAdjYu*Sqr(g1) + 1.6*mHu2*traceYuAdjYu*Sqr(g1) + 12*mHd2*AbsSqr(
      Lambdax)*Sqr(g2) + 12*mHu2*AbsSqr(Lambdax)*Sqr(g2) + 12*mT2*AbsSqr(
      Lambdax)*Sqr(g2) + 12*AbsSqr(TLambdax)*Sqr(g2) - 12*MassWB*Conj(TLambdax)
      *Lambdax*Sqr(g2) + 0.04*Conj(MassB)*Sqr(g1)*(-40*traceAdjYuTYu + 80*MassB
      *traceYuAdjYu + 621*MassB*Sqr(g1) + 90*MassB*Sqr(g2) + 45*MassWB*Sqr(g2))
      + 32*traceconjTYuTpTYu*Sqr(g3) - 32*MassG*traceconjTYuTpYu*Sqr(g3) + 32*
      tracemq2AdjYuYu*Sqr(g3) + 32*tracemu2YuAdjYu*Sqr(g3) + 32*mHu2*
      traceYuAdjYu*Sqr(g3) + 64*traceYuAdjYu*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjYuTYu*Conj(MassG)*Sqr(g3) - 15*mHd2*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 15*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 15*mT2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + 0.6*Conj(MassWB)*Sqr(g2)*(3*(MassB + 2*MassWB)*
      Sqr(g1) + 115*MassWB*Sqr(g2) + 20*Conj(Lambdax)*(2*MassWB*Lambdax -
      TLambdax)) - 9*traceconjTYdTpYd*Conj(Lambdax)*TLambdax - 3*
      traceconjTYeTpYe*Conj(Lambdax)*TLambdax));


   return beta_mHu2;
}

/**
 * Calculates the three-loop beta function of mHu2.
 *
 * @return three-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mHu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
