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

// File generated at Tue 10 Oct 2017 21:35:13

#include "TMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mHd2.
 *
 * @return 1-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mHd2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHd2;

   beta_mHd2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 6*
      traceconjTYdTpTYd + 2*traceconjTYeTpTYe + 6*tracemd2YdAdjYd + 2*
      traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*tracemq2AdjYdYd + 6*mHd2*
      traceYdAdjYd + 2*mHd2*traceYeAdjYe + 3*mHd2*AbsSqr(Lambdax) + 3*mHu2*
      AbsSqr(Lambdax) + 3*mT2*AbsSqr(Lambdax) + 3*AbsSqr(TLambdax) - 1.2*AbsSqr
      (MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2)));


   return beta_mHd2;
}

/**
 * Calculates the 2-loop beta function of mHd2.
 *
 * @return 2-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mHd2_2_loop(const Soft_traces& soft_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHd2;

   beta_mHd2 = Re(twoLoop*(-3.0983866769659336*g1*Tr31 - 36*
      tracemd2YdAdjYdYdAdjYd - 6*tracemd2YdAdjYuYuAdjYd - 12*
      traceme2YeAdjYeYeAdjYe - 12*traceml2AdjYeYeAdjYeYe - 36*
      tracemq2AdjYdYdAdjYdYd - 6*tracemq2AdjYdYdAdjYuYu - 6*
      tracemq2AdjYuYuAdjYdYd - 6*tracemu2YuAdjYdYdAdjYu - 36*
      traceYdAdjTYdTYdAdjYd - 6*traceYdAdjTYuTYuAdjYd - 36*
      traceYdAdjYdTYdAdjTYd - 36*mHd2*traceYdAdjYdYdAdjYd - 6*
      traceYdAdjYuTYuAdjTYd - 6*mHd2*traceYdAdjYuYuAdjYd - 6*mHu2*
      traceYdAdjYuYuAdjYd - 12*traceYeAdjTYeTYeAdjYe - 12*traceYeAdjYeTYeAdjTYe
      - 12*mHd2*traceYeAdjYeYeAdjYe - 6*traceYuAdjTYdTYdAdjYu - 6*
      traceYuAdjYdTYdAdjTYu - 9*traceconjTYuTpTYu*AbsSqr(Lambdax) - 9*
      tracemq2AdjYuYu*AbsSqr(Lambdax) - 9*tracemu2YuAdjYu*AbsSqr(Lambdax) - 9*
      mHd2*traceYuAdjYu*AbsSqr(Lambdax) - 18*mHu2*traceYuAdjYu*AbsSqr(Lambdax)
      - 9*mT2*traceYuAdjYu*AbsSqr(Lambdax) - 9*traceYuAdjYu*AbsSqr(TLambdax) -
      30*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 9*traceAdjYuTYu*Conj(TLambdax)*
      Lambdax + 6*Tr22*Quad(g2) + 1.2*Tr2U111*Sqr(g1) - 0.8*traceconjTYdTpTYd*
      Sqr(g1) + 0.8*MassB*traceconjTYdTpYd*Sqr(g1) + 2.4*traceconjTYeTpTYe*Sqr(
      g1) - 2.4*MassB*traceconjTYeTpYe*Sqr(g1) - 0.8*tracemd2YdAdjYd*Sqr(g1) +
      2.4*traceme2YeAdjYe*Sqr(g1) + 2.4*traceml2AdjYeYe*Sqr(g1) - 0.8*
      tracemq2AdjYdYd*Sqr(g1) - 0.8*mHd2*traceYdAdjYd*Sqr(g1) + 2.4*mHd2*
      traceYeAdjYe*Sqr(g1) + 12*mHd2*AbsSqr(Lambdax)*Sqr(g2) + 12*mHu2*AbsSqr(
      Lambdax)*Sqr(g2) + 12*mT2*AbsSqr(Lambdax)*Sqr(g2) + 12*AbsSqr(TLambdax)*
      Sqr(g2) - 12*MassWB*Conj(TLambdax)*Lambdax*Sqr(g2) + 0.04*Conj(MassB)*Sqr
      (g1)*(621*MassB*Sqr(g1) + 5*(4*(traceAdjYdTYd - 3*traceAdjYeTYe - 2*MassB
      *traceYdAdjYd + 6*MassB*traceYeAdjYe) + 9*(2*MassB + MassWB)*Sqr(g2))) +
      32*traceconjTYdTpTYd*Sqr(g3) - 32*MassG*traceconjTYdTpYd*Sqr(g3) + 32*
      tracemd2YdAdjYd*Sqr(g3) + 32*tracemq2AdjYdYd*Sqr(g3) + 32*mHd2*
      traceYdAdjYd*Sqr(g3) + 64*traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjYdTYd*Conj(MassG)*Sqr(g3) - 15*mHd2*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 15*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 15*mT2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + 0.6*Conj(MassWB)*Sqr(g2)*(3*(MassB + 2*MassWB)*
      Sqr(g1) + 115*MassWB*Sqr(g2) + 20*Conj(Lambdax)*(2*MassWB*Lambdax -
      TLambdax)) - 9*traceconjTYuTpYu*Conj(Lambdax)*TLambdax));


   return beta_mHd2;
}

/**
 * Calculates the 3-loop beta function of mHd2.
 *
 * @return 3-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mHd2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy