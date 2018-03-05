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

// File generated at Mon 5 Mar 2018 17:55:31

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
 * Calculates the 1-loop beta function of mT2.
 *
 * @return 1-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mT2_1_loop(const Soft_traces& soft_traces) const
{


   double beta_mT2;

   beta_mT2 = Re(2*oneOver16PiSqr*((mHd2 + mHu2 + mT2)*AbsSqr(Lambdax) +
      AbsSqr(TLambdax) - 8*AbsSqr(MassWB)*Sqr(g2)));


   return beta_mT2;
}

/**
 * Calculates the 2-loop beta function of mT2.
 *
 * @return 2-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mT2_2_loop(const Soft_traces& soft_traces) const
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
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mT2;

   beta_mT2 = Re(0.4*twoLoop*(40*Tr22*Quad(g2) - 30*(mHd2 + mHu2 + mT2)*
      Sqr(Conj(Lambdax))*Sqr(Lambdax) + Conj(Lambdax)*(Lambdax*(3*(mHd2 + mHu2
      + mT2)*Sqr(g1) - 5*(3*traceconjTYdTpTYd + traceconjTYeTpTYe + 3*
      traceconjTYuTpTYu + 3*tracemd2YdAdjYd + traceme2YeAdjYe + traceml2AdjYeYe
      + 3*tracemq2AdjYdYd + 3*tracemq2AdjYuYu + 3*tracemu2YuAdjYu + 6*mHd2*
      traceYdAdjYd + 3*mHu2*traceYdAdjYd + 3*mT2*traceYdAdjYd + 2*mHd2*
      traceYeAdjYe + mHu2*traceYeAdjYe + mT2*traceYeAdjYe + 3*mHd2*traceYuAdjYu
      + 6*mHu2*traceYuAdjYu + 3*mT2*traceYuAdjYu + (mHd2 + mHu2 + mT2)*Sqr(g2)
      )) + 3*Conj(MassB)*Sqr(g1)*(2*MassB*Lambdax - TLambdax) - 5*(3*
      traceconjTYdTpYd + traceconjTYeTpYe + 3*traceconjTYuTpYu + 12*Conj(
      TLambdax)*Lambdax)*TLambdax) - Conj(TLambdax)*(Lambdax*(3*MassB*Sqr(g1) +
      5*(3*traceAdjYdTYd + traceAdjYeTYe + 3*traceAdjYuTYu - MassWB*Sqr(g2)))
      + (-3*Sqr(g1) + 5*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + Sqr(
      g2)))*TLambdax) + 5*Conj(MassWB)*Sqr(g2)*(152*MassWB*Sqr(g2) + Conj(
      Lambdax)*(-2*MassWB*Lambdax + TLambdax))));


   return beta_mT2;
}

/**
 * Calculates the 3-loop beta function of mT2.
 *
 * @return 3-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mT2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

/**
 * Calculates the 4-loop beta function of mT2.
 *
 * @return 4-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mT2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

} // namespace flexiblesusy
