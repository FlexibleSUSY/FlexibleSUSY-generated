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

// File generated at Tue 8 Mar 2016 16:13:43

#include "MRSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHd2;

   beta_mHd2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 6*
      tracemd2YdAdjYd + 2*traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*
      tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 2*mHd2*traceYeAdjYe + 2*(mHd2 +
      mRd2 + mS2)*AbsSqr(LamSD) + 3*(mHd2 + mRd2 + mT2)*AbsSqr(LamTD)));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
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

   beta_mHd2 = Re(twoLoop*(-2*AbsSqr(LamSD)*(2*(mHd2 + mHu2 + mRd2 + mRu2
      + 2*mS2)*AbsSqr(LamSU) + 3*(2*mHd2 + 2*mRd2 + mS2 + mT2)*AbsSqr(LamTD))
      + 3*AbsSqr(LamTD)*(-((mHd2 + mHu2 + mRd2 + mRu2 + 2*mT2)*AbsSqr(LamTU)) +
      4*(mHd2 + mRd2 + mT2)*Sqr(g2)) + 0.4*(15*Power(g2,4)*Tr22 -
      7.745966692414834*g1*Tr31 - 90*tracemd2YdAdjYdYdAdjYd - 15*
      tracemd2YdAdjYuYuAdjYd - 30*traceme2YeAdjYeYeAdjYe - 30*
      traceml2AdjYeYeAdjYeYe - 90*tracemq2AdjYdYdAdjYdYd - 15*
      tracemq2AdjYdYdAdjYuYu - 15*tracemq2AdjYuYuAdjYdYd - 15*
      tracemu2YuAdjYdYdAdjYu - 90*mHd2*traceYdAdjYdYdAdjYd - 15*mHd2*
      traceYdAdjYuYuAdjYd - 15*mHu2*traceYdAdjYuYuAdjYd - 30*mHd2*
      traceYeAdjYeYeAdjYe + 3*Tr2U111*Sqr(g1) - 2*tracemd2YdAdjYd*Sqr(g1) + 6*
      traceme2YeAdjYe*Sqr(g1) + 6*traceml2AdjYeYe*Sqr(g1) - 2*tracemq2AdjYdYd*
      Sqr(g1) + 6*mHd2*traceYeAdjYe*Sqr(g1) - 2*mHd2*traceYdAdjYd*(Sqr(g1) - 40
      *Sqr(g3)) + 80*tracemd2YdAdjYd*Sqr(g3) + 80*tracemq2AdjYdYd*Sqr(g3)) - 12
      *(mHd2 + mRd2 + mS2)*Sqr(LamSD)*Sqr(Conj(LamSD)) - 15*(mHd2 + mRd2 + mT2)
      *Sqr(LamTD)*Sqr(Conj(LamTD))));


   return beta_mHd2;
}

/**
 * Calculates the three-loop beta function of mHd2.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mHd2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
