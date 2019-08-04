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

// File generated at Sun 4 Aug 2019 19:25:18

#include "MRSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mS2.
 *
 * @return 1-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_1_loop(const Soft_traces& soft_traces) const
{


   double beta_mS2;

   beta_mS2 = Re(4*oneOver16PiSqr*(mHd2*AbsSqr(LamSD) + mRd2*AbsSqr(LamSD) +
      mS2*AbsSqr(LamSD) + mHu2*AbsSqr(LamSU) + mRu2*AbsSqr(LamSU) + mS2*AbsSqr(
      LamSU)));


   return beta_mS2;
}

/**
 * Calculates the 2-loop beta function of mS2.
 *
 * @return 2-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;


   double beta_mS2;

   beta_mS2 = Re(0.8*twoLoop*(-15*tracemd2YdAdjYd*AbsSqr(LamSD) - 5*
      traceme2YeAdjYe*AbsSqr(LamSD) - 5*traceml2AdjYeYe*AbsSqr(LamSD) - 15*
      tracemq2AdjYdYd*AbsSqr(LamSD) - 30*mHd2*traceYdAdjYd*AbsSqr(LamSD) - 15*
      mRd2*traceYdAdjYd*AbsSqr(LamSD) - 15*mS2*traceYdAdjYd*AbsSqr(LamSD) - 10*
      mHd2*traceYeAdjYe*AbsSqr(LamSD) - 5*mRd2*traceYeAdjYe*AbsSqr(LamSD) - 5*
      mS2*traceYeAdjYe*AbsSqr(LamSD) - 15*tracemq2AdjYuYu*AbsSqr(LamSU) - 15*
      tracemu2YuAdjYu*AbsSqr(LamSU) - 30*mHu2*traceYuAdjYu*AbsSqr(LamSU) - 15*
      mRu2*traceYuAdjYu*AbsSqr(LamSU) - 15*mS2*traceYuAdjYu*AbsSqr(LamSU) - 30*
      mHd2*AbsSqr(LamSD)*AbsSqr(LamTD) - 30*mRd2*AbsSqr(LamSD)*AbsSqr(LamTD) -
      15*mS2*AbsSqr(LamSD)*AbsSqr(LamTD) - 15*mT2*AbsSqr(LamSD)*AbsSqr(LamTD) -
      30*mHu2*AbsSqr(LamSU)*AbsSqr(LamTU) - 30*mRu2*AbsSqr(LamSU)*AbsSqr(LamTU)
      - 15*mS2*AbsSqr(LamSU)*AbsSqr(LamTU) - 15*mT2*AbsSqr(LamSU)*AbsSqr(LamTU)
      + 3*mHd2*AbsSqr(LamSD)*Sqr(g1) + 3*mRd2*AbsSqr(LamSD)*Sqr(g1) + 3*mS2*
      AbsSqr(LamSD)*Sqr(g1) + 3*mHu2*AbsSqr(LamSU)*Sqr(g1) + 3*mRu2*AbsSqr(
      LamSU)*Sqr(g1) + 3*mS2*AbsSqr(LamSU)*Sqr(g1) + 15*mHd2*AbsSqr(LamSD)*Sqr(
      g2) + 15*mRd2*AbsSqr(LamSD)*Sqr(g2) + 15*mS2*AbsSqr(LamSD)*Sqr(g2) + 15*
      mHu2*AbsSqr(LamSU)*Sqr(g2) + 15*mRu2*AbsSqr(LamSU)*Sqr(g2) + 15*mS2*
      AbsSqr(LamSU)*Sqr(g2) - 20*mHd2*Sqr(LamSD)*Sqr(Conj(LamSD)) - 20*mRd2*Sqr
      (LamSD)*Sqr(Conj(LamSD)) - 20*mS2*Sqr(LamSD)*Sqr(Conj(LamSD)) - 20*mHu2*
      Sqr(LamSU)*Sqr(Conj(LamSU)) - 20*mRu2*Sqr(LamSU)*Sqr(Conj(LamSU)) - 20*
      mS2*Sqr(LamSU)*Sqr(Conj(LamSU))));


   return beta_mS2;
}

/**
 * Calculates the 3-loop beta function of mS2.
 *
 * @return 3-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mS2;

   beta_mS2 = 0;


   return beta_mS2;
}

/**
 * Calculates the 4-loop beta function of mS2.
 *
 * @return 4-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mS2;

   beta_mS2 = 0;


   return beta_mS2;
}

/**
 * Calculates the 5-loop beta function of mS2.
 *
 * @return 5-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mS2;

   beta_mS2 = 0;


   return beta_mS2;
}

} // namespace flexiblesusy
