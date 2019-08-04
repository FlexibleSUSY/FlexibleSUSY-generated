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

// File generated at Sun 4 Aug 2019 17:36:09

#include "MRSSMEFTHiggs_soft_parameters.hpp"
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
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_1_loop(const Soft_traces& soft_traces) const
{


   double beta_mT2;

   beta_mT2 = Re(2*oneOver16PiSqr*(mHd2*AbsSqr(LamTD) + mRd2*AbsSqr(LamTD) +
      mT2*AbsSqr(LamTD) + mHu2*AbsSqr(LamTU) + mRu2*AbsSqr(LamTU) + mT2*AbsSqr(
      LamTU)));


   return beta_mT2;
}

/**
 * Calculates the 2-loop beta function of mT2.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_2_loop(const Soft_traces& soft_traces) const
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
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mT2;

   beta_mT2 = Re(0.4*twoLoop*(-15*tracemd2YdAdjYd*AbsSqr(LamTD) - 5*
      traceme2YeAdjYe*AbsSqr(LamTD) - 5*traceml2AdjYeYe*AbsSqr(LamTD) - 15*
      tracemq2AdjYdYd*AbsSqr(LamTD) - 30*mHd2*traceYdAdjYd*AbsSqr(LamTD) - 15*
      mRd2*traceYdAdjYd*AbsSqr(LamTD) - 15*mT2*traceYdAdjYd*AbsSqr(LamTD) - 10*
      mHd2*traceYeAdjYe*AbsSqr(LamTD) - 5*mRd2*traceYeAdjYe*AbsSqr(LamTD) - 5*
      mT2*traceYeAdjYe*AbsSqr(LamTD) - 20*mHd2*AbsSqr(LamSD)*AbsSqr(LamTD) - 20
      *mRd2*AbsSqr(LamSD)*AbsSqr(LamTD) - 10*mS2*AbsSqr(LamSD)*AbsSqr(LamTD) -
      10*mT2*AbsSqr(LamSD)*AbsSqr(LamTD) - 15*tracemq2AdjYuYu*AbsSqr(LamTU) -
      15*tracemu2YuAdjYu*AbsSqr(LamTU) - 30*mHu2*traceYuAdjYu*AbsSqr(LamTU) -
      15*mRu2*traceYuAdjYu*AbsSqr(LamTU) - 15*mT2*traceYuAdjYu*AbsSqr(LamTU) -
      20*mHu2*AbsSqr(LamSU)*AbsSqr(LamTU) - 20*mRu2*AbsSqr(LamSU)*AbsSqr(LamTU)
      - 10*mS2*AbsSqr(LamSU)*AbsSqr(LamTU) - 10*mT2*AbsSqr(LamSU)*AbsSqr(LamTU)
      + 40*Tr22*Quad(g2) + 3*mHd2*AbsSqr(LamTD)*Sqr(g1) + 3*mRd2*AbsSqr(LamTD)*
      Sqr(g1) + 3*mT2*AbsSqr(LamTD)*Sqr(g1) + 3*mHu2*AbsSqr(LamTU)*Sqr(g1) + 3*
      mRu2*AbsSqr(LamTU)*Sqr(g1) + 3*mT2*AbsSqr(LamTU)*Sqr(g1) - 5*mHd2*AbsSqr(
      LamTD)*Sqr(g2) - 5*mRd2*AbsSqr(LamTD)*Sqr(g2) - 5*mT2*AbsSqr(LamTD)*Sqr(
      g2) - 5*mHu2*AbsSqr(LamTU)*Sqr(g2) - 5*mRu2*AbsSqr(LamTU)*Sqr(g2) - 5*mT2
      *AbsSqr(LamTU)*Sqr(g2) - 30*mHd2*Sqr(LamTD)*Sqr(Conj(LamTD)) - 30*mRd2*
      Sqr(LamTD)*Sqr(Conj(LamTD)) - 30*mT2*Sqr(LamTD)*Sqr(Conj(LamTD)) - 30*
      mHu2*Sqr(LamTU)*Sqr(Conj(LamTU)) - 30*mRu2*Sqr(LamTU)*Sqr(Conj(LamTU)) -
      30*mT2*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_mT2;
}

/**
 * Calculates the 3-loop beta function of mT2.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_3_loop(const Soft_traces& soft_traces) const
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
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

/**
 * Calculates the 5-loop beta function of mT2.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

} // namespace flexiblesusy
