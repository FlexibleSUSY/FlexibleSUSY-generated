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
 * Calculates the 1-loop beta function of mRu2.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mRu2;

   beta_mRu2 = Re(0.2*(-3.872983346207417*g1*Tr11 + 10*mHu2*AbsSqr(LamSU) + 10*
      mRu2*AbsSqr(LamSU) + 10*mS2*AbsSqr(LamSU) + 15*mHu2*AbsSqr(LamTU) + 15*
      mRu2*AbsSqr(LamTU) + 15*mT2*AbsSqr(LamTU)));


   return oneLoop * beta_mRu2;
}

/**
 * Calculates the 2-loop beta function of mRu2.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mRu2;

   beta_mRu2 = Re(0.2*(-15.491933384829668*g1*Tr31 - 30*tracemq2AdjYuYu*AbsSqr(
      LamSU) - 30*tracemu2YuAdjYu*AbsSqr(LamSU) - 60*mHu2*traceYuAdjYu*AbsSqr(
      LamSU) - 30*mRu2*traceYuAdjYu*AbsSqr(LamSU) - 30*mS2*traceYuAdjYu*AbsSqr(
      LamSU) - 20*mHd2*AbsSqr(LamSD)*AbsSqr(LamSU) - 20*mHu2*AbsSqr(LamSD)*
      AbsSqr(LamSU) - 20*mRd2*AbsSqr(LamSD)*AbsSqr(LamSU) - 20*mRu2*AbsSqr(
      LamSD)*AbsSqr(LamSU) - 40*mS2*AbsSqr(LamSD)*AbsSqr(LamSU) - 45*
      tracemq2AdjYuYu*AbsSqr(LamTU) - 45*tracemu2YuAdjYu*AbsSqr(LamTU) - 90*
      mHu2*traceYuAdjYu*AbsSqr(LamTU) - 45*mRu2*traceYuAdjYu*AbsSqr(LamTU) - 45
      *mT2*traceYuAdjYu*AbsSqr(LamTU) - 60*mHu2*AbsSqr(LamSU)*AbsSqr(LamTU) -
      60*mRu2*AbsSqr(LamSU)*AbsSqr(LamTU) - 30*mS2*AbsSqr(LamSU)*AbsSqr(LamTU)
      - 30*mT2*AbsSqr(LamSU)*AbsSqr(LamTU) - 15*mHd2*AbsSqr(LamTD)*AbsSqr(LamTU
      ) - 15*mHu2*AbsSqr(LamTD)*AbsSqr(LamTU) - 15*mRd2*AbsSqr(LamTD)*AbsSqr(
      LamTU) - 15*mRu2*AbsSqr(LamTD)*AbsSqr(LamTU) - 30*mT2*AbsSqr(LamTD)*
      AbsSqr(LamTU) + 30*Tr22*Quad(g2) + 6*Tr2U111*Sqr(g1) + 60*mHu2*AbsSqr(
      LamTU)*Sqr(g2) + 60*mRu2*AbsSqr(LamTU)*Sqr(g2) + 60*mT2*AbsSqr(LamTU)*Sqr
      (g2) - 60*mHu2*Sqr(LamSU)*Sqr(Conj(LamSU)) - 60*mRu2*Sqr(LamSU)*Sqr(Conj(
      LamSU)) - 60*mS2*Sqr(LamSU)*Sqr(Conj(LamSU)) - 75*mHu2*Sqr(LamTU)*Sqr(
      Conj(LamTU)) - 75*mRu2*Sqr(LamTU)*Sqr(Conj(LamTU)) - 75*mT2*Sqr(LamTU)*
      Sqr(Conj(LamTU))));


   return twoLoop * beta_mRu2;
}

/**
 * Calculates the 3-loop beta function of mRu2.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRu2;

   beta_mRu2 = 0;


   return threeLoop * beta_mRu2;
}

/**
 * Calculates the 4-loop beta function of mRu2.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRu2;

   beta_mRu2 = 0;


   return fourLoop * beta_mRu2;
}

/**
 * Calculates the 5-loop beta function of mRu2.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRu2;

   beta_mRu2 = 0;


   return fiveLoop * beta_mRu2;
}

} // namespace flexiblesusy
