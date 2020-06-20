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
 * Calculates the 1-loop beta function of mRd2.
 *
 * @return 1-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mRd2;

   beta_mRd2 = Re(0.2*(3.872983346207417*g1*Tr11 + 10*mHd2*AbsSqr(LamSD) + 10*
      mRd2*AbsSqr(LamSD) + 10*mS2*AbsSqr(LamSD) + 15*mHd2*AbsSqr(LamTD) + 15*
      mRd2*AbsSqr(LamTD) + 15*mT2*AbsSqr(LamTD)));


   return oneLoop * beta_mRd2;
}

/**
 * Calculates the 2-loop beta function of mRd2.
 *
 * @return 2-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mRd2;

   beta_mRd2 = Re(0.2*(15.491933384829668*g1*Tr31 - 30*tracemd2YdAdjYd*AbsSqr(
      LamSD) - 10*traceme2YeAdjYe*AbsSqr(LamSD) - 10*traceml2AdjYeYe*AbsSqr(
      LamSD) - 30*tracemq2AdjYdYd*AbsSqr(LamSD) - 60*mHd2*traceYdAdjYd*AbsSqr(
      LamSD) - 30*mRd2*traceYdAdjYd*AbsSqr(LamSD) - 30*mS2*traceYdAdjYd*AbsSqr(
      LamSD) - 20*mHd2*traceYeAdjYe*AbsSqr(LamSD) - 10*mRd2*traceYeAdjYe*AbsSqr
      (LamSD) - 10*mS2*traceYeAdjYe*AbsSqr(LamSD) - 20*mHd2*AbsSqr(LamSD)*
      AbsSqr(LamSU) - 20*mHu2*AbsSqr(LamSD)*AbsSqr(LamSU) - 20*mRd2*AbsSqr(
      LamSD)*AbsSqr(LamSU) - 20*mRu2*AbsSqr(LamSD)*AbsSqr(LamSU) - 40*mS2*
      AbsSqr(LamSD)*AbsSqr(LamSU) - 45*tracemd2YdAdjYd*AbsSqr(LamTD) - 15*
      traceme2YeAdjYe*AbsSqr(LamTD) - 15*traceml2AdjYeYe*AbsSqr(LamTD) - 45*
      tracemq2AdjYdYd*AbsSqr(LamTD) - 90*mHd2*traceYdAdjYd*AbsSqr(LamTD) - 45*
      mRd2*traceYdAdjYd*AbsSqr(LamTD) - 45*mT2*traceYdAdjYd*AbsSqr(LamTD) - 30*
      mHd2*traceYeAdjYe*AbsSqr(LamTD) - 15*mRd2*traceYeAdjYe*AbsSqr(LamTD) - 15
      *mT2*traceYeAdjYe*AbsSqr(LamTD) - 60*mHd2*AbsSqr(LamSD)*AbsSqr(LamTD) -
      60*mRd2*AbsSqr(LamSD)*AbsSqr(LamTD) - 30*mS2*AbsSqr(LamSD)*AbsSqr(LamTD)
      - 30*mT2*AbsSqr(LamSD)*AbsSqr(LamTD) - 15*mHd2*AbsSqr(LamTD)*AbsSqr(LamTU
      ) - 15*mHu2*AbsSqr(LamTD)*AbsSqr(LamTU) - 15*mRd2*AbsSqr(LamTD)*AbsSqr(
      LamTU) - 15*mRu2*AbsSqr(LamTD)*AbsSqr(LamTU) - 30*mT2*AbsSqr(LamTD)*
      AbsSqr(LamTU) + 30*Tr22*Quad(g2) + 6*Tr2U111*Sqr(g1) + 60*mHd2*AbsSqr(
      LamTD)*Sqr(g2) + 60*mRd2*AbsSqr(LamTD)*Sqr(g2) + 60*mT2*AbsSqr(LamTD)*Sqr
      (g2) - 60*mHd2*Sqr(LamSD)*Sqr(Conj(LamSD)) - 60*mRd2*Sqr(LamSD)*Sqr(Conj(
      LamSD)) - 60*mS2*Sqr(LamSD)*Sqr(Conj(LamSD)) - 75*mHd2*Sqr(LamTD)*Sqr(
      Conj(LamTD)) - 75*mRd2*Sqr(LamTD)*Sqr(Conj(LamTD)) - 75*mT2*Sqr(LamTD)*
      Sqr(Conj(LamTD))));


   return twoLoop * beta_mRd2;
}

/**
 * Calculates the 3-loop beta function of mRd2.
 *
 * @return 3-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRd2;

   beta_mRd2 = 0;


   return threeLoop * beta_mRd2;
}

/**
 * Calculates the 4-loop beta function of mRd2.
 *
 * @return 4-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRd2;

   beta_mRd2 = 0;


   return fourLoop * beta_mRd2;
}

/**
 * Calculates the 5-loop beta function of mRd2.
 *
 * @return 5-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRd2;

   beta_mRd2 = 0;


   return fiveLoop * beta_mRd2;
}

} // namespace flexiblesusy
