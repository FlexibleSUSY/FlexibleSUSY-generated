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
 * Calculates the 1-loop beta function of mHd2.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mHd2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHd2;

   beta_mHd2 = Re(0.2*(-3.872983346207417*g1*Tr11 + 30*tracemd2YdAdjYd + 10*
      traceme2YeAdjYe + 10*traceml2AdjYeYe + 30*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*AbsSqr(LamSD) + 10*mRd2*
      AbsSqr(LamSD) + 10*mS2*AbsSqr(LamSD) + 15*mHd2*AbsSqr(LamTD) + 15*mRd2*
      AbsSqr(LamTD) + 15*mT2*AbsSqr(LamTD)));


   return oneLoop * beta_mHd2;
}

/**
 * Calculates the 2-loop beta function of mHd2.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mHd2_2_loop(const Soft_traces& soft_traces) const
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
   const double tracemd2YdAdjYdYdAdjYd = TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd = TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe = TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe = TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd = TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu = TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd = TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu = TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHd2;

   beta_mHd2 = Re(0.2*(-15.491933384829668*g1*Tr31 - 180*tracemd2YdAdjYdYdAdjYd
       - 30*tracemd2YdAdjYuYuAdjYd - 60*traceme2YeAdjYeYeAdjYe - 60*
      traceml2AdjYeYeAdjYeYe - 180*tracemq2AdjYdYdAdjYdYd - 30*
      tracemq2AdjYdYdAdjYuYu - 30*tracemq2AdjYuYuAdjYdYd - 30*
      tracemu2YuAdjYdYdAdjYu - 180*mHd2*traceYdAdjYdYdAdjYd - 30*mHd2*
      traceYdAdjYuYuAdjYd - 30*mHu2*traceYdAdjYuYuAdjYd - 60*mHd2*
      traceYeAdjYeYeAdjYe - 20*mHd2*AbsSqr(LamSD)*AbsSqr(LamSU) - 20*mHu2*
      AbsSqr(LamSD)*AbsSqr(LamSU) - 20*mRd2*AbsSqr(LamSD)*AbsSqr(LamSU) - 20*
      mRu2*AbsSqr(LamSD)*AbsSqr(LamSU) - 40*mS2*AbsSqr(LamSD)*AbsSqr(LamSU) -
      60*mHd2*AbsSqr(LamSD)*AbsSqr(LamTD) - 60*mRd2*AbsSqr(LamSD)*AbsSqr(LamTD)
      - 30*mS2*AbsSqr(LamSD)*AbsSqr(LamTD) - 30*mT2*AbsSqr(LamSD)*AbsSqr(LamTD)
      - 15*mHd2*AbsSqr(LamTD)*AbsSqr(LamTU) - 15*mHu2*AbsSqr(LamTD)*AbsSqr(
      LamTU) - 15*mRd2*AbsSqr(LamTD)*AbsSqr(LamTU) - 15*mRu2*AbsSqr(LamTD)*
      AbsSqr(LamTU) - 30*mT2*AbsSqr(LamTD)*AbsSqr(LamTU) + 30*Tr22*Quad(g2) + 6
      *Tr2U111*Sqr(g1) - 4*tracemd2YdAdjYd*Sqr(g1) + 12*traceme2YeAdjYe*Sqr(g1)
      + 12*traceml2AdjYeYe*Sqr(g1) - 4*tracemq2AdjYdYd*Sqr(g1) - 4*mHd2*
      traceYdAdjYd*Sqr(g1) + 12*mHd2*traceYeAdjYe*Sqr(g1) + 60*mHd2*AbsSqr(
      LamTD)*Sqr(g2) + 60*mRd2*AbsSqr(LamTD)*Sqr(g2) + 60*mT2*AbsSqr(LamTD)*Sqr
      (g2) + 160*tracemd2YdAdjYd*Sqr(g3) + 160*tracemq2AdjYdYd*Sqr(g3) + 160*
      mHd2*traceYdAdjYd*Sqr(g3) - 60*mHd2*Sqr(LamSD)*Sqr(Conj(LamSD)) - 60*mRd2
      *Sqr(LamSD)*Sqr(Conj(LamSD)) - 60*mS2*Sqr(LamSD)*Sqr(Conj(LamSD)) - 75*
      mHd2*Sqr(LamTD)*Sqr(Conj(LamTD)) - 75*mRd2*Sqr(LamTD)*Sqr(Conj(LamTD)) -
      75*mT2*Sqr(LamTD)*Sqr(Conj(LamTD))));


   return twoLoop * beta_mHd2;
}

/**
 * Calculates the 3-loop beta function of mHd2.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mHd2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return threeLoop * beta_mHd2;
}

/**
 * Calculates the 4-loop beta function of mHd2.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mHd2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return fourLoop * beta_mHd2;
}

/**
 * Calculates the 5-loop beta function of mHd2.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mHd2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return fiveLoop * beta_mHd2;
}

} // namespace flexiblesusy
