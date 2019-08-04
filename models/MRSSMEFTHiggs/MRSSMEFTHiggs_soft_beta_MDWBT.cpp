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

// File generated at Sun 4 Aug 2019 17:36:13

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
 * Calculates the 1-loop beta function of MDWBT.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDWBT_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MDWBT;

   beta_MDWBT = Re(MDWBT*oneOver16PiSqr*(AbsSqr(LamTD) + AbsSqr(LamTU)));


   return beta_MDWBT;
}

/**
 * Calculates the 2-loop beta function of MDWBT.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDWBT_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MDWBT;

   beta_MDWBT = Re(0.2*MDWBT*twoLoop*(-15*traceYdAdjYd*AbsSqr(LamTD) - 5*
      traceYeAdjYe*AbsSqr(LamTD) - 10*AbsSqr(LamSD)*AbsSqr(LamTD) - 15*
      traceYuAdjYu*AbsSqr(LamTU) - 10*AbsSqr(LamSU)*AbsSqr(LamTU) + 440*Quad(g2
      ) + 3*AbsSqr(LamTD)*Sqr(g1) + 3*AbsSqr(LamTU)*Sqr(g1) - 30*traceYdAdjYd*
      Sqr(g2) - 10*traceYeAdjYe*Sqr(g2) - 30*traceYuAdjYu*Sqr(g2) - 10*AbsSqr(
      LamSD)*Sqr(g2) - 10*AbsSqr(LamSU)*Sqr(g2) - 40*AbsSqr(LamTD)*Sqr(g2) - 40
      *AbsSqr(LamTU)*Sqr(g2) + 12*Sqr(g1)*Sqr(g2) + 120*Sqr(g2)*Sqr(g3) - 15*
      Sqr(LamTD)*Sqr(Conj(LamTD)) - 15*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_MDWBT;
}

/**
 * Calculates the 3-loop beta function of MDWBT.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDWBT_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDWBT;

   beta_MDWBT = 0;


   return beta_MDWBT;
}

/**
 * Calculates the 4-loop beta function of MDWBT.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDWBT_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDWBT;

   beta_MDWBT = 0;


   return beta_MDWBT;
}

/**
 * Calculates the 5-loop beta function of MDWBT.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_MDWBT_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDWBT;

   beta_MDWBT = 0;


   return beta_MDWBT;
}

} // namespace flexiblesusy
