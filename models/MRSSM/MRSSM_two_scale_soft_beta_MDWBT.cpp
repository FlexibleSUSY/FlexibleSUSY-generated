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

// File generated at Fri 8 Jan 2016 12:05:10

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
 * Calculates the one-loop beta function of MDWBT.
 *
 * @return one-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_MDWBT_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MDWBT;

   beta_MDWBT = Re(MDWBT*oneOver16PiSqr*(AbsSqr(LamTD) + AbsSqr(LamTU)));


   return beta_MDWBT;
}

/**
 * Calculates the two-loop beta function of MDWBT.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_MDWBT_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MDWBT;

   beta_MDWBT = Re(0.2*MDWBT*twoLoop*(440*Power(g2,4) - 15*traceYdAdjYd*
      AbsSqr(LamTD) - 5*traceYeAdjYe*AbsSqr(LamTD) - 15*traceYuAdjYu*AbsSqr(
      LamTU) + 3*AbsSqr(LamTD)*Sqr(g1) + 3*AbsSqr(LamTU)*Sqr(g1) - 30*
      traceYdAdjYd*Sqr(g2) - 10*traceYeAdjYe*Sqr(g2) - 30*traceYuAdjYu*Sqr(g2)
      - 40*AbsSqr(LamTD)*Sqr(g2) - 40*AbsSqr(LamTU)*Sqr(g2) + 12*Sqr(g1)*Sqr(g2
      ) - 10*AbsSqr(LamSD)*(AbsSqr(LamTD) + Sqr(g2)) - 10*AbsSqr(LamSU)*(AbsSqr
      (LamTU) + Sqr(g2)) + 120*Sqr(g2)*Sqr(g3) - 15*Sqr(LamTD)*Sqr(Conj(LamTD))
      - 15*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_MDWBT;
}

/**
 * Calculates the three-loop beta function of MDWBT.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_MDWBT_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MDWBT;

   beta_MDWBT = 0;


   return beta_MDWBT;
}

} // namespace flexiblesusy
