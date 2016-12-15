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

// File generated at Thu 15 Dec 2016 12:39:30

#include "MRSSMtower_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of BMuU.
 *
 * @return one-loop beta function
 */
double MRSSMtower_soft_parameters::calc_beta_BMuU_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_BMuU;

   beta_BMuU = Re(oneOver16PiSqr*(4*LamSU*BMuD*Conj(LamSD) + BMuU*(3*
      traceYuAdjYu + 6*AbsSqr(LamSU) + 3*AbsSqr(LamTU) - 0.6*Sqr(g1) - 3*Sqr(g2
      ))));


   return beta_BMuU;
}

/**
 * Calculates the two-loop beta function of BMuU.
 *
 * @return two-loop beta function
 */
double MRSSMtower_soft_parameters::calc_beta_BMuU_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_BMuU;

   beta_BMuU = Re(0.1*twoLoop*(-8*LamSU*BMuD*Conj(LamSD)*(15*traceYdAdjYd
      + 5*traceYeAdjYe + 10*AbsSqr(LamSD) + 15*AbsSqr(LamTD) - 9*Sqr(g1) - 45*
      Sqr(g2)) + BMuU*(45*Power(g1,4) + 165*Power(g2,4) - 30*
      traceYdAdjYuYuAdjYd - 90*traceYuAdjYuYuAdjYu + 8*traceYuAdjYu*Sqr(g1) +
      18*Sqr(g1)*Sqr(g2) + 15*AbsSqr(LamTU)*(-3*traceYuAdjYu - 2*AbsSqr(LamTD)
      + 8*Sqr(g2)) + 2*AbsSqr(LamSU)*(-75*traceYuAdjYu - 20*AbsSqr(LamSD) - 90*
      AbsSqr(LamTU) + 36*Sqr(g1) + 180*Sqr(g2)) + 160*traceYuAdjYu*Sqr(g3) -
      140*Sqr(LamSU)*Sqr(Conj(LamSU)) - 75*Sqr(LamTU)*Sqr(Conj(LamTU)))));


   return beta_BMuU;
}

/**
 * Calculates the three-loop beta function of BMuU.
 *
 * @return three-loop beta function
 */
double MRSSMtower_soft_parameters::calc_beta_BMuU_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuU;

   beta_BMuU = 0;


   return beta_BMuU;
}

} // namespace flexiblesusy
