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

// File generated at Tue 10 Oct 2017 21:35:14

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
 * Calculates the 1-loop beta function of BMuD.
 *
 * @return 1-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_BMuD_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_BMuD;

   beta_BMuD = Re(oneOver16PiSqr*(4*LamSD*BMuU*Conj(LamSU) + BMuD*(3*
      traceYdAdjYd + traceYeAdjYe + 6*AbsSqr(LamSD) + 3*AbsSqr(LamTD) - 0.6*Sqr
      (g1) - 3*Sqr(g2))));


   return beta_BMuD;
}

/**
 * Calculates the 2-loop beta function of BMuD.
 *
 * @return 2-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_BMuD_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_BMuD;

   beta_BMuD = Re(0.1*twoLoop*(-8*LamSD*BMuU*Conj(LamSU)*(15*traceYuAdjYu
      + 10*AbsSqr(LamSU) + 15*AbsSqr(LamTU) - 9*Sqr(g1) - 45*Sqr(g2)) + BMuD*(
      -90*traceYdAdjYdYdAdjYd - 30*traceYdAdjYuYuAdjYd - 30*traceYeAdjYeYeAdjYe
      + 45*Quad(g1) + 165*Quad(g2) - 4*traceYdAdjYd*Sqr(g1) + 12*traceYeAdjYe*
      Sqr(g1) + 18*Sqr(g1)*Sqr(g2) + 15*AbsSqr(LamTD)*(-3*traceYdAdjYd -
      traceYeAdjYe - 2*AbsSqr(LamTU) + 8*Sqr(g2)) + 2*AbsSqr(LamSD)*(-75*
      traceYdAdjYd - 25*traceYeAdjYe - 20*AbsSqr(LamSU) - 90*AbsSqr(LamTD) + 36
      *Sqr(g1) + 180*Sqr(g2)) + 160*traceYdAdjYd*Sqr(g3) - 140*Sqr(LamSD)*Sqr(
      Conj(LamSD)) - 75*Sqr(LamTD)*Sqr(Conj(LamTD)))));


   return beta_BMuD;
}

/**
 * Calculates the 3-loop beta function of BMuD.
 *
 * @return 3-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_BMuD_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuD;

   beta_BMuD = 0;


   return beta_BMuD;
}

} // namespace flexiblesusy