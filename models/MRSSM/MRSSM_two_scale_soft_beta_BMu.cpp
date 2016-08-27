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

// File generated at Sat 27 Aug 2016 12:04:30

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
 * Calculates the one-loop beta function of BMu.
 *
 * @return one-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_BMu_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_BMu;

   beta_BMu = Re(0.1*oneOver16PiSqr*BMu*(30*traceYdAdjYd + 10*
      traceYeAdjYe + 30*traceYuAdjYu + 10*AbsSqr(LamSD) + 10*AbsSqr(LamSU) + 15
      *AbsSqr(LamTD) + 15*AbsSqr(LamTU) - 6*Sqr(g1) - 30*Sqr(g2)));


   return beta_BMu;
}

/**
 * Calculates the two-loop beta function of BMu.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_BMu_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_BMu;

   beta_BMu = Re(0.05*twoLoop*BMu*(90*Power(g1,4) + 330*Power(g2,4) - 180
      *traceYdAdjYdYdAdjYd - 120*traceYdAdjYuYuAdjYd - 60*traceYeAdjYeYeAdjYe -
      180*traceYuAdjYuYuAdjYu - 20*AbsSqr(LamSD)*(4*AbsSqr(LamSU) + 3*AbsSqr(
      LamTD)) - 60*AbsSqr(LamSU)*AbsSqr(LamTU) - 60*AbsSqr(LamTD)*AbsSqr(LamTU)
      - 8*traceYdAdjYd*Sqr(g1) + 24*traceYeAdjYe*Sqr(g1) + 16*traceYuAdjYu*Sqr
      (g1) + 120*AbsSqr(LamTD)*Sqr(g2) + 120*AbsSqr(LamTU)*Sqr(g2) + 36*Sqr(g1)
      *Sqr(g2) + 320*traceYdAdjYd*Sqr(g3) + 320*traceYuAdjYu*Sqr(g3) - 60*Sqr(
      LamSD)*Sqr(Conj(LamSD)) - 60*Sqr(LamSU)*Sqr(Conj(LamSU)) - 75*Sqr(LamTD)*
      Sqr(Conj(LamTD)) - 75*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_BMu;
}

/**
 * Calculates the three-loop beta function of BMu.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_BMu_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMu;

   beta_BMu = 0;


   return beta_BMu;
}

} // namespace flexiblesusy
