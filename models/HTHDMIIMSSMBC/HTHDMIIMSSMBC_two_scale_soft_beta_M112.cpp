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

// File generated at Wed 12 Apr 2017 11:00:08

#include "HTHDMIIMSSMBC_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of M112.
 *
 * @return one-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M112_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_M112;

   beta_M112 = Re(oneOver16PiSqr*(2*(6*Lambda1*M112 - 6*Lambda6*M122 + 2*
      Lambda3*M222 + Lambda4*M222 + 3*M112*traceYdAdjYd + M112*traceYeAdjYe) -
      0.9*M112*Sqr(g1) - 4.5*M112*Sqr(g2)));


   return beta_M112;
}

/**
 * Calculates the two-loop beta function of M112.
 *
 * @return two-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M112_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_M112;

   beta_M112 = Re(twoLoop*(-2*Lambda3*Lambda4*M112 + 72*Lambda1*Lambda6*
      M122 + 24*Lambda3*Lambda6*M122 + 24*Lambda4*Lambda6*M122 + 24*Lambda5*
      Lambda6*M122 + 12*Lambda3*Lambda7*M122 + 12*Lambda4*Lambda7*M122 + 12*
      Lambda5*Lambda7*M122 - 8*Lambda3*Lambda4*M222 - 72*Lambda1*M112*
      traceYdAdjYd + 36*Lambda6*M122*traceYdAdjYd - 13.5*M112*
      traceYdAdjYdYdAdjYd - 4.5*M112*traceYdAdjYuYuAdjYd - 24*Lambda1*M112*
      traceYeAdjYe + 12*Lambda6*M122*traceYeAdjYe - 4.5*M112*
      traceYeAdjYeYeAdjYe + 36*Lambda6*M122*traceYuAdjYu - 24*Lambda3*M222*
      traceYuAdjYu - 12*Lambda4*M222*traceYuAdjYu + 0.75*(96*Lambda1*M112 - 96*
      Lambda6*M122 + 32*Lambda3*M222 + 16*Lambda4*M222 + 15*M112*traceYdAdjYd +
      5*M112*traceYeAdjYe)*Sqr(g2) + 0.025*Sqr(g1)*(576*Lambda1*M112 - 576*
      Lambda6*M122 + 192*Lambda3*M222 + 96*Lambda4*M222 + 50*M112*traceYdAdjYd
      + 150*M112*traceYeAdjYe + 45*M112*Sqr(g2)) + 40*M112*traceYdAdjYd*Sqr(g3)
      - 60*M112*Sqr(Lambda1) - 2*M112*Sqr(Lambda3) - 8*M222*Sqr(Lambda3) - 2*
      M112*Sqr(Lambda4) - 8*M222*Sqr(Lambda4) - 3*M112*Sqr(Lambda5) - 12*M222*
      Sqr(Lambda5) - 27*M112*Sqr(Lambda6) - 18*M222*Sqr(Lambda6) + 3*M112*Sqr(
      Lambda7) - 18*M222*Sqr(Lambda7) + 0.0075*Power(g1,4)*(619*M112 + 120*M222
      - 288*Sqr(Mu)) + Power(g2,4)*(-5.1875*M112 + 7.5*M222 - 18*Sqr(Mu))));


   return beta_M112;
}

/**
 * Calculates the three-loop beta function of M112.
 *
 * @return three-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M112_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

} // namespace flexiblesusy
