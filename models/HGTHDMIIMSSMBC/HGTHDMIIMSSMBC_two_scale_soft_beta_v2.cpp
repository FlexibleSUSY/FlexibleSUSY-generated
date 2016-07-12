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

// File generated at Tue 12 Jul 2016 10:37:47

#include "HGTHDMIIMSSMBC_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of v2.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v2;

   beta_v2 = Re(0.1*oneOver16PiSqr*v2*(-30*traceYuAdjYu + 6*Sqr(g1) + 30*
      Sqr(g2) - 15*Sqr(g2u) - 5*Sqr(g2up)));


   return beta_v2;
}

/**
 * Calculates the two-loop beta function of v2.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v2;

   beta_v2 = Re(twoLoop*(-3*Lambda1*Lambda6*v1 - 1.5*Lambda3*Lambda6*v1 -
      1.5*Lambda4*Lambda6*v1 - 1.5*Lambda5*Lambda6*v1 - 3*Lambda2*Lambda7*v1 -
      1.5*Lambda3*Lambda7*v1 - 1.5*Lambda4*Lambda7*v1 - 1.5*Lambda5*Lambda7*v1
      - 1.80375*Power(g1,4)*v2 + 6.28125*Power(g2,4)*v2 + 2.8125*Power(g2u,4)*
      v2 + 0.5625*Power(g2up,4)*v2 - Lambda3*Lambda4*v2 + 2.25*
      traceYdAdjYuYuAdjYd*v2 + 6.75*traceYuAdjYuYuAdjYu*v2 + 0.1125*v2*Sqr(g1)*
      Sqr(g2) - 0.7875*v2*Sqr(g1)*Sqr(g2u) + 1.125*v2*Sqr(g1d)*Sqr(g2u) -
      11.4375*v2*Sqr(g2)*Sqr(g2u) - 0.2625*v2*Sqr(g1)*Sqr(g2up) + 0.375*v2*Sqr(
      g1dp)*Sqr(g2up) - 1.3125*v2*Sqr(g2)*Sqr(g2up) + 1.125*v2*Sqr(g2u)*Sqr(
      g2up) - 0.025*traceYuAdjYu*v2*(103*Sqr(g1) + 315*Sqr(g2) + 800*Sqr(g3)) -
      6*v2*Sqr(Lambda2) - v2*Sqr(Lambda3) - v2*Sqr(Lambda4) - 1.5*v2*Sqr(
      Lambda5) - 1.5*v2*Sqr(Lambda6) - 4.5*v2*Sqr(Lambda7)));


   return beta_v2;
}

/**
 * Calculates the three-loop beta function of v2.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v2;

   beta_v2 = 0;


   return beta_v2;
}

} // namespace flexiblesusy
