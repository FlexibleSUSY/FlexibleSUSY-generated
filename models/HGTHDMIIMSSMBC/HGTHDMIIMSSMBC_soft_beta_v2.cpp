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

// File generated at Tue 10 Oct 2017 20:54:06

#include "HGTHDMIIMSSMBC_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of v2.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v2;

   beta_v2 = Re(0.1*oneOver16PiSqr*v2*(6*Sqr(g1) + 5*(-6*traceYuAdjYu + 6
      *Sqr(g2) - 3*Sqr(g2u) - Sqr(g2up))));


   return beta_v2;
}

/**
 * Calculates the 2-loop beta function of v2.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v2;

   beta_v2 = Re(twoLoop*(-3*Lambda1*Lambda6*v1 - 1.5*Lambda4*Lambda6*v1 -
      1.5*Lambda5*Lambda6*v1 - 3*Lambda2*Lambda7*v1 - 1.5*Lambda4*Lambda7*v1 -
      1.5*Lambda5*Lambda7*v1 + 2.25*traceYdAdjYuYuAdjYd*v2 + 6.75*
      traceYuAdjYuYuAdjYu*v2 - 0.5*Lambda3*(3*Lambda6*v1 + 3*Lambda7*v1 + 2*
      Lambda4*v2) - 1.75875*v2*Quad(g1) + 7.40625*v2*Quad(g2) + 2.8125*v2*Quad(
      g2u) + 0.5625*v2*Quad(g2up) - 3.025*traceYuAdjYu*v2*Sqr(g1) - 10.125*
      traceYuAdjYu*v2*Sqr(g2) + 0.5625*v2*Sqr(g1)*Sqr(g2) - 1.0125*v2*Sqr(g1)*
      Sqr(g2u) + 1.125*v2*Sqr(g1d)*Sqr(g2u) - 12.5625*v2*Sqr(g2)*Sqr(g2u) -
      0.3375*v2*Sqr(g1)*Sqr(g2up) + 0.375*v2*Sqr(g1dp)*Sqr(g2up) - 1.6875*v2*
      Sqr(g2)*Sqr(g2up) + 1.125*v2*Sqr(g2u)*Sqr(g2up) - 20*traceYuAdjYu*v2*Sqr(
      g3) - 6*v2*Sqr(Lambda2) - v2*Sqr(Lambda3) - v2*Sqr(Lambda4) - 1.5*v2*Sqr(
      Lambda5) - 1.5*v2*Sqr(Lambda6) - 4.5*v2*Sqr(Lambda7)));


   return beta_v2;
}

/**
 * Calculates the 3-loop beta function of v2.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v2;

   beta_v2 = 0;


   return beta_v2;
}

} // namespace flexiblesusy
