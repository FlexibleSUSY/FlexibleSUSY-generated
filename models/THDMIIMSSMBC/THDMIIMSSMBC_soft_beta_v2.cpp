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

// File generated at Sun 26 Aug 2018 14:08:54

#include "THDMIIMSSMBC_soft_parameters.hpp"
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
double THDMIIMSSMBC_soft_parameters::calc_beta_v2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v2;

   beta_v2 = Re(0.6*oneOver16PiSqr*v2*(-5*traceYuAdjYu + Sqr(g1) + 5*Sqr(g2)));


   return beta_v2;
}

/**
 * Calculates the 2-loop beta function of v2.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v2;

   beta_v2 = Re(twoLoop*(-1.5*Lambda1*Lambda6*v1 - 0.75*Lambda3*Lambda6*v1 -
      0.75*Lambda4*Lambda6*v1 - 0.75*Lambda5*Lambda6*v1 - 1.5*Lambda2*Lambda7*
      v1 - 0.75*Lambda3*Lambda7*v1 - 0.75*Lambda4*Lambda7*v1 - 0.75*Lambda5*
      Lambda7*v1 - Lambda3*Lambda4*v2 + 2.25*traceYdAdjYuYuAdjYd*v2 + 6.75*
      traceYuAdjYuYuAdjYu*v2 - 4.5*v2*AbsSqr(Lambda7) - 0.75*(2*Lambda1*v1 +
      Lambda3*v1 + Lambda4*v1 + 2*Lambda6*v2)*Conj(Lambda6) - 1.5*Lambda2*v1*
      Conj(Lambda7) - 0.75*Lambda3*v1*Conj(Lambda7) - 0.75*Lambda4*v1*Conj(
      Lambda7) - 0.75*Conj(Lambda5)*(2*Lambda5*v2 + v1*Conj(Lambda6) + v1*Conj(
      Lambda7)) - 1.60875*v2*Quad(g1) + 11.15625*v2*Quad(g2) - 3.025*
      traceYuAdjYu*v2*Sqr(g1) - 10.125*traceYuAdjYu*v2*Sqr(g2) + 0.5625*v2*Sqr(
      g1)*Sqr(g2) - 20*traceYuAdjYu*v2*Sqr(g3) - 6*v2*Sqr(Lambda2) - v2*Sqr(
      Lambda3) - v2*Sqr(Lambda4)));


   return beta_v2;
}

/**
 * Calculates the 3-loop beta function of v2.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v2;

   beta_v2 = 0;


   return beta_v2;
}

/**
 * Calculates the 4-loop beta function of v2.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v2;

   beta_v2 = 0;


   return beta_v2;
}

} // namespace flexiblesusy
