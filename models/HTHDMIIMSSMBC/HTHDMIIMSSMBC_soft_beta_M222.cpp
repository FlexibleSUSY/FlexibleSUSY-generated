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

// File generated at Fri 20 Oct 2017 08:36:29

#include "HTHDMIIMSSMBC_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of M222.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M222_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M222;

   beta_M222 = Re(oneOver16PiSqr*(4*Lambda3*M112 + 2*Lambda4*M112 - 12*
      Lambda7*M122 + 12*Lambda2*M222 + 6*M222*traceYuAdjYu - 0.9*M222*Sqr(g1) -
      4.5*M222*Sqr(g2)));


   return beta_M222;
}

/**
 * Calculates the 2-loop beta function of M222.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M222_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M222;

   beta_M222 = Re(twoLoop*(-8*Lambda3*Lambda4*M112 + 12*Lambda3*Lambda6*
      M122 + 12*Lambda4*Lambda6*M122 + 12*Lambda5*Lambda6*M122 + 72*Lambda2*
      Lambda7*M122 + 24*Lambda3*Lambda7*M122 + 24*Lambda4*Lambda7*M122 + 24*
      Lambda5*Lambda7*M122 - 2*Lambda3*Lambda4*M222 - 24*Lambda3*M112*
      traceYdAdjYd - 12*Lambda4*M112*traceYdAdjYd + 36*Lambda7*M122*
      traceYdAdjYd - 4.5*M222*traceYdAdjYuYuAdjYd - 8*Lambda3*M112*traceYeAdjYe
      - 4*Lambda4*M112*traceYeAdjYe + 12*Lambda7*M122*traceYeAdjYe + 36*
      Lambda7*M122*traceYuAdjYu - 72*Lambda2*M222*traceYuAdjYu - 13.5*M222*
      traceYuAdjYuYuAdjYu + (24*Lambda3*M112 + 12*Lambda4*M112 - 72*Lambda7*
      M122 + 72*Lambda2*M222 + 11.25*M222*traceYuAdjYu)*Sqr(g2) + 0.025*Sqr(g1)
      *(192*Lambda3*M112 + 96*Lambda4*M112 - 576*Lambda7*M122 + 576*Lambda2*
      M222 + 170*M222*traceYuAdjYu + 45*M222*Sqr(g2)) + 40*M222*traceYuAdjYu*
      Sqr(g3) - 60*M222*Sqr(Lambda2) - 8*M112*Sqr(Lambda3) - 2*M222*Sqr(Lambda3
      ) - 8*M112*Sqr(Lambda4) - 2*M222*Sqr(Lambda4) - 12*M112*Sqr(Lambda5) - 3*
      M222*Sqr(Lambda5) - 18*M112*Sqr(Lambda6) + 3*M222*Sqr(Lambda6) - 18*M112*
      Sqr(Lambda7) - 27*M222*Sqr(Lambda7) + 0.0075*Quad(g1)*(120*M112 + 619*
      M222 - 288*Sqr(Mu)) + Quad(g2)*(7.5*M112 - 5.1875*M222 - 18*Sqr(Mu))));


   return beta_M222;
}

/**
 * Calculates the 3-loop beta function of M222.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M222_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

} // namespace flexiblesusy
