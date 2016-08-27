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

// File generated at Sat 27 Aug 2016 11:41:46

#include "THDMIIMSSMBC_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of M222.
 *
 * @return one-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M222_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M222;

   beta_M222 = Re(oneOver16PiSqr*(4*Lambda3*M112 + 2*Lambda4*M112 - 12*
      Lambda7*M122 + 12*Lambda2*M222 + 6*M222*traceYuAdjYu - 0.9*M222*Sqr(g1) -
      4.5*M222*Sqr(g2)));


   return beta_M222;
}

/**
 * Calculates the two-loop beta function of M222.
 *
 * @return two-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M222_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M222;

   beta_M222 = Re(twoLoop*(0.9*Power(g1,4)*M112 + 7.5*Power(g2,4)*M112 -
      8*Lambda3*Lambda4*M112 + 12*Lambda3*Lambda6*M122 + 12*Lambda4*Lambda6*
      M122 + 12*Lambda5*Lambda6*M122 + 72*Lambda2*Lambda7*M122 + 24*Lambda3*
      Lambda7*M122 + 24*Lambda4*Lambda7*M122 + 24*Lambda5*Lambda7*M122 + 4.3425
      *Power(g1,4)*M222 - 7.6875*Power(g2,4)*M222 - 2*Lambda3*Lambda4*M222 - 12
      *(2*Lambda3*M112 + Lambda4*M112 - 3*Lambda7*M122)*traceYdAdjYd - 4.5*M222
      *traceYdAdjYuYuAdjYd - 4*(2*Lambda3*M112 + Lambda4*M112 - 3*Lambda7*M122)
      *traceYeAdjYe + 36*Lambda7*M122*traceYuAdjYu - 72*Lambda2*M222*
      traceYuAdjYu - 13.5*M222*traceYuAdjYuYuAdjYu + 4.8*Lambda3*M112*Sqr(g1) +
      2.4*Lambda4*M112*Sqr(g1) - 14.4*Lambda7*M122*Sqr(g1) + 14.4*Lambda2*M222
      *Sqr(g1) + 4.25*M222*traceYuAdjYu*Sqr(g1) + 24*Lambda3*M112*Sqr(g2) + 12*
      Lambda4*M112*Sqr(g2) - 72*Lambda7*M122*Sqr(g2) + 72*Lambda2*M222*Sqr(g2)
      + 11.25*M222*traceYuAdjYu*Sqr(g2) + 1.125*M222*Sqr(g1)*Sqr(g2) + 40*M222*
      traceYuAdjYu*Sqr(g3) - 60*M222*Sqr(Lambda2) - 8*M112*Sqr(Lambda3) - 2*
      M222*Sqr(Lambda3) - 8*M112*Sqr(Lambda4) - 2*M222*Sqr(Lambda4) - 12*M112*
      Sqr(Lambda5) - 3*M222*Sqr(Lambda5) - 18*M112*Sqr(Lambda6) + 3*M222*Sqr(
      Lambda6) - 18*M112*Sqr(Lambda7) - 27*M222*Sqr(Lambda7)));


   return beta_M222;
}

/**
 * Calculates the three-loop beta function of M222.
 *
 * @return three-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M222_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

} // namespace flexiblesusy
