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

// File generated at Tue 10 Oct 2017 21:12:37

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
 * Calculates the 1-loop beta function of M112.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M112_1_loop(const Soft_traces& soft_traces) const
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
 * Calculates the 2-loop beta function of M112.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M112_2_loop(const Soft_traces& soft_traces) const
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
      traceYuAdjYu - 12*Lambda4*M222*traceYuAdjYu + 0.0225*(193*M112 + 40*M222)
      *Quad(g1) - 0.1875*(41*M112 - 40*M222)*Quad(g2) + 0.75*(96*Lambda1*M112 -
      96*Lambda6*M122 + 32*Lambda3*M222 + 16*Lambda4*M222 + 15*M112*
      traceYdAdjYd + 5*M112*traceYeAdjYe)*Sqr(g2) + 0.025*Sqr(g1)*(576*Lambda1*
      M112 - 576*Lambda6*M122 + 192*Lambda3*M222 + 96*Lambda4*M222 + 50*M112*
      traceYdAdjYd + 150*M112*traceYeAdjYe + 45*M112*Sqr(g2)) + 40*M112*
      traceYdAdjYd*Sqr(g3) - 60*M112*Sqr(Lambda1) - 2*M112*Sqr(Lambda3) - 8*
      M222*Sqr(Lambda3) - 2*M112*Sqr(Lambda4) - 8*M222*Sqr(Lambda4) - 3*M112*
      Sqr(Lambda5) - 12*M222*Sqr(Lambda5) - 27*M112*Sqr(Lambda6) - 18*M222*Sqr(
      Lambda6) + 3*M112*Sqr(Lambda7) - 18*M222*Sqr(Lambda7)));


   return beta_M112;
}

/**
 * Calculates the 3-loop beta function of M112.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M112_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

} // namespace flexiblesusy
