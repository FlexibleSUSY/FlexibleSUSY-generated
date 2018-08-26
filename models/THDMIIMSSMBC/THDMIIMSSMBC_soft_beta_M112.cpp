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

// File generated at Sun 26 Aug 2018 14:08:52

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

   beta_M112 = Re(oneOver16PiSqr*(12*Lambda1*M112 + 4*Lambda3*M222 + 2*Lambda4*
      M222 + 6*M112*traceYdAdjYd + 2*M112*traceYeAdjYe - 6*M122*Conj(Lambda6) -
      6*Lambda6*Conj(M122) - 0.9*M112*Sqr(g1) - 4.5*M112*Sqr(g2)));


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

   beta_M112 = Re(twoLoop*(-2*Lambda3*Lambda4*M112 + 12*Lambda5*Lambda6*M122 +
      6*Lambda5*Lambda7*M122 - 8*Lambda3*Lambda4*M222 - 72*Lambda1*M112*
      traceYdAdjYd - 13.5*M112*traceYdAdjYdYdAdjYd - 4.5*M112*
      traceYdAdjYuYuAdjYd - 24*Lambda1*M112*traceYeAdjYe - 4.5*M112*
      traceYeAdjYeYeAdjYe - 24*Lambda3*M222*traceYuAdjYu - 12*Lambda4*M222*
      traceYuAdjYu + 3*M112*AbsSqr(Lambda7) - 18*M222*AbsSqr(Lambda7) + 6*
      Lambda3*M122*Conj(Lambda7) + 6*Lambda4*M122*Conj(Lambda7) + 36*Lambda1*
      Lambda6*Conj(M122) + 12*Lambda3*Lambda6*Conj(M122) + 12*Lambda4*Lambda6*
      Conj(M122) + 6*Lambda3*Lambda7*Conj(M122) + 6*Lambda4*Lambda7*Conj(M122)
      + 18*Lambda6*traceYdAdjYd*Conj(M122) + 6*Lambda6*traceYeAdjYe*Conj(M122)
      + 18*Lambda6*traceYuAdjYu*Conj(M122) - 3*Conj(Lambda5)*(Lambda5*M112 + 4*
      Lambda5*M222 - 4*Conj(Lambda6)*Conj(M122) - 2*Conj(Lambda7)*Conj(M122)) +
      4.3425*M112*Quad(g1) + 0.9*M222*Quad(g1) - 7.6875*M112*Quad(g2) + 7.5*
      M222*Quad(g2) + 14.4*Lambda1*M112*Sqr(g1) + 4.8*Lambda3*M222*Sqr(g1) +
      2.4*Lambda4*M222*Sqr(g1) + 1.25*M112*traceYdAdjYd*Sqr(g1) + 3.75*M112*
      traceYeAdjYe*Sqr(g1) - 7.2*Lambda6*Conj(M122)*Sqr(g1) + 72*Lambda1*M112*
      Sqr(g2) + 24*Lambda3*M222*Sqr(g2) + 12*Lambda4*M222*Sqr(g2) + 11.25*M112*
      traceYdAdjYd*Sqr(g2) + 3.75*M112*traceYeAdjYe*Sqr(g2) - 36*Lambda6*Conj(
      M122)*Sqr(g2) + 1.125*M112*Sqr(g1)*Sqr(g2) - 0.6*Conj(Lambda6)*(15*
      Lambda6*(3*M112 + 2*M222) + 2*M122*(6*Sqr(g1) + 5*(-6*Lambda1 - 2*Lambda3
       - 2*Lambda4 - 3*traceYdAdjYd - traceYeAdjYe - 3*traceYuAdjYu + 6*Sqr(g2)
      ))) + 40*M112*traceYdAdjYd*Sqr(g3) - 60*M112*Sqr(Lambda1) - 2*M112*Sqr(
      Lambda3) - 8*M222*Sqr(Lambda3) - 2*M112*Sqr(Lambda4) - 8*M222*Sqr(Lambda4
      )));


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

/**
 * Calculates the 4-loop beta function of M112.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M112_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

} // namespace flexiblesusy
