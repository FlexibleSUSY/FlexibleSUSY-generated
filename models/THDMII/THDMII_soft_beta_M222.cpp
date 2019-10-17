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

// File generated at Wed 16 Oct 2019 21:28:00

#include "THDMII_soft_parameters.hpp"
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
double THDMII_soft_parameters::calc_beta_M222_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M222;

   beta_M222 = Re(0.1*oneOver16PiSqr*(40*Lambda3*M112 + 20*Lambda4*M112 + 120*
      Lambda2*M222 + 60*M222*traceYuAdjYu - 60*M122*Conj(Lambda7) - 60*Lambda7*
      Conj(M122) - 9*M222*Sqr(g1) - 45*M222*Sqr(g2)));


   return beta_M222;
}

/**
 * Calculates the 2-loop beta function of M222.
 *
 * @return 2-loop beta function
 */
double THDMII_soft_parameters::calc_beta_M222_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M222;

   beta_M222 = Re(0.0025*twoLoop*(-3200*Lambda3*Lambda4*M112 + 1800*Lambda5*
      Lambda6*M122 + 4200*Lambda5*Lambda7*M122 - 800*Lambda3*Lambda4*M222 -
      9600*Lambda3*M112*traceYdAdjYd - 4800*Lambda4*M112*traceYdAdjYd - 1800*
      M222*traceYdAdjYuYuAdjYd - 3200*Lambda3*M112*traceYeAdjYe - 1600*Lambda4*
      M112*traceYeAdjYe - 28800*Lambda2*M222*traceYuAdjYu - 5400*M222*
      traceYuAdjYuYuAdjYu - 4800*M112*AbsSqr(Lambda5) - 1200*M222*AbsSqr(
      Lambda5) - 7200*M112*AbsSqr(Lambda6) + 1200*M222*AbsSqr(Lambda6) - 7200*
      M112*AbsSqr(Lambda7) - 10800*M222*AbsSqr(Lambda7) - 1200*Lambda1*M122*
      Conj(Lambda6) + 1800*Lambda3*M122*Conj(Lambda6) + 1800*Lambda4*M122*Conj(
      Lambda6) + 13200*Lambda2*M122*Conj(Lambda7) + 4200*Lambda3*M122*Conj(
      Lambda7) + 4200*Lambda4*M122*Conj(Lambda7) + 7200*M122*traceYdAdjYd*Conj(
      Lambda7) + 2400*M122*traceYeAdjYe*Conj(Lambda7) + 7200*M122*traceYuAdjYu*
      Conj(Lambda7) - 1200*Lambda1*Lambda6*Conj(M122) + 1800*Lambda3*Lambda6*
      Conj(M122) + 1800*Lambda4*Lambda6*Conj(M122) + 13200*Lambda2*Lambda7*Conj
      (M122) + 4200*Lambda3*Lambda7*Conj(M122) + 4200*Lambda4*Lambda7*Conj(M122
      ) + 7200*Lambda7*traceYdAdjYd*Conj(M122) + 2400*Lambda7*traceYeAdjYe*Conj
      (M122) + 7200*Lambda7*traceYuAdjYu*Conj(M122) + 1800*Conj(Lambda5)*Conj(
      Lambda6)*Conj(M122) + 4200*Conj(Lambda5)*Conj(Lambda7)*Conj(M122) + 360*
      M112*Quad(g1) + 1737*M222*Quad(g1) + 3000*M112*Quad(g2) - 3075*M222*Quad(
      g2) + 1920*Lambda3*M112*Sqr(g1) + 960*Lambda4*M112*Sqr(g1) + 5760*Lambda2
      *M222*Sqr(g1) + 1700*M222*traceYuAdjYu*Sqr(g1) - 2880*M122*Conj(Lambda7)*
      Sqr(g1) - 2880*Lambda7*Conj(M122)*Sqr(g1) + 9600*Lambda3*M112*Sqr(g2) +
      4800*Lambda4*M112*Sqr(g2) + 28800*Lambda2*M222*Sqr(g2) + 4500*M222*
      traceYuAdjYu*Sqr(g2) - 14400*M122*Conj(Lambda7)*Sqr(g2) - 14400*Lambda7*
      Conj(M122)*Sqr(g2) + 450*M222*Sqr(g1)*Sqr(g2) + 16000*M222*traceYuAdjYu*
      Sqr(g3) - 24000*M222*Sqr(Lambda2) - 3200*M112*Sqr(Lambda3) - 800*M222*Sqr
      (Lambda3) - 3200*M112*Sqr(Lambda4) - 800*M222*Sqr(Lambda4)));


   return beta_M222;
}

/**
 * Calculates the 3-loop beta function of M222.
 *
 * @return 3-loop beta function
 */
double THDMII_soft_parameters::calc_beta_M222_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

/**
 * Calculates the 4-loop beta function of M222.
 *
 * @return 4-loop beta function
 */
double THDMII_soft_parameters::calc_beta_M222_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

/**
 * Calculates the 5-loop beta function of M222.
 *
 * @return 5-loop beta function
 */
double THDMII_soft_parameters::calc_beta_M222_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return beta_M222;
}

} // namespace flexiblesusy
