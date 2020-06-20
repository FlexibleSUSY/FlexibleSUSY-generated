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
 * Calculates the 1-loop beta function of M122.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M122_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M122;

   beta_M122 = Re(0.1*(-60*Lambda6*M112 + 20*Lambda3*M122 + 40*Lambda4*M122 -
      60*Lambda7*M222 + 30*M122*traceYdAdjYd + 10*M122*traceYeAdjYe + 30*M122*
      traceYuAdjYu + 60*Conj(Lambda5)*Conj(M122) - 9*M122*Sqr(g1) - 45*M122*Sqr
      (g2)));


   return oneLoop * beta_M122;
}

/**
 * Calculates the 2-loop beta function of M122.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M122_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M122;

   beta_M122 = Re(0.0025*(13200*Lambda1*Lambda6*M112 + 4200*Lambda3*Lambda6*
      M112 + 4200*Lambda4*Lambda6*M112 - 1200*Lambda2*Lambda7*M112 + 1800*
      Lambda3*Lambda7*M112 + 1800*Lambda4*Lambda7*M112 - 4800*Lambda1*Lambda3*
      M122 - 4800*Lambda2*Lambda3*M122 - 4800*Lambda1*Lambda4*M122 - 4800*
      Lambda2*Lambda4*M122 - 2400*Lambda3*Lambda4*M122 - 1200*Lambda1*Lambda6*
      M222 + 1800*Lambda3*Lambda6*M222 + 1800*Lambda4*Lambda6*M222 + 13200*
      Lambda2*Lambda7*M222 + 4200*Lambda3*Lambda7*M222 + 4200*Lambda4*Lambda7*
      M222 + 14400*Lambda6*M112*traceYdAdjYd - 2400*Lambda3*M122*traceYdAdjYd -
      4800*Lambda4*M122*traceYdAdjYd - 2700*M122*traceYdAdjYdYdAdjYd - 6600*
      M122*traceYdAdjYuYuAdjYd + 4800*Lambda6*M112*traceYeAdjYe - 800*Lambda3*
      M122*traceYeAdjYe - 1600*Lambda4*M122*traceYeAdjYe - 900*M122*
      traceYeAdjYeYeAdjYe - 2400*Lambda3*M122*traceYuAdjYu - 4800*Lambda4*M122*
      traceYuAdjYu + 14400*Lambda7*M222*traceYuAdjYu - 2700*M122*
      traceYuAdjYuYuAdjYu + 1200*M122*AbsSqr(Lambda5) - 4800*Lambda7*M122*Conj(
      Lambda6) + 4200*M112*Conj(Lambda5)*Conj(Lambda6) + 1800*M222*Conj(Lambda5
      )*Conj(Lambda6) - 4800*Lambda6*M122*Conj(Lambda7) + 1800*M112*Conj(
      Lambda5)*Conj(Lambda7) + 4200*M222*Conj(Lambda5)*Conj(Lambda7) - 4800*
      Lambda6*Lambda7*Conj(M122) - 4800*Lambda1*Conj(Lambda5)*Conj(M122) - 4800
      *Lambda2*Conj(Lambda5)*Conj(M122) - 4800*Lambda3*Conj(Lambda5)*Conj(M122)
      - 4800*Lambda4*Conj(Lambda5)*Conj(M122) - 7200*traceYdAdjYd*Conj(Lambda5)
      *Conj(M122) - 2400*traceYeAdjYe*Conj(Lambda5)*Conj(M122) - 7200*
      traceYuAdjYu*Conj(Lambda5)*Conj(M122) + 1377*M122*Quad(g1) - 6075*M122*
      Quad(g2) - 2880*Lambda6*M112*Sqr(g1) + 960*Lambda3*M122*Sqr(g1) + 1920*
      Lambda4*M122*Sqr(g1) - 2880*Lambda7*M222*Sqr(g1) + 250*M122*traceYdAdjYd*
      Sqr(g1) + 750*M122*traceYeAdjYe*Sqr(g1) + 850*M122*traceYuAdjYu*Sqr(g1) +
      2880*Conj(Lambda5)*Conj(M122)*Sqr(g1) - 14400*Lambda6*M112*Sqr(g2) + 4800
      *Lambda3*M122*Sqr(g2) + 9600*Lambda4*M122*Sqr(g2) - 14400*Lambda7*M222*
      Sqr(g2) + 2250*M122*traceYdAdjYd*Sqr(g2) + 750*M122*traceYeAdjYe*Sqr(g2)
      + 2250*M122*traceYuAdjYu*Sqr(g2) + 14400*Conj(Lambda5)*Conj(M122)*Sqr(g2)
      + 450*M122*Sqr(g1)*Sqr(g2) + 8000*M122*traceYdAdjYd*Sqr(g3) + 8000*M122*
      traceYuAdjYu*Sqr(g3) + 2400*M122*Sqr(Lambda1) + 2400*M122*Sqr(Lambda2) -
      4800*Conj(M122)*Sqr(Lambda6) - 4800*Conj(M122)*Sqr(Lambda7)));


   return twoLoop * beta_M122;
}

/**
 * Calculates the 3-loop beta function of M122.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M122_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return threeLoop * beta_M122;
}

/**
 * Calculates the 4-loop beta function of M122.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M122_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return fourLoop * beta_M122;
}

/**
 * Calculates the 5-loop beta function of M122.
 *
 * @return 5-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_M122_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return fiveLoop * beta_M122;
}

} // namespace flexiblesusy
