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
 * Calculates the 1-loop beta function of v2.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_v2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v2;

   beta_v2 = Re(0.6*v2*(-5*traceYuAdjYu + Sqr(g1) + 5*Sqr(g2)));


   return oneLoop * beta_v2;
}

/**
 * Calculates the 2-loop beta function of v2.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_v2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v2;

   beta_v2 = Re(0.00125*(-1200*Lambda1*Lambda6*v1 - 600*Lambda3*Lambda6*v1 -
      600*Lambda4*Lambda6*v1 - 600*Lambda5*Lambda6*v1 - 1200*Lambda2*Lambda7*v1
       - 600*Lambda3*Lambda7*v1 - 600*Lambda4*Lambda7*v1 - 600*Lambda5*Lambda7*
      v1 - 800*Lambda3*Lambda4*v2 + 1800*traceYdAdjYuYuAdjYd*v2 + 5400*
      traceYuAdjYuYuAdjYu*v2 - 1200*v2*AbsSqr(Lambda5) - 1200*v2*AbsSqr(Lambda6
      ) - 3600*v2*AbsSqr(Lambda7) - 1200*Lambda1*v1*Conj(Lambda6) - 600*Lambda3
      *v1*Conj(Lambda6) - 600*Lambda4*v1*Conj(Lambda6) - 600*v1*Conj(Lambda5)*
      Conj(Lambda6) - 1200*Lambda2*v1*Conj(Lambda7) - 600*Lambda3*v1*Conj(
      Lambda7) - 600*Lambda4*v1*Conj(Lambda7) - 600*v1*Conj(Lambda5)*Conj(
      Lambda7) - 1407*v2*Quad(g1) + 7925*v2*Quad(g2) - 2420*traceYuAdjYu*v2*Sqr
      (g1) - 8100*traceYuAdjYu*v2*Sqr(g2) + 450*v2*Sqr(g1)*Sqr(g2) - 16000*
      traceYuAdjYu*v2*Sqr(g3) - 4800*v2*Sqr(Lambda2) - 800*v2*Sqr(Lambda3) -
      800*v2*Sqr(Lambda4)));


   return twoLoop * beta_v2;
}

/**
 * Calculates the 3-loop beta function of v2.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_v2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v2;

   beta_v2 = 0;


   return threeLoop * beta_v2;
}

/**
 * Calculates the 4-loop beta function of v2.
 *
 * @return 4-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_v2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v2;

   beta_v2 = 0;


   return fourLoop * beta_v2;
}

/**
 * Calculates the 5-loop beta function of v2.
 *
 * @return 5-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_v2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v2;

   beta_v2 = 0;


   return fiveLoop * beta_v2;
}

} // namespace flexiblesusy
