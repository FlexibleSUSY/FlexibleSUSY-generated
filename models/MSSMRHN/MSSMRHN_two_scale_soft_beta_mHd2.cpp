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

// File generated at Wed 12 Apr 2017 13:08:26

#include "MSSMRHN_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHd2;

   beta_mHd2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 6*
      traceconjTYdTpTYd + 2*traceconjTYeTpTYe + 6*tracemd2YdAdjYd + 2*
      traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*tracemq2AdjYdYd + 6*mHd2*
      traceYdAdjYd + 2*mHd2*traceYeAdjYe - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*AbsSqr
      (MassWB)*Sqr(g2)));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYeAdjYvTYvAdjTYe =
      TRACE_STRUCT.traceYeAdjYvTYvAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYeAdjTYvTYvAdjYe =
      TRACE_STRUCT.traceYeAdjTYvTYvAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYvAdjYeTYeAdjTYv =
      TRACE_STRUCT.traceYvAdjYeTYeAdjTYv;
   const double traceYvAdjTYeTYeAdjYv =
      TRACE_STRUCT.traceYvAdjTYeTYeAdjYv;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceme2YeAdjYvYvAdjYe =
      TRACE_STRUCT.traceme2YeAdjYvYvAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double traceml2AdjYeYeAdjYvYv =
      TRACE_STRUCT.traceml2AdjYeYeAdjYvYv;
   const double traceml2AdjYvYvAdjYeYe =
      TRACE_STRUCT.traceml2AdjYvYvAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemv2YvAdjYeYeAdjYv =
      TRACE_STRUCT.tracemv2YvAdjYeYeAdjYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHd2;

   beta_mHd2 = Re(0.04*twoLoop*(Conj(MassB)*Sqr(g1)*(621*MassB*Sqr(g1) +
      5*(4*(traceAdjYdTYd - 3*traceAdjYeTYe - 2*MassB*traceYdAdjYd + 6*MassB*
      traceYeAdjYe) + 9*(2*MassB + MassWB)*Sqr(g2))) + 5*(3*Conj(MassWB)*Sqr(g2
      )*(3*(MassB + 2*MassWB)*Sqr(g1) + 55*MassWB*Sqr(g2)) - 160*(traceAdjYdTYd
      - 2*MassG*traceYdAdjYd)*Conj(MassG)*Sqr(g3) + 2*(15*Power(g2,4)*Tr22 -
      7.745966692414834*g1*Tr31 + (3*Tr2U111 - 2*(traceconjTYdTpTYd - MassB*
      traceconjTYdTpYd - 3*traceconjTYeTpTYe + 3*MassB*traceconjTYeTpYe +
      tracemd2YdAdjYd - 3*traceme2YeAdjYe - 3*traceml2AdjYeYe + tracemq2AdjYdYd
      + mHd2*traceYdAdjYd - 3*mHd2*traceYeAdjYe))*Sqr(g1) - 5*(18*
      tracemd2YdAdjYdYdAdjYd + 3*tracemd2YdAdjYuYuAdjYd + 6*
      traceme2YeAdjYeYeAdjYe + traceme2YeAdjYvYvAdjYe + 6*
      traceml2AdjYeYeAdjYeYe + traceml2AdjYeYeAdjYvYv + traceml2AdjYvYvAdjYeYe
      + 18*tracemq2AdjYdYdAdjYdYd + 3*tracemq2AdjYdYdAdjYuYu + 3*
      tracemq2AdjYuYuAdjYdYd + 3*tracemu2YuAdjYdYdAdjYu +
      tracemv2YvAdjYeYeAdjYv + 18*traceYdAdjTYdTYdAdjYd + 3*
      traceYdAdjTYuTYuAdjYd + 18*traceYdAdjYdTYdAdjTYd + 18*mHd2*
      traceYdAdjYdYdAdjYd + 3*traceYdAdjYuTYuAdjTYd + 3*mHd2*
      traceYdAdjYuYuAdjYd + 3*mHu2*traceYdAdjYuYuAdjYd + 6*
      traceYeAdjTYeTYeAdjYe + traceYeAdjTYvTYvAdjYe + 6*traceYeAdjYeTYeAdjTYe +
      6*mHd2*traceYeAdjYeYeAdjYe + traceYeAdjYvTYvAdjTYe + mHd2*
      traceYeAdjYvYvAdjYe + mHu2*traceYeAdjYvYvAdjYe + 3*traceYuAdjTYdTYdAdjYu
      + 3*traceYuAdjYdTYdAdjTYu + traceYvAdjTYeTYeAdjYv + traceYvAdjYeTYeAdjTYv
      - 16*(traceconjTYdTpTYd - MassG*traceconjTYdTpYd + tracemd2YdAdjYd +
      tracemq2AdjYdYd + mHd2*traceYdAdjYd)*Sqr(g3))))));


   return beta_mHd2;
}

/**
 * Calculates the three-loop beta function of mHd2.
 *
 * @return three-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHd2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
