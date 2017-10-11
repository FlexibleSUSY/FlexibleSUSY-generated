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

// File generated at Tue 10 Oct 2017 22:49:09

#include "MSSMRHN_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mHu2.
 *
 * @return 1-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHu2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double traceml2AdjYvYv = TRACE_STRUCT.traceml2AdjYvYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double tracemv2YvAdjYv = TRACE_STRUCT.tracemv2YvAdjYv;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHu2;

   beta_mHu2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 6*
      traceconjTYuTpTYu + 2*traceconjTYvTpTYv + 2*traceml2AdjYvYv + 6*
      tracemq2AdjYuYu + 6*tracemu2YuAdjYu + 2*tracemv2YvAdjYv + 6*mHu2*
      traceYuAdjYu + 2*mHu2*traceYvAdjYv - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*AbsSqr
      (MassWB)*Sqr(g2)));


   return beta_mHu2;
}

/**
 * Calculates the 2-loop beta function of mHu2.
 *
 * @return 2-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYeAdjYvTYvAdjTYe =
      TRACE_STRUCT.traceYeAdjYvTYvAdjTYe;
   const double traceYeAdjTYvTYvAdjYe =
      TRACE_STRUCT.traceYeAdjTYvTYvAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu =
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu =
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double traceYvAdjYeTYeAdjTYv =
      TRACE_STRUCT.traceYvAdjYeTYeAdjTYv;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTYvAdjTYv =
      TRACE_STRUCT.traceYvAdjYvTYvAdjTYv;
   const double traceYvAdjTYeTYeAdjYv =
      TRACE_STRUCT.traceYvAdjTYeTYeAdjYv;
   const double traceYvAdjTYvTYvAdjYv =
      TRACE_STRUCT.traceYvAdjTYvTYvAdjYv;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYvYvAdjYe =
      TRACE_STRUCT.traceme2YeAdjYvYvAdjYe;
   const double traceml2AdjYeYeAdjYvYv =
      TRACE_STRUCT.traceml2AdjYeYeAdjYvYv;
   const double traceml2AdjYvYvAdjYeYe =
      TRACE_STRUCT.traceml2AdjYvYvAdjYeYe;
   const double traceml2AdjYvYvAdjYvYv =
      TRACE_STRUCT.traceml2AdjYvYvAdjYvYv;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double tracemv2YvAdjYeYeAdjYv =
      TRACE_STRUCT.tracemv2YvAdjYeYeAdjYv;
   const double tracemv2YvAdjYvYvAdjYv =
      TRACE_STRUCT.tracemv2YvAdjYvYvAdjYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHu2;

   beta_mHu2 = Re(twoLoop*(3.0983866769659336*g1*Tr31 - 6*
      tracemd2YdAdjYuYuAdjYd - 2*traceme2YeAdjYvYvAdjYe - 2*
      traceml2AdjYeYeAdjYvYv - 2*traceml2AdjYvYvAdjYeYe - 12*
      traceml2AdjYvYvAdjYvYv - 6*tracemq2AdjYdYdAdjYuYu - 6*
      tracemq2AdjYuYuAdjYdYd - 36*tracemq2AdjYuYuAdjYuYu - 6*
      tracemu2YuAdjYdYdAdjYu - 36*tracemu2YuAdjYuYuAdjYu - 2*
      tracemv2YvAdjYeYeAdjYv - 12*tracemv2YvAdjYvYvAdjYv - 6*
      traceYdAdjTYuTYuAdjYd - 6*traceYdAdjYuTYuAdjTYd - 6*mHd2*
      traceYdAdjYuYuAdjYd - 6*mHu2*traceYdAdjYuYuAdjYd - 2*
      traceYeAdjTYvTYvAdjYe - 2*traceYeAdjYvTYvAdjTYe - 2*mHd2*
      traceYeAdjYvYvAdjYe - 2*mHu2*traceYeAdjYvYvAdjYe - 6*
      traceYuAdjTYdTYdAdjYu - 36*traceYuAdjTYuTYuAdjYu - 6*
      traceYuAdjYdTYdAdjTYu - 36*traceYuAdjYuTYuAdjTYu - 36*mHu2*
      traceYuAdjYuYuAdjYu - 2*traceYvAdjTYeTYeAdjYv - 12*traceYvAdjTYvTYvAdjYv
      - 2*traceYvAdjYeTYeAdjTYv - 12*traceYvAdjYvTYvAdjTYv - 12*mHu2*
      traceYvAdjYvYvAdjYv + 6*Tr22*Quad(g2) + 1.2*Tr2U111*Sqr(g1) + 1.6*
      traceconjTYuTpTYu*Sqr(g1) - 1.6*MassB*traceconjTYuTpYu*Sqr(g1) + 1.6*
      tracemq2AdjYuYu*Sqr(g1) + 1.6*tracemu2YuAdjYu*Sqr(g1) + 1.6*mHu2*
      traceYuAdjYu*Sqr(g1) + 0.6*Conj(MassWB)*Sqr(g2)*(3*(MassB + 2*MassWB)*Sqr
      (g1) + 55*MassWB*Sqr(g2)) + 0.04*Conj(MassB)*Sqr(g1)*(621*MassB*Sqr(g1) +
      5*(-8*(traceAdjYuTYu - 2*MassB*traceYuAdjYu) + 9*(2*MassB + MassWB)*Sqr(
      g2))) + 32*traceconjTYuTpTYu*Sqr(g3) - 32*MassG*traceconjTYuTpYu*Sqr(g3)
      + 32*tracemq2AdjYuYu*Sqr(g3) + 32*tracemu2YuAdjYu*Sqr(g3) + 32*mHu2*
      traceYuAdjYu*Sqr(g3) + 64*traceYuAdjYu*AbsSqr(MassG)*Sqr(g3) - 32*
      traceAdjYuTYu*Conj(MassG)*Sqr(g3)));


   return beta_mHu2;
}

/**
 * Calculates the 3-loop beta function of mHu2.
 *
 * @return 3-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy