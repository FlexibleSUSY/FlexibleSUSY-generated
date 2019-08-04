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

// File generated at Sun 4 Aug 2019 19:58:01

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

   beta_mHu2 = Re(0.2*oneOver16PiSqr*(3.872983346207417*g1*Tr11 + 30*
      traceconjTYuTpTYu + 10*traceconjTYvTpTYv + 10*traceml2AdjYvYv + 30*
      tracemq2AdjYuYu + 30*tracemu2YuAdjYu + 10*tracemv2YvAdjYv + 30*mHu2*
      traceYuAdjYu + 10*mHu2*traceYvAdjYv - 6*AbsSqr(MassB)*Sqr(g1) - 30*AbsSqr
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
   const double traceYdAdjYuTYuAdjTYd = TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd = TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYvYvAdjYe = TRACE_STRUCT.traceYeAdjYvYvAdjYe;
   const double traceYeAdjYvTYvAdjTYe = TRACE_STRUCT.traceYeAdjYvTYvAdjTYe;
   const double traceYeAdjTYvTYvAdjYe = TRACE_STRUCT.traceYeAdjTYvTYvAdjYe;
   const double traceYuAdjYdTYdAdjTYu = TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu = TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu = TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu = TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double traceYvAdjYeTYeAdjTYv = TRACE_STRUCT.traceYvAdjYeTYeAdjTYv;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTYvAdjTYv = TRACE_STRUCT.traceYvAdjYvTYvAdjTYv;
   const double traceYvAdjTYeTYeAdjYv = TRACE_STRUCT.traceYvAdjTYeTYeAdjYv;
   const double traceYvAdjTYvTYvAdjYv = TRACE_STRUCT.traceYvAdjTYvTYvAdjYv;
   const double tracemd2YdAdjYuYuAdjYd = TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYvYvAdjYe = TRACE_STRUCT.traceme2YeAdjYvYvAdjYe;
   const double traceml2AdjYeYeAdjYvYv = TRACE_STRUCT.traceml2AdjYeYeAdjYvYv;
   const double traceml2AdjYvYvAdjYeYe = TRACE_STRUCT.traceml2AdjYvYvAdjYeYe;
   const double traceml2AdjYvYvAdjYvYv = TRACE_STRUCT.traceml2AdjYvYvAdjYvYv;
   const double tracemq2AdjYdYdAdjYuYu = TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd = TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu = TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu = TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu = TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double tracemv2YvAdjYeYeAdjYv = TRACE_STRUCT.tracemv2YvAdjYeYeAdjYv;
   const double tracemv2YvAdjYvYvAdjYv = TRACE_STRUCT.tracemv2YvAdjYvYvAdjYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHu2;

   beta_mHu2 = Re(0.04*twoLoop*(77.45966692414834*g1*Tr31 - 150*
      tracemd2YdAdjYuYuAdjYd - 50*traceme2YeAdjYvYvAdjYe - 50*
      traceml2AdjYeYeAdjYvYv - 50*traceml2AdjYvYvAdjYeYe - 300*
      traceml2AdjYvYvAdjYvYv - 150*tracemq2AdjYdYdAdjYuYu - 150*
      tracemq2AdjYuYuAdjYdYd - 900*tracemq2AdjYuYuAdjYuYu - 150*
      tracemu2YuAdjYdYdAdjYu - 900*tracemu2YuAdjYuYuAdjYu - 50*
      tracemv2YvAdjYeYeAdjYv - 300*tracemv2YvAdjYvYvAdjYv - 150*
      traceYdAdjTYuTYuAdjYd - 150*traceYdAdjYuTYuAdjTYd - 150*mHd2*
      traceYdAdjYuYuAdjYd - 150*mHu2*traceYdAdjYuYuAdjYd - 50*
      traceYeAdjTYvTYvAdjYe - 50*traceYeAdjYvTYvAdjTYe - 50*mHd2*
      traceYeAdjYvYvAdjYe - 50*mHu2*traceYeAdjYvYvAdjYe - 150*
      traceYuAdjTYdTYdAdjYu - 900*traceYuAdjTYuTYuAdjYu - 150*
      traceYuAdjYdTYdAdjTYu - 900*traceYuAdjYuTYuAdjTYu - 900*mHu2*
      traceYuAdjYuYuAdjYu - 50*traceYvAdjTYeTYeAdjYv - 300*
      traceYvAdjTYvTYvAdjYv - 50*traceYvAdjYeTYeAdjTYv - 300*
      traceYvAdjYvTYvAdjTYv - 300*mHu2*traceYvAdjYvYvAdjYv + 621*AbsSqr(MassB)*
      Quad(g1) + 150*Tr22*Quad(g2) + 825*AbsSqr(MassWB)*Quad(g2) + 30*Tr2U111*
      Sqr(g1) + 40*traceconjTYuTpTYu*Sqr(g1) - 40*MassB*traceconjTYuTpYu*Sqr(g1
      ) + 40*tracemq2AdjYuYu*Sqr(g1) + 40*tracemu2YuAdjYu*Sqr(g1) + 40*mHu2*
      traceYuAdjYu*Sqr(g1) + 80*traceYuAdjYu*AbsSqr(MassB)*Sqr(g1) - 40*
      traceAdjYuTYu*Conj(MassB)*Sqr(g1) + 90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2) + 90
      *AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2) +
      45*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 800*traceconjTYuTpTYu*Sqr(g3) -
      800*MassG*traceconjTYuTpYu*Sqr(g3) + 800*tracemq2AdjYuYu*Sqr(g3) + 800*
      tracemu2YuAdjYu*Sqr(g3) + 800*mHu2*traceYuAdjYu*Sqr(g3) + 1600*
      traceYuAdjYu*AbsSqr(MassG)*Sqr(g3) - 800*traceAdjYuTYu*Conj(MassG)*Sqr(g3
      )));


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

/**
 * Calculates the 4-loop beta function of mHu2.
 *
 * @return 4-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

/**
 * Calculates the 5-loop beta function of mHu2.
 *
 * @return 5-loop beta function
 */
double MSSMRHN_soft_parameters::calc_beta_mHu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
