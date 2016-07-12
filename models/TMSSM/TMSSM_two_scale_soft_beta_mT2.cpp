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

// File generated at Tue 12 Jul 2016 10:45:09

#include "TMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mT2.
 *
 * @return one-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mT2_one_loop(const Soft_traces& soft_traces) const
{


   double beta_mT2;

   beta_mT2 = Re(2*oneOver16PiSqr*((mHd2 + mHu2 + mT2)*AbsSqr(Lambdax) +
      AbsSqr(TLambdax) - 8*AbsSqr(MassWB)*Sqr(g2)));


   return beta_mT2;
}

/**
 * Calculates the two-loop beta function of mT2.
 *
 * @return two-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mT2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mT2;

   beta_mT2 = Re(twoLoop*(-12*(mHd2 + mHu2 + mT2)*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) + Conj(Lambdax)*(-6*traceconjTYdTpTYd*Lambdax - 2*
      traceconjTYeTpTYe*Lambdax - 6*traceconjTYuTpTYu*Lambdax - 6*
      tracemd2YdAdjYd*Lambdax - 2*traceme2YeAdjYe*Lambdax - 2*traceml2AdjYeYe*
      Lambdax - 6*tracemq2AdjYdYd*Lambdax - 6*tracemq2AdjYuYu*Lambdax - 6*
      tracemu2YuAdjYu*Lambdax - 12*mHd2*traceYdAdjYd*Lambdax - 6*mHu2*
      traceYdAdjYd*Lambdax - 6*mT2*traceYdAdjYd*Lambdax - 4*mHd2*traceYeAdjYe*
      Lambdax - 2*mHu2*traceYeAdjYe*Lambdax - 2*mT2*traceYeAdjYe*Lambdax - 6*
      mHd2*traceYuAdjYu*Lambdax - 12*mHu2*traceYuAdjYu*Lambdax - 6*mT2*
      traceYuAdjYu*Lambdax - 24*AbsSqr(TLambdax)*Lambdax + 1.2*mHd2*Lambdax*Sqr
      (g1) + 1.2*mHu2*Lambdax*Sqr(g1) + 1.2*mT2*Lambdax*Sqr(g1) - 2*mHd2*
      Lambdax*Sqr(g2) - 2*mHu2*Lambdax*Sqr(g2) - 2*mT2*Lambdax*Sqr(g2) + 1.2*
      Conj(MassB)*Sqr(g1)*(2*MassB*Lambdax - TLambdax) - 6*traceconjTYdTpYd*
      TLambdax - 2*traceconjTYeTpYe*TLambdax - 6*traceconjTYuTpYu*TLambdax) + 2
      *Conj(MassWB)*Sqr(g2)*(152*MassWB*Sqr(g2) + Conj(Lambdax)*(-2*MassWB*
      Lambdax + TLambdax)) - 0.4*(-40*Power(g2,4)*Tr22 + Conj(TLambdax)*(
      Lambdax*(15*traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu + 3*MassB*
      Sqr(g1) - 5*MassWB*Sqr(g2)) + (15*traceYdAdjYd + 5*traceYeAdjYe + 15*
      traceYuAdjYu - 3*Sqr(g1) + 5*Sqr(g2))*TLambdax))));


   return beta_mT2;
}

/**
 * Calculates the three-loop beta function of mT2.
 *
 * @return three-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_mT2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

} // namespace flexiblesusy
