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

// File generated at Thu 15 Dec 2016 12:42:54

#include "E6SSMtower_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;


   double beta_TLambdax;

   beta_TLambdax = Re(oneOver16PiSqr*(0.2*Lambdax*(6*MassB*Sqr(g1) + 10*(
      3*traceAdjKappaTKappa + 2*traceAdjLambda12TLambda12 + 3*traceAdjYdTYd +
      traceAdjYeTYe + 3*traceAdjYuTYu + 3*MassWB*Sqr(g2)) + 19*MassBp*Sqr(gN))
      + (3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + 12*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(
      g2) - 1.9*Sqr(gN))*TLambdax));


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;


   double beta_TLambdax;

   beta_TLambdax = Re(0.005*twoLoop*(-4*Lambdax*(1188*Power(g1,4)*MassB +
      3933*Power(gN,4)*MassBp + 3300*Power(g2,4)*MassWB + 1200*
      traceKappaAdjKappaTKappaAdjKappa + 800*
      traceLambda12AdjLambda12TLambda12AdjLambda12 + 1800*traceYdAdjYdTYdAdjYd
      + 600*traceYdAdjYuTYuAdjYd + 600*traceYeAdjYeTYeAdjYe + 600*
      traceYuAdjYdTYdAdjYu + 1800*traceYuAdjYuTYuAdjYu - 1600*
      traceAdjKappaTKappa*Sqr(g3) - 1600*traceAdjYdTYd*Sqr(g3) - 1600*
      traceAdjYuTYu*Sqr(g3) + 1600*MassG*traceKappaAdjKappa*Sqr(g3) + 1600*
      MassG*traceYdAdjYd*Sqr(g3) + 1600*MassG*traceYuAdjYu*Sqr(g3) + 180*
      traceAdjKappaTKappa*Sqr(gN) + 120*traceAdjLambda12TLambda12*Sqr(gN) + 60*
      traceAdjYdTYd*Sqr(gN) + 20*traceAdjYeTYe*Sqr(gN) + 30*traceAdjYuTYu*Sqr(
      gN) - 180*MassBp*traceKappaAdjKappa*Sqr(gN) - 120*MassBp*
      traceLambda12AdjLambda12*Sqr(gN) - 60*MassBp*traceYdAdjYd*Sqr(gN) - 20*
      MassBp*traceYeAdjYe*Sqr(gN) - 30*MassBp*traceYuAdjYu*Sqr(gN) + Sqr(g1)*(
      40*(-2*traceAdjKappaTKappa - 3*traceAdjLambda12TLambda12 + traceAdjYdTYd
      - 3*traceAdjYeTYe - 2*traceAdjYuTYu + 2*MassB*traceKappaAdjKappa + 3*
      MassB*traceLambda12AdjLambda12 - MassB*traceYdAdjYd + 3*MassB*
      traceYeAdjYe + 2*MassB*traceYuAdjYu) + 180*(MassB + MassWB)*Sqr(g2) + 27*
      (MassB + MassBp)*Sqr(gN)) + 15*Sqr(g2)*(40*(-traceAdjLambda12TLambda12 +
      MassWB*traceLambda12AdjLambda12) + 13*(MassBp + MassWB)*Sqr(gN))) + (1188
      *Power(g1,4) + 3300*Power(g2,4) + 3933*Power(gN,4) - 1200*
      traceKappaAdjKappaKappaAdjKappa - 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 1800*traceYdAdjYdYdAdjYd -
      1200*traceYdAdjYuYuAdjYd - 600*traceYeAdjYeYeAdjYe - 1800*
      traceYuAdjYuYuAdjYu + 3200*traceKappaAdjKappa*Sqr(g3) + 3200*traceYdAdjYd
      *Sqr(g3) + 3200*traceYuAdjYu*Sqr(g3) - 360*traceKappaAdjKappa*Sqr(gN) -
      240*traceLambda12AdjLambda12*Sqr(gN) - 120*traceYdAdjYd*Sqr(gN) - 40*
      traceYeAdjYe*Sqr(gN) - 60*traceYuAdjYu*Sqr(gN) + 30*Sqr(g2)*(40*
      traceLambda12AdjLambda12 + 13*Sqr(gN)) + 2*Sqr(g1)*(40*(2*
      traceKappaAdjKappa + 3*traceLambda12AdjLambda12 - traceYdAdjYd + 3*
      traceYeAdjYe + 2*traceYuAdjYu) + 180*Sqr(g2) + 27*Sqr(gN)))*TLambdax -
      10000*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TLambdax - 20*AbsSqr(Lambdax)*(2*
      Lambdax*(12*MassB*Sqr(g1) + 10*(6*traceAdjKappaTKappa + 4*
      traceAdjLambda12TLambda12 + 9*traceAdjYdTYd + 3*traceAdjYeTYe + 9*
      traceAdjYuTYu + 6*MassWB*Sqr(g2)) + 13*MassBp*Sqr(gN)) - 3*(-60*
      traceKappaAdjKappa - 40*traceLambda12AdjLambda12 - 90*traceYdAdjYd - 30*
      traceYeAdjYe - 90*traceYuAdjYu + 12*Sqr(g1) + 60*Sqr(g2) + 13*Sqr(gN))*
      TLambdax)));


   return beta_TLambdax;
}

/**
 * Calculates the three-loop beta function of TLambdax.
 *
 * @return three-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_TLambdax_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

} // namespace flexiblesusy
