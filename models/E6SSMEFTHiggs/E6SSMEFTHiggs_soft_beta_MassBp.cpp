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

// File generated at Mon 5 Mar 2018 16:36:19

#include "E6SSMEFTHiggs_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of MassBp.
 *
 * @return 1-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_MassBp_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassBp;

   beta_MassBp = Re(18.8*MassBp*oneOver16PiSqr*Sqr(gN));


   return beta_MassBp;
}

/**
 * Calculates the 2-loop beta function of MassBp.
 *
 * @return 2-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_MassBp_2_loop(const Soft_traces& soft_traces) const
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


   double beta_MassBp;

   beta_MassBp = Re(0.04*twoLoop*Sqr(gN)*(285*traceAdjKappaTKappa + 190*
      traceAdjLambda12TLambda12 + 210*traceAdjYdTYd + 70*traceAdjYeTYe + 90*
      traceAdjYuTYu - 285*MassBp*traceKappaAdjKappa - 190*MassBp*
      traceLambda12AdjLambda12 - 210*MassBp*traceYdAdjYd - 70*MassBp*
      traceYeAdjYe - 90*MassBp*traceYuAdjYu + 162*(MassB + MassBp)*Sqr(g1) +
      510*MassBp*Sqr(g2) + 510*MassWB*Sqr(g2) + 1200*MassBp*Sqr(g3) + 1200*
      MassG*Sqr(g3) + 916*MassBp*Sqr(gN) - 190*Conj(Lambdax)*(MassBp*Lambdax -
      TLambdax)));


   return beta_MassBp;
}

/**
 * Calculates the 3-loop beta function of MassBp.
 *
 * @return 3-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_MassBp_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassBp;

   beta_MassBp = 0;


   return beta_MassBp;
}

/**
 * Calculates the 4-loop beta function of MassBp.
 *
 * @return 4-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_MassBp_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassBp;

   beta_MassBp = 0;


   return beta_MassBp;
}

} // namespace flexiblesusy
