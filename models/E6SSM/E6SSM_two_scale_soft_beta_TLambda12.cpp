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

// File generated at Thu 15 Dec 2016 12:50:57

#include "E6SSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of TLambda12.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSM_soft_parameters::calc_beta_TLambda12_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = (oneOver16PiSqr*(0.2*Lambda12*(30*traceAdjKappaTKappa
      + 20*traceAdjLambda12TLambda12 + 6*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) +
      19*MassBp*Sqr(gN)) + (3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12 -
      0.6*Sqr(g1) - 3*Sqr(g2) - 1.9*Sqr(gN))*TLambda12 + 2*Conj(Lambdax)*(2*
      Lambda12*TLambdax + Lambdax*TLambda12) + 3*(Lambda12*(Lambda12).adjoint()
      *TLambda12) + 3*(TLambda12*(Lambda12).adjoint()*Lambda12))).real();


   return beta_TLambda12;
}

/**
 * Calculates the two-loop beta function of TLambda12.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSM_soft_parameters::calc_beta_TLambda12_two_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = (twoLoop*(-4*Lambdax*Sqr(Conj(Lambdax))*(4*Lambda12*
      TLambdax + Lambdax*TLambda12) + 0.005*(-4*Lambda12*(1188*Power(g1,4)*
      MassB + 3933*Power(gN,4)*MassBp + 3300*Power(g2,4)*MassWB + 1200*
      traceKappaAdjKappaTKappaAdjKappa + 800*
      traceLambda12AdjLambda12TLambda12AdjLambda12 - 1600*traceAdjKappaTKappa*
      Sqr(g3) + 1600*MassG*traceKappaAdjKappa*Sqr(g3) + 180*traceAdjKappaTKappa
      *Sqr(gN) + 120*traceAdjLambda12TLambda12*Sqr(gN) - 180*MassBp*
      traceKappaAdjKappa*Sqr(gN) - 120*MassBp*traceLambda12AdjLambda12*Sqr(gN)
      + Sqr(g1)*(40*(-2*traceAdjKappaTKappa - 3*traceAdjLambda12TLambda12 + 2*
      MassB*traceKappaAdjKappa + 3*MassB*traceLambda12AdjLambda12) + 180*(MassB
      + MassWB)*Sqr(g2) + 27*(MassB + MassBp)*Sqr(gN)) + 15*Sqr(g2)*(40*(
      -traceAdjLambda12TLambda12 + MassWB*traceLambda12AdjLambda12) + 13*(
      MassBp + MassWB)*Sqr(gN))) + (1188*Power(g1,4) + 3300*Power(g2,4) + 3933*
      Power(gN,4) - 1200*traceKappaAdjKappaKappaAdjKappa - 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 3200*traceKappaAdjKappa*Sqr
      (g3) - 360*traceKappaAdjKappa*Sqr(gN) - 240*traceLambda12AdjLambda12*Sqr(
      gN) + 30*Sqr(g2)*(40*traceLambda12AdjLambda12 + 13*Sqr(gN)) + 2*Sqr(g1)*(
      80*traceKappaAdjKappa + 120*traceLambda12AdjLambda12 + 180*Sqr(g2) + 27*
      Sqr(gN)))*TLambda12) - 0.4*Conj(Lambdax)*(2*Lambda12*(15*traceYdAdjYd + 5
      *traceYeAdjYe + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2) + 3*Sqr(gN))*
      TLambdax + Lambdax*(2*Lambda12*(3*MassB*Sqr(g1) + 5*(3*traceAdjYdTYd +
      traceAdjYeTYe + 3*traceAdjYuTYu + 3*MassWB*Sqr(g2)) - 3*MassBp*Sqr(gN)) +
      (15*traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr
      (g2) + 3*Sqr(gN))*TLambda12)) + (-12*traceAdjKappaTKappa - 8*
      traceAdjLambda12TLambda12 - 5*MassBp*Sqr(gN) - 8*Conj(Lambdax)*TLambdax)*
      (Lambda12*(Lambda12).adjoint()*Lambda12) + (-9*traceKappaAdjKappa - 6*
      traceLambda12AdjLambda12 - 6*AbsSqr(Lambdax) + 3.5*Sqr(gN))*(Lambda12*(
      Lambda12).adjoint()*TLambda12) + (-9*traceKappaAdjKappa - 6*
      traceLambda12AdjLambda12 - 6*AbsSqr(Lambdax) + 4*Sqr(gN))*(TLambda12*(
      Lambda12).adjoint()*Lambda12) - 3*(Lambda12*(Lambda12).adjoint()*Lambda12
      *(Lambda12).adjoint()*TLambda12) - 4*(Lambda12*(Lambda12).adjoint()*
      TLambda12*(Lambda12).adjoint()*Lambda12) - 3*(TLambda12*(Lambda12)
      .adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12))).real();


   return beta_TLambda12;
}

/**
 * Calculates the three-loop beta function of TLambda12.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSM_soft_parameters::calc_beta_TLambda12_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = ZEROMATRIX(2,2);


   return beta_TLambda12;
}

} // namespace flexiblesusy
