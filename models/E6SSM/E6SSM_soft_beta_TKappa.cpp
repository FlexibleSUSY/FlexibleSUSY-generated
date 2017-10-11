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

// File generated at Tue 10 Oct 2017 21:48:11

#include "E6SSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of TKappa.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = (oneOver16PiSqr*((3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) - 0.26666666666666666*Sqr(g1
      ) - 5.333333333333333*Sqr(g3) - 1.9*Sqr(gN))*TKappa + 0.06666666666666667
      *Kappa*(90*traceAdjKappaTKappa + 60*traceAdjLambda12TLambda12 + 8*MassB*
      Sqr(g1) + 160*MassG*Sqr(g3) + 57*MassBp*Sqr(gN) + 60*Conj(Lambdax)*
      TLambdax) + 3*(Kappa*(Kappa).adjoint()*TKappa) + 3*(TKappa*(Kappa)
      .adjoint()*Kappa))).real();


   return beta_TKappa;
}

/**
 * Calculates the 2-loop beta function of TKappa.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_2_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = (twoLoop*(-0.0022222222222222222*Kappa*(4672*MassB*Quad(
      g1) + 25600*MassG*Quad(g3) + 240*Sqr(g3)*(60*(-traceAdjKappaTKappa +
      MassG*traceKappaAdjKappa) + 13*(MassBp + MassG)*Sqr(gN)) + 9*(3933*MassBp
      *Quad(gN) + 200*(6*traceKappaAdjKappaTKappaAdjKappa + 4*
      traceLambda12AdjLambda12TLambda12AdjLambda12 - 3*(
      traceAdjLambda12TLambda12 - MassWB*traceLambda12AdjLambda12)*Sqr(g2)) +
      60*(3*traceAdjKappaTKappa + 2*traceAdjLambda12TLambda12 - 3*MassBp*
      traceKappaAdjKappa - 2*MassBp*traceLambda12AdjLambda12)*Sqr(gN)) + 4*Sqr(
      g1)*(320*(MassB + MassG)*Sqr(g3) + 3*(30*(-2*traceAdjKappaTKappa - 3*
      traceAdjLambda12TLambda12 + 2*MassB*traceKappaAdjKappa + 3*MassB*
      traceLambda12AdjLambda12) + 19*(MassB + MassBp)*Sqr(gN)))) +
      0.0005555555555555556*(4672*Quad(g1) + 25600*Quad(g3) + 480*Sqr(g3)*(60*
      traceKappaAdjKappa + 13*Sqr(gN)) + 8*Sqr(g1)*(180*traceKappaAdjKappa +
      270*traceLambda12AdjLambda12 + 320*Sqr(g3) + 57*Sqr(gN)) + 9*(3933*Quad(
      gN) + 400*(-3*traceKappaAdjKappaKappaAdjKappa - 2*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 3*traceLambda12AdjLambda12*
      Sqr(g2)) - 120*(3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12)*Sqr(gN
      )))*TKappa - 4*Lambdax*Sqr(Conj(Lambdax))*(Lambdax*TKappa + 4*Kappa*
      TLambdax) - 0.4*Conj(Lambdax)*(Lambdax*(15*traceYdAdjYd + 5*traceYeAdjYe
      + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2) + 3*Sqr(gN))*TKappa + 2*Kappa*
      (Lambdax*(3*MassB*Sqr(g1) + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 3*
      traceAdjYuTYu + 3*MassWB*Sqr(g2)) - 3*MassBp*Sqr(gN)) + (15*traceYdAdjYd
      + 5*traceYeAdjYe + 15*traceYuAdjYu - 3*Sqr(g1) - 15*Sqr(g2) + 3*Sqr(gN))*
      TLambdax)) + (-12*traceAdjKappaTKappa - 8*traceAdjLambda12TLambda12 - 5*
      MassBp*Sqr(gN) - 8*Conj(Lambdax)*TLambdax)*(Kappa*(Kappa).adjoint()*Kappa
      ) + (-9*traceKappaAdjKappa - 6*traceLambda12AdjLambda12 - 6*AbsSqr(
      Lambdax) + 3.5*Sqr(gN))*(Kappa*(Kappa).adjoint()*TKappa) + (-9*
      traceKappaAdjKappa - 6*traceLambda12AdjLambda12 - 6*AbsSqr(Lambdax) + 4*
      Sqr(gN))*(TKappa*(Kappa).adjoint()*Kappa) - 3*(Kappa*(Kappa).adjoint()*
      Kappa*(Kappa).adjoint()*TKappa) - 4*(Kappa*(Kappa).adjoint()*TKappa*(
      Kappa).adjoint()*Kappa) - 3*(TKappa*(Kappa).adjoint()*Kappa*(Kappa)
      .adjoint()*Kappa))).real();


   return beta_TKappa;
}

/**
 * Calculates the 3-loop beta function of TKappa.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = ZEROMATRIX(3,3);


   return beta_TKappa;
}

} // namespace flexiblesusy
