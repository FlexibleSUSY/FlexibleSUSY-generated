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

// File generated at Mon 5 Mar 2018 17:57:57

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
 * Calculates the 1-loop beta function of TYe.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYe_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*((3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 1.8*Sqr(g1) - 3*Sqr(g2) - 0.7*Sqr(gN))*TYe + 0.2*Ye*(30*
      traceAdjYdTYd + 10*traceAdjYeTYe + 18*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) +
      7*MassBp*Sqr(gN) + 10*Conj(Lambdax)*TLambdax) + 4*(Ye*Ye.adjoint()*TYe)
      + 5*(TYe*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the 2-loop beta function of TYe.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYe_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(0.025*(-4*Ye*(360*traceYdAdjYdTYdAdjYd + 60*
      traceYdAdjYuTYuAdjYd + 120*traceYeAdjYeTYeAdjYe + 60*traceYuAdjYdTYdAdjYu
      + 756*MassB*Quad(g1) + 660*MassWB*Quad(g2) + 273*MassBp*Quad(gN) - 320*
      traceAdjYdTYd*Sqr(g3) + 320*MassG*traceYdAdjYd*Sqr(g3) + 12*traceAdjYdTYd
      *Sqr(gN) + 4*traceAdjYeTYe*Sqr(gN) - 12*MassBp*traceYdAdjYd*Sqr(gN) - 4*
      MassBp*traceYeAdjYe*Sqr(gN) + 39*(MassBp + MassWB)*Sqr(g2)*Sqr(gN) + Sqr(
      g1)*(8*(traceAdjYdTYd - 3*traceAdjYeTYe - MassB*traceYdAdjYd + 3*MassB*
      traceYeAdjYe) + 36*(MassB + MassWB)*Sqr(g2) + 3*(MassB + MassBp)*Sqr(gN))
      ) + (-360*traceYdAdjYdYdAdjYd - 120*traceYdAdjYuYuAdjYd - 120*
      traceYeAdjYeYeAdjYe + 756*Quad(g1) + 660*Quad(g2) + 273*Quad(gN) + 640*
      traceYdAdjYd*Sqr(g3) - 24*traceYdAdjYd*Sqr(gN) - 8*traceYeAdjYe*Sqr(gN) +
      78*Sqr(g2)*Sqr(gN) + 2*Sqr(g1)*(-8*traceYdAdjYd + 24*traceYeAdjYe + 36*
      Sqr(g2) + 3*Sqr(gN)))*TYe - 120*Lambdax*Sqr(Conj(Lambdax))*(Lambdax*TYe +
      4*Ye*TLambdax) - 40*Conj(Lambdax)*(Lambdax*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 3*traceYuAdjYu - Sqr(gN))*TYe + 2*Ye*(Lambdax*
      (3*traceAdjKappaTKappa + 2*traceAdjLambda12TLambda12 + 3*traceAdjYuTYu +
      MassBp*Sqr(gN)) + (3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*
      traceYuAdjYu - Sqr(gN))*TLambdax))) - 3*(6*traceAdjYdTYd + 2*
      traceAdjYeTYe + 4*MassWB*Sqr(g2) + MassBp*Sqr(gN) + 2*Conj(Lambdax)*
      TLambdax)*(Ye*Ye.adjoint()*Ye) + (-12*traceYdAdjYd - 4*traceYeAdjYe - 4*
      AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*Sqr(g2) + 1.8*Sqr(gN))*(Ye*Ye.adjoint()
      *TYe) + (-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) - 1.2*Sqr(
      g1) + 12*Sqr(g2) + 2.7*Sqr(gN))*(TYe*Ye.adjoint()*Ye) - 6*(Ye*Ye.adjoint(
      )*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 6*(TYe
      *Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the 3-loop beta function of TYe.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYe_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 4-loop beta function of TYe.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYe_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
