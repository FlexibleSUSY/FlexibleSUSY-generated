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

// File generated at Wed 12 Apr 2017 11:50:57

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
 * Calculates the one-loop beta function of me2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (oneOver16PiSqr*(4*mHd2*(Ye*Ye.adjoint()) + 4*(TYe*(TYe)
      .adjoint()) + 2*(me2*Ye*Ye.adjoint()) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*
      Ye.adjoint()*me2) + 0.1*(15.491933384829668*g1*Tr11 + 3.1622776601683795*
      gN*Tr14 - 48*AbsSqr(MassB)*Sqr(g1) - 2*AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX
      (3))).real();


   return beta_me2;
}

/**
 * Calculates the two-loop beta function of me2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (twoLoop*(-0.8*(15*traceconjTYdTpTYd + 5*traceconjTYeTpTYe
      + 15*tracemd2YdAdjYd + 5*traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*
      tracemq2AdjYdYd + 30*mHd2*traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*
      AbsSqr(Lambdax) + 5*mHu2*AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*
      AbsSqr(TLambdax) + 3*mHd2*Sqr(g1) + 6*AbsSqr(MassB)*Sqr(g1) - 15*mHd2*Sqr
      (g2) - 30*AbsSqr(MassWB)*Sqr(g2) - 3*mHd2*Sqr(gN) - 6*AbsSqr(MassBp)*Sqr(
      gN))*(Ye*Ye.adjoint()) + 0.8*(-15*traceAdjYdTYd - 5*traceAdjYeTYe + 3*
      MassB*Sqr(g1) - 15*MassWB*Sqr(g2) - 3*MassBp*Sqr(gN) - 5*Conj(Lambdax)*
      TLambdax)*(Ye*(TYe).adjoint()) - 0.8*(-3*Conj(MassB)*Sqr(g1) + 5*(3*
      traceconjTYdTpYd + traceconjTYeTpYe + Conj(TLambdax)*Lambdax + 3*Conj(
      MassWB)*Sqr(g2)) + 3*Conj(MassBp)*Sqr(gN))*(TYe*Ye.adjoint()) - 0.8*(15*
      traceYdAdjYd + 5*traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*Sqr(g1) - 15*Sqr(g2
      ) - 3*Sqr(gN))*(TYe*(TYe).adjoint()) + (-6*traceYdAdjYd - 2*traceYeAdjYe
      - 2*AbsSqr(Lambdax) - 1.2*Sqr(g1) + 6*Sqr(g2) + 1.2*Sqr(gN))*(me2*Ye*
      Ye.adjoint()) - 0.8*(15*traceYdAdjYd + 5*traceYeAdjYe + 5*AbsSqr(Lambdax)
      + 3*Sqr(g1) - 15*Sqr(g2) - 3*Sqr(gN))*(Ye*ml2*Ye.adjoint()) + (-6*
      traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr(Lambdax) - 1.2*Sqr(g1) + 6*Sqr(
      g2) + 1.2*Sqr(gN))*(Ye*Ye.adjoint()*me2) - 8*mHd2*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()) - 4*(Ye*Ye.adjoint()*TYe*(TYe).adjoint()) - 4*(Ye*(TYe)
      .adjoint()*TYe*Ye.adjoint()) - 4*(TYe*Ye.adjoint()*Ye*(TYe).adjoint()) -
      4*(TYe*(TYe).adjoint()*Ye*Ye.adjoint()) - 2*(me2*Ye*Ye.adjoint()*Ye*
      Ye.adjoint()) - 4*(Ye*ml2*Ye.adjoint()*Ye*Ye.adjoint()) - 4*(Ye*
      Ye.adjoint()*me2*Ye*Ye.adjoint()) - 4*(Ye*Ye.adjoint()*Ye*ml2*Ye.adjoint(
      )) - 2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*me2) + 0.01*(20*(2*g1*(
      2.449489742783178*gN*(Tr2U114 + Tr2U141) + 15.491933384829668*Tr31) + gN*
      (gN*Tr2U144 + 6.324555320336759*Tr34) + 24*Tr2U111*Sqr(g1)) + 3*Conj(
      MassBp)*Sqr(gN)*(-8*(MassB + 2*MassBp)*Sqr(g1) + 189*MassBp*Sqr(gN)) + 24
      *Conj(MassB)*Sqr(g1)*(648*MassB*Sqr(g1) - (2*MassB + MassBp)*Sqr(gN)))*
      UNITMATRIX(3))).real();


   return beta_me2;
}

/**
 * Calculates the three-loop beta function of me2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

} // namespace flexiblesusy
