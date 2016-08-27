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

// File generated at Sat 27 Aug 2016 12:44:30

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
 * Calculates the one-loop beta function of ml2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_ml2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + 2*((TYe)
      .adjoint()*TYe) + ml2*Ye.adjoint()*Ye + 2*(Ye.adjoint()*me2*Ye) +
      Ye.adjoint()*Ye*ml2 - 0.7745966692414834*g1*Tr11*UNITMATRIX(3) +
      0.6324555320336759*gN*Tr14*UNITMATRIX(3) - 1.2*AbsSqr(MassB)*Sqr(g1)*
      UNITMATRIX(3) - 6*AbsSqr(MassWB)*Sqr(g2)*UNITMATRIX(3) - 0.8*AbsSqr(
      MassBp)*Sqr(gN)*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the two-loop beta function of ml2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_ml2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (twoLoop*((-6*traceconjTYdTpTYd - 2*traceconjTYeTpTYe - 6*
      tracemd2YdAdjYd - 2*traceme2YeAdjYe - 2*traceml2AdjYeYe - 6*
      tracemq2AdjYdYd - 12*mHd2*traceYdAdjYd - 4*mHd2*traceYeAdjYe - 4*mHd2*
      AbsSqr(Lambdax) - 2*mHu2*AbsSqr(Lambdax) - 2*ms2*AbsSqr(Lambdax) - 2*
      AbsSqr(TLambdax) + 2.4*mHd2*Sqr(g1) + 4.8*AbsSqr(MassB)*Sqr(g1) + 0.6*
      mHd2*Sqr(gN) + 1.2*AbsSqr(MassBp)*Sqr(gN))*(Ye.adjoint()*Ye) + (-6*
      traceconjTYdTpYd - 2*traceconjTYeTpYe - 2*Conj(TLambdax)*Lambdax - 2.4*
      Conj(MassB)*Sqr(g1) - 0.6*Conj(MassBp)*Sqr(gN))*(Ye.adjoint()*TYe) + (-6*
      traceAdjYdTYd - 2*traceAdjYeTYe - 2.4*MassB*Sqr(g1) - 0.6*MassBp*Sqr(gN)
      - 2*Conj(Lambdax)*TLambdax)*((TYe).adjoint()*Ye) + (-6*traceYdAdjYd - 2*
      traceYeAdjYe - 2*AbsSqr(Lambdax) + 2.4*Sqr(g1) + 0.6*Sqr(gN))*((TYe)
      .adjoint()*TYe) + (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 1.2
      *Sqr(g1) + 0.3*Sqr(gN))*(ml2*Ye.adjoint()*Ye) + (-6*traceYdAdjYd - 2*
      traceYeAdjYe - 2*AbsSqr(Lambdax) + 2.4*Sqr(g1) + 0.6*Sqr(gN))*(Ye.adjoint
      ()*me2*Ye) + (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 1.2*Sqr(
      g1) + 0.3*Sqr(gN))*(Ye.adjoint()*Ye*ml2) - 8*mHd2*(Ye.adjoint()*Ye*
      Ye.adjoint()*Ye) - 4*(Ye.adjoint()*Ye*(TYe).adjoint()*TYe) - 4*(
      Ye.adjoint()*TYe*(TYe).adjoint()*Ye) - 4*((TYe).adjoint()*Ye*Ye.adjoint()
      *TYe) - 4*((TYe).adjoint()*TYe*Ye.adjoint()*Ye) - 2*(ml2*Ye.adjoint()*Ye*
      Ye.adjoint()*Ye) - 4*(Ye.adjoint()*me2*Ye*Ye.adjoint()*Ye) - 4*(
      Ye.adjoint()*Ye*ml2*Ye.adjoint()*Ye) - 4*(Ye.adjoint()*Ye*Ye.adjoint()*
      me2*Ye) - 2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye*ml2) + 6*Power(g2,4)*Tr22*
      UNITMATRIX(3) - 0.9797958971132712*g1*gN*Tr2U114*UNITMATRIX(3) -
      0.9797958971132712*g1*gN*Tr2U141*UNITMATRIX(3) - 3.0983866769659336*g1*
      Tr31*UNITMATRIX(3) + 2.5298221281347035*gN*Tr34*UNITMATRIX(3) + 87*Power(
      g2,4)*AbsSqr(MassWB)*UNITMATRIX(3) + 1.2*Tr2U111*Sqr(g1)*UNITMATRIX(3) +
      3.6*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(3) + 1.8*MassB*Conj(MassWB)
      *Sqr(g1)*Sqr(g2)*UNITMATRIX(3) + 0.8*Tr2U144*Sqr(gN)*UNITMATRIX(3) + 2.4*
      AbsSqr(MassWB)*Sqr(g2)*Sqr(gN)*UNITMATRIX(3) + 1.2*MassBp*Conj(MassWB)*
      Sqr(g2)*Sqr(gN)*UNITMATRIX(3) + 0.24*Conj(MassBp)*Sqr(gN)*(3*(MassB + 2*
      MassBp)*Sqr(g1) + 5*(2*MassBp + MassWB)*Sqr(g2) + 96*MassBp*Sqr(gN))*
      UNITMATRIX(3) + 0.36*Conj(MassB)*Sqr(g1)*(99*MassB*Sqr(g1) + 5*(2*MassB +
      MassWB)*Sqr(g2) + 2*(2*MassB + MassBp)*Sqr(gN))*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the three-loop beta function of ml2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_ml2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
