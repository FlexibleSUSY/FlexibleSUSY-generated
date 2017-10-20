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

// File generated at Fri 20 Oct 2017 08:52:15

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
 * Calculates the 1-loop beta function of md2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_md2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd)
      .adjoint()) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*
      Yd.adjoint()*md2) + 0.06666666666666667*(7.745966692414834*g1*Tr11 +
      9.486832980505138*gN*Tr14 - 8*AbsSqr(MassB)*Sqr(g1) - 160*AbsSqr(MassG)*
      Sqr(g3) - 12*AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the 2-loop beta function of md2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_md2_2_loop(const Soft_traces& soft_traces) const
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
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (twoLoop*((-12*traceconjTYdTpTYd - 4*traceconjTYeTpTYe - 12
      *tracemd2YdAdjYd - 4*traceme2YeAdjYe - 4*traceml2AdjYeYe - 12*
      tracemq2AdjYdYd - 24*mHd2*traceYdAdjYd - 8*mHd2*traceYeAdjYe - 8*mHd2*
      AbsSqr(Lambdax) - 4*mHu2*AbsSqr(Lambdax) - 4*ms2*AbsSqr(Lambdax) - 4*
      AbsSqr(TLambdax) + 0.8*mHd2*Sqr(g1) + 1.6*AbsSqr(MassB)*Sqr(g1) + 12*mHd2
      *Sqr(g2) + 24*AbsSqr(MassWB)*Sqr(g2) + 1.2*mHd2*Sqr(gN) + 2.4*AbsSqr(
      MassBp)*Sqr(gN))*(Yd*Yd.adjoint()) - 0.4*(30*traceAdjYdTYd + 10*
      traceAdjYeTYe + 2*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 3*MassBp*Sqr(gN) +
      10*Conj(Lambdax)*TLambdax)*(Yd*(TYd).adjoint()) - 0.4*(2*Conj(MassB)*Sqr(
      g1) + 10*(3*traceconjTYdTpYd + traceconjTYeTpYe + Conj(TLambdax)*Lambdax
      + 3*Conj(MassWB)*Sqr(g2)) + 3*Conj(MassBp)*Sqr(gN))*(TYd*Yd.adjoint()) +
      (-12*traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(Lambdax) + 0.8*Sqr(g1) + 12
      *Sqr(g2) + 1.2*Sqr(gN))*(TYd*(TYd).adjoint()) + (-6*traceYdAdjYd - 2*
      traceYeAdjYe - 2*AbsSqr(Lambdax) + 0.4*Sqr(g1) + 6*Sqr(g2) + 0.6*Sqr(gN))
      *(md2*Yd*Yd.adjoint()) + (-12*traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(
      Lambdax) + 0.8*Sqr(g1) + 12*Sqr(g2) + 1.2*Sqr(gN))*(Yd*mq2*Yd.adjoint())
      + (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr(Lambdax) + 0.4*Sqr(g1) + 6
      *Sqr(g2) + 0.6*Sqr(gN))*(Yd*Yd.adjoint()*md2) - 8*mHd2*(Yd*Yd.adjoint()*
      Yd*Yd.adjoint()) - 4*(Yd*Yd.adjoint()*TYd*(TYd).adjoint()) - 4*(mHd2 +
      mHu2)*(Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*TYu*(TYd)
      .adjoint()) - 4*(Yd*(TYd).adjoint()*TYd*Yd.adjoint()) - 4*(Yd*(TYu)
      .adjoint()*TYu*Yd.adjoint()) - 4*(TYd*Yd.adjoint()*Yd*(TYd).adjoint()) -
      4*(TYd*Yu.adjoint()*Yu*(TYd).adjoint()) - 4*(TYd*(TYd).adjoint()*Yd*
      Yd.adjoint()) - 4*(TYd*(TYu).adjoint()*Yu*Yd.adjoint()) - 2*(md2*Yd*
      Yd.adjoint()*Yd*Yd.adjoint()) - 2*(md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()) -
      4*(Yd*mq2*Yd.adjoint()*Yd*Yd.adjoint()) - 4*(Yd*mq2*Yu.adjoint()*Yu*
      Yd.adjoint()) - 4*(Yd*Yd.adjoint()*md2*Yd*Yd.adjoint()) - 4*(Yd*
      Yd.adjoint()*Yd*mq2*Yd.adjoint()) - 2*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*
      md2) - 4*(Yd*Yu.adjoint()*mu2*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*Yu*
      mq2*Yd.adjoint()) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*md2) +
      0.017777777777777778*(4*Conj(MassB)*Sqr(g1)*(219*MassB*Sqr(g1) + 20*(2*
      MassB + MassG)*Sqr(g3) - 3*(2*MassB + MassBp)*Sqr(gN)) + 12*Conj(MassBp)*
      Sqr(gN)*(-((MassB + 2*MassBp)*Sqr(g1)) + 2*(5*(2*MassBp + MassG)*Sqr(g3)
      + 54*MassBp*Sqr(gN))) + 5*(3*(g1*(2.449489742783178*gN*(Tr2U114 + Tr2U141
      ) + 7.745966692414834*Tr31) + 3*gN*(gN*Tr2U144 + 3.1622776601683795*Tr34)
      + 40*Tr23*Quad(g3) + 2*Tr2U111*Sqr(g1)) + 8*Conj(MassG)*Sqr(g3)*(2*(
      MassB + 2*MassG)*Sqr(g1) + 75*MassG*Sqr(g3) + 3*(MassBp + 2*MassG)*Sqr(gN
      ))))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the 3-loop beta function of md2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_md2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
