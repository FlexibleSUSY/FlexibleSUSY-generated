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

// File generated at Mon 9 May 2016 12:51:52

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
 * Calculates the one-loop beta function of mq2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_mq2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = (oneOver16PiSqr*(2*mHd2*(Yd.adjoint()*Yd) + 2*mHu2*(
      Yu.adjoint()*Yu) + 2*((TYd).adjoint()*TYd) + 2*((TYu).adjoint()*TYu) +
      mq2*Yd.adjoint()*Yd + mq2*Yu.adjoint()*Yu + 2*(Yd.adjoint()*md2*Yd) +
      Yd.adjoint()*Yd*mq2 + 2*(Yu.adjoint()*mu2*Yu) + Yu.adjoint()*Yu*mq2 +
      0.2581988897471611*g1*Tr11*UNITMATRIX(3) + 0.31622776601683794*gN*Tr14*
      UNITMATRIX(3) - 0.13333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) -
      6*AbsSqr(MassWB)*Sqr(g2)*UNITMATRIX(3) - 10.666666666666666*AbsSqr(MassG
      )*Sqr(g3)*UNITMATRIX(3) - 0.2*AbsSqr(MassBp)*Sqr(gN)*UNITMATRIX(3))).real
      ();


   return beta_mq2;
}

/**
 * Calculates the two-loop beta function of mq2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_mq2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mq2;

   const Eigen::Matrix<double,3,3> beta_mq2_1 = (-0.2*twoLoop*(2*(15*
      traceconjTYdTpTYd + 5*traceconjTYeTpTYe + 15*tracemd2YdAdjYd + 5*
      traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*AbsSqr(Lambdax) + 5*mHu2*
      AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*AbsSqr(TLambdax) - 2*mHd2*Sqr
      (g1) - 4*AbsSqr(MassB)*Sqr(g1) - 3*mHd2*Sqr(gN) - 6*AbsSqr(MassBp)*Sqr(gN
      ))*(Yd.adjoint()*Yd) + 2*(5*(3*traceconjTYdTpYd + traceconjTYeTpYe + Conj
      (TLambdax)*Lambdax) + 2*Conj(MassB)*Sqr(g1) + 3*Conj(MassBp)*Sqr(gN))*(
      Yd.adjoint()*TYd) + 30*traceconjTYuTpTYu*(Yu.adjoint()*Yu) + 30*
      tracemq2AdjYuYu*(Yu.adjoint()*Yu) + 30*tracemu2YuAdjYu*(Yu.adjoint()*Yu)
      + 60*mHu2*traceYuAdjYu*(Yu.adjoint()*Yu) + 10*mHd2*AbsSqr(Lambdax)*(
      Yu.adjoint()*Yu) + 20*mHu2*AbsSqr(Lambdax)*(Yu.adjoint()*Yu) + 10*ms2*
      AbsSqr(Lambdax)*(Yu.adjoint()*Yu) + 10*AbsSqr(TLambdax)*(Yu.adjoint()*Yu)
      - 8*mHu2*Sqr(g1)*(Yu.adjoint()*Yu) - 16*AbsSqr(MassB)*Sqr(g1)*(
      Yu.adjoint()*Yu) - 2*mHu2*Sqr(gN)*(Yu.adjoint()*Yu) - 4*AbsSqr(MassBp)*
      Sqr(gN)*(Yu.adjoint()*Yu) + 30*traceconjTYuTpYu*(Yu.adjoint()*TYu) + 10*
      Conj(TLambdax)*Lambdax*(Yu.adjoint()*TYu) + 8*Conj(MassB)*Sqr(g1)*(
      Yu.adjoint()*TYu) + 2*Conj(MassBp)*Sqr(gN)*(Yu.adjoint()*TYu) + 30*
      traceAdjYdTYd*((TYd).adjoint()*Yd) + 10*traceAdjYeTYe*((TYd).adjoint()*Yd
      ) + 4*MassB*Sqr(g1)*((TYd).adjoint()*Yd) + 6*MassBp*Sqr(gN)*((TYd)
      .adjoint()*Yd) + 30*traceYdAdjYd*((TYd).adjoint()*TYd) + 10*traceYeAdjYe*
      ((TYd).adjoint()*TYd) + 10*AbsSqr(Lambdax)*((TYd).adjoint()*TYd) - 4*Sqr(
      g1)*((TYd).adjoint()*TYd) - 6*Sqr(gN)*((TYd).adjoint()*TYd) + 30*
      traceAdjYuTYu*((TYu).adjoint()*Yu) + 8*MassB*Sqr(g1)*((TYu).adjoint()*Yu)
      + 2*MassBp*Sqr(gN)*((TYu).adjoint()*Yu) + 30*traceYuAdjYu*((TYu).adjoint
      ()*TYu) + 10*AbsSqr(Lambdax)*((TYu).adjoint()*TYu) - 8*Sqr(g1)*((TYu)
      .adjoint()*TYu) - 2*Sqr(gN)*((TYu).adjoint()*TYu) + 15*traceYdAdjYd*(mq2*
      Yd.adjoint()*Yd) + 5*traceYeAdjYe*(mq2*Yd.adjoint()*Yd) + 5*AbsSqr(
      Lambdax)*(mq2*Yd.adjoint()*Yd) - 2*Sqr(g1)*(mq2*Yd.adjoint()*Yd) - 3*Sqr(
      gN)*(mq2*Yd.adjoint()*Yd) + 15*traceYuAdjYu*(mq2*Yu.adjoint()*Yu) + 5*
      AbsSqr(Lambdax)*(mq2*Yu.adjoint()*Yu) - 4*Sqr(g1)*(mq2*Yu.adjoint()*Yu) -
      Sqr(gN)*(mq2*Yu.adjoint()*Yu) + 30*traceYdAdjYd*(Yd.adjoint()*md2*Yd) +
      10*traceYeAdjYe*(Yd.adjoint()*md2*Yd) + 10*AbsSqr(Lambdax)*(Yd.adjoint()*
      md2*Yd) - 4*Sqr(g1)*(Yd.adjoint()*md2*Yd) - 6*Sqr(gN)*(Yd.adjoint()*md2*
      Yd) + 15*traceYdAdjYd*(Yd.adjoint()*Yd*mq2) + 5*traceYeAdjYe*(Yd.adjoint(
      )*Yd*mq2) + 5*AbsSqr(Lambdax)*(Yd.adjoint()*Yd*mq2) - 2*Sqr(g1)*(
      Yd.adjoint()*Yd*mq2) - 3*Sqr(gN)*(Yd.adjoint()*Yd*mq2) + 30*traceYuAdjYu*
      (Yu.adjoint()*mu2*Yu) + 10*AbsSqr(Lambdax)*(Yu.adjoint()*mu2*Yu) - 8*Sqr(
      g1)*(Yu.adjoint()*mu2*Yu) - 2*Sqr(gN)*(Yu.adjoint()*mu2*Yu) + 15*
      traceYuAdjYu*(Yu.adjoint()*Yu*mq2) + 5*AbsSqr(Lambdax)*(Yu.adjoint()*Yu*
      mq2) - 4*Sqr(g1)*(Yu.adjoint()*Yu*mq2) - Sqr(gN)*(Yu.adjoint()*Yu*mq2) +
      40*mHd2*(Yd.adjoint()*Yd*Yd.adjoint()*Yd) + 20*(Yd.adjoint()*Yd*(TYd)
      .adjoint()*TYd) + 20*(Yd.adjoint()*TYd*(TYd).adjoint()*Yd) + 40*mHu2*(
      Yu.adjoint()*Yu*Yu.adjoint()*Yu) + 20*(Yu.adjoint()*Yu*(TYu).adjoint()*
      TYu) + 20*(Yu.adjoint()*TYu*(TYu).adjoint()*Yu) + 20*((TYd).adjoint()*Yd*
      Yd.adjoint()*TYd) + 20*((TYd).adjoint()*TYd*Yd.adjoint()*Yd) + 20*((TYu)
      .adjoint()*Yu*Yu.adjoint()*TYu) + 20*((TYu).adjoint()*TYu*Yu.adjoint()*Yu
      ) + 10*(mq2*Yd.adjoint()*Yd*Yd.adjoint()*Yd) + 10*(mq2*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu) + 20*(Yd.adjoint()*md2*Yd*Yd.adjoint()*Yd) + 20*(
      Yd.adjoint()*Yd*mq2*Yd.adjoint()*Yd) + 20*(Yd.adjoint()*Yd*Yd.adjoint()*
      md2*Yd) + 10*(Yd.adjoint()*Yd*Yd.adjoint()*Yd*mq2) + 20*(Yu.adjoint()*mu2
      *Yu*Yu.adjoint()*Yu) + 20*(Yu.adjoint()*Yu*mq2*Yu.adjoint()*Yu) + 20*(
      Yu.adjoint()*Yu*Yu.adjoint()*mu2*Yu) + 10*(Yu.adjoint()*Yu*Yu.adjoint()*
      Yu*mq2))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_mq2_2 = (-0.0011111111111111111*
      twoLoop*UNITMATRIX(3)*(1800*Conj(Lambdax)*TLambdax*((TYd).adjoint()*Yd +
      (TYu).adjoint()*Yu) - (2*Conj(MassB)*Sqr(g1)*(1734*MassB*Sqr(g1) + 90*(2*
      MassB + MassWB)*Sqr(g2) + 320*MassB*Sqr(g3) + 160*MassG*Sqr(g3) - 66*
      MassB*Sqr(gN) - 33*MassBp*Sqr(gN)) + 3*Conj(MassBp)*Sqr(gN)*(-22*(MassB +
      2*MassBp)*Sqr(g1) + 90*(2*MassBp + MassWB)*Sqr(g2) + 320*MassBp*Sqr(g3)
      + 160*MassG*Sqr(g3) + 1701*MassBp*Sqr(gN)) + 10*(16*Conj(MassG)*Sqr(g3)*(
      2*(MassB + 2*MassG)*Sqr(g1) + 3*(10*(3*(2*MassG + MassWB)*Sqr(g2) + 10*
      MassG*Sqr(g3)) + (MassBp + 2*MassG)*Sqr(gN))) + 3*(3*Conj(MassWB)*Sqr(g2)
      *(2*(MassB + 2*MassWB)*Sqr(g1) + 10*(87*MassWB*Sqr(g2) + 16*(MassG + 2*
      MassWB)*Sqr(g3)) + 3*(MassBp + 2*MassWB)*Sqr(gN)) + 2*(90*Power(g2,4)*
      Tr22 + 160*Power(g3,4)*Tr23 + 2.449489742783178*g1*gN*Tr2U114 +
      2.449489742783178*g1*gN*Tr2U141 + 15.491933384829668*g1*Tr31 +
      18.973665961010276*gN*Tr34 + 2*Tr2U111*Sqr(g1) + 3*Tr2U144*Sqr(gN)))))*
      UNITMATRIX(3))).real();

   beta_mq2 = beta_mq2_1 + beta_mq2_2;


   return beta_mq2;
}

/**
 * Calculates the three-loop beta function of mq2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_mq2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

} // namespace flexiblesusy
