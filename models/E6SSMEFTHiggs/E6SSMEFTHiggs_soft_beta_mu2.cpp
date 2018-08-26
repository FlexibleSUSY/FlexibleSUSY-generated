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

// File generated at Sun 26 Aug 2018 13:58:59

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
 * Calculates the 1-loop beta function of mu2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mu2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (oneOver16PiSqr*(4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu).adjoint(
      )) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*Yu.adjoint
      ()*mu2) - 0.03333333333333333*(30.983866769659336*g1*Tr11 -
      9.486832980505138*gN*Tr14 + 64*AbsSqr(MassB)*Sqr(g1) + 320*AbsSqr(MassG)*
      Sqr(g3) + 6*AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the 2-loop beta function of mu2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (twoLoop*(-0.8*(15*traceconjTYuTpTYu + 15*tracemq2AdjYuYu + 15*
      tracemu2YuAdjYu + 30*mHu2*traceYuAdjYu + 5*mHd2*AbsSqr(Lambdax) + 10*mHu2
      *AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*AbsSqr(TLambdax) + mHu2*Sqr(
      g1) + 2*AbsSqr(MassB)*Sqr(g1) - 15*mHu2*Sqr(g2) - 30*AbsSqr(MassWB)*Sqr(
      g2) - mHu2*Sqr(gN) - 2*AbsSqr(MassBp)*Sqr(gN))*(Yu*Yu.adjoint()) + 0.8*(-
      15*traceAdjYuTYu + MassB*Sqr(g1) - 15*MassWB*Sqr(g2) - MassBp*Sqr(gN) - 5
      *Conj(Lambdax)*TLambdax)*(Yu*(TYu).adjoint()) - 0.8*(-(Conj(MassB)*Sqr(g1
      )) + 5*(3*traceconjTYuTpYu + Conj(TLambdax)*Lambdax + 3*Conj(MassWB)*Sqr(
      g2)) + Conj(MassBp)*Sqr(gN))*(TYu*Yu.adjoint()) - 0.8*(15*traceYuAdjYu +
      5*AbsSqr(Lambdax) + Sqr(g1) - 15*Sqr(g2) - Sqr(gN))*(TYu*(TYu).adjoint())
      - 0.4*(15*traceYuAdjYu + 5*AbsSqr(Lambdax) + Sqr(g1) - 15*Sqr(g2) - Sqr(
      gN))*(mu2*Yu*Yu.adjoint()) - 0.8*(15*traceYuAdjYu + 5*AbsSqr(Lambdax) +
      Sqr(g1) - 15*Sqr(g2) - Sqr(gN))*(Yu*mq2*Yu.adjoint()) - 0.4*(15*
      traceYuAdjYu + 5*AbsSqr(Lambdax) + Sqr(g1) - 15*Sqr(g2) - Sqr(gN))*(Yu*Yu
      .adjoint()*mu2) - 4*(mHd2 + mHu2)*(Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 4*(
      Yu*Yd.adjoint()*TYd*(TYu).adjoint()) - 8*mHu2*(Yu*Yu.adjoint()*Yu*Yu.
      adjoint()) - 4*(Yu*Yu.adjoint()*TYu*(TYu).adjoint()) - 4*(Yu*(TYd).
      adjoint()*TYd*Yu.adjoint()) - 4*(Yu*(TYu).adjoint()*TYu*Yu.adjoint()) - 4
      *(TYu*Yd.adjoint()*Yd*(TYu).adjoint()) - 4*(TYu*Yu.adjoint()*Yu*(TYu).
      adjoint()) - 4*(TYu*(TYd).adjoint()*Yd*Yu.adjoint()) - 4*(TYu*(TYu).
      adjoint()*Yu*Yu.adjoint()) - 2*(mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 2*
      (mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*mq2*Yd.adjoint()*Yd*Yu.
      adjoint()) - 4*(Yu*mq2*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*Yd.adjoint()
      *md2*Yd*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*Yd*mq2*Yu.adjoint()) - 2*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*mu2) - 4*(Yu*Yu.adjoint()*mu2*Yu*Yu.adjoint(
      )) - 4*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint()) - 2*(Yu*Yu.adjoint()*Yu*Yu.
      adjoint()*mu2) + 0.0011111111111111111*(3*Conj(MassBp)*Sqr(gN)*(128*(
      MassB + 2*MassBp)*Sqr(g1) + 160*(2*MassBp + MassG)*Sqr(g3) + 1701*MassBp*
      Sqr(gN)) + 128*Conj(MassB)*Sqr(g1)*(456*MassB*Sqr(g1) + 40*(2*MassB +
      MassG)*Sqr(g3) + 3*(2*MassB + MassBp)*Sqr(gN)) + 20*(3*(-4*g1*(
      2.449489742783178*gN*(Tr2U114 + Tr2U141) + 15.491933384829668*Tr31) + 3*
      gN*(gN*Tr2U144 + 6.324555320336759*Tr34) + 160*Tr23*Quad(g3) + 32*Tr2U111
      *Sqr(g1)) + 8*Conj(MassG)*Sqr(g3)*(32*(MassB + 2*MassG)*Sqr(g1) + 300*
      MassG*Sqr(g3) + 3*(MassBp + 2*MassG)*Sqr(gN))))*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the 3-loop beta function of mu2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

/**
 * Calculates the 4-loop beta function of mu2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

} // namespace flexiblesusy
