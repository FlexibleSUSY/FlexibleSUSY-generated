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

// File generated at Sun 26 Aug 2018 13:58:25

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
 * Calculates the 1-loop beta function of TYu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYu_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (oneOver16PiSqr*(0.03333333333333333*(-((-90*traceYuAdjYu - 30*
      AbsSqr(Lambdax) + 26*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3) + 9*Sqr(gN))*TYu)
      + 2*Yu*(90*traceAdjYuTYu + 26*MassB*Sqr(g1) + 90*MassWB*Sqr(g2) + 160*
      MassG*Sqr(g3) + 9*MassBp*Sqr(gN) + 30*Conj(Lambdax)*TLambdax)) + 2*(Yu*Yd
      .adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) + TYu*Yd.adjoint()*Yd + 5*(TYu*
      Yu.adjoint()*Yu))).real();


   return beta_TYu;
}

/**
 * Calculates the 2-loop beta function of TYu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYu_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 = TRACE_STRUCT.
      traceAdjLambda12TLambda12;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (twoLoop*(-0.0022222222222222222*Yu*(2700*traceYdAdjYuTYuAdjYd +
      2700*traceYuAdjYdTYdAdjYu + 16200*traceYuAdjYuTYuAdjYu + 15652*MassB*Quad
      (g1) + 29700*MassWB*Quad(g2) + 25600*MassG*Quad(g3) + 5157*MassBp*Quad(gN
      ) - 14400*traceAdjYuTYu*Sqr(g3) + 14400*MassG*traceYuAdjYu*Sqr(g3) + 270*
      traceAdjYuTYu*Sqr(gN) - 270*MassBp*traceYuAdjYu*Sqr(gN) + 480*MassBp*Sqr(
      g3)*Sqr(gN) + 480*MassG*Sqr(g3)*Sqr(gN) + Sqr(g1)*(-720*traceAdjYuTYu +
      720*MassB*traceYuAdjYu + 900*(MassB + MassWB)*Sqr(g2) + 2720*(MassB +
      MassG)*Sqr(g3) + 483*MassB*Sqr(gN) + 483*MassBp*Sqr(gN)) + 225*Sqr(g2)*(
      32*(MassG + MassWB)*Sqr(g3) + 3*(MassBp + MassWB)*Sqr(gN))) + (-3*
      traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu + 8.695555555555556*Quad(g1)
      + 16.5*Quad(g2) + 14.222222222222221*Quad(g3) + 2.865*Quad(gN) + 16*
      traceYuAdjYu*Sqr(g3) + Sqr(g1)*(0.8*traceYuAdjYu + Sqr(g2) +
      3.022222222222222*Sqr(g3) + 0.5366666666666666*Sqr(gN)) + Sqr(g2)*(8*Sqr(
      g3) + 0.75*Sqr(gN)) - 0.3*traceYuAdjYu*Sqr(gN) + 0.5333333333333333*Sqr(
      g3)*Sqr(gN))*TYu - 3*Lambdax*Sqr(Conj(Lambdax))*(Lambdax*TYu + 4*Yu*
      TLambdax) - 0.5*Conj(Lambdax)*(Lambdax*(6*traceKappaAdjKappa + 4*
      traceLambda12AdjLambda12 + 6*traceYdAdjYd + 2*traceYeAdjYe - 3*Sqr(gN))*
      TYu + 2*Yu*(Lambdax*(6*traceAdjKappaTKappa + 4*traceAdjLambda12TLambda12
      + 6*traceAdjYdTYd + 2*traceAdjYeTYe + 3*MassBp*Sqr(gN)) + (6*
      traceKappaAdjKappa + 4*traceLambda12AdjLambda12 + 6*traceYdAdjYd + 2*
      traceYeAdjYe - 3*Sqr(gN))*TLambdax)) - 0.4*(15*traceAdjYdTYd + 5*
      traceAdjYeTYe + 2*MassB*Sqr(g1) + 3*MassBp*Sqr(gN) + 5*Conj(Lambdax)*
      TLambdax)*(Yu*Yd.adjoint()*Yd) + (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*
      AbsSqr(Lambdax) + 0.8*Sqr(g1) + 1.2*Sqr(gN))*(Yu*Yd.adjoint()*TYd) - 0.4*
      (45*traceAdjYuTYu + 2*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 3*MassBp*Sqr(gN
      ) + 15*Conj(Lambdax)*TLambdax)*(Yu*Yu.adjoint()*Yu) + (-12*traceYuAdjYu -
      4*AbsSqr(Lambdax) + 1.2*Sqr(g1) + 6*Sqr(g2) + 0.8*Sqr(gN))*(Yu*Yu.adjoint
      ()*TYu) + (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 0.4*Sqr(g1)
      + 0.6*Sqr(gN))*(TYu*Yd.adjoint()*Yd) + (-15*traceYuAdjYu - 5*AbsSqr(
      Lambdax) + 12*Sqr(g2) + Sqr(gN))*(TYu*Yu.adjoint()*Yu) - 4*(Yu*Yd.adjoint
      ()*Yd*Yd.adjoint()*TYd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*(Yu
      *Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 4*(Yu*Yd.adjoint()*TYd*Yu.adjoint()*
      Yu) - 6*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 8*(Yu*Yu.adjoint()*TYu*Yu
      .adjoint()*Yu) - 2*(TYu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYu*Yd.
      adjoint()*Yd*Yu.adjoint()*Yu) - 6*(TYu*Yu.adjoint()*Yu*Yu.adjoint()*Yu)))
      .real();


   return beta_TYu;
}

/**
 * Calculates the 3-loop beta function of TYu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

/**
 * Calculates the 4-loop beta function of TYu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return beta_TYu;
}

} // namespace flexiblesusy
