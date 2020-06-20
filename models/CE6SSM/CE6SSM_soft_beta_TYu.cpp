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


#include "CE6SSM_soft_parameters.hpp"
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
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYu_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (0.03333333333333333*(180*traceAdjYuTYu*Yu + 52*MassB*Yu*Sqr(g1)
      + 180*MassWB*Yu*Sqr(g2) + 320*MassG*Yu*Sqr(g3) + 18*MassBp*Yu*Sqr(gN) +
      90*traceYuAdjYu*TYu + 30*AbsSqr(Lambdax)*TYu - 26*Sqr(g1)*TYu - 90*Sqr(g2
      )*TYu - 160*Sqr(g3)*TYu - 9*Sqr(gN)*TYu + 60*Yu*Conj(Lambdax)*TLambdax) +
      2*(Yu*Yd.adjoint()*TYd) + 4*(Yu*Yu.adjoint()*TYu) + TYu*Yd.adjoint()*Yd +
      5*(TYu*Yu.adjoint()*Yu)).real();


   return oneLoop * beta_TYu;
}

/**
 * Calculates the 2-loop beta function of TYu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYu_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 = TRACE_STRUCT.
      traceAdjLambda12TLambda12;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = (0.0005555555555555556*(-10800*traceYdAdjYuTYuAdjYd*Yu - 10800*
      traceYuAdjYdTYdAdjYu*Yu - 64800*traceYuAdjYuTYuAdjYu*Yu - 10800*
      traceAdjKappaTKappa*Yu*AbsSqr(Lambdax) - 7200*traceAdjLambda12TLambda12*
      Yu*AbsSqr(Lambdax) - 10800*traceAdjYdTYd*Yu*AbsSqr(Lambdax) - 3600*
      traceAdjYeTYe*Yu*AbsSqr(Lambdax) - 62608*MassB*Yu*Quad(g1) - 118800*
      MassWB*Yu*Quad(g2) - 102400*MassG*Yu*Quad(g3) - 20628*MassBp*Yu*Quad(gN)
      + 2880*traceAdjYuTYu*Yu*Sqr(g1) - 2880*MassB*traceYuAdjYu*Yu*Sqr(g1) -
      3600*MassB*Yu*Sqr(g1)*Sqr(g2) - 3600*MassWB*Yu*Sqr(g1)*Sqr(g2) + 57600*
      traceAdjYuTYu*Yu*Sqr(g3) - 57600*MassG*traceYuAdjYu*Yu*Sqr(g3) - 10880*
      MassB*Yu*Sqr(g1)*Sqr(g3) - 10880*MassG*Yu*Sqr(g1)*Sqr(g3) - 28800*MassG*
      Yu*Sqr(g2)*Sqr(g3) - 28800*MassWB*Yu*Sqr(g2)*Sqr(g3) - 1080*traceAdjYuTYu
      *Yu*Sqr(gN) + 1080*MassBp*traceYuAdjYu*Yu*Sqr(gN) - 5400*MassBp*Yu*AbsSqr
      (Lambdax)*Sqr(gN) - 1932*MassB*Yu*Sqr(g1)*Sqr(gN) - 1932*MassBp*Yu*Sqr(g1
      )*Sqr(gN) - 2700*MassBp*Yu*Sqr(g2)*Sqr(gN) - 2700*MassWB*Yu*Sqr(g2)*Sqr(
      gN) - 1920*MassBp*Yu*Sqr(g3)*Sqr(gN) - 1920*MassG*Yu*Sqr(g3)*Sqr(gN) -
      5400*traceYdAdjYuYuAdjYd*TYu - 16200*traceYuAdjYuYuAdjYu*TYu - 5400*
      traceKappaAdjKappa*AbsSqr(Lambdax)*TYu - 3600*traceLambda12AdjLambda12*
      AbsSqr(Lambdax)*TYu - 5400*traceYdAdjYd*AbsSqr(Lambdax)*TYu - 1800*
      traceYeAdjYe*AbsSqr(Lambdax)*TYu + 15652*Quad(g1)*TYu + 29700*Quad(g2)*
      TYu + 25600*Quad(g3)*TYu + 5157*Quad(gN)*TYu + 1440*traceYuAdjYu*Sqr(g1)*
      TYu + 1800*Sqr(g1)*Sqr(g2)*TYu + 28800*traceYuAdjYu*Sqr(g3)*TYu + 5440*
      Sqr(g1)*Sqr(g3)*TYu + 14400*Sqr(g2)*Sqr(g3)*TYu - 540*traceYuAdjYu*Sqr(gN
      )*TYu + 2700*AbsSqr(Lambdax)*Sqr(gN)*TYu + 966*Sqr(g1)*Sqr(gN)*TYu + 1350
      *Sqr(g2)*Sqr(gN)*TYu + 960*Sqr(g3)*Sqr(gN)*TYu - 5400*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)*TYu - 10800*traceKappaAdjKappa*Yu*Conj(Lambdax)*TLambdax -
      7200*traceLambda12AdjLambda12*Yu*Conj(Lambdax)*TLambdax - 10800*
      traceYdAdjYd*Yu*Conj(Lambdax)*TLambdax - 3600*traceYeAdjYe*Yu*Conj(
      Lambdax)*TLambdax + 5400*Yu*Conj(Lambdax)*Sqr(gN)*TLambdax - 21600*Yu*
      Lambdax*Sqr(Conj(Lambdax))*TLambdax) - 0.4*(15*traceAdjYdTYd + 5*
      traceAdjYeTYe + 2*MassB*Sqr(g1) + 3*MassBp*Sqr(gN) + 5*Conj(Lambdax)*
      TLambdax)*(Yu*Yd.adjoint()*Yd) + 0.4*(-15*traceYdAdjYd - 5*traceYeAdjYe -
      5*AbsSqr(Lambdax) + 2*Sqr(g1) + 3*Sqr(gN))*(Yu*Yd.adjoint()*TYd) - 0.4*(
      45*traceAdjYuTYu + 2*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 3*MassBp*Sqr(gN)
      + 15*Conj(Lambdax)*TLambdax)*(Yu*Yu.adjoint()*Yu) + 0.4*(-30*traceYuAdjYu
       - 10*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(g2) + 2*Sqr(gN))*(Yu*Yu.
      adjoint()*TYu) + 0.2*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(
      Lambdax) + 2*Sqr(g1) + 3*Sqr(gN))*(TYu*Yd.adjoint()*Yd) + (-15*
      traceYuAdjYu - 5*AbsSqr(Lambdax) + 12*Sqr(g2) + Sqr(gN))*(TYu*Yu.adjoint(
      )*Yu) - 4*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 2*(Yu*Yd.adjoint()*Yd*
      Yu.adjoint()*TYu) - 4*(Yu*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 4*(Yu*Yd.
      adjoint()*TYd*Yu.adjoint()*Yu) - 6*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu)
      - 8*(Yu*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 2*(TYu*Yd.adjoint()*Yd*Yd.
      adjoint()*Yd) - 4*(TYu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 6*(TYu*Yu.
      adjoint()*Yu*Yu.adjoint()*Yu)).real();


   return twoLoop * beta_TYu;
}

/**
 * Calculates the 3-loop beta function of TYu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return threeLoop * beta_TYu;
}

/**
 * Calculates the 4-loop beta function of TYu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return fourLoop * beta_TYu;
}

/**
 * Calculates the 5-loop beta function of TYu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_TYu_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYu;

   beta_TYu = ZEROMATRIX(3,3);


   return fiveLoop * beta_TYu;
}

} // namespace flexiblesusy
