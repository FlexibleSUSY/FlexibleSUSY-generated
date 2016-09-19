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

// File generated at Mon 19 Sep 2016 10:14:46

#include "MSSMNoFV_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMNoFV_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*(Yd*(6*traceAdjYdTYd + 2*traceAdjYeTYe +
      0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) + 10.666666666666666*
      MassG*Sqr(g3)) + 3*traceYdAdjYd*TYd + traceYeAdjYe*TYd -
      0.4666666666666667*Sqr(g1)*TYd - 3*Sqr(g2)*TYd - 5.333333333333333*Sqr(g3
      )*TYd + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*
      Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu)).real();


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMNoFV_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(-0.044444444444444446*Yd*(287*Power(g1,4)*MassB -
      160*Power(g3,4)*MassG + 675*Power(g2,4)*MassWB + 810*
      traceYdAdjYdTYdAdjYd + 135*traceYdAdjYuTYuAdjYd + 270*
      traceYeAdjYeTYeAdjYe + 135*traceYuAdjYdTYdAdjYu + 18*traceAdjYdTYd*Sqr(g1
      ) - 54*traceAdjYeTYe*Sqr(g1) + 54*MassB*traceYeAdjYe*Sqr(g1) + 45*MassB*
      Sqr(g1)*Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(g2) - 720*traceAdjYdTYd*Sqr(g3) +
      40*MassB*Sqr(g1)*Sqr(g3) + 40*MassG*Sqr(g1)*Sqr(g3) + 360*MassG*Sqr(g2)*
      Sqr(g3) + 360*MassWB*Sqr(g2)*Sqr(g3) - 18*traceYdAdjYd*(MassB*Sqr(g1) -
      40*MassG*Sqr(g3))) + 3.188888888888889*Power(g1,4)*TYd + 7.5*Power(g2,4)*
      TYd - 1.7777777777777777*Power(g3,4)*TYd - 9*traceYdAdjYdYdAdjYd*TYd - 3*
      traceYdAdjYuYuAdjYd*TYd - 3*traceYeAdjYeYeAdjYe*TYd - 0.4*traceYdAdjYd*
      Sqr(g1)*TYd + 1.2*traceYeAdjYe*Sqr(g1)*TYd + Sqr(g1)*Sqr(g2)*TYd + 16*
      traceYdAdjYd*Sqr(g3)*TYd + 0.8888888888888888*Sqr(g1)*Sqr(g3)*TYd + 8*Sqr
      (g2)*Sqr(g3)*TYd - 0.4*(45*traceAdjYdTYd + 15*traceAdjYeTYe + 4*MassB*Sqr
      (g1) + 30*MassWB*Sqr(g2))*(Yd*Yd.adjoint()*Yd) + (-12*traceYdAdjYd - 4*
      traceYeAdjYe + 1.2*Sqr(g1) + 6*Sqr(g2))*(Yd*Yd.adjoint()*TYd) + (-6*
      traceAdjYuTYu - 1.6*MassB*Sqr(g1))*(Yd*Yu.adjoint()*Yu) + (-6*
      traceYuAdjYu + 1.6*Sqr(g1))*(Yd*Yu.adjoint()*TYu) + (-15*traceYdAdjYd - 5
      *traceYeAdjYe + 1.2*Sqr(g1) + 12*Sqr(g2))*(TYd*Yd.adjoint()*Yd) + (-3*
      traceYuAdjYu + 0.8*Sqr(g1))*(TYd*Yu.adjoint()*Yu) - 6*(Yd*Yd.adjoint()*Yd
      *Yd.adjoint()*TYd) - 8*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*TYd) - 4*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*
      TYu) - 4*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*(Yd*Yu.adjoint()*TYu*
      Yu.adjoint()*Yu) - 6*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(TYd*
      Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu
      ))).real();


   return beta_TYd;
}

/**
 * Calculates the three-loop beta function of TYd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSMNoFV_soft_parameters::calc_beta_TYd_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceTYdAdjYd = TRACE_STRUCT.traceTYdAdjYd;
   const double traceTYeAdjYe = TRACE_STRUCT.traceTYeAdjYe;
   const double traceTYuAdjYu = TRACE_STRUCT.traceTYuAdjYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYdTYdAdjYdYd = TRACE_STRUCT.traceAdjYdTYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYeTYeAdjYeYe = TRACE_STRUCT.traceAdjYeTYeAdjYeYe;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;
   const double traceAdjYuTYuAdjYdYd = TRACE_STRUCT.traceAdjYuTYuAdjYdYd;
   const double traceAdjYuTYuAdjYuYu = TRACE_STRUCT.traceAdjYuTYuAdjYuYu;
   const double traceTYdAdjYuYuAdjYd = TRACE_STRUCT.traceTYdAdjYuYuAdjYd;
   const double traceAdjYdYdAdjYdYdAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYdYdAdjYdYd;
   const double traceAdjYdYdAdjYdTYdAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYdTYdAdjYdYd;
   const double traceAdjYeYeAdjYeYeAdjYeYe =
      TRACE_STRUCT.traceAdjYeYeAdjYeYeAdjYeYe;
   const double traceAdjYeYeAdjYeTYeAdjYeYe =
      TRACE_STRUCT.traceAdjYeYeAdjYeTYeAdjYeYe;
   const double traceAdjYuYuAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuYuAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuTYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuYuAdjYuTYuAdjYdYd;
   const double traceAdjYuTYuAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuTYuAdjYuYuAdjYdYd;
   const double traceTYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceTYdAdjYuYuAdjYuYuAdjYd;


   Eigen::Matrix<double,3,3> beta_TYd;

   const Eigen::Matrix<double,3,3> beta_TYd_1 = ((-0.00044444444444444447
      *threeLoop*Yd*(389302*Power(g1,6)*MassB + 5*Power(g1,4)*(-15*(1890*MassB*
      traceAdjYdYd + 3798*MassB*traceAdjYeYe + 728*MassB*traceAdjYuYu - 945*
      traceTYdAdjYd - 1899*traceTYeAdjYe - 364*traceTYuAdjYu) + 1962*(2*MassB +
      MassWB)*Sqr(g2) + 15568*(2*MassB + MassG)*Sqr(g3)) + 50*Sqr(g1)*(765*
      Power(g2,4)*(MassB + 2*MassWB) - 9*Sqr(g2)*(3*(MassB*traceAdjYdYd +
      MassWB*traceAdjYdYd + 27*MassB*traceAdjYeYe + 27*MassWB*traceAdjYeYe -
      traceTYdAdjYd - 27*traceTYeAdjYe) + 16*(MassB + MassG + MassWB)*Sqr(g3))
      + 2*(1060*Power(g3,4)*(MassB + 2*MassG) + 27*(-10*traceAdjYdTYdAdjYdYd +
      5*MassB*traceAdjYdYdAdjYdYd - 30*traceAdjYeTYeAdjYeYe + 15*MassB*
      traceAdjYeYeAdjYeYe + 4*traceAdjYuTYuAdjYdYd - 4*MassB*
      traceAdjYuYuAdjYdYd + 4*traceTYdAdjYuYuAdjYd) - 852*(MassB*traceAdjYdYd +
      MassG*traceAdjYdYd - traceTYdAdjYd)*Sqr(g3))) + 125*(18630*Power(g2,6)*
      MassWB + 45*Power(g2,4)*(-6*MassWB*(21*traceAdjYdYd + 7*traceAdjYeYe + 12
      *traceAdjYuYu) + 63*traceTYdAdjYd + 21*traceTYeAdjYe + 36*traceTYuAdjYu +
      112*(MassG + 2*MassWB)*Sqr(g3)) + 36*Sqr(g2)*(68*Power(g3,4)*(2*MassG +
      MassWB) + 3*(-6*traceAdjYdTYdAdjYdYd + MassWB*(3*traceAdjYdYdAdjYdYd +
      traceAdjYeYeAdjYeYe + 6*traceAdjYuYuAdjYdYd) - 2*(traceAdjYeTYeAdjYeYe +
      3*(traceAdjYuTYuAdjYdYd + traceTYdAdjYuYuAdjYd))) - 132*(MassG*
      traceAdjYdYd + MassWB*traceAdjYdYd - traceTYdAdjYd)*Sqr(g3)) + 4*(5440*
      Power(g3,6)*MassG - 480*Power(g3,4)*(2*MassG*(2*traceAdjYdYd +
      traceAdjYuYu) - 2*traceTYdAdjYd - traceTYuAdjYu) - 27*(3*
      traceAdjYdYdAdjYdTYdAdjYdYd + 12*traceAdjYdYd*traceAdjYeTYeAdjYeYe + 4*
      traceAdjYeTYeAdjYeYe*traceAdjYeYe + 12*traceAdjYdTYdAdjYdYd*(3*
      traceAdjYdYd + traceAdjYeYe) + traceAdjYeYeAdjYeTYeAdjYeYe + 3*
      traceAdjYuTYuAdjYuYuAdjYdYd + 6*traceAdjYuTYuAdjYdYd*traceAdjYuYu + 3*
      traceAdjYuYuAdjYuTYuAdjYdYd + 18*traceAdjYdYdAdjYdYd*traceTYdAdjYd + 6*
      traceAdjYeYeAdjYeYe*traceTYdAdjYd + 6*traceAdjYuYu*traceTYdAdjYuYuAdjYd +
      3*traceTYdAdjYuYuAdjYuYuAdjYd + 6*traceAdjYdYdAdjYdYd*traceTYeAdjYe + 2*
      traceAdjYeYeAdjYeYe*traceTYeAdjYe + 6*traceAdjYuYuAdjYdYd*traceTYuAdjYu)
      - 216*(6*traceAdjYdTYdAdjYdYd + traceAdjYuTYuAdjYdYd - MassG*(3*
      traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd) + traceTYdAdjYuYuAdjYd)*Sqr(g3
      )))) + 0.016*threeLoop*(4179*Power(g1,6)*MassB + 25*(81*Power(g2,4)*(
      MassB + 2*MassWB) + 4*(44*Power(g3,4)*(MassB + 2*MassG) + 27*
      traceAdjYdTYdAdjYdYd))*Sqr(g1) - 375*(640*Power(g3,6)*MassG + 315*Power(
      g2,6)*MassWB - 12*(4*Power(g3,4)*(2*MassG + MassWB) + 3*
      traceAdjYdTYdAdjYdYd)*Sqr(g2) - 72*Power(g2,4)*(MassG + 2*MassWB)*Sqr(g3)
      + 96*traceAdjYdTYdAdjYdYd*Sqr(g3)) + 35*Power(g1,4)*(27*(2*MassB +
      MassWB)*Sqr(g2) + 11*(MassB*traceAdjYdYd + 8*(2*MassB + MassG)*Sqr(g3))))
      *(Yd*1.2020569031595942))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_2 = ((-0.04*threeLoop*(Power(
      g1,4)*(162*MassB*traceAdjYeYe + 77*traceTYdAdjYd - 81*traceTYeAdjYe) - 10
      *Sqr(g1)*(9*(5*MassB*traceAdjYdYd + 5*MassWB*traceAdjYdYd - 9*MassB*
      traceAdjYeYe - 9*MassWB*traceAdjYeYe - 5*traceTYdAdjYd + 9*traceTYeAdjYe)
      *Sqr(g2) - 2*(3*(18*traceAdjYeTYeAdjYeYe + MassB*(9*traceAdjYdYdAdjYdYd -
      9*traceAdjYeYeAdjYeYe + 7*traceAdjYuYuAdjYdYd) - 7*(traceAdjYuTYuAdjYdYd
      + traceTYdAdjYuYuAdjYd)) + 56*(MassB*traceAdjYdYd + MassG*traceAdjYdYd -
      traceTYdAdjYd)*Sqr(g3))) - 25*(63*Power(g2,4)*(2*MassWB*(3*traceAdjYdYd
      + traceAdjYeYe) - 3*traceTYdAdjYd - traceTYeAdjYe) - 36*Sqr(g2)*(3*MassWB
      *traceAdjYdYdAdjYdYd - 2*traceAdjYeTYeAdjYeYe + MassWB*
      traceAdjYeYeAdjYeYe + 8*(MassG*traceAdjYdYd + MassWB*traceAdjYdYd -
      traceTYdAdjYd)*Sqr(g3)) + 4*(9*(3*traceAdjYdYdAdjYdTYdAdjYdYd +
      traceAdjYeYeAdjYeTYeAdjYeYe) + 8*Power(g3,4)*(2*MassG*traceAdjYdYd -
      traceTYdAdjYd) - 24*(traceAdjYuTYuAdjYdYd - MassG*(3*traceAdjYdYdAdjYdYd
      + traceAdjYuYuAdjYdYd) + traceTYdAdjYuYuAdjYd)*Sqr(g3))))*(Yd*
      1.2020569031595942) - 0.004*threeLoop*(2786*Power(g1,6) + 5*Power(g1,4)*(
      77*traceAdjYdYd - 81*traceAdjYeYe + 378*Sqr(g2) + 1232*Sqr(g3)) + 50*Sqr(
      g1)*(81*Power(g2,4) + 9*(5*traceAdjYdYd - 9*traceAdjYeYe)*Sqr(g2) + 2*(88
      *Power(g3,4) - 27*traceAdjYdYdAdjYdYd + 27*traceAdjYeYeAdjYeYe - 21*
      traceAdjYuYuAdjYdYd - 56*traceAdjYdYd*Sqr(g3))) - 125*(630*Power(g2,6) -
      9*Power(g2,4)*(7*(3*traceAdjYdYd + traceAdjYeYe) + 48*Sqr(g3)) - 36*Sqr(
      g2)*(8*Power(g3,4) - 3*traceAdjYdYdAdjYdYd - traceAdjYeYeAdjYeYe - 8*
      traceAdjYdYd*Sqr(g3)) + 4*(320*Power(g3,6) - 8*Power(g3,4)*traceAdjYdYd +
      3*(3*traceAdjYdYdAdjYdYdAdjYdYd + traceAdjYeYeAdjYeYeAdjYeYe) - 24*(3*
      traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd)*Sqr(g3))))*(TYd*
      1.2020569031595942) + 0.013333333333333334*threeLoop*(5269*Power(g1,4)*
      MassB + 5*Sqr(g1)*(381*(MassB + MassWB)*Sqr(g2) + 2*(-153*MassB*
      traceAdjYdYd + 69*MassB*traceAdjYeYe + 153*traceTYdAdjYd - 69*
      traceTYeAdjYe + 148*(MassB + MassG)*Sqr(g3))) + 75*(219*Power(g2,4)*
      MassWB + 2*Sqr(g2)*(-15*(3*MassWB*traceAdjYdYd + MassWB*traceAdjYeYe - 3*
      traceTYdAdjYd - traceTYeAdjYe) + 92*(MassG + MassWB)*Sqr(g3)) - 4*(8*
      Power(g3,4)*MassG + 3*(-18*traceAdjYdTYdAdjYdYd - 6*traceAdjYeTYeAdjYeYe
      - 3*traceAdjYuTYuAdjYdYd + 9*traceAdjYdYd*traceTYdAdjYd + 3*traceAdjYeYe*
      traceTYdAdjYd - 3*traceTYdAdjYuYuAdjYd + 3*traceAdjYdYd*traceTYeAdjYe +
      traceAdjYeYe*traceTYeAdjYe) - 12*(MassG*traceAdjYdYd - MassG*traceAdjYeYe
      - traceTYdAdjYd + traceTYeAdjYe)*Sqr(g3))))*(Yd*Yd.adjoint()*Yd) -
      0.013333333333333334*threeLoop*(1792*Power(g1,4) + 5*Sqr(g1)*(-201*
      traceAdjYdYd + 252*Sqr(g2) + 224*Sqr(g3)) + 25*(198*Power(g2,4) - 32*
      Power(g3,4) + 9*Sqr(g2)*(-21*traceAdjYdYd + 32*Sqr(g3))))*(Yd*Yd.adjoint(
      )*TYd))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_3 = ((-0.2*threeLoop*(-360*
      traceAdjYdYdAdjYdYd + 120*traceAdjYdYd*traceAdjYeYe - 120*
      traceAdjYeYeAdjYeYe - 120*traceAdjYuYuAdjYdYd + 31*traceAdjYeYe*Sqr(g1) -
      105*traceAdjYeYe*Sqr(g2) + 160*(traceAdjYdYd - traceAdjYeYe)*Sqr(g3) +
      180*Sqr(traceAdjYdYd) + 20*Sqr(traceAdjYeYe))*(Yd*Yd.adjoint()*TYd) +
      0.013333333333333334*threeLoop*(3767*Power(g1,4)*MassB + 15*Sqr(g1)*(59*(
      MassB + MassWB)*Sqr(g2) + 4*(-5*MassB*traceAdjYuYu + 5*traceTYuAdjYu + 34
      *(MassB + MassG)*Sqr(g3))) + 25*(135*Power(g2,4)*MassWB + 12*Sqr(g2)*(-9*
      MassWB*traceAdjYuYu + 9*traceTYuAdjYu + 2*(MassG + MassWB)*Sqr(g3)) - 4*(
      8*Power(g3,4)*MassG - 9*(traceAdjYuTYuAdjYdYd + 6*traceAdjYuTYuAdjYuYu +
      traceTYdAdjYuYuAdjYd - 3*traceAdjYuYu*traceTYuAdjYu) - 12*(MassG*
      traceAdjYuYu - traceTYuAdjYu)*Sqr(g3))))*(Yd*Yu.adjoint()*Yu) -
      0.0033333333333333335*threeLoop*(7534*Power(g1,4) + 6750*Power(g2,4) -
      1600*Power(g3,4) - 3600*traceAdjYuYuAdjYdYd - 10800*traceAdjYuYuAdjYuYu -
      1200*traceAdjYuYu*Sqr(g1) - 10800*traceAdjYuYu*Sqr(g2) + 3540*Sqr(g1)*
      Sqr(g2) + 4800*traceAdjYuYu*Sqr(g3) + 8160*Sqr(g1)*Sqr(g3) + 2400*Sqr(g2)
      *Sqr(g3) + 5400*Sqr(traceAdjYuYu))*(Yd*Yu.adjoint()*TYu) -
      0.0033333333333333335*threeLoop*(8639*Power(g1,4) + 29475*Power(g2,4) -
      4000*Power(g3,4) - 27000*traceAdjYdYdAdjYdYd + 9000*traceAdjYdYd*
      traceAdjYeYe - 9000*traceAdjYeYeAdjYeYe - 9000*traceAdjYuYuAdjYdYd - 5160
      *traceAdjYdYd*Sqr(g1) + 2280*traceAdjYeYe*Sqr(g1) - 21600*traceAdjYdYd*
      Sqr(g2) - 7200*traceAdjYeYe*Sqr(g2) + 6390*Sqr(g1)*Sqr(g2) + 12000*
      traceAdjYdYd*Sqr(g3) - 12000*traceAdjYeYe*Sqr(g3) + 4400*Sqr(g1)*Sqr(g3)
      + 54000*Sqr(g2)*Sqr(g3) + 13500*Sqr(traceAdjYdYd) + 1500*Sqr(traceAdjYeYe
      ))*(TYd*Yd.adjoint()*Yd) - 0.0033333333333333335*threeLoop*(3767*Power(g1
      ,4) + 3375*Power(g2,4) - 800*Power(g3,4) - 1800*traceAdjYuYuAdjYdYd -
      5400*traceAdjYuYuAdjYuYu - 600*traceAdjYuYu*Sqr(g1) - 5400*traceAdjYuYu*
      Sqr(g2) + 1770*Sqr(g1)*Sqr(g2) + 2400*traceAdjYuYu*Sqr(g3) + 4080*Sqr(g1)
      *Sqr(g3) + 1200*Sqr(g2)*Sqr(g3) + 2700*Sqr(traceAdjYuYu))*(TYd*Yu.adjoint
      ()*Yu) - 0.0033333333333333335*threeLoop*(168*Power(g1,4)*MassB - 326400*
      Power(g3,4)*MassG - 48600*Power(g2,4)*MassWB - 6480*MassB*traceAdjYdYd*
      Sqr(g1) + 5040*MassB*traceAdjYeYe*Sqr(g1) + 6480*traceTYdAdjYd*Sqr(g1) -
      5040*traceTYeAdjYe*Sqr(g1) - 32400*MassWB*traceAdjYdYd*Sqr(g2) - 10800*
      MassWB*traceAdjYeYe*Sqr(g2) + 32400*traceTYdAdjYd*Sqr(g2) + 10800*
      traceTYeAdjYe*Sqr(g2) + 6120*MassB*Sqr(g1)*Sqr(g2) + 6120*MassWB*Sqr(g1)*
      Sqr(g2) + 86400*MassG*traceAdjYdYd*Sqr(g3) - 86400*traceTYdAdjYd*Sqr(g3)
      + 11520*MassB*Sqr(g1)*Sqr(g3) + 11520*MassG*Sqr(g1)*Sqr(g3) + 57600*MassG
      *Sqr(g2)*Sqr(g3) + 57600*MassWB*Sqr(g2)*Sqr(g3))*(Yd*Yd.adjoint()*Yd*
      1.2020569031595942) - 0.0033333333333333335*threeLoop*(-112*Power(g1,4) +
      21600*Power(g2,4) + 108800*Power(g3,4) + 4680*traceAdjYdYd*Sqr(g1) -
      3240*traceAdjYeYe*Sqr(g1) + 16200*traceAdjYdYd*Sqr(g2) - 3600*Sqr(g1)*Sqr
      (g2) - 57600*traceAdjYdYd*Sqr(g3) - 8320*Sqr(g1)*Sqr(g3) - 28800*Sqr(g2)*
      Sqr(g3))*(Yd*Yd.adjoint()*TYd*1.2020569031595942))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_4 = ((-18*threeLoop*
      traceAdjYeYe*Sqr(g2)*(Yd*Yd.adjoint()*TYd*1.2020569031595942) -
      0.02666666666666667*threeLoop*(143*Power(g1,4)*MassB - 25*(544*Power(g3,4
      )*MassG + 189*Power(g2,4)*MassWB - 144*(MassG*traceAdjYuYu -
      traceTYuAdjYu)*Sqr(g3)) + 5*Sqr(g1)*(135*(MassB + MassWB)*Sqr(g2) + 8*(-9
      *MassB*traceAdjYuYu + 9*traceTYuAdjYu + 16*(MassB + MassG)*Sqr(g3))))*(Yd
      *Yu.adjoint()*Yu*1.2020569031595942) + 0.006666666666666667*threeLoop*(
      286*Power(g1,4) - 9450*Power(g2,4) - 27200*Power(g3,4) - 1440*
      traceAdjYuYu*Sqr(g1) + 2700*Sqr(g1)*Sqr(g2) + 14400*traceAdjYuYu*Sqr(g3)
      + 2560*Sqr(g1)*Sqr(g3))*(Yd*Yu.adjoint()*TYu*1.2020569031595942) +
      0.006666666666666667*threeLoop*(7*Power(g1,4) - 7425*Power(g2,4) - 68000*
      Power(g3,4) - 2520*traceAdjYdYd*Sqr(g1) + 2160*traceAdjYeYe*Sqr(g1) -
      16200*traceAdjYdYd*Sqr(g2) - 5400*traceAdjYeYe*Sqr(g2) + 2790*Sqr(g1)*Sqr
      (g2) + 36000*traceAdjYdYd*Sqr(g3) + 4480*Sqr(g1)*Sqr(g3) + 28800*Sqr(g2)*
      Sqr(g3))*(TYd*Yd.adjoint()*Yd*1.2020569031595942) + 0.006666666666666667*
      threeLoop*(143*Power(g1,4) - 4725*Power(g2,4) - 13600*Power(g3,4) - 720*
      traceAdjYuYu*Sqr(g1) + 1350*Sqr(g1)*Sqr(g2) + 7200*traceAdjYuYu*Sqr(g3) +
      1280*Sqr(g1)*Sqr(g3))*(TYd*Yu.adjoint()*Yu*1.2020569031595942) +
      0.006666666666666667*threeLoop*(3600*traceTYdAdjYd + 1200*traceTYeAdjYe -
      40*MassB*Sqr(g1) - 1800*MassWB*Sqr(g2) - 12800*MassG*Sqr(g3))*(Yd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) + 0.006666666666666667*threeLoop*(2700*
      traceAdjYdYd + 900*traceAdjYeYe + 90*Sqr(g1) + 450*Sqr(g2) + 9600*Sqr(g3)
      )*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) + 0.006666666666666667*threeLoop*
      (3600*traceAdjYdYd + 1200*traceAdjYeYe + 40*Sqr(g1) + 1800*Sqr(g2) +
      12800*Sqr(g3))*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) +
      0.006666666666666667*threeLoop*(-1800*traceTYdAdjYd - 600*traceTYeAdjYe +
      3600*traceTYuAdjYu + 580*MassB*Sqr(g1) - 2700*MassWB*Sqr(g2) - 6400*
      MassG*Sqr(g3))*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) +
      0.006666666666666667*threeLoop*(-900*traceAdjYdYd - 300*traceAdjYeYe +
      1800*traceAdjYuYu - 290*Sqr(g1) + 1350*Sqr(g2) + 3200*Sqr(g3))*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*TYd) + 0.006666666666666667*threeLoop*(1800*
      traceTYuAdjYu - 1100*MassB*Sqr(g1) + 900*MassWB*Sqr(g2) - 6400*MassG*Sqr(
      g3))*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu) + 0.006666666666666667*
      threeLoop*(1800*traceAdjYuYu + 1100*Sqr(g1) - 900*Sqr(g2) + 6400*Sqr(g3))
      *(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu) + 0.006666666666666667*threeLoop*(
      -1800*traceAdjYdYd - 600*traceAdjYeYe + 3600*traceAdjYuYu - 580*Sqr(g1) +
      2700*Sqr(g2) + 6400*Sqr(g3))*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) +
      0.006666666666666667*threeLoop*(1800*traceAdjYuYu + 1100*Sqr(g1) - 900*
      Sqr(g2) + 6400*Sqr(g3))*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) +
      0.006666666666666667*threeLoop*(2700*traceAdjYdYd + 900*traceAdjYeYe - 30
      *Sqr(g1) + 2250*Sqr(g2) + 9600*Sqr(g3))*(TYd*Yd.adjoint()*Yd*Yd.adjoint()
      *Yd) + 0.006666666666666667*threeLoop*(-1800*traceAdjYdYd - 600*
      traceAdjYeYe + 3600*traceAdjYuYu - 580*Sqr(g1) + 2700*Sqr(g2) + 6400*Sqr(
      g3))*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) + 0.006666666666666667*
      threeLoop*(900*traceAdjYuYu + 550*Sqr(g1) - 450*Sqr(g2) + 3200*Sqr(g3))*(
      TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu) + 0.006666666666666667*threeLoop*(
      -180*Sqr(g1) + 2700*Sqr(g2))*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd*
      1.2020569031595942) - 7.2*MassB*threeLoop*Sqr(g1)*(Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd*1.2020569031595942))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_5 = ((threeLoop*(
      28.837185185185184*Power(g1,6)*TYd + 172.5*Power(g2,6)*TYd +
      201.4814814814815*Power(g3,6)*TYd - 15.75*Power(g1,4)*traceAdjYdYd*TYd -
      78.75*Power(g2,4)*traceAdjYdYd*TYd - 106.66666666666667*Power(g3,4)*
      traceAdjYdYd*TYd + 54*traceAdjYdYd*traceAdjYdYdAdjYdYd*TYd + 3*
      traceAdjYdYdAdjYdYdAdjYdYd*TYd - 31.65*Power(g1,4)*traceAdjYeYe*TYd -
      26.25*Power(g2,4)*traceAdjYeYe*TYd + 18*traceAdjYdYdAdjYdYd*traceAdjYeYe*
      TYd + 18*traceAdjYdYd*traceAdjYeYeAdjYeYe*TYd + 6*traceAdjYeYe*
      traceAdjYeYeAdjYeYe*TYd + traceAdjYeYeAdjYeYeAdjYeYe*TYd -
      6.066666666666666*Power(g1,4)*traceAdjYuYu*TYd - 45*Power(g2,4)*
      traceAdjYuYu*TYd - 53.333333333333336*Power(g3,4)*traceAdjYuYu*TYd + 18*
      traceAdjYuYu*traceAdjYuYuAdjYdYd*TYd + 9*traceAdjYuYuAdjYuYuAdjYdYd*TYd +
      8.5*Power(g2,4)*Sqr(g1)*TYd + 23.555555555555557*Power(g3,4)*Sqr(g1)*TYd
      + 3*traceAdjYdYdAdjYdYd*Sqr(g1)*TYd + 9*traceAdjYeYeAdjYeYe*Sqr(g1)*TYd
      - 2.4*traceAdjYuYuAdjYdYd*Sqr(g1)*TYd + 2.18*Power(g1,4)*Sqr(g2)*TYd + 68
      *Power(g3,4)*Sqr(g2)*TYd + 9*traceAdjYdYdAdjYdYd*Sqr(g2)*TYd + 3*
      traceAdjYeYeAdjYeYe*Sqr(g2)*TYd + 18*traceAdjYuYuAdjYdYd*Sqr(g2)*TYd -
      0.3*traceAdjYdYd*Sqr(g1)*Sqr(g2)*TYd - 8.1*traceAdjYeYe*Sqr(g1)*Sqr(g2)*
      TYd + 17.297777777777778*Power(g1,4)*Sqr(g3)*TYd + 140*Power(g2,4)*Sqr(g3
      )*TYd + 72*traceAdjYdYdAdjYdYd*Sqr(g3)*TYd + 24*traceAdjYuYuAdjYdYd*Sqr(
      g3)*TYd - 18.933333333333334*traceAdjYdYd*Sqr(g1)*Sqr(g3)*TYd - 132*
      traceAdjYdYd*Sqr(g2)*Sqr(g3)*TYd - 1.6*Sqr(g1)*Sqr(g2)*Sqr(g3)*TYd) + 36*
      MassWB*threeLoop*Sqr(g2)*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*
      1.2020569031595942) + 3.6*threeLoop*(Sqr(g1) - 5*Sqr(g2))*(Yd*Yu.adjoint(
      )*Yu*Yd.adjoint()*TYd*1.2020569031595942) + threeLoop*(12*MassB*Sqr(g1) -
      36*MassWB*Sqr(g2))*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*
      1.2020569031595942) + threeLoop*(-12*Sqr(g1) + 36*Sqr(g2))*(Yd*Yu.adjoint
      ()*Yu*Yu.adjoint()*TYu*1.2020569031595942) + threeLoop*(7.2*Sqr(g1) - 36*
      Sqr(g2))*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd*1.2020569031595942) +
      threeLoop*(-12*Sqr(g1) + 36*Sqr(g2))*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu
      *1.2020569031595942) + threeLoop*(1.2*Sqr(g1) - 18*Sqr(g2))*(TYd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd*1.2020569031595942) + threeLoop*(7.2*Sqr(
      g1) - 36*Sqr(g2))*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*1.2020569031595942
      ) + threeLoop*(-6*Sqr(g1) + 18*Sqr(g2))*(TYd*Yu.adjoint()*Yu*Yu.adjoint()
      *Yu*1.2020569031595942) + 6*threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd
      *Yd.adjoint()*TYd) + 12*threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd*
      Yd.adjoint()*Yd) + 6*threeLoop*(Yd*Yd.adjoint()*Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*TYd) + 4*threeLoop*(Yd*Yd.adjoint()*Yd*Yu.adjoint()*TYu*
      Yd.adjoint()*Yd) + 12*threeLoop*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) + 4*threeLoop*(Yd*Yd.adjoint()*TYd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) - 2*threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*
      Yd.adjoint()*TYd) + 8*threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*
      Yu.adjoint()*TYu) - 4*threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd*
      Yd.adjoint()*Yd) + 8*threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd*
      Yu.adjoint()*Yu) + 6*threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*
      Yd.adjoint()*TYd) + 12*threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu*
      Yd.adjoint()*Yd) - 4*threeLoop*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) + 8*threeLoop*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu) + 12*threeLoop*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) + 12*threeLoop*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) - 4*threeLoop*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) + 4*threeLoop*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu) + 12*threeLoop*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) + 24*threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*
      Yd.adjoint()*TYd*1.2020569031595942) + 36*threeLoop*(Yd*Yd.adjoint()*Yd*
      Yd.adjoint()*TYd*Yd.adjoint()*Yd*1.2020569031595942) + 36*threeLoop*(Yd*
      Yd.adjoint()*TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*1.2020569031595942) + 12
      *threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu*
      1.2020569031595942) + 12*threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu*
      Yu.adjoint()*Yu*1.2020569031595942) + 12*threeLoop*(Yd*Yu.adjoint()*TYu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu*1.2020569031595942) + 30*threeLoop*(TYd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*1.2020569031595942) + 6*
      threeLoop*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu*
      1.2020569031595942))*UNITMATRIX(3)).real();

   beta_TYd = beta_TYd_1 + beta_TYd_2 + beta_TYd_3 + beta_TYd_4 +
      beta_TYd_5;


   return beta_TYd;
}

} // namespace flexiblesusy
