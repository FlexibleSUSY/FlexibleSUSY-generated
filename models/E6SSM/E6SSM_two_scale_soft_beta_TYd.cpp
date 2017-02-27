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

// File generated at Mon 27 Feb 2017 13:33:22

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
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*((3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr
      (g3) - 0.7*Sqr(gN))*TYd + 0.06666666666666667*Yd*(90*traceAdjYdTYd + 30*
      traceAdjYeTYe + 14*MassB*Sqr(g1) + 90*MassWB*Sqr(g2) + 160*MassG*Sqr(g3)
      + 21*MassBp*Sqr(gN) + 30*Conj(Lambdax)*TLambdax) + 4*(Yd*Yd.adjoint()*TYd
      ) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*
      Yu)).real();


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYd;

   const Eigen::Matrix<double,3,3> beta_TYd_1 = ((-0.002777777777777778*
      twoLoop*(4*Yd*(1652*Power(g1,4)*MassB + 2457*Power(gN,4)*MassBp + 5120*
      Power(g3,4)*MassG + 5940*Power(g2,4)*MassWB + 3240*traceYdAdjYdTYdAdjYd +
      540*traceYdAdjYuTYuAdjYd + 1080*traceYeAdjYeTYeAdjYe + 540*
      traceYuAdjYdTYdAdjYu - 2880*traceAdjYdTYd*Sqr(g3) + 2880*MassG*
      traceYdAdjYd*Sqr(g3) + 108*traceAdjYdTYd*Sqr(gN) + 36*traceAdjYeTYe*Sqr(
      gN) - 108*MassBp*traceYdAdjYd*Sqr(gN) - 36*MassBp*traceYeAdjYe*Sqr(gN) +
      240*MassBp*Sqr(g3)*Sqr(gN) + 240*MassG*Sqr(g3)*Sqr(gN) + 90*Sqr(g2)*(16*(
      MassG + MassWB)*Sqr(g3) + 3*(MassBp + MassWB)*Sqr(gN)) + 2*Sqr(g1)*(90*(
      MassB + MassWB)*Sqr(g2) + 80*(MassB + MassG)*Sqr(g3) - 3*(12*(
      -traceAdjYdTYd + 3*traceAdjYeTYe + MassB*traceYdAdjYd - 3*MassB*
      traceYeAdjYe) + 7*(MassB + MassBp)*Sqr(gN)))) - (1652*Power(g1,4) + 5940*
      Power(g2,4) + 5120*Power(g3,4) + 2457*Power(gN,4) - 3240*
      traceYdAdjYdYdAdjYd - 1080*traceYdAdjYuYuAdjYd - 1080*traceYeAdjYeYeAdjYe
      + 5760*traceYdAdjYd*Sqr(g3) + 4*Sqr(g1)*(-36*traceYdAdjYd + 108*
      traceYeAdjYe + 90*Sqr(g2) + 80*Sqr(g3) - 21*Sqr(gN)) - 216*traceYdAdjYd*
      Sqr(gN) - 72*traceYeAdjYe*Sqr(gN) + 480*Sqr(g3)*Sqr(gN) + 180*Sqr(g2)*(16
      *Sqr(g3) + 3*Sqr(gN)))*TYd + 1080*Lambdax*Sqr(Conj(Lambdax))*(Lambdax*TYd
      + 4*Yd*TLambdax) + 360*Conj(Lambdax)*(Lambdax*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 3*traceYuAdjYu - Sqr(gN))*TYd + 2*Yd*(Lambdax*
      (3*traceAdjKappaTKappa + 2*traceAdjLambda12TLambda12 + 3*traceAdjYuTYu +
      MassBp*Sqr(gN)) + (3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*
      traceYuAdjYu - Sqr(gN))*TLambdax))) - 0.4*twoLoop*(45*traceAdjYdTYd + 15*
      traceAdjYeTYe + 4*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 6*MassBp*Sqr(gN) +
      15*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) + 0.2*twoLoop*(-60*
      traceYdAdjYd - 20*traceYeAdjYe - 20*AbsSqr(Lambdax) + 6*Sqr(g1) + 30*Sqr(
      g2) + 9*Sqr(gN))*(Yd*Yd.adjoint()*TYd) - 0.4*twoLoop*(15*traceAdjYuTYu +
      4*MassB*Sqr(g1) + MassBp*Sqr(gN))*(Yd*Yu.adjoint()*Yu) + 0.4*twoLoop*(-15
      *traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + Sqr(gN))*(Yd*Yu.adjoint()
      *TYu) + 0.2*twoLoop*(-75*traceYdAdjYd - 25*traceYeAdjYe - 25*AbsSqr(
      Lambdax) + 6*Sqr(g1) + 60*Sqr(g2) + 9*Sqr(gN))*(TYd*Yd.adjoint()*Yd) +
      0.2*twoLoop*(-15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + Sqr(gN))*
      (TYd*Yu.adjoint()*Yu) - 6*twoLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) -
      8*twoLoop*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*twoLoop*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*TYd) - 4*twoLoop*(Yd*Yu.adjoint()*Yu*
      Yu.adjoint()*TYu) - 4*twoLoop*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*
      twoLoop*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*twoLoop*(TYd*Yd.adjoint
      ()*Yd*Yd.adjoint()*Yd) - 4*twoLoop*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd)
      - 2*twoLoop*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_2 = (-2*twoLoop*Conj(Lambdax)
      *TLambdax*(Yd*Yu.adjoint()*Yu)*UNITMATRIX(3)).real();

   beta_TYd = beta_TYd_1 + beta_TYd_2;


   return beta_TYd;
}

/**
 * Calculates the three-loop beta function of TYd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

} // namespace flexiblesusy
