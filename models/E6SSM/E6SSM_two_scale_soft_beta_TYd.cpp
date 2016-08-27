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

// File generated at Sat 27 Aug 2016 12:44:15

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

   beta_TYd = (oneOver16PiSqr*(3*traceYdAdjYd*TYd + traceYeAdjYe*TYd +
      AbsSqr(Lambdax)*TYd - 0.4666666666666667*Sqr(g1)*TYd - 3*Sqr(g2)*TYd -
      5.333333333333333*Sqr(g3)*TYd - 0.7*Sqr(gN)*TYd + Yd*(6*traceAdjYdTYd + 2
      *traceAdjYeTYe + 0.9333333333333333*MassB*Sqr(g1) + 6*MassWB*Sqr(g2) +
      10.666666666666666*MassG*Sqr(g3) + 1.4*MassBp*Sqr(gN) + 2*Conj(Lambdax)*
      TLambdax) + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*
      Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu)).real();


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
      twoLoop*(6608*Power(g1,4)*MassB*Yd + 9828*Power(gN,4)*MassBp*Yd + 20480*
      Power(g3,4)*MassG*Yd + 23760*Power(g2,4)*MassWB*Yd + 12960*
      traceYdAdjYdTYdAdjYd*Yd + 2160*traceYdAdjYuTYuAdjYd*Yd + 4320*
      traceYeAdjYeTYeAdjYe*Yd + 2160*traceYuAdjYdTYdAdjYu*Yd + 2160*
      traceAdjKappaTKappa*Yd*AbsSqr(Lambdax) + 1440*traceAdjLambda12TLambda12*
      Yd*AbsSqr(Lambdax) + 2160*traceAdjYuTYu*Yd*AbsSqr(Lambdax) + 288*
      traceAdjYdTYd*Yd*Sqr(g1) - 864*traceAdjYeTYe*Yd*Sqr(g1) - 288*MassB*
      traceYdAdjYd*Yd*Sqr(g1) + 864*MassB*traceYeAdjYe*Yd*Sqr(g1) + 720*MassB*
      Yd*Sqr(g1)*Sqr(g2) + 720*MassWB*Yd*Sqr(g1)*Sqr(g2) - 11520*traceAdjYdTYd*
      Yd*Sqr(g3) + 11520*MassG*traceYdAdjYd*Yd*Sqr(g3) + 640*MassB*Yd*Sqr(g1)*
      Sqr(g3) + 640*MassG*Yd*Sqr(g1)*Sqr(g3) + 5760*MassG*Yd*Sqr(g2)*Sqr(g3) +
      5760*MassWB*Yd*Sqr(g2)*Sqr(g3) + 432*traceAdjYdTYd*Yd*Sqr(gN) + 144*
      traceAdjYeTYe*Yd*Sqr(gN) - 432*MassBp*traceYdAdjYd*Yd*Sqr(gN) - 144*
      MassBp*traceYeAdjYe*Yd*Sqr(gN) + 720*MassBp*Yd*AbsSqr(Lambdax)*Sqr(gN) -
      168*MassB*Yd*Sqr(g1)*Sqr(gN) - 168*MassBp*Yd*Sqr(g1)*Sqr(gN) + 1080*
      MassBp*Yd*Sqr(g2)*Sqr(gN) + 1080*MassWB*Yd*Sqr(g2)*Sqr(gN) + 960*MassBp*
      Yd*Sqr(g3)*Sqr(gN) + 960*MassG*Yd*Sqr(g3)*Sqr(gN) - 1652*Power(g1,4)*TYd
      - 5940*Power(g2,4)*TYd - 5120*Power(g3,4)*TYd - 2457*Power(gN,4)*TYd +
      3240*traceYdAdjYdYdAdjYd*TYd + 1080*traceYdAdjYuYuAdjYd*TYd + 1080*
      traceYeAdjYeYeAdjYe*TYd + 1080*traceKappaAdjKappa*AbsSqr(Lambdax)*TYd +
      720*traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYd + 1080*traceYuAdjYu*
      AbsSqr(Lambdax)*TYd + 144*traceYdAdjYd*Sqr(g1)*TYd - 432*traceYeAdjYe*Sqr
      (g1)*TYd - 360*Sqr(g1)*Sqr(g2)*TYd - 5760*traceYdAdjYd*Sqr(g3)*TYd - 320*
      Sqr(g1)*Sqr(g3)*TYd - 2880*Sqr(g2)*Sqr(g3)*TYd + 216*traceYdAdjYd*Sqr(gN)
      *TYd + 72*traceYeAdjYe*Sqr(gN)*TYd - 360*AbsSqr(Lambdax)*Sqr(gN)*TYd + 84
      *Sqr(g1)*Sqr(gN)*TYd - 540*Sqr(g2)*Sqr(gN)*TYd - 480*Sqr(g3)*Sqr(gN)*TYd
      + 2160*traceKappaAdjKappa*Yd*Conj(Lambdax)*TLambdax + 1440*
      traceLambda12AdjLambda12*Yd*Conj(Lambdax)*TLambdax + 2160*traceYuAdjYu*Yd
      *Conj(Lambdax)*TLambdax - 720*Yd*Conj(Lambdax)*Sqr(gN)*TLambdax + 1080*
      Lambdax*Sqr(Conj(Lambdax))*(Lambdax*TYd + 4*Yd*TLambdax)) -
      0.002777777777777778*twoLoop*(144*(4*MassB*Sqr(g1) + 15*(3*traceAdjYdTYd
      + traceAdjYeTYe + 2*MassWB*Sqr(g2)) + 6*MassBp*Sqr(gN)) + 2160*Conj(
      Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) - 0.002777777777777778*twoLoop*(
      4320*traceYdAdjYd + 1440*traceYeAdjYe + 1440*AbsSqr(Lambdax) - 432*Sqr(g1
      ) - 2160*Sqr(g2) - 648*Sqr(gN))*(Yd*Yd.adjoint()*TYd) -
      0.002777777777777778*twoLoop*(2160*traceAdjYuTYu + 576*MassB*Sqr(g1) +
      144*MassBp*Sqr(gN))*(Yd*Yu.adjoint()*Yu) - 0.002777777777777778*twoLoop*(
      2160*traceYuAdjYu + 720*AbsSqr(Lambdax) - 576*Sqr(g1) - 144*Sqr(gN))*(Yd*
      Yu.adjoint()*TYu) - 0.002777777777777778*twoLoop*(5400*traceYdAdjYd +
      1800*traceYeAdjYe + 1800*AbsSqr(Lambdax) - 432*Sqr(g1) - 4320*Sqr(g2) -
      648*Sqr(gN))*(TYd*Yd.adjoint()*Yd) - 0.002777777777777778*twoLoop*(1080*
      traceYuAdjYu + 360*AbsSqr(Lambdax) - 288*Sqr(g1) - 72*Sqr(gN))*(TYd*
      Yu.adjoint()*Yu) - 6*twoLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 8*
      twoLoop*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*twoLoop*(Yd*Yu.adjoint(
      )*Yu*Yd.adjoint()*TYd) - 4*twoLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu)
      - 4*twoLoop*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*twoLoop*(Yd*
      Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*twoLoop*(TYd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) - 4*twoLoop*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*
      twoLoop*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))*UNITMATRIX(3)).real();
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
