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

// File generated at Tue 22 Jan 2019 17:02:00

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
 * Calculates the 1-loop beta function of TYd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*(0.03333333333333333*(180*traceAdjYdTYd*Yd + 60*
      traceAdjYeTYe*Yd + 28*MassB*Yd*Sqr(g1) + 180*MassWB*Yd*Sqr(g2) + 320*
      MassG*Yd*Sqr(g3) + 42*MassBp*Yd*Sqr(gN) + 90*traceYdAdjYd*TYd + 30*
      traceYeAdjYe*TYd + 30*AbsSqr(Lambdax)*TYd - 14*Sqr(g1)*TYd - 90*Sqr(g2)*
      TYd - 160*Sqr(g3)*TYd - 21*Sqr(gN)*TYd + 60*Yd*Conj(Lambdax)*TLambdax) +
      4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*Yu.adjoint()*TYu) + 5*(TYd*Yd.adjoint()*
      Yd) + TYd*Yu.adjoint()*Yu)).real();


   return beta_TYd;
}

/**
 * Calculates the 2-loop beta function of TYd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_2_loop(const Soft_traces& soft_traces) const
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
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   const Eigen::Matrix<double,3,3> beta_TYd_1 = ((-0.002777777777777778*twoLoop
      *(12960*traceYdAdjYdTYdAdjYd*Yd + 2160*traceYdAdjYuTYuAdjYd*Yd + 4320*
      traceYeAdjYeTYeAdjYe*Yd + 2160*traceYuAdjYdTYdAdjYu*Yd + 2160*
      traceAdjKappaTKappa*Yd*AbsSqr(Lambdax) + 1440*traceAdjLambda12TLambda12*
      Yd*AbsSqr(Lambdax) + 2160*traceAdjYuTYu*Yd*AbsSqr(Lambdax) + 6608*MassB*
      Yd*Quad(g1) + 23760*MassWB*Yd*Quad(g2) + 20480*MassG*Yd*Quad(g3) + 9828*
      MassBp*Yd*Quad(gN) + 288*traceAdjYdTYd*Yd*Sqr(g1) - 864*traceAdjYeTYe*Yd*
      Sqr(g1) - 288*MassB*traceYdAdjYd*Yd*Sqr(g1) + 864*MassB*traceYeAdjYe*Yd*
      Sqr(g1) + 720*MassB*Yd*Sqr(g1)*Sqr(g2) + 720*MassWB*Yd*Sqr(g1)*Sqr(g2) -
      11520*traceAdjYdTYd*Yd*Sqr(g3) + 11520*MassG*traceYdAdjYd*Yd*Sqr(g3) +
      640*MassB*Yd*Sqr(g1)*Sqr(g3) + 640*MassG*Yd*Sqr(g1)*Sqr(g3) + 5760*MassG*
      Yd*Sqr(g2)*Sqr(g3) + 5760*MassWB*Yd*Sqr(g2)*Sqr(g3) + 432*traceAdjYdTYd*
      Yd*Sqr(gN) + 144*traceAdjYeTYe*Yd*Sqr(gN) - 432*MassBp*traceYdAdjYd*Yd*
      Sqr(gN) - 144*MassBp*traceYeAdjYe*Yd*Sqr(gN) + 720*MassBp*Yd*AbsSqr(
      Lambdax)*Sqr(gN) - 168*MassB*Yd*Sqr(g1)*Sqr(gN) - 168*MassBp*Yd*Sqr(g1)*
      Sqr(gN) + 1080*MassBp*Yd*Sqr(g2)*Sqr(gN) + 1080*MassWB*Yd*Sqr(g2)*Sqr(gN)
      + 960*MassBp*Yd*Sqr(g3)*Sqr(gN) + 960*MassG*Yd*Sqr(g3)*Sqr(gN) + 3240*
      traceYdAdjYdYdAdjYd*TYd + 1080*traceYdAdjYuYuAdjYd*TYd + 1080*
      traceYeAdjYeYeAdjYe*TYd + 1080*traceKappaAdjKappa*AbsSqr(Lambdax)*TYd +
      720*traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYd + 1080*traceYuAdjYu*
      AbsSqr(Lambdax)*TYd - 1652*Quad(g1)*TYd - 5940*Quad(g2)*TYd - 5120*Quad(
      g3)*TYd - 2457*Quad(gN)*TYd + 144*traceYdAdjYd*Sqr(g1)*TYd - 432*
      traceYeAdjYe*Sqr(g1)*TYd - 360*Sqr(g1)*Sqr(g2)*TYd - 5760*traceYdAdjYd*
      Sqr(g3)*TYd - 320*Sqr(g1)*Sqr(g3)*TYd - 2880*Sqr(g2)*Sqr(g3)*TYd + 216*
      traceYdAdjYd*Sqr(gN)*TYd + 72*traceYeAdjYe*Sqr(gN)*TYd - 360*AbsSqr(
      Lambdax)*Sqr(gN)*TYd + 84*Sqr(g1)*Sqr(gN)*TYd - 540*Sqr(g2)*Sqr(gN)*TYd -
      480*Sqr(g3)*Sqr(gN)*TYd + 1080*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TYd + 2160
      *traceKappaAdjKappa*Yd*Conj(Lambdax)*TLambdax + 1440*
      traceLambda12AdjLambda12*Yd*Conj(Lambdax)*TLambdax + 2160*traceYuAdjYu*Yd
      *Conj(Lambdax)*TLambdax - 720*Yd*Conj(Lambdax)*Sqr(gN)*TLambdax + 4320*Yd
      *Lambdax*Sqr(Conj(Lambdax))*TLambdax) - 0.4*twoLoop*(45*traceAdjYdTYd +
      15*traceAdjYeTYe + 4*MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 6*MassBp*Sqr(gN)
      + 15*Conj(Lambdax)*TLambdax)*(Yd*Yd.adjoint()*Yd) + 0.2*twoLoop*(-60*
      traceYdAdjYd - 20*traceYeAdjYe - 20*AbsSqr(Lambdax) + 6*Sqr(g1) + 30*Sqr(
      g2) + 9*Sqr(gN))*(Yd*Yd.adjoint()*TYd) - 0.4*twoLoop*(15*traceAdjYuTYu +
      4*MassB*Sqr(g1) + MassBp*Sqr(gN))*(Yd*Yu.adjoint()*Yu) + 0.4*twoLoop*(-15
      *traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + Sqr(gN))*(Yd*Yu.adjoint()
      *TYu) + 0.2*twoLoop*(-75*traceYdAdjYd - 25*traceYeAdjYe - 25*AbsSqr(
      Lambdax) + 6*Sqr(g1) + 60*Sqr(g2) + 9*Sqr(gN))*(TYd*Yd.adjoint()*Yd) +
      0.2*twoLoop*(-15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + Sqr(gN))*
      (TYd*Yu.adjoint()*Yu) - 6*twoLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) -
      8*twoLoop*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) - 2*twoLoop*(Yd*Yu.
      adjoint()*Yu*Yd.adjoint()*TYd) - 4*twoLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint
      ()*TYu) - 4*twoLoop*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*twoLoop*(Yd
      *Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*twoLoop*(TYd*Yd.adjoint()*Yd*Yd.
      adjoint()*Yd) - 4*twoLoop*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*
      twoLoop*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_2 = (-2*twoLoop*Conj(Lambdax)*
      TLambdax*(Yd*Yu.adjoint()*Yu)*UNITMATRIX(3)).real();

   beta_TYd = beta_TYd_1 + beta_TYd_2;


   return beta_TYd;
}

/**
 * Calculates the 3-loop beta function of TYd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

/**
 * Calculates the 4-loop beta function of TYd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

/**
 * Calculates the 5-loop beta function of TYd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TYd_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = ZEROMATRIX(3,3);


   return beta_TYd;
}

} // namespace flexiblesusy
