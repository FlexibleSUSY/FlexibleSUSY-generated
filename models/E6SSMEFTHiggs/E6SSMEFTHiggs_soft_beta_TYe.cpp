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

// File generated at Wed 16 Oct 2019 19:09:42

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
 * Calculates the 1-loop beta function of TYe.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYe_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(0.1*(60*traceAdjYdTYd*Ye + 20*traceAdjYeTYe*Ye +
      36*MassB*Ye*Sqr(g1) + 60*MassWB*Ye*Sqr(g2) + 14*MassBp*Ye*Sqr(gN) + 30*
      traceYdAdjYd*TYe + 10*traceYeAdjYe*TYe + 10*AbsSqr(Lambdax)*TYe - 18*Sqr(
      g1)*TYe - 30*Sqr(g2)*TYe - 7*Sqr(gN)*TYe + 20*Ye*Conj(Lambdax)*TLambdax)
      + 4*(Ye*Ye.adjoint()*TYe) + 5*(TYe*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the 2-loop beta function of TYe.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYe_2_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(0.025*(-1440*traceYdAdjYdTYdAdjYd*Ye - 240*
      traceYdAdjYuTYuAdjYd*Ye - 480*traceYeAdjYeTYeAdjYe*Ye - 240*
      traceYuAdjYdTYdAdjYu*Ye - 240*traceAdjKappaTKappa*Ye*AbsSqr(Lambdax) -
      160*traceAdjLambda12TLambda12*Ye*AbsSqr(Lambdax) - 240*traceAdjYuTYu*Ye*
      AbsSqr(Lambdax) - 3024*MassB*Ye*Quad(g1) - 2640*MassWB*Ye*Quad(g2) - 1092
      *MassBp*Ye*Quad(gN) - 32*traceAdjYdTYd*Ye*Sqr(g1) + 96*traceAdjYeTYe*Ye*
      Sqr(g1) + 32*MassB*traceYdAdjYd*Ye*Sqr(g1) - 96*MassB*traceYeAdjYe*Ye*Sqr
      (g1) - 144*MassB*Ye*Sqr(g1)*Sqr(g2) - 144*MassWB*Ye*Sqr(g1)*Sqr(g2) +
      1280*traceAdjYdTYd*Ye*Sqr(g3) - 1280*MassG*traceYdAdjYd*Ye*Sqr(g3) - 48*
      traceAdjYdTYd*Ye*Sqr(gN) - 16*traceAdjYeTYe*Ye*Sqr(gN) + 48*MassBp*
      traceYdAdjYd*Ye*Sqr(gN) + 16*MassBp*traceYeAdjYe*Ye*Sqr(gN) - 80*MassBp*
      Ye*AbsSqr(Lambdax)*Sqr(gN) - 12*MassB*Ye*Sqr(g1)*Sqr(gN) - 12*MassBp*Ye*
      Sqr(g1)*Sqr(gN) - 156*MassBp*Ye*Sqr(g2)*Sqr(gN) - 156*MassWB*Ye*Sqr(g2)*
      Sqr(gN) - 360*traceYdAdjYdYdAdjYd*TYe - 120*traceYdAdjYuYuAdjYd*TYe - 120
      *traceYeAdjYeYeAdjYe*TYe - 120*traceKappaAdjKappa*AbsSqr(Lambdax)*TYe -
      80*traceLambda12AdjLambda12*AbsSqr(Lambdax)*TYe - 120*traceYuAdjYu*AbsSqr
      (Lambdax)*TYe + 756*Quad(g1)*TYe + 660*Quad(g2)*TYe + 273*Quad(gN)*TYe -
      16*traceYdAdjYd*Sqr(g1)*TYe + 48*traceYeAdjYe*Sqr(g1)*TYe + 72*Sqr(g1)*
      Sqr(g2)*TYe + 640*traceYdAdjYd*Sqr(g3)*TYe - 24*traceYdAdjYd*Sqr(gN)*TYe
      - 8*traceYeAdjYe*Sqr(gN)*TYe + 40*AbsSqr(Lambdax)*Sqr(gN)*TYe + 6*Sqr(g1)
      *Sqr(gN)*TYe + 78*Sqr(g2)*Sqr(gN)*TYe - 120*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)*TYe - 240*traceKappaAdjKappa*Ye*Conj(Lambdax)*TLambdax - 160*
      traceLambda12AdjLambda12*Ye*Conj(Lambdax)*TLambdax - 240*traceYuAdjYu*Ye*
      Conj(Lambdax)*TLambdax + 80*Ye*Conj(Lambdax)*Sqr(gN)*TLambdax - 480*Ye*
      Lambdax*Sqr(Conj(Lambdax))*TLambdax) - 3*(6*traceAdjYdTYd + 2*
      traceAdjYeTYe + 4*MassWB*Sqr(g2) + MassBp*Sqr(gN) + 2*Conj(Lambdax)*
      TLambdax)*(Ye*Ye.adjoint()*Ye) + 0.2*(-60*traceYdAdjYd - 20*traceYeAdjYe
      - 20*AbsSqr(Lambdax) + 6*Sqr(g1) + 30*Sqr(g2) + 9*Sqr(gN))*(Ye*Ye.adjoint
      ()*TYe) + 0.1*(-150*traceYdAdjYd - 50*traceYeAdjYe - 50*AbsSqr(Lambdax) -
      12*Sqr(g1) + 120*Sqr(g2) + 27*Sqr(gN))*(TYe*Ye.adjoint()*Ye) - 6*(Ye*Ye.
      adjoint()*Ye*Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye)
      - 6*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the 3-loop beta function of TYe.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYe_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 4-loop beta function of TYe.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYe_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

/**
 * Calculates the 5-loop beta function of TYe.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_TYe_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = ZEROMATRIX(3,3);


   return beta_TYe;
}

} // namespace flexiblesusy
