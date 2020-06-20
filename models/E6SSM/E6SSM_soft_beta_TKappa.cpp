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
 * Calculates the 1-loop beta function of TKappa.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 = TRACE_STRUCT.
      traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = (0.03333333333333333*(180*traceAdjKappaTKappa*Kappa + 120*
      traceAdjLambda12TLambda12*Kappa + 16*MassB*Kappa*Sqr(g1) + 320*MassG*
      Kappa*Sqr(g3) + 114*MassBp*Kappa*Sqr(gN) + 90*traceKappaAdjKappa*TKappa +
      60*traceLambda12AdjLambda12*TKappa + 60*AbsSqr(Lambdax)*TKappa - 8*Sqr(g1
      )*TKappa - 160*Sqr(g3)*TKappa - 57*Sqr(gN)*TKappa + 120*Conj(Lambdax)*
      Kappa*TLambdax) + 3*(Kappa*(Kappa).adjoint()*TKappa) + 3*(TKappa*(Kappa).
      adjoint()*Kappa)).real();


   return oneLoop * beta_TKappa;
}

/**
 * Calculates the 2-loop beta function of TKappa.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_2_loop(const Soft_traces& soft_traces) const
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
   const double traceKappaAdjKappaTKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12TLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12Lambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = (0.0005555555555555556*(-43200*
      traceKappaAdjKappaTKappaAdjKappa*Kappa - 28800*
      traceLambda12AdjLambda12TLambda12AdjLambda12*Kappa - 21600*traceAdjYdTYd*
      AbsSqr(Lambdax)*Kappa - 7200*traceAdjYeTYe*AbsSqr(Lambdax)*Kappa - 21600*
      traceAdjYuTYu*AbsSqr(Lambdax)*Kappa - 18688*MassB*Kappa*Quad(g1) - 102400
      *MassG*Kappa*Quad(g3) - 141588*MassBp*Kappa*Quad(gN) + 2880*
      traceAdjKappaTKappa*Kappa*Sqr(g1) + 4320*traceAdjLambda12TLambda12*Kappa*
      Sqr(g1) - 2880*MassB*traceKappaAdjKappa*Kappa*Sqr(g1) - 4320*MassB*
      traceLambda12AdjLambda12*Kappa*Sqr(g1) - 4320*MassB*AbsSqr(Lambdax)*Kappa
      *Sqr(g1) + 21600*traceAdjLambda12TLambda12*Kappa*Sqr(g2) - 21600*MassWB*
      traceLambda12AdjLambda12*Kappa*Sqr(g2) - 21600*MassWB*AbsSqr(Lambdax)*
      Kappa*Sqr(g2) + 57600*traceAdjKappaTKappa*Kappa*Sqr(g3) - 57600*MassG*
      traceKappaAdjKappa*Kappa*Sqr(g3) - 5120*MassB*Kappa*Sqr(g1)*Sqr(g3) -
      5120*MassG*Kappa*Sqr(g1)*Sqr(g3) - 6480*traceAdjKappaTKappa*Kappa*Sqr(gN)
      - 4320*traceAdjLambda12TLambda12*Kappa*Sqr(gN) + 6480*MassBp*
      traceKappaAdjKappa*Kappa*Sqr(gN) + 4320*MassBp*traceLambda12AdjLambda12*
      Kappa*Sqr(gN) + 4320*MassBp*AbsSqr(Lambdax)*Kappa*Sqr(gN) - 912*MassB*
      Kappa*Sqr(g1)*Sqr(gN) - 912*MassBp*Kappa*Sqr(g1)*Sqr(gN) - 12480*MassBp*
      Kappa*Sqr(g3)*Sqr(gN) - 12480*MassG*Kappa*Sqr(g3)*Sqr(gN) - 10800*
      traceKappaAdjKappaKappaAdjKappa*TKappa - 7200*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TKappa - 10800*traceYdAdjYd*
      AbsSqr(Lambdax)*TKappa - 3600*traceYeAdjYe*AbsSqr(Lambdax)*TKappa - 10800
      *traceYuAdjYu*AbsSqr(Lambdax)*TKappa + 4672*Quad(g1)*TKappa + 25600*Quad(
      g3)*TKappa + 35397*Quad(gN)*TKappa + 1440*traceKappaAdjKappa*Sqr(g1)*
      TKappa + 2160*traceLambda12AdjLambda12*Sqr(g1)*TKappa + 2160*AbsSqr(
      Lambdax)*Sqr(g1)*TKappa + 10800*traceLambda12AdjLambda12*Sqr(g2)*TKappa +
      10800*AbsSqr(Lambdax)*Sqr(g2)*TKappa + 28800*traceKappaAdjKappa*Sqr(g3)*
      TKappa + 2560*Sqr(g1)*Sqr(g3)*TKappa - 3240*traceKappaAdjKappa*Sqr(gN)*
      TKappa - 2160*traceLambda12AdjLambda12*Sqr(gN)*TKappa - 2160*AbsSqr(
      Lambdax)*Sqr(gN)*TKappa + 456*Sqr(g1)*Sqr(gN)*TKappa + 6240*Sqr(g3)*Sqr(
      gN)*TKappa - 7200*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TKappa - 21600*
      traceYdAdjYd*Conj(Lambdax)*Kappa*TLambdax - 7200*traceYeAdjYe*Conj(
      Lambdax)*Kappa*TLambdax - 21600*traceYuAdjYu*Conj(Lambdax)*Kappa*TLambdax
       + 4320*Conj(Lambdax)*Kappa*Sqr(g1)*TLambdax + 21600*Conj(Lambdax)*Kappa*
      Sqr(g2)*TLambdax - 4320*Conj(Lambdax)*Kappa*Sqr(gN)*TLambdax - 28800*
      Kappa*Lambdax*Sqr(Conj(Lambdax))*TLambdax) + (-12*traceAdjKappaTKappa - 8
      *traceAdjLambda12TLambda12 - 5*MassBp*Sqr(gN) - 8*Conj(Lambdax)*TLambdax)
      *(Kappa*(Kappa).adjoint()*Kappa) + 0.5*(-18*traceKappaAdjKappa - 12*
      traceLambda12AdjLambda12 - 12*AbsSqr(Lambdax) + 7*Sqr(gN))*(Kappa*(Kappa)
      .adjoint()*TKappa) + (-9*traceKappaAdjKappa - 6*traceLambda12AdjLambda12
      - 6*AbsSqr(Lambdax) + 4*Sqr(gN))*(TKappa*(Kappa).adjoint()*Kappa) - 3*(
      Kappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint()*TKappa) - 4*(Kappa*(Kappa
      ).adjoint()*TKappa*(Kappa).adjoint()*Kappa) - 3*(TKappa*(Kappa).adjoint()
      *Kappa*(Kappa).adjoint()*Kappa)).real();


   return twoLoop * beta_TKappa;
}

/**
 * Calculates the 3-loop beta function of TKappa.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = ZEROMATRIX(3,3);


   return threeLoop * beta_TKappa;
}

/**
 * Calculates the 4-loop beta function of TKappa.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = ZEROMATRIX(3,3);


   return fourLoop * beta_TKappa;
}

/**
 * Calculates the 5-loop beta function of TKappa.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_TKappa_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_TKappa;

   beta_TKappa = ZEROMATRIX(3,3);


   return fiveLoop * beta_TKappa;
}

} // namespace flexiblesusy
