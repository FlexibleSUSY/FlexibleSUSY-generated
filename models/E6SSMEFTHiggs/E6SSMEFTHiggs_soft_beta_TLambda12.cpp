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
 * Calculates the 1-loop beta function of TLambda12.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_TLambda12_1_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 = TRACE_STRUCT.
      traceAdjLambda12TLambda12;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = (0.1*(60*traceAdjKappaTKappa*Lambda12 + 40*
      traceAdjLambda12TLambda12*Lambda12 + 12*MassB*Lambda12*Sqr(g1) + 60*
      MassWB*Lambda12*Sqr(g2) + 38*MassBp*Lambda12*Sqr(gN) + 40*Conj(Lambdax)*
      Lambda12*TLambdax + 30*traceKappaAdjKappa*TLambda12 + 20*
      traceLambda12AdjLambda12*TLambda12 + 20*AbsSqr(Lambdax)*TLambda12 - 6*Sqr
      (g1)*TLambda12 - 30*Sqr(g2)*TLambda12 - 19*Sqr(gN)*TLambda12) + 3*(
      Lambda12*(Lambda12).adjoint()*TLambda12) + 3*(TLambda12*(Lambda12).
      adjoint()*Lambda12)).real();


   return oneLoop * beta_TLambda12;
}

/**
 * Calculates the 2-loop beta function of TLambda12.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_TLambda12_2_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = (0.005*(-4800*traceKappaAdjKappaTKappaAdjKappa*Lambda12 -
      3200*traceLambda12AdjLambda12TLambda12AdjLambda12*Lambda12 - 2400*
      traceAdjYdTYd*AbsSqr(Lambdax)*Lambda12 - 800*traceAdjYeTYe*AbsSqr(Lambdax
      )*Lambda12 - 2400*traceAdjYuTYu*AbsSqr(Lambdax)*Lambda12 - 4752*MassB*
      Lambda12*Quad(g1) - 13200*MassWB*Lambda12*Quad(g2) - 15732*MassBp*
      Lambda12*Quad(gN) + 320*traceAdjKappaTKappa*Lambda12*Sqr(g1) + 480*
      traceAdjLambda12TLambda12*Lambda12*Sqr(g1) - 320*MassB*traceKappaAdjKappa
      *Lambda12*Sqr(g1) - 480*MassB*traceLambda12AdjLambda12*Lambda12*Sqr(g1) -
      480*MassB*AbsSqr(Lambdax)*Lambda12*Sqr(g1) + 2400*
      traceAdjLambda12TLambda12*Lambda12*Sqr(g2) - 2400*MassWB*
      traceLambda12AdjLambda12*Lambda12*Sqr(g2) - 2400*MassWB*AbsSqr(Lambdax)*
      Lambda12*Sqr(g2) - 720*MassB*Lambda12*Sqr(g1)*Sqr(g2) - 720*MassWB*
      Lambda12*Sqr(g1)*Sqr(g2) + 6400*traceAdjKappaTKappa*Lambda12*Sqr(g3) -
      6400*MassG*traceKappaAdjKappa*Lambda12*Sqr(g3) - 720*traceAdjKappaTKappa*
      Lambda12*Sqr(gN) - 480*traceAdjLambda12TLambda12*Lambda12*Sqr(gN) + 720*
      MassBp*traceKappaAdjKappa*Lambda12*Sqr(gN) + 480*MassBp*
      traceLambda12AdjLambda12*Lambda12*Sqr(gN) + 480*MassBp*AbsSqr(Lambdax)*
      Lambda12*Sqr(gN) - 108*MassB*Lambda12*Sqr(g1)*Sqr(gN) - 108*MassBp*
      Lambda12*Sqr(g1)*Sqr(gN) - 780*MassBp*Lambda12*Sqr(g2)*Sqr(gN) - 780*
      MassWB*Lambda12*Sqr(g2)*Sqr(gN) - 2400*traceYdAdjYd*Conj(Lambdax)*
      Lambda12*TLambdax - 800*traceYeAdjYe*Conj(Lambdax)*Lambda12*TLambdax -
      2400*traceYuAdjYu*Conj(Lambdax)*Lambda12*TLambdax + 480*Conj(Lambdax)*
      Lambda12*Sqr(g1)*TLambdax + 2400*Conj(Lambdax)*Lambda12*Sqr(g2)*TLambdax
      - 480*Conj(Lambdax)*Lambda12*Sqr(gN)*TLambdax - 3200*Lambdax*Lambda12*Sqr
      (Conj(Lambdax))*TLambdax - 1200*traceKappaAdjKappaKappaAdjKappa*TLambda12
       - 800*traceLambda12AdjLambda12Lambda12AdjLambda12*TLambda12 - 1200*
      traceYdAdjYd*AbsSqr(Lambdax)*TLambda12 - 400*traceYeAdjYe*AbsSqr(Lambdax)
      *TLambda12 - 1200*traceYuAdjYu*AbsSqr(Lambdax)*TLambda12 + 1188*Quad(g1)*
      TLambda12 + 3300*Quad(g2)*TLambda12 + 3933*Quad(gN)*TLambda12 + 160*
      traceKappaAdjKappa*Sqr(g1)*TLambda12 + 240*traceLambda12AdjLambda12*Sqr(
      g1)*TLambda12 + 240*AbsSqr(Lambdax)*Sqr(g1)*TLambda12 + 1200*
      traceLambda12AdjLambda12*Sqr(g2)*TLambda12 + 1200*AbsSqr(Lambdax)*Sqr(g2)
      *TLambda12 + 360*Sqr(g1)*Sqr(g2)*TLambda12 + 3200*traceKappaAdjKappa*Sqr(
      g3)*TLambda12 - 360*traceKappaAdjKappa*Sqr(gN)*TLambda12 - 240*
      traceLambda12AdjLambda12*Sqr(gN)*TLambda12 - 240*AbsSqr(Lambdax)*Sqr(gN)*
      TLambda12 + 54*Sqr(g1)*Sqr(gN)*TLambda12 + 390*Sqr(g2)*Sqr(gN)*TLambda12
      - 800*Sqr(Conj(Lambdax))*Sqr(Lambdax)*TLambda12) + (-12*
      traceAdjKappaTKappa - 8*traceAdjLambda12TLambda12 - 5*MassBp*Sqr(gN) - 8*
      Conj(Lambdax)*TLambdax)*(Lambda12*(Lambda12).adjoint()*Lambda12) + 0.5*(-
      18*traceKappaAdjKappa - 12*traceLambda12AdjLambda12 - 12*AbsSqr(Lambdax)
      + 7*Sqr(gN))*(Lambda12*(Lambda12).adjoint()*TLambda12) + (-9*
      traceKappaAdjKappa - 6*traceLambda12AdjLambda12 - 6*AbsSqr(Lambdax) + 4*
      Sqr(gN))*(TLambda12*(Lambda12).adjoint()*Lambda12) - 3*(Lambda12*(
      Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*TLambda12) - 4*(
      Lambda12*(Lambda12).adjoint()*TLambda12*(Lambda12).adjoint()*Lambda12) -
      3*(TLambda12*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12)
      ).real();


   return twoLoop * beta_TLambda12;
}

/**
 * Calculates the 3-loop beta function of TLambda12.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_TLambda12_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = ZEROMATRIX(2,2);


   return threeLoop * beta_TLambda12;
}

/**
 * Calculates the 4-loop beta function of TLambda12.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_TLambda12_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = ZEROMATRIX(2,2);


   return fourLoop * beta_TLambda12;
}

/**
 * Calculates the 5-loop beta function of TLambda12.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_TLambda12_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,2,2> beta_TLambda12;

   beta_TLambda12 = ZEROMATRIX(2,2);


   return fiveLoop * beta_TLambda12;
}

} // namespace flexiblesusy
