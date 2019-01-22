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

// File generated at Tue 22 Jan 2019 14:42:31

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
 * Calculates the 1-loop beta function of TLambdax.
 *
 * @return 1-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_TLambdax_1_loop(const Soft_traces& soft_traces) const
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


   double beta_TLambdax;

   beta_TLambdax = Re(0.1*oneOver16PiSqr*(60*traceAdjKappaTKappa*Lambdax + 40*
      traceAdjLambda12TLambda12*Lambdax + 60*traceAdjYdTYd*Lambdax + 20*
      traceAdjYeTYe*Lambdax + 60*traceAdjYuTYu*Lambdax + 12*MassB*Lambdax*Sqr(
      g1) + 60*MassWB*Lambdax*Sqr(g2) + 38*MassBp*Lambdax*Sqr(gN) + 30*
      traceKappaAdjKappa*TLambdax + 20*traceLambda12AdjLambda12*TLambdax + 30*
      traceYdAdjYd*TLambdax + 10*traceYeAdjYe*TLambdax + 30*traceYuAdjYu*
      TLambdax + 120*AbsSqr(Lambdax)*TLambdax - 6*Sqr(g1)*TLambdax - 30*Sqr(g2)
      *TLambdax - 19*Sqr(gN)*TLambdax));


   return beta_TLambdax;
}

/**
 * Calculates the 2-loop beta function of TLambdax.
 *
 * @return 2-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_TLambdax_2_loop(const Soft_traces& soft_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12TLambda12AdjLambda12;


   double beta_TLambdax;

   beta_TLambdax = Re(0.005*twoLoop*(-4800*traceKappaAdjKappaTKappaAdjKappa*
      Lambdax - 3200*traceLambda12AdjLambda12TLambda12AdjLambda12*Lambdax -
      7200*traceYdAdjYdTYdAdjYd*Lambdax - 2400*traceYdAdjYuTYuAdjYd*Lambdax -
      2400*traceYeAdjYeTYeAdjYe*Lambdax - 2400*traceYuAdjYdTYdAdjYu*Lambdax -
      7200*traceYuAdjYuTYuAdjYu*Lambdax - 4752*MassB*Lambdax*Quad(g1) - 13200*
      MassWB*Lambdax*Quad(g2) - 15732*MassBp*Lambdax*Quad(gN) + 320*
      traceAdjKappaTKappa*Lambdax*Sqr(g1) + 480*traceAdjLambda12TLambda12*
      Lambdax*Sqr(g1) - 160*traceAdjYdTYd*Lambdax*Sqr(g1) + 480*traceAdjYeTYe*
      Lambdax*Sqr(g1) + 320*traceAdjYuTYu*Lambdax*Sqr(g1) - 320*MassB*
      traceKappaAdjKappa*Lambdax*Sqr(g1) - 480*MassB*traceLambda12AdjLambda12*
      Lambdax*Sqr(g1) + 160*MassB*traceYdAdjYd*Lambdax*Sqr(g1) - 480*MassB*
      traceYeAdjYe*Lambdax*Sqr(g1) - 320*MassB*traceYuAdjYu*Lambdax*Sqr(g1) +
      2400*traceAdjLambda12TLambda12*Lambdax*Sqr(g2) - 2400*MassWB*
      traceLambda12AdjLambda12*Lambdax*Sqr(g2) - 720*MassB*Lambdax*Sqr(g1)*Sqr(
      g2) - 720*MassWB*Lambdax*Sqr(g1)*Sqr(g2) + 6400*traceAdjKappaTKappa*
      Lambdax*Sqr(g3) + 6400*traceAdjYdTYd*Lambdax*Sqr(g3) + 6400*traceAdjYuTYu
      *Lambdax*Sqr(g3) - 6400*MassG*traceKappaAdjKappa*Lambdax*Sqr(g3) - 6400*
      MassG*traceYdAdjYd*Lambdax*Sqr(g3) - 6400*MassG*traceYuAdjYu*Lambdax*Sqr(
      g3) - 720*traceAdjKappaTKappa*Lambdax*Sqr(gN) - 480*
      traceAdjLambda12TLambda12*Lambdax*Sqr(gN) - 240*traceAdjYdTYd*Lambdax*Sqr
      (gN) - 80*traceAdjYeTYe*Lambdax*Sqr(gN) - 120*traceAdjYuTYu*Lambdax*Sqr(
      gN) + 720*MassBp*traceKappaAdjKappa*Lambdax*Sqr(gN) + 480*MassBp*
      traceLambda12AdjLambda12*Lambdax*Sqr(gN) + 240*MassBp*traceYdAdjYd*
      Lambdax*Sqr(gN) + 80*MassBp*traceYeAdjYe*Lambdax*Sqr(gN) + 120*MassBp*
      traceYuAdjYu*Lambdax*Sqr(gN) - 108*MassB*Lambdax*Sqr(g1)*Sqr(gN) - 108*
      MassBp*Lambdax*Sqr(g1)*Sqr(gN) - 780*MassBp*Lambdax*Sqr(g2)*Sqr(gN) - 780
      *MassWB*Lambdax*Sqr(g2)*Sqr(gN) - 2400*traceAdjKappaTKappa*Conj(Lambdax)*
      Sqr(Lambdax) - 1600*traceAdjLambda12TLambda12*Conj(Lambdax)*Sqr(Lambdax)
      - 3600*traceAdjYdTYd*Conj(Lambdax)*Sqr(Lambdax) - 1200*traceAdjYeTYe*Conj
      (Lambdax)*Sqr(Lambdax) - 3600*traceAdjYuTYu*Conj(Lambdax)*Sqr(Lambdax) -
      480*MassB*Conj(Lambdax)*Sqr(g1)*Sqr(Lambdax) - 2400*MassWB*Conj(Lambdax)*
      Sqr(g2)*Sqr(Lambdax) - 520*MassBp*Conj(Lambdax)*Sqr(gN)*Sqr(Lambdax) -
      1200*traceKappaAdjKappaKappaAdjKappa*TLambdax - 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12*TLambdax - 1800*
      traceYdAdjYdYdAdjYd*TLambdax - 1200*traceYdAdjYuYuAdjYd*TLambdax - 600*
      traceYeAdjYeYeAdjYe*TLambdax - 1800*traceYuAdjYuYuAdjYu*TLambdax - 3600*
      traceKappaAdjKappa*AbsSqr(Lambdax)*TLambdax - 2400*
      traceLambda12AdjLambda12*AbsSqr(Lambdax)*TLambdax - 5400*traceYdAdjYd*
      AbsSqr(Lambdax)*TLambdax - 1800*traceYeAdjYe*AbsSqr(Lambdax)*TLambdax -
      5400*traceYuAdjYu*AbsSqr(Lambdax)*TLambdax + 1188*Quad(g1)*TLambdax +
      3300*Quad(g2)*TLambdax + 3933*Quad(gN)*TLambdax + 160*traceKappaAdjKappa*
      Sqr(g1)*TLambdax + 240*traceLambda12AdjLambda12*Sqr(g1)*TLambdax - 80*
      traceYdAdjYd*Sqr(g1)*TLambdax + 240*traceYeAdjYe*Sqr(g1)*TLambdax + 160*
      traceYuAdjYu*Sqr(g1)*TLambdax + 720*AbsSqr(Lambdax)*Sqr(g1)*TLambdax +
      1200*traceLambda12AdjLambda12*Sqr(g2)*TLambdax + 3600*AbsSqr(Lambdax)*Sqr
      (g2)*TLambdax + 360*Sqr(g1)*Sqr(g2)*TLambdax + 3200*traceKappaAdjKappa*
      Sqr(g3)*TLambdax + 3200*traceYdAdjYd*Sqr(g3)*TLambdax + 3200*traceYuAdjYu
      *Sqr(g3)*TLambdax - 360*traceKappaAdjKappa*Sqr(gN)*TLambdax - 240*
      traceLambda12AdjLambda12*Sqr(gN)*TLambdax - 120*traceYdAdjYd*Sqr(gN)*
      TLambdax - 40*traceYeAdjYe*Sqr(gN)*TLambdax - 60*traceYuAdjYu*Sqr(gN)*
      TLambdax + 780*AbsSqr(Lambdax)*Sqr(gN)*TLambdax + 54*Sqr(g1)*Sqr(gN)*
      TLambdax + 390*Sqr(g2)*Sqr(gN)*TLambdax - 10000*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)*TLambdax));


   return beta_TLambdax;
}

/**
 * Calculates the 3-loop beta function of TLambdax.
 *
 * @return 3-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_TLambdax_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

/**
 * Calculates the 4-loop beta function of TLambdax.
 *
 * @return 4-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_TLambdax_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

/**
 * Calculates the 5-loop beta function of TLambdax.
 *
 * @return 5-loop beta function
 */
double E6SSMEFTHiggs_soft_parameters::calc_beta_TLambdax_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

} // namespace flexiblesusy
