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

// File generated at Tue 22 Jan 2019 17:18:21

#include "NUTSMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of BMu.
 *
 * @return 1-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_BMu;

   beta_BMu = Re(0.2*oneOver16PiSqr*(15*traceYdAdjYd*BMu + 5*traceYeAdjYe*BMu +
      15*traceYuAdjYu*BMu + 30*AbsSqr(Lambdax)*BMu + 10*BMS*Conj(Kappa)*Lambdax
       + 30*traceAdjYdTYd*Mu + 10*traceAdjYeTYe*Mu + 30*traceAdjYuTYu*Mu - 3*
      BMu*Sqr(g1) + 6*MassB*Mu*Sqr(g1) - 15*BMu*Sqr(g2) + 30*MassWB*Mu*Sqr(g2)
      + 20*Conj(Lambdax)*Mu*TLambdax));


   return beta_BMu;
}

/**
 * Calculates the 2-loop beta function of BMu.
 *
 * @return 2-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;


   double beta_BMu;

   beta_BMu = Re(0.02*twoLoop*(-450*traceYdAdjYdYdAdjYd*BMu - 300*
      traceYdAdjYuYuAdjYd*BMu - 150*traceYeAdjYeYeAdjYe*BMu - 450*
      traceYuAdjYuYuAdjYu*BMu - 750*traceYdAdjYd*AbsSqr(Lambdax)*BMu - 250*
      traceYeAdjYe*AbsSqr(Lambdax)*BMu - 750*traceYuAdjYu*AbsSqr(Lambdax)*BMu -
      200*AbsSqr(Kappa)*AbsSqr(Lambdax)*BMu - 1800*traceYdAdjYdTYdAdjYd*Mu -
      600*traceYdAdjYuTYuAdjYd*Mu - 600*traceYeAdjYeTYeAdjYe*Mu - 600*
      traceYuAdjYdTYdAdjYu*Mu - 1800*traceYuAdjYuTYuAdjYu*Mu - 900*
      traceAdjYdTYd*AbsSqr(Lambdax)*Mu - 300*traceAdjYeTYe*AbsSqr(Lambdax)*Mu -
      900*traceAdjYuTYu*AbsSqr(Lambdax)*Mu + 207*BMu*Quad(g1) - 828*MassB*Mu*
      Quad(g1) + 375*BMu*Quad(g2) - 1500*MassWB*Mu*Quad(g2) - 20*traceYdAdjYd*
      BMu*Sqr(g1) + 60*traceYeAdjYe*BMu*Sqr(g1) + 40*traceYuAdjYu*BMu*Sqr(g1) +
      360*AbsSqr(Lambdax)*BMu*Sqr(g1) - 40*traceAdjYdTYd*Mu*Sqr(g1) + 120*
      traceAdjYeTYe*Mu*Sqr(g1) + 80*traceAdjYuTYu*Mu*Sqr(g1) + 40*MassB*
      traceYdAdjYd*Mu*Sqr(g1) - 120*MassB*traceYeAdjYe*Mu*Sqr(g1) - 80*MassB*
      traceYuAdjYu*Mu*Sqr(g1) - 360*MassB*AbsSqr(Lambdax)*Mu*Sqr(g1) + 1800*
      AbsSqr(Lambdax)*BMu*Sqr(g2) - 1800*MassWB*AbsSqr(Lambdax)*Mu*Sqr(g2) + 90
      *BMu*Sqr(g1)*Sqr(g2) - 180*MassB*Mu*Sqr(g1)*Sqr(g2) - 180*MassWB*Mu*Sqr(
      g1)*Sqr(g2) + 800*traceYdAdjYd*BMu*Sqr(g3) + 800*traceYuAdjYu*BMu*Sqr(g3)
      + 1600*traceAdjYdTYd*Mu*Sqr(g3) + 1600*traceAdjYuTYu*Mu*Sqr(g3) - 1600*
      MassG*traceYdAdjYd*Mu*Sqr(g3) - 1600*MassG*traceYuAdjYu*Mu*Sqr(g3) - 400*
      BMS*Kappa*Lambdax*Sqr(Conj(Kappa)) - 400*BMS*Conj(Kappa)*Conj(Lambdax)*
      Sqr(Lambdax) - 700*BMu*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 400*AbsSqr(
      Lambdax)*Conj(Kappa)*Mu*TKappa - 400*MS*Lambdax*Sqr(Conj(Kappa))*TKappa -
      400*MS*AbsSqr(Lambdax)*Conj(Kappa)*TLambdax - 300*traceYdAdjYd*Conj(
      Lambdax)*Mu*TLambdax - 100*traceYeAdjYe*Conj(Lambdax)*Mu*TLambdax - 300*
      traceYuAdjYu*Conj(Lambdax)*Mu*TLambdax - 400*AbsSqr(Kappa)*Conj(Lambdax)*
      Mu*TLambdax - 1600*Lambdax*Mu*Sqr(Conj(Lambdax))*TLambdax));


   return beta_BMu;
}

/**
 * Calculates the 3-loop beta function of BMu.
 *
 * @return 3-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMu;

   beta_BMu = 0;


   return beta_BMu;
}

/**
 * Calculates the 4-loop beta function of BMu.
 *
 * @return 4-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMu;

   beta_BMu = 0;


   return beta_BMu;
}

/**
 * Calculates the 5-loop beta function of BMu.
 *
 * @return 5-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMu;

   beta_BMu = 0;


   return beta_BMu;
}

} // namespace flexiblesusy
