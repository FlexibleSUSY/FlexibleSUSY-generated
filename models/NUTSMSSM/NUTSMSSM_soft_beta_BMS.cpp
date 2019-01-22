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

// File generated at Tue 22 Jan 2019 17:18:22

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
 * Calculates the 1-loop beta function of BMS.
 *
 * @return 1-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMS_1_loop(const Soft_traces& soft_traces) const
{


   double beta_BMS;

   beta_BMS = Re(4*oneOver16PiSqr*(2*AbsSqr(Kappa)*BMS + AbsSqr(Lambdax)*BMS +
      2*BMu*Conj(Lambdax)*Kappa + 2*MS*Conj(Kappa)*TKappa + 2*MS*Conj(Lambdax)*
      TLambdax));


   return beta_BMS;
}

/**
 * Calculates the 2-loop beta function of BMS.
 *
 * @return 2-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMS_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_BMS;

   beta_BMS = Re(-0.8*twoLoop*(30*MS*traceAdjYdTYd*AbsSqr(Lambdax) + 10*MS*
      traceAdjYeTYe*AbsSqr(Lambdax) + 30*MS*traceAdjYuTYu*AbsSqr(Lambdax) + 15*
      traceYdAdjYd*AbsSqr(Lambdax)*BMS + 5*traceYeAdjYe*AbsSqr(Lambdax)*BMS +
      15*traceYuAdjYu*AbsSqr(Lambdax)*BMS + 40*AbsSqr(Kappa)*AbsSqr(Lambdax)*
      BMS + 30*traceYdAdjYd*BMu*Conj(Lambdax)*Kappa + 10*traceYeAdjYe*BMu*Conj(
      Lambdax)*Kappa + 30*traceYuAdjYu*BMu*Conj(Lambdax)*Kappa + 30*
      traceAdjYdTYd*Conj(Lambdax)*Kappa*Mu + 10*traceAdjYeTYe*Conj(Lambdax)*
      Kappa*Mu + 30*traceAdjYuTYu*Conj(Lambdax)*Kappa*Mu + 6*MassB*MS*AbsSqr(
      Lambdax)*Sqr(g1) - 3*AbsSqr(Lambdax)*BMS*Sqr(g1) - 18*BMu*Conj(Lambdax)*
      Kappa*Sqr(g1) + 18*MassB*Conj(Lambdax)*Kappa*Mu*Sqr(g1) + 30*MassWB*MS*
      AbsSqr(Lambdax)*Sqr(g2) - 15*AbsSqr(Lambdax)*BMS*Sqr(g2) - 90*BMu*Conj(
      Lambdax)*Kappa*Sqr(g2) + 90*MassWB*Conj(Lambdax)*Kappa*Mu*Sqr(g2) + 20*
      BMu*Kappa*Lambdax*Sqr(Conj(Lambdax)) + 40*BMS*Sqr(Conj(Kappa))*Sqr(Kappa)
      + 10*BMS*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 40*MS*AbsSqr(Lambdax)*Conj(
      Kappa)*TKappa + 100*MS*Kappa*Sqr(Conj(Kappa))*TKappa + 30*MS*traceYdAdjYd
      *Conj(Lambdax)*TLambdax + 10*MS*traceYeAdjYe*Conj(Lambdax)*TLambdax + 30*
      MS*traceYuAdjYu*Conj(Lambdax)*TLambdax + 60*MS*AbsSqr(Kappa)*Conj(Lambdax
      )*TLambdax - 6*MS*Conj(Lambdax)*Sqr(g1)*TLambdax - 30*MS*Conj(Lambdax)*
      Sqr(g2)*TLambdax + 40*MS*Lambdax*Sqr(Conj(Lambdax))*TLambdax + 20*Kappa*
      Mu*Sqr(Conj(Lambdax))*TLambdax));


   return beta_BMS;
}

/**
 * Calculates the 3-loop beta function of BMS.
 *
 * @return 3-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMS_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMS;

   beta_BMS = 0;


   return beta_BMS;
}

/**
 * Calculates the 4-loop beta function of BMS.
 *
 * @return 4-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMS_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMS;

   beta_BMS = 0;


   return beta_BMS;
}

/**
 * Calculates the 5-loop beta function of BMS.
 *
 * @return 5-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMS_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMS;

   beta_BMS = 0;


   return beta_BMS;
}

} // namespace flexiblesusy
