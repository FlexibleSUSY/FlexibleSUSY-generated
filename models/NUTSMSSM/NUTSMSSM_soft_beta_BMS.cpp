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

// File generated at Mon 5 Mar 2018 18:25:02

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

   beta_BMS = Re(oneOver16PiSqr*(4*(2*AbsSqr(Kappa) + AbsSqr(Lambdax))*
      BMS + 8*(BMu*Conj(Lambdax)*Kappa + MS*Conj(Kappa)*TKappa + MS*Conj(
      Lambdax)*TLambdax)));


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

   beta_BMS = Re(-0.8*twoLoop*(BMS*(40*AbsSqr(Kappa)*AbsSqr(Lambdax) +
      AbsSqr(Lambdax)*(10*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu - 3*Sqr(g2))) + 40*Sqr(Conj(Kappa))*Sqr(
      Kappa)) + 2*(50*MS*Kappa*Sqr(Conj(Kappa))*TKappa + 10*Sqr(Conj(Lambdax))*
      (BMu*Kappa*Lambdax + (2*MS*Lambdax + Kappa*Mu)*TLambdax) + Conj(Lambdax)*
      (15*MS*traceAdjYdTYd*Lambdax + 5*MS*traceAdjYeTYe*Lambdax + 15*MS*
      traceAdjYuTYu*Lambdax + 15*traceAdjYdTYd*Kappa*Mu + 5*traceAdjYeTYe*Kappa
      *Mu + 15*traceAdjYuTYu*Kappa*Mu + 3*MassB*MS*Lambdax*Sqr(g1) + 9*MassB*
      Kappa*Mu*Sqr(g1) + BMu*Kappa*(-9*Sqr(g1) + 5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu - 9*Sqr(g2))) + 15*MassWB*MS*Lambdax*Sqr(g2
      ) + 45*MassWB*Kappa*Mu*Sqr(g2) + 15*MS*traceYdAdjYd*TLambdax + 5*MS*
      traceYeAdjYe*TLambdax + 15*MS*traceYuAdjYu*TLambdax - 3*MS*Sqr(g1)*
      TLambdax - 15*MS*Sqr(g2)*TLambdax + 10*MS*Conj(Kappa)*(2*Lambdax*TKappa +
      3*Kappa*TLambdax)))));


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

} // namespace flexiblesusy
