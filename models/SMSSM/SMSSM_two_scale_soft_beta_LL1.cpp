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

// File generated at Mon 27 Feb 2017 13:34:57

#include "SMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of LL1.
 *
 * @return one-loop beta function
 */
double SMSSM_soft_parameters::calc_beta_LL1_one_loop(const Soft_traces& soft_traces) const
{


   double beta_LL1;

   beta_LL1 = Re(2*oneOver16PiSqr*(MS*BMS*Conj(Kappa) + 2*MS*BMu*Conj(
      Lambdax) + 2*ms2*Conj(MS)*Kappa + 2*mHd2*Conj(Mu)*Lambdax + 2*mHu2*Conj(
      Mu)*Lambdax + AbsSqr(Kappa)*LL1 + AbsSqr(Lambdax)*LL1 + Conj(BMS)*TKappa
      + 2*L1*Conj(Kappa)*TKappa + 2*Conj(BMu)*TLambdax + 2*L1*Conj(Lambdax)*
      TLambdax));


   return beta_LL1;
}

/**
 * Calculates the two-loop beta function of LL1.
 *
 * @return two-loop beta function
 */
double SMSSM_soft_parameters::calc_beta_LL1_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceYdAdjYdconjmd2 = TRACE_STRUCT.traceYdAdjYdconjmd2;
   const double traceYdconjmq2AdjYd = TRACE_STRUCT.traceYdconjmq2AdjYd;
   const double traceYeAdjYeconjme2 = TRACE_STRUCT.traceYeAdjYeconjme2;
   const double traceYeconjml2AdjYe = TRACE_STRUCT.traceYeconjml2AdjYe;
   const double traceYuAdjYuconjmu2 = TRACE_STRUCT.traceYuAdjYuconjmu2;
   const double traceYuconjmq2AdjYu = TRACE_STRUCT.traceYuconjmq2AdjYu;


   double beta_LL1;

   beta_LL1 = Re(-0.4*twoLoop*(30*L1*traceAdjYdTYd*AbsSqr(Lambdax) + 10*
      L1*traceAdjYeTYe*AbsSqr(Lambdax) + 30*L1*traceAdjYuTYu*AbsSqr(Lambdax) +
      20*MS*(AbsSqr(Kappa) + AbsSqr(Lambdax))*BMS*Conj(Kappa) + 30*MS*
      traceYdAdjYd*BMu*Conj(Lambdax) + 10*MS*traceYeAdjYe*BMu*Conj(Lambdax) +
      30*MS*traceYuAdjYu*BMu*Conj(Lambdax) + 30*traceAdjYdTYd*Conj(BMu)*Lambdax
      + 10*traceAdjYeTYe*Conj(BMu)*Lambdax + 30*traceAdjYuTYu*Conj(BMu)*
      Lambdax + 30*traceconjTYdTpTYd*Conj(Mu)*Lambdax + 10*traceconjTYeTpTYe*
      Conj(Mu)*Lambdax + 30*traceconjTYuTpTYu*Conj(Mu)*Lambdax + 60*mHd2*
      traceYdAdjYd*Conj(Mu)*Lambdax + 30*mHu2*traceYdAdjYd*Conj(Mu)*Lambdax +
      30*traceYdAdjYdconjmd2*Conj(Mu)*Lambdax + 30*traceYdconjmq2AdjYd*Conj(Mu)
      *Lambdax + 20*mHd2*traceYeAdjYe*Conj(Mu)*Lambdax + 10*mHu2*traceYeAdjYe*
      Conj(Mu)*Lambdax + 10*traceYeAdjYeconjme2*Conj(Mu)*Lambdax + 10*
      traceYeconjml2AdjYe*Conj(Mu)*Lambdax + 30*mHd2*traceYuAdjYu*Conj(Mu)*
      Lambdax + 60*mHu2*traceYuAdjYu*Conj(Mu)*Lambdax + 30*traceYuAdjYuconjmu2*
      Conj(Mu)*Lambdax + 30*traceYuconjmq2AdjYu*Conj(Mu)*Lambdax + 40*AbsSqr(
      TLambdax)*Conj(Mu)*Lambdax + 30*MS*traceAdjYdTYd*Conj(Lambdax)*Mu + 10*MS
      *traceAdjYeTYe*Conj(Lambdax)*Mu + 30*MS*traceAdjYuTYu*Conj(Lambdax)*Mu +
      15*traceYdAdjYd*AbsSqr(Lambdax)*LL1 + 5*traceYeAdjYe*AbsSqr(Lambdax)*LL1
      + 15*traceYuAdjYu*AbsSqr(Lambdax)*LL1 + 20*AbsSqr(Kappa)*AbsSqr(Lambdax)*
      LL1 + 6*L1*MassB*AbsSqr(Lambdax)*Sqr(g1) - 6*MS*BMu*Conj(Lambdax)*Sqr(g1)
      + 6*MassB*Conj(BMu)*Lambdax*Sqr(g1) - 6*mHd2*Conj(Mu)*Lambdax*Sqr(g1) -
      6*mHu2*Conj(Mu)*Lambdax*Sqr(g1) - 12*AbsSqr(MassB)*Conj(Mu)*Lambdax*Sqr(
      g1) + 6*MassB*MS*Conj(Lambdax)*Mu*Sqr(g1) - 3*AbsSqr(Lambdax)*LL1*Sqr(g1)
      + 30*L1*MassWB*AbsSqr(Lambdax)*Sqr(g2) - 30*MS*BMu*Conj(Lambdax)*Sqr(g2)
      + 30*MassWB*Conj(BMu)*Lambdax*Sqr(g2) - 30*mHd2*Conj(Mu)*Lambdax*Sqr(g2)
      - 30*mHu2*Conj(Mu)*Lambdax*Sqr(g2) - 60*AbsSqr(MassWB)*Conj(Mu)*Lambdax*
      Sqr(g2) + 30*MassWB*MS*Conj(Lambdax)*Mu*Sqr(g2) - 15*AbsSqr(Lambdax)*LL1*
      Sqr(g2) + 20*MS*BMu*Lambdax*Sqr(Conj(Lambdax)) + 20*LL1*Sqr(Conj(Kappa))*
      Sqr(Kappa) + 40*mHd2*Conj(Lambdax)*Conj(Mu)*Sqr(Lambdax) + 40*mHu2*Conj(
      Lambdax)*Conj(Mu)*Sqr(Lambdax) + 20*ms2*Conj(Lambdax)*Conj(Mu)*Sqr(
      Lambdax) + 10*LL1*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 40*AbsSqr(Kappa)*Conj
      (BMS)*TKappa + 20*AbsSqr(Lambdax)*Conj(BMS)*TKappa + 40*L1*AbsSqr(Lambdax
      )*Conj(Kappa)*TKappa + 80*L1*Kappa*Sqr(Conj(Kappa))*TKappa + 20*Sqr(MS)*
      Sqr(Conj(Kappa))*TKappa + 20*Conj(MS)*((mHd2 + mHu2 + 3*ms2)*AbsSqr(
      Lambdax)*Kappa + 2*AbsSqr(TKappa)*Kappa + AbsSqr(TLambdax)*Kappa + 5*ms2*
      Conj(Kappa)*Sqr(Kappa) + Conj(TLambdax)*Lambdax*TKappa) + 30*traceYdAdjYd
      *Conj(BMu)*TLambdax + 10*traceYeAdjYe*Conj(BMu)*TLambdax + 30*
      traceYuAdjYu*Conj(BMu)*TLambdax + 40*AbsSqr(Lambdax)*Conj(BMu)*TLambdax +
      30*L1*traceYdAdjYd*Conj(Lambdax)*TLambdax + 10*L1*traceYeAdjYe*Conj(
      Lambdax)*TLambdax + 30*L1*traceYuAdjYu*Conj(Lambdax)*TLambdax + 40*L1*
      AbsSqr(Kappa)*Conj(Lambdax)*TLambdax + 30*traceconjTYdTpYd*Conj(Mu)*
      TLambdax + 10*traceconjTYeTpYe*Conj(Mu)*TLambdax + 30*traceconjTYuTpYu*
      Conj(Mu)*TLambdax + 20*Conj(BMS)*Conj(Lambdax)*Kappa*TLambdax - 6*Conj(
      BMu)*Sqr(g1)*TLambdax - 6*L1*Conj(Lambdax)*Sqr(g1)*TLambdax + 6*MassB*
      Conj(Mu)*Sqr(g1)*TLambdax - 30*Conj(BMu)*Sqr(g2)*TLambdax - 30*L1*Conj(
      Lambdax)*Sqr(g2)*TLambdax + 30*MassWB*Conj(Mu)*Sqr(g2)*TLambdax + 20*Conj
      (Kappa)*Conj(Lambdax)*Sqr(MS)*TLambdax + 40*L1*Lambdax*Sqr(Conj(Lambdax))
      *TLambdax + 20*MS*Mu*Sqr(Conj(Lambdax))*TLambdax));


   return beta_LL1;
}

/**
 * Calculates the three-loop beta function of LL1.
 *
 * @return three-loop beta function
 */
double SMSSM_soft_parameters::calc_beta_LL1_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LL1;

   beta_LL1 = 0;


   return beta_LL1;
}

} // namespace flexiblesusy
