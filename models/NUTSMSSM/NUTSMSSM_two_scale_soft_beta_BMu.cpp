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

// File generated at Fri 8 Jan 2016 15:15:35

#include "NUTSMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of BMu.
 *
 * @return one-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_BMu;

   beta_BMu = Re(oneOver16PiSqr*(2*BMS*Conj(Kappa)*Lambdax + 6*
      traceAdjYdTYd*Mu + 2*traceAdjYeTYe*Mu + 6*traceAdjYuTYu*Mu + 1.2*MassB*Mu
      *Sqr(g1) + BMu*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + 6*AbsSqr
      (Lambdax) - 0.6*Sqr(g1) - 3*Sqr(g2)) + 6*MassWB*Mu*Sqr(g2) + 4*Conj(
      Lambdax)*Mu*TLambdax));


   return beta_BMu;
}

/**
 * Calculates the two-loop beta function of BMu.
 *
 * @return two-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;


   double beta_BMu;

   beta_BMu = Re(0.02*twoLoop*(BMu*(207*Power(g1,4) + 375*Power(g2,4) -
      450*traceYdAdjYdYdAdjYd - 300*traceYdAdjYuYuAdjYd - 150*
      traceYeAdjYeYeAdjYe - 450*traceYuAdjYuYuAdjYu + 60*traceYeAdjYe*Sqr(g1) +
      40*traceYuAdjYu*Sqr(g1) + 90*Sqr(g1)*Sqr(g2) + 10*AbsSqr(Lambdax)*(-75*
      traceYdAdjYd - 25*traceYeAdjYe - 75*traceYuAdjYu - 20*AbsSqr(Kappa) + 36*
      Sqr(g1) + 180*Sqr(g2)) - 20*traceYdAdjYd*(Sqr(g1) - 40*Sqr(g3)) + 800*
      traceYuAdjYu*Sqr(g3) - 700*Sqr(Conj(Lambdax))*Sqr(Lambdax)) - 4*(100*(
      AbsSqr(Kappa) + AbsSqr(Lambdax))*BMS*Conj(Kappa)*Lambdax + 207*Power(g1,4
      )*MassB*Mu + 375*Power(g2,4)*MassWB*Mu + 450*traceYdAdjYdTYdAdjYd*Mu +
      150*traceYdAdjYuTYuAdjYd*Mu + 150*traceYeAdjYeTYeAdjYe*Mu + 150*
      traceYuAdjYdTYdAdjYu*Mu + 450*traceYuAdjYuTYuAdjYu*Mu + 10*traceAdjYdTYd*
      Mu*Sqr(g1) - 30*traceAdjYeTYe*Mu*Sqr(g1) - 20*traceAdjYuTYu*Mu*Sqr(g1) -
      10*MassB*traceYdAdjYd*Mu*Sqr(g1) + 30*MassB*traceYeAdjYe*Mu*Sqr(g1) + 20*
      MassB*traceYuAdjYu*Mu*Sqr(g1) + 45*MassB*Mu*Sqr(g1)*Sqr(g2) + 45*MassWB*
      Mu*Sqr(g1)*Sqr(g2) - 400*traceAdjYdTYd*Mu*Sqr(g3) - 400*traceAdjYuTYu*Mu*
      Sqr(g3) + 400*MassG*traceYdAdjYd*Mu*Sqr(g3) + 400*MassG*traceYuAdjYu*Mu*
      Sqr(g3) + 100*MS*Lambdax*Sqr(Conj(Kappa))*TKappa + 400*Lambdax*Mu*Sqr(
      Conj(Lambdax))*TLambdax + 5*Conj(Lambdax)*(Mu*(3*Lambdax*(15*
      traceAdjYdTYd + 5*traceAdjYeTYe + 15*traceAdjYuTYu + 6*MassB*Sqr(g1) + 30
      *MassWB*Sqr(g2)) + 5*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)*
      TLambdax) + 20*Conj(Kappa)*(Lambdax*Mu*TKappa + (MS*Lambdax + Kappa*Mu)*
      TLambdax)))));


   return beta_BMu;
}

/**
 * Calculates the three-loop beta function of BMu.
 *
 * @return three-loop beta function
 */
double NUTSMSSM_soft_parameters::calc_beta_BMu_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMu;

   beta_BMu = 0;


   return beta_BMu;
}

} // namespace flexiblesusy