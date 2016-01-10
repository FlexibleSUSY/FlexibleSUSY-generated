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

// File generated at Sun 10 Jan 2016 15:30:59

#include "TMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of TLambdax.
 *
 * @return one-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_TLambdax_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_TLambdax;

   beta_TLambdax = Re(oneOver16PiSqr*(0.4*Lambdax*(15*traceAdjYdTYd + 5*
      traceAdjYeTYe + 15*traceAdjYuTYu + 3*MassB*Sqr(g1) + 35*MassWB*Sqr(g2)) +
      (3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + 12*AbsSqr(Lambdax) -
      0.6*Sqr(g1) - 7*Sqr(g2))*TLambdax));


   return beta_TLambdax;
}

/**
 * Calculates the two-loop beta function of TLambdax.
 *
 * @return two-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_TLambdax_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYuAdjYuTYuAdjYu = TRACE_STRUCT.traceYuAdjYuTYuAdjYu;


   double beta_TLambdax;

   beta_TLambdax = Re(twoLoop*(-0.08*Lambdax*(207*Power(g1,4)*MassB +
      2075*Power(g2,4)*MassWB + 450*traceYdAdjYdTYdAdjYd + 150*
      traceYdAdjYuTYuAdjYd + 150*traceYeAdjYeTYeAdjYe + 150*
      traceYuAdjYdTYdAdjYu + 450*traceYuAdjYuTYuAdjYu + 10*traceAdjYdTYd*Sqr(g1
      ) - 30*traceAdjYeTYe*Sqr(g1) - 20*traceAdjYuTYu*Sqr(g1) + 30*MassB*
      traceYeAdjYe*Sqr(g1) + 20*MassB*traceYuAdjYu*Sqr(g1) + 45*MassB*Sqr(g1)*
      Sqr(g2) + 45*MassWB*Sqr(g1)*Sqr(g2) - 400*traceAdjYdTYd*Sqr(g3) - 400*
      traceAdjYuTYu*Sqr(g3) + 400*MassG*traceYuAdjYu*Sqr(g3) - 10*traceYdAdjYd*
      (MassB*Sqr(g1) - 40*MassG*Sqr(g3))) + (4.14*Power(g1,4) + 41.5*Power(g2,4
      ) - 9*traceYdAdjYdYdAdjYd - 6*traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe
      - 9*traceYuAdjYuYuAdjYu + 1.2*traceYeAdjYe*Sqr(g1) + 0.8*traceYuAdjYu*
      Sqr(g1) + 1.8*Sqr(g1)*Sqr(g2) - 0.4*traceYdAdjYd*(Sqr(g1) - 40*Sqr(g3)) +
      16*traceYuAdjYu*Sqr(g3))*TLambdax - 52.5*Sqr(Conj(Lambdax))*Sqr(Lambdax)
      *TLambdax - 0.1*AbsSqr(Lambdax)*(2*Lambdax*(75*traceAdjYdTYd + 25*
      traceAdjYeTYe + 75*traceAdjYuTYu + 6*MassB*Sqr(g1) + 110*MassWB*Sqr(g2))
      - 3*(-75*traceYdAdjYd - 25*traceYeAdjYe - 75*traceYuAdjYu + 6*Sqr(g1) +
      110*Sqr(g2))*TLambdax)));


   return beta_TLambdax;
}

/**
 * Calculates the three-loop beta function of TLambdax.
 *
 * @return three-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_TLambdax_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_TLambdax;

   beta_TLambdax = 0;


   return beta_TLambdax;
}

} // namespace flexiblesusy
