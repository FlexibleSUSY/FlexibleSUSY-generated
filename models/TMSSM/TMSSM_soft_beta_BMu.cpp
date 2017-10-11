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

// File generated at Tue 10 Oct 2017 21:35:00

#include "TMSSM_soft_parameters.hpp"
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
double TMSSM_soft_parameters::calc_beta_BMu_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_BMu;

   beta_BMu = Re(oneOver16PiSqr*(6*traceAdjYdTYd*Mu + 2*traceAdjYeTYe*Mu
      + 6*traceAdjYuTYu*Mu + 1.2*MassB*Mu*Sqr(g1) + BMu*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + 3*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(
      g2)) + 6*MassWB*Mu*Sqr(g2) + 6*Conj(Lambdax)*Mu*TLambdax));


   return beta_BMu;
}

/**
 * Calculates the 2-loop beta function of BMu.
 *
 * @return 2-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMu_2_loop(const Soft_traces& soft_traces) const
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

   beta_BMu = Re(0.02*twoLoop*(BMu*(207*Quad(g1) + 75*AbsSqr(Lambdax)*(-3
      *traceYdAdjYd - traceYeAdjYe - 3*traceYuAdjYu + 8*Sqr(g2)) + 10*Sqr(g1)*(
      -2*traceYdAdjYd + 6*traceYeAdjYe + 4*traceYuAdjYu + 9*Sqr(g2)) + 25*(27*
      Quad(g2) + 2*(-3*(3*traceYdAdjYdYdAdjYd + 2*traceYdAdjYuYuAdjYd +
      traceYeAdjYeYeAdjYe + 3*traceYuAdjYuYuAdjYu) + 16*(traceYdAdjYd +
      traceYuAdjYu)*Sqr(g3))) - 375*Sqr(Conj(Lambdax))*Sqr(Lambdax)) - 2*Mu*(
      414*MassB*Quad(g1) + 10*Sqr(g1)*(2*(traceAdjYdTYd - 3*traceAdjYeTYe - 2*
      traceAdjYuTYu - MassB*traceYdAdjYd + 3*MassB*traceYeAdjYe + 2*MassB*
      traceYuAdjYu) + 9*(MassB + MassWB)*Sqr(g2)) + 50*(27*MassWB*Quad(g2) + 2*
      (3*(3*traceYdAdjYdTYdAdjYd + traceYdAdjYuTYuAdjYd + traceYeAdjYeTYeAdjYe
      + traceYuAdjYdTYdAdjYu + 3*traceYuAdjYuTYuAdjYu) - 8*(traceAdjYdTYd +
      traceAdjYuTYu - MassG*(traceYdAdjYd + traceYuAdjYu))*Sqr(g3))) + 750*
      Lambdax*Sqr(Conj(Lambdax))*TLambdax + 75*Conj(Lambdax)*(Lambdax*(3*
      traceAdjYdTYd + traceAdjYeTYe + 3*traceAdjYuTYu + 8*MassWB*Sqr(g2)) + (3*
      traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu - 8*Sqr(g2))*TLambdax))));


   return beta_BMu;
}

/**
 * Calculates the 3-loop beta function of BMu.
 *
 * @return 3-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMu;

   beta_BMu = 0;


   return beta_BMu;
}

} // namespace flexiblesusy