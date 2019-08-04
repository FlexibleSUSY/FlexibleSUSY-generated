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

// File generated at Sun 4 Aug 2019 19:19:07

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
 * Calculates the 1-loop beta function of BMT.
 *
 * @return 1-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_1_loop(const Soft_traces& soft_traces) const
{


   double beta_BMT;

   beta_BMT = Re(2*oneOver16PiSqr*(AbsSqr(Lambdax)*BMT + 8*MassWB*MT*Sqr(g2) -
      4*BMT*Sqr(g2) + 2*MT*Conj(Lambdax)*TLambdax));


   return beta_BMT;
}

/**
 * Calculates the 2-loop beta function of BMT.
 *
 * @return 2-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_BMT;

   beta_BMT = Re(-0.4*twoLoop*(30*MT*traceAdjYdTYd*AbsSqr(Lambdax) + 10*MT*
      traceAdjYeTYe*AbsSqr(Lambdax) + 30*MT*traceAdjYuTYu*AbsSqr(Lambdax) + 15*
      traceYdAdjYd*AbsSqr(Lambdax)*BMT + 5*traceYeAdjYe*AbsSqr(Lambdax)*BMT +
      15*traceYuAdjYu*AbsSqr(Lambdax)*BMT + 560*MassWB*MT*Quad(g2) - 140*BMT*
      Quad(g2) + 6*MassB*MT*AbsSqr(Lambdax)*Sqr(g1) - 3*AbsSqr(Lambdax)*BMT*Sqr
      (g1) - 10*MassWB*MT*AbsSqr(Lambdax)*Sqr(g2) + 5*AbsSqr(Lambdax)*BMT*Sqr(
      g2) + 15*BMT*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 30*MT*traceYdAdjYd*Conj(
      Lambdax)*TLambdax + 10*MT*traceYeAdjYe*Conj(Lambdax)*TLambdax + 30*MT*
      traceYuAdjYu*Conj(Lambdax)*TLambdax - 6*MT*Conj(Lambdax)*Sqr(g1)*TLambdax
       + 10*MT*Conj(Lambdax)*Sqr(g2)*TLambdax + 60*MT*Lambdax*Sqr(Conj(Lambdax)
      )*TLambdax));


   return beta_BMT;
}

/**
 * Calculates the 3-loop beta function of BMT.
 *
 * @return 3-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMT;

   beta_BMT = 0;


   return beta_BMT;
}

/**
 * Calculates the 4-loop beta function of BMT.
 *
 * @return 4-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMT;

   beta_BMT = 0;


   return beta_BMT;
}

/**
 * Calculates the 5-loop beta function of BMT.
 *
 * @return 5-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMT;

   beta_BMT = 0;


   return beta_BMT;
}

} // namespace flexiblesusy
