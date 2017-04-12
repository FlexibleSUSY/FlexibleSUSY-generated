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

// File generated at Wed 12 Apr 2017 11:36:29

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
 * Calculates the one-loop beta function of BMT.
 *
 * @return one-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_one_loop(const Soft_traces& soft_traces) const
{


   double beta_BMT;

   beta_BMT = Re(2*oneOver16PiSqr*(AbsSqr(Lambdax)*BMT + 8*MassWB*MT*Sqr(
      g2) - 4*BMT*Sqr(g2) + 2*MT*Conj(Lambdax)*TLambdax));


   return beta_BMT;
}

/**
 * Calculates the two-loop beta function of BMT.
 *
 * @return two-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_BMT;

   beta_BMT = Re(-0.4*twoLoop*(BMT*(-140*Power(g2,4) + AbsSqr(Lambdax)*(
      -3*Sqr(g1) + 5*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + Sqr(g2))
      ) + 15*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 2*MT*(280*Power(g2,4)*MassWB +
      30*Lambdax*Sqr(Conj(Lambdax))*TLambdax + Conj(Lambdax)*(Lambdax*(3*MassB*
      Sqr(g1) + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 3*traceAdjYuTYu - MassWB*
      Sqr(g2))) + (-3*Sqr(g1) + 5*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + Sqr(g2)))*TLambdax))));


   return beta_BMT;
}

/**
 * Calculates the three-loop beta function of BMT.
 *
 * @return three-loop beta function
 */
double TMSSM_soft_parameters::calc_beta_BMT_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMT;

   beta_BMT = 0;


   return beta_BMT;
}

} // namespace flexiblesusy
