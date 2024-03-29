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


#include "CE6SSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of MassWB.
 *
 * @return 1-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_MassWB_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassWB;

   beta_MassWB = Re(8*MassWB*Sqr(g2));


   return oneLoop * beta_MassWB;
}

/**
 * Calculates the 2-loop beta function of MassWB.
 *
 * @return 2-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_MassWB_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjLambda12TLambda12 = TRACE_STRUCT.
      traceAdjLambda12TLambda12;


   double beta_MassWB;

   beta_MassWB = Re(0.4*Sqr(g2)*(10*traceAdjLambda12TLambda12 + 30*
      traceAdjYdTYd + 10*traceAdjYeTYe + 30*traceAdjYuTYu - 10*MassWB*
      traceLambda12AdjLambda12 - 30*MassWB*traceYdAdjYd - 10*MassWB*
      traceYeAdjYe - 30*MassWB*traceYuAdjYu - 10*MassWB*AbsSqr(Lambdax) + 18*
      MassB*Sqr(g1) + 18*MassWB*Sqr(g1) + 460*MassWB*Sqr(g2) + 120*MassG*Sqr(g3
      ) + 120*MassWB*Sqr(g3) + 17*MassBp*Sqr(gN) + 17*MassWB*Sqr(gN) + 10*Conj(
      Lambdax)*TLambdax));


   return twoLoop * beta_MassWB;
}

/**
 * Calculates the 3-loop beta function of MassWB.
 *
 * @return 3-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_MassWB_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return threeLoop * beta_MassWB;
}

/**
 * Calculates the 4-loop beta function of MassWB.
 *
 * @return 4-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_MassWB_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return fourLoop * beta_MassWB;
}

/**
 * Calculates the 5-loop beta function of MassWB.
 *
 * @return 5-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_MassWB_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return fiveLoop * beta_MassWB;
}

} // namespace flexiblesusy
