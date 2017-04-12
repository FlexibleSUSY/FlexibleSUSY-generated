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

// File generated at Wed 12 Apr 2017 12:27:56

#include "UMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of MassU.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassU_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);


   double beta_MassU;

   beta_MassU = Re(2*MassU*oneOver16PiSqr*Sqr(gp)*(9*Sqr(Qd) + 3*Sqr(Qe)
      + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu)
      + 3*Sqr(Qv)));


   return beta_MassU;
}

/**
 * Calculates the two-loop beta function of MassU.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassU_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;


   double beta_MassU;

   beta_MassU = Re(0.8*twoLoop*Sqr(gp)*(3*(MassB + MassU)*Sqr(g1)*(2*Sqr(
      Qd) + 6*Sqr(Qe) + Sqr(QHd) + Sqr(QHu) + 3*Sqr(Ql) + Sqr(Qq) + 8*Sqr(Qu))
      + 5*(4*MassU*(9*Power(Qd,4) + 3*Power(Qe,4) + 2*Power(QHd,4) + 2*Power(
      QHu,4) + 6*Power(Ql,4) + 18*Power(Qq,4) + Power(Qs,4) + 9*Power(Qu,4) + 3
      *Power(Qv,4))*Sqr(gp) + 6*traceAdjYdTYd*Sqr(Qd) - 6*MassU*traceYdAdjYd*
      Sqr(Qd) + 2*traceAdjYeTYe*Sqr(Qe) - 2*MassU*traceYeAdjYe*Sqr(Qe) + 6*
      traceAdjYdTYd*Sqr(QHd) + 2*traceAdjYeTYe*Sqr(QHd) - 6*MassU*traceYdAdjYd*
      Sqr(QHd) - 2*MassU*traceYeAdjYe*Sqr(QHd) + 3*MassU*Sqr(g2)*Sqr(QHd) + 3*
      MassWB*Sqr(g2)*Sqr(QHd) + 6*traceAdjYuTYu*Sqr(QHu) + 2*traceAdjYvTYv*Sqr(
      QHu) - 6*MassU*traceYuAdjYu*Sqr(QHu) - 2*MassU*traceYvAdjYv*Sqr(QHu) + 3*
      MassU*Sqr(g2)*Sqr(QHu) + 3*MassWB*Sqr(g2)*Sqr(QHu) + 2*traceAdjYeTYe*Sqr(
      Ql) + 2*traceAdjYvTYv*Sqr(Ql) - 2*MassU*traceYeAdjYe*Sqr(Ql) - 2*MassU*
      traceYvAdjYv*Sqr(Ql) + 9*MassU*Sqr(g2)*Sqr(Ql) + 9*MassWB*Sqr(g2)*Sqr(Ql)
      + 6*traceAdjYdTYd*Sqr(Qq) + 6*traceAdjYuTYu*Sqr(Qq) - 6*MassU*
      traceYdAdjYd*Sqr(Qq) - 6*MassU*traceYuAdjYu*Sqr(Qq) + 27*MassU*Sqr(g2)*
      Sqr(Qq) + 27*MassWB*Sqr(g2)*Sqr(Qq) + 6*traceAdjYuTYu*Sqr(Qu) - 6*MassU*
      traceYuAdjYu*Sqr(Qu) + 24*(MassG + MassU)*Sqr(g3)*(Sqr(Qd) + 2*Sqr(Qq) +
      Sqr(Qu)) + 2*traceAdjYvTYv*Sqr(Qv) - 2*MassU*traceYvAdjYv*Sqr(Qv)) - 10*
      Conj(Lambdax)*(Sqr(QHd) + Sqr(QHu) + Sqr(Qs))*(MassU*Lambdax - TLambdax))
      );


   return beta_MassU;
}

/**
 * Calculates the three-loop beta function of MassU.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_MassU_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassU;

   beta_MassU = 0;


   return beta_MassU;
}

} // namespace flexiblesusy
